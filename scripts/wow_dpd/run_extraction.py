#!/usr/bin/env python3
"""Relax the closed W/O/W particle, then exchange external DCM until the sampled in-particle target."""
import argparse
import csv
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import time
from types import SimpleNamespace

import numpy as np
from scipy.spatial import QhullError

from engine import frame, quote, static_commands
from particle_series import measure as particle_composition
from solvent_sink import eligible_ids, exchange_ids

HERE = Path(__file__).resolve().parent
STREAMS = ('sink.csv', 'particle.csv', 'trajectory.lammpstrj')
STATE_VERSION = 1


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--initial', type=Path, default=HERE/'initial')
    p.add_argument('--interactions', type=Path, default=HERE/'interactions.json')
    p.add_argument('--out', type=Path, help='Persistent directory, resumed automatically if it has state.json')
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--seed', type=int, default=48727)
    p.add_argument('--dt', type=float, default=0.01)
    p.add_argument('--chunk', type=int, default=100, help='Integration steps per sink decision')
    p.add_argument('--sample-every', type=int, default=1000, help='Steps between particle measurements and checkpoints')
    p.add_argument('--hold-steps', type=int, default=1000, help='Initial steps with a closed sink')
    p.add_argument('--max-steps', type=int, default=2000000, help='Cumulative extraction-step cap, excluding closed hold')
    p.add_argument('--target', type=float, default=0.05, help='DCM/(DCM+PLGA) mass fraction within convex envelope')
    p.add_argument('--evap-rate', type=float, default=0.02, help='First-order exchange rate per reduced time unit')
    p.add_argument('--time-budget', type=float, help='Pause after this many wall seconds; allow checkpoint grace in allocation')
    p.add_argument('--check-input', action='store_true', help='Check LAMMPS forces, sink, and particle composition at step zero')
    a = p.parse_args()
    if not a.check_input and not a.out:
        p.error('--out is required for production')
    if min(a.threads, a.chunk, a.sample_every, a.max_steps) < 1 or a.hold_steps < 0:
        p.error('Invalid step or thread count')
    if not (0 < a.seed < 899999999 and 0 < a.target < 1 and math.isfinite(a.dt) and a.dt > 0
            and math.isfinite(a.evap_rate) and a.evap_rate >= 0):
        p.error('Invalid seed, target, timestep or rate')
    if a.time_budget is not None and (not math.isfinite(a.time_budget) or a.time_budget <= 0):
        p.error('Wall budget must be a positive number of seconds')
    if a.sample_every % a.chunk:
        p.error('--sample-every must be a multiple of --chunk, so measurements coincide with exchange decisions')
    if max(int(os.environ.get(k, '1')) for k in ('OMPI_COMM_WORLD_SIZE','PMI_SIZE','MV2_COMM_WORLD_SIZE','SLURM_NTASKS')) > 1:
        p.error('Start one Python process in one MPI task. Use OpenMP CPU threads within that task.')
    return a


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for data in iter(lambda: stream.read(1024*1024), b''):
            h.update(data)
    return h.hexdigest()


def model_config(a):
    paths = [a.initial/n for n in ('system.data','geometry.npz','layout.json')] + [a.interactions]
    return dict(input_sha256=[digest(p) for p in paths], **{k:getattr(a,k) for k in (
        'threads','seed','dt','chunk','sample_every','hold_steps','target','evap_rate')})


def sync_dir(folder):
    fd = os.open(folder, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def atomic_json(path, value):
    tmp = path.with_name(path.name+'.tmp')
    with tmp.open('w') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)
    sync_dir(path.parent)


def initial_state(config):
    return dict(version=STATE_VERSION, config=config, step=0, extraction_step=0,
                last_sample=None, last_sample_step=None, exported_DCM=0, checkpoint=None,
                offsets={name:0 for name in STREAMS}, status='initializing')


def recover(folder, config):
    statefile = folder/'state.json'
    if not statefile.exists():
        if any((folder/n).exists() for n in (*STREAMS,'restart.A','restart.B')):
            raise RuntimeError('No state.json for existing run output. Use a new directory or restore the manifest.')
        state = initial_state(config)
        atomic_json(statefile, state)
        return state
    state = json.loads(statefile.read_text())
    if state['version'] != STATE_VERSION or config != state['config']:
        raise ValueError('Input files or physical/sampling settings differ from the saved run in state.json')
    checkpoint = state['checkpoint']
    if checkpoint and (checkpoint not in ('restart.A','restart.B') or not (folder/checkpoint).is_file()):
        raise RuntimeError('Committed restart is missing')
    # A previous job may have appended rows or a partial dump after the last
    # atomic commit. Validate all stream lengths before truncating any tail.
    for name in STREAMS:
        path = folder/name
        size = path.stat().st_size if path.exists() else 0
        if size < state['offsets'][name]:
            raise RuntimeError(f'Output is shorter than the committed offset: {name}')
    for name in STREAMS:
        path = folder/name
        if path.exists() and path.stat().st_size > state['offsets'][name]:
            with path.open('r+b') as f:
                f.truncate(state['offsets'][name])
    return state


def commit(lmp, folder, state, status):
    # On interruption the old manifest and old alternate restart still agree.
    dest = 'restart.B' if state['checkpoint'] == 'restart.A' else 'restart.A'
    restart = folder/dest
    lmp.command(f'write_restart {quote(restart)}')
    with restart.open('rb') as f:
        os.fsync(f.fileno())
    offsets = {}
    for name in STREAMS:
        path = folder/name
        with path.open('rb') as f:
            os.fsync(f.fileno())
            offsets[name] = os.fstat(f.fileno()).st_size
    new = dict(state, checkpoint=dest, offsets=offsets, status=status)
    atomic_json(folder/'state.json', new)
    state.update(new)


def append_csv(path, row):
    fresh = not path.exists() or path.stat().st_size == 0
    with path.open('a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(row))
        if fresh:
            writer.writeheader()
        writer.writerow(row)


def stock(frame_types):
    return (int(np.count_nonzero(frame_types == 2)), int(np.count_nonzero(frame_types == 1)))


def sink_step(lmp, a, state, nsteps, phase, initial_inventory, job_start):
    ids, types, xyz, length = frame(lmp)
    sink_error = ''
    try:
        eligible, diag = eligible_ids(ids, types, xyz, length)
    except ValueError as exc:
        if str(exc) != 'Unable to center the droplet on its polymer distribution' or phase == 'initial':
            raise
        # Without an identifiable exterior the flux is undefined. Exchange
        # nothing in this interval and retry classification at the next one.
        sink_error = str(exc)
        eligible = np.empty(0, dtype=ids.dtype)
        diag = dict.fromkeys(('polymer_center_x_rc','polymer_center_y_rc','polymer_center_z_rc',
                              'polymer_max_radius_rc','outside_radial_guard','eligible'))
    rng = np.random.default_rng(np.random.SeedSequence([a.seed,state['step']]))
    selected = exchange_ids(eligible, a.evap_rate if phase == 'extraction' else 0.0, nsteps*a.dt, rng)
    if len(selected):
        mask = np.isin(ids, selected)
        if mask.sum() != len(selected) or (types[mask] != 2).any():
            raise RuntimeError('Sink attempted to exchange an invalid bead')
        types[mask] = 1
    state['exported_DCM'] += len(selected)
    n_d, n_w = stock(types)
    if (n_d+state['exported_DCM'], n_w-state['exported_DCM']) != initial_inventory:
        raise RuntimeError('Global water/DCM exchange ledger disagrees with bead counts')
    return dict(step=state['step'], phase=phase, extraction_step=state['extraction_step'],
                interval_steps=nsteps, job_wall_seconds=round(time.monotonic()-job_start, 3),
                exported_this_interval=len(selected), exported_DCM=state['exported_DCM'],
                remaining_DCM=n_d, water_beads=n_w, sink_error=sink_error, **diag)


def sample_particle(lmp, step, layout, orig, pockets):
    ids, types, xyz, length = frame(lmp)
    if len(ids) != len(orig) or np.min(ids) != 1 or np.max(ids) != len(ids):
        raise RuntimeError('LAMMPS bead IDs disagree with the reference geometry')
    order = np.argsort(ids)
    if not np.all(np.diff(ids[order]) == 1):
        raise RuntimeError('LAMMPS bead IDs are not unique and contiguous')
    rc = layout['rc_nm']
    return particle_composition(step, types[order], xyz[order]*rc, length*rc, orig, pockets,
                                layout, 0.5, 0.7, 1.4, 2.1)


def new_engine(a, state, interactions, log):
    from lammps import lammps
    switches = ['-screen','none','-log',str(log)]
    if a.threads > 1:
        switches += ['-sf','omp','-pk','omp',str(a.threads)]
    lmp = lammps(cmdargs=switches)
    try:
        if a.threads > 1 and not lmp.has_style('pair', 'dpd/omp'):
            raise RuntimeError('LAMMPS needs OPENMP and dpd/omp for the requested CPU threads')
        restart = state['checkpoint']
        cfg = SimpleNamespace(data=a.initial/'system.data',restart=a.out/restart if restart else None,
                              seed=a.seed, dt=a.dt)
        if restart:
            # PairDPD does not save RNG state in a binary restart. Generate a
            # fresh stream per absolute checkpoint instead of recycling step-zero noise.
            if state['step']:
                seed32 = np.random.SeedSequence([a.seed,state['step'],719]).generate_state(1)[0]
                cfg.seed = 1 + int(seed32) % 899999998
        lmp.commands_list(static_commands(cfg, interactions))
        if int(lmp.extract_global('ntimestep')) != state['step']:
            raise RuntimeError('Binary restart step disagrees with state.json')
        if not restart:
            lmp.command(f'velocity all create 1.0 {a.seed+1} mom yes rot no dist gaussian')
        lmp.command(f'thermo {a.chunk}')
        lmp.command('run 0 post no')
        return lmp
    except BaseException:
        lmp.close()
        raise


def observe(lmp, a, state, layout, orig, pockets, phase):
    try:
        row = sample_particle(lmp, state['step'], layout, orig, pockets)
        err = ''
    except (ValueError, QhullError) as exc:
        if not state.get('particle_fields'):
            raise   # The assembled starting geometry must be measurable.
        # A transient channel or fragmented carrier can defeat a single-particle
        # boundary. Keep its frame and ledger; a missing measurement cannot stop extraction.
        err = ' '.join(str(exc).split())[:500]
        row = dict.fromkeys(state['particle_fields'])
        row['step'] = state['step']
    state['particle_fields'] = list(row)
    append_csv(a.out/'particle.csv', dict(phase=phase, extraction_step=state['extraction_step'],
                                          analysis_error=err, **row))
    lmp.command(f'write_dump all custom {quote(a.out/"trajectory.lammpstrj")} id mol type x y z modify sort id append yes')
    state['last_sample_step'] = state['step']
    if err:
        state['last_sample'] = None
        print(f"Step {state['step']} ({phase}): particle boundary unavailable: {err}", flush=True)
        return
    state['last_sample'] = {k:float(row[k]) for k in (
        'dcm_wt_fraction_particle','probe_dcm_wt_fraction_particle','phi_water_all_sites',
        'probe_phi_water_all_sites','phi_water_additive_carrier','polymer_sites_inside_fraction')}
    print(f"Step {state['step']} ({phase}): local DCM {100*row['dcm_wt_fraction_particle']:.3f} wt%, "
          f"internal water {row['phi_water_all_sites']:.4f} by bead volume", flush=True)


def run(a, request, job_start, state, layout, orig, pockets, interactions):
    cap = a.hold_steps+max(a.max_steps, state.get('max_steps', 0))
    if not a.check_input and (state['status'] == 'target_reached' or (state['checkpoint'] and state['step'] >= cap)):
        print(f"Stored status {state['status']}, absolute step {state['step']}. "
              "For a capped incomplete run, increase --max-steps to continue.", flush=True)
        return
    log = 'none' if a.check_input else (a.out/f'log.lammps.{time.time_ns()}').resolve()
    lmp = new_engine(a, state, interactions, log)
    try:
        initial_inventory = stock(orig)
        if a.check_input:
            row = sample_particle(lmp, 0, layout, orig, pockets)
            diag = sink_step(lmp, a, state, 0, 'initial', initial_inventory, job_start)
            print(f"Static check OK: {len(orig)} beads; step 0; eligible exterior DCM {diag['eligible']}; "
                  f"isolated W1 domains {row['isolated_occupied_wet_domains']}; "
                  f"local DCM {100*row['dcm_wt_fraction_particle']:.3f} wt%.")
            return
        if not state['checkpoint']:
            diag = sink_step(lmp, a, state, 0, 'initial', initial_inventory, job_start)
            append_csv(a.out/'sink.csv', diag)
            observe(lmp, a, state, layout, orig, pockets, 'initial')
            already_at_target = a.hold_steps == 0 and state['last_sample']['dcm_wt_fraction_particle'] <= a.target
            commit(lmp, a.out, state, 'target_reached' if already_at_target else 'running')
            if already_at_target:
                print("The initial open-bath frame already meets the requested mass fraction.", flush=True)
                return
        while True:
            if request or (a.time_budget is not None and time.monotonic()-job_start >= a.time_budget):
                commit(lmp, a.out, state, 'paused')
                print(f"Paused at step {state['step']}. Resubmit with the same output directory.", flush=True)
                return
            old = state['step']
            in_hold = old < a.hold_steps
            phase = 'closed_hold' if in_hold else 'extraction'
            relative = old if in_hold else old-a.hold_steps
            boundary = a.hold_steps if in_hold else cap
            steps = min(boundary-old, a.chunk-relative%a.chunk)
            lmp.command(f'run {steps} post no')
            state['step'] = old+steps
            if int(lmp.extract_global('ntimestep')) != state['step']:
                raise RuntimeError('LAMMPS step does not match integration request')
            state['extraction_step'] = max(0,state['step']-a.hold_steps)
            diag = sink_step(lmp, a, state, steps, phase, initial_inventory, job_start)
            append_csv(a.out/'sink.csv', diag)
            if (relative+steps)%a.sample_every == 0 or state['step'] == boundary:
                observe(lmp, a, state, layout, orig, pockets, phase)
                sample = state['last_sample']
                reached = sample is not None and sample['dcm_wt_fraction_particle'] <= a.target
                status = ('target_reached' if state['step'] >= a.hold_steps and reached else
                          'max_steps' if state['step'] == cap else 'running')
                commit(lmp, a.out, state, status)
                if status != 'running':
                    print(f"Stored status {status} at step {state['step']}.", flush=True)
                    return
    finally:
        lmp.close()


def main():
    job_start = time.monotonic()
    a = arguments()
    request = []
    for sig in (signal.SIGTERM, signal.SIGUSR1):
        signal.signal(sig, lambda sig, frame: request.append(sig))
    layout = json.loads((a.initial/'layout.json').read_text())
    interactions = json.loads(a.interactions.read_text())
    with np.load(a.initial/'geometry.npz', allow_pickle=False) as g:
        orig, pockets = g['types'], g['pocket']
    if a.check_input:
        run(a,request,job_start,initial_state({}),layout,orig,pockets,interactions)
        return
    a.out.mkdir(parents=True, exist_ok=True)
    with (a.out/'.writer.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit('Another extraction job is writing this run directory.')
        state = recover(a.out, model_config(a))
        if a.max_steps > state.get('max_steps', 0):
            state['max_steps'] = a.max_steps
            atomic_json(a.out/'state.json', state)
        run(a,request,job_start,state,layout,orig,pockets,interactions)


if __name__ == '__main__':
    main()
