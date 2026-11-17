#!/usr/bin/env python3
"""Prepare a short run from the supplied, audited single-semaglutide asset.

No downloads, online parameter generators, or changes to the existing environment.
The prebuilt box contains neutralizing Na+ plus approximately 0.150 M NaCl.
Its topology is a documented initial model, not a claim of validated energetics.
"""
import argparse
import hashlib
import json
import math
import shutil
from pathlib import Path


def positive(value):
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise argparse.ArgumentTypeError('must be finite and positive')
    return value


def write_mdp(folder, stage, args, restraints):
    common = dict(cutoff_scheme='Verlet', nstlist=20, pbc='xyz',
                  coulombtype='PME', rcoulomb=1.0, rvdw=1.0,
                  vdwtype='Cut-off', DispCorr='EnerPres',
                  fourierspacing=0.12, pme_order=4, ewald_rtol=1e-5,
                  constraints='h-bonds', constraint_algorithm='lincs',
                  lincs_order=4, lincs_iter=1)
    if stage == 'em':
        common.update(integrator='steep', nsteps=20000, emtol=1000,
                      emstep=0.01, nstenergy=100, nstlog=100)
    else:
        duration_ps = {'nvt': args.nvt_ps, 'npt': args.npt_ps,
                       'production': args.production_ns*1000}[stage]
        nsteps = max(1, int(round(duration_ps/.002)))
        common.update(integrator='md', dt=.002, nsteps=nsteps,
                      continuation='no' if stage == 'nvt' else 'yes',
                      tcoupl='v-rescale', tc_grps='System', tau_t=.5,
                      ref_t=args.temperature, nstxout=0, nstvout=0,
                      nstfout=0, nstxout_compressed=1000,
                      compressed_x_precision=1000, nstenergy=1000,
                      nstlog=1000, comm_mode='Linear', nstcomm=100)
        if stage == 'nvt':
            common.update(gen_vel='yes', gen_temp=args.temperature,
                          gen_seed=args.seed, pcoupl='no')
        else:
            common.update(gen_vel='no', pcoupl='c-rescale',
                          pcoupltype='isotropic', tau_p=2.0, ref_p=1.0,
                          compressibility=4.5e-5, refcoord_scaling='com')
        if stage in ('nvt', 'npt') and restraints:
            common['define'] = '-DSEMA_POSRES'
    header = ('; Single semaglutide: short numerical/structure screen.\n'
              '; Fixed protonation state; NOT constant-pH MD.\n'
              '; Production has no positional restraints.\n')
    (folder/(stage+'.mdp')).write_text(header + ''.join(
        '%-25s = %s\n' % (k.replace('_','-'),v) for k,v in common.items()))


def prepare(args):
    assets = args.assets.resolve()
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError('Output is not empty: %s. Start a fresh run.' % output)
    checks = json.loads((assets/'checksums.json').read_text())
    if not checks:
        raise ValueError('Empty asset checksum list')
    for relative, expected in checks.items():
        path = (assets/relative).resolve()
        if assets not in path.parents or not path.is_file():
            raise ValueError('Invalid/missing asset: '+relative)
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError('Asset checksum mismatch: '+relative)
    manifest = json.loads((assets/'manifest.json').read_text())
    if manifest.get('chemical_audit', {}).get('status') != 'pass':
        raise ValueError('Supplied molecule has not passed the chemical identity audit')
    if not manifest.get('parameter_audit', {}).get('coverage_complete'):
        raise ValueError('Incomplete topology parameters; refusing simulation')
    output.mkdir(parents=True, exist_ok=True)
    for relative in checks:
        if relative == 'manifest.json':
            continue
        destination = output/relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(assets/relative, destination)
    manifest.update(temperature_K=args.temperature, production_ns=args.production_ns,
                    nvt_ps=args.nvt_ps, npt_ps=args.npt_ps, seed=args.seed,
                    cutoff_nm=1.0, source_asset_checksums=checks)
    manifest['requested'] = dict(production_ns=args.production_ns,
                                nvt_ps=args.nvt_ps, npt_ps=args.npt_ps,
                                temperature_K=args.temperature, seed=args.seed)
    manifest['scope'] = ('One fixed-protonation semaglutide in OPC water and ions. '
                         'Short implementation screen; parameter energetics and '
                         'experimental agreement remain unvalidated.')
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    mdp = output/'mdp'
    mdp.mkdir(exist_ok=True)
    topology = (output/'topol.top').read_text()
    restraints = 'SEMA_POSRES' in topology or bool(manifest.get('position_restraints'))
    for stage in ('em','nvt','npt','production'):
        write_mdp(mdp,stage,args,restraints)
    print(json.dumps({
        'solute': 'semaglutide', 'peptide_residues':31,
        'solute_atom_count':len(manifest['atoms']),
        'total_simulation_sites':manifest.get('system_atom_count'),
        'water_molecules':manifest.get('water_molecules'),
        'ions':manifest.get('ions'), 'box_volume_nm3':manifest.get('box_volume_nm3'),
        'force_field':manifest.get('provenance',{}).get('force_field'),
        'production_ns':args.production_ns, 'nvt_ps':args.nvt_ps, 'npt_ps':args.npt_ps,
        'physical_force_field_validation':'NOT established by this short run'
    },indent=2),flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--assets',required=True,type=Path)
    parser.add_argument('--gmx',required=True,help='GROMACS executable, recorded by the launcher')
    parser.add_argument('--production-ns',type=positive,default=1.)
    parser.add_argument('--nvt-ps',type=positive,default=50.)
    parser.add_argument('--npt-ps',type=positive,default=200.)
    parser.add_argument('--temperature',type=positive,default=298.15)
    parser.add_argument('--seed',type=int,default=1729)
    args = parser.parse_args()
    if not 0 < args.seed < 2147483647:
        parser.error('seed must be in 1..2147483646')
    prepare(args)


if __name__ == '__main__':
    main()
