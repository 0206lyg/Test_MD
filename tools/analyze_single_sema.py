#!/usr/bin/env python3
"""Inspect a short single-semaglutide run; never equate it to force-field validation.

Requires Python >= 3.8 and NumPy, both in Test_MD's preparation environment.
GROMACS exports whole solute coordinates, so Rg and internal distances do not
depend on periodic wrapping. No trajectory is fitted before the image-distance
check. The dihedral metric is explicitly a helical-region proxy, not DSSP.
"""
import argparse
import csv
import itertools
import json
import math
import re
import subprocess
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np


def read_json(path, default=None):
    return json.loads(path.read_text()) if path.is_file() else default


def gromacs_export(run, analysis, gmx):
    """Export the latest available trajectory; preserve diagnostic failures."""
    messages = []
    stage = next((s for s in ('production', 'npt', 'nvt')
                  if (run / (s + '.xtc')).is_file()
                  and (run / (s + '.tpr')).is_file()), None)
    if stage is None:
        return None, ['No MD trajectory is available yet.']
    commands = [
        ('trjconv', [gmx, 'trjconv', '-f', str(run / (stage + '.xtc')),
                     '-s', str(run / (stage + '.tpr')),
                     '-n', str(run / 'input/index.ndx'),
                     '-o', str(analysis / 'frames.gro'), '-pbc', 'whole',
                     '-dt', '10'], 'Semaglutide\n'),
        ('solute_xtc', [gmx, 'trjconv', '-f', str(run / (stage + '.xtc')),
                        '-s', str(run / (stage + '.tpr')),
                        '-n', str(run / 'input/index.ndx'),
                        '-o', str(analysis / 'semaglutide.xtc'),
                        '-pbc', 'whole'], 'Semaglutide\n'),
        ('solute_gro', [gmx, 'trjconv', '-f', str(run / (stage + '.xtc')),
                        '-s', str(run / (stage + '.tpr')),
                        '-n', str(run / 'input/index.ndx'),
                        '-o', str(analysis / 'semaglutide.gro'),
                        '-dump', '0', '-pbc', 'whole'], 'Semaglutide\n'),
    ]
    if (run / (stage + '.edr')).is_file():
        # NVT has no barostat density term in some GROMACS outputs.
        selections = ('Temperature\nPressure\nPotential\n0\n' if stage == 'nvt'
                      else 'Temperature\nPressure\nDensity\nPotential\n0\n')
        commands.append(('energy', [gmx, 'energy', '-f', str(run / (stage + '.edr')),
                                   '-o', str(analysis / 'thermodynamics.xvg')], selections))
    for label, command, selection in commands:
        with (analysis / (label + '.log')).open('w') as log:
            log.write('COMMAND: ' + ' '.join(command) + '\n')
            log.flush()
            try:
                result = subprocess.run(command, input=selection, text=True,
                                        stdout=log, stderr=subprocess.STDOUT, cwd=run)
            except OSError as error:
                messages.append('%s could not run: %s' % (label, error))
                continue
        if result.returncode:
            messages.append('%s failed; inspect analysis/%s.log.' % (label, label))
    return stage, messages


def gro_frames(path):
    """Yield time_ps, Nx3 coordinates_nm, row-vector triclinic box_nm."""
    with path.open() as handle:
        frame = 0
        while True:
            title = handle.readline()
            if not title:
                break
            line = handle.readline()
            if not line:
                raise ValueError('Truncated GRO atom count')
            count = int(line.strip())
            xyz = np.empty((count, 3), dtype=float)
            for i in range(count):
                line = handle.readline()
                if len(line) < 44:
                    raise ValueError('Truncated GRO atom record')
                xyz[i] = [float(line[20:28]), float(line[28:36]), float(line[36:44])]
            fields = [float(v) for v in handle.readline().split()]
            if len(fields) == 3:
                box = np.diag(fields)
            elif len(fields) == 9:
                box = np.array([[fields[0], fields[3], fields[4]],
                                [fields[5], fields[1], fields[6]],
                                [fields[7], fields[8], fields[2]]])
            else:
                raise ValueError('GRO box must contain 3 or 9 numbers')
            match = re.search(r'\bt\s*=\s*([-+0-9.eE]+)', title)
            time = float(match.group(1)) if match else float(frame * 10)
            if not np.isfinite(xyz).all() or not np.isfinite(box).all():
                raise ValueError('Non-finite GRO coordinate or box')
            if np.linalg.det(box) <= 0:
                raise ValueError('Invalid/non-positive GRO box volume')
            yield time, xyz, box
            frame += 1


def kabsch_rmsd(xyz, reference):
    x = np.asarray(xyz, dtype=float) - np.mean(xyz, axis=0)
    y = np.asarray(reference, dtype=float) - np.mean(reference, axis=0)
    u, _, vt = np.linalg.svd(x.T @ y)
    correction = np.diag([1., 1., np.linalg.det(u @ vt)])
    rotation = u @ correction @ vt
    return float(np.sqrt(np.mean(np.sum((x @ rotation - y) ** 2, axis=1))))


def radius_of_gyration(xyz, masses):
    center = np.average(xyz, weights=masses, axis=0)
    return float(np.sqrt(np.average(np.sum((xyz - center) ** 2, axis=1), weights=masses)))


def dihedral(a, b, c, d):
    b0, b1, b2 = -(b-a), c-b, d-c
    norm = np.linalg.norm(b1)
    if norm < 1e-12:
        return float('nan')
    b1 = b1 / norm
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    if min(np.linalg.norm(v), np.linalg.norm(w)) < 1e-12:
        return float('nan')
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def minimum_image_copy_distance(xyz, box):
    """Minimum solute-to-copy distance for a reduced GROMACS triclinic cell.

    Uses all heavy-atom pairs and the 27 lattice points surrounding each pair's
    fractional nearest-integer position. Excludes only the ORIGINAL molecule
    (zero translation), not i==j in another periodic copy. This also catches a
    solute longer than one cell: translations are not limited to +/-1 globally.
    """
    inverse = np.linalg.inv(box)
    minimum = float('inf')
    offsets = np.array(list(itertools.product((-1, 0, 1), repeat=3)), dtype=int)
    for start in range(0, len(xyz), 48):
        delta = (xyz[start:start+48, None, :] - xyz[None, :, :]).reshape(-1, 3)
        nearest = np.rint(delta @ inverse).astype(int)
        for offset in offsets:
            lattice = nearest + offset
            keep = np.any(lattice != 0, axis=1)
            if not keep.any():
                continue
            distances = delta[keep] - lattice[keep] @ box
            d2 = np.einsum('ij,ij->i', distances, distances)
            minimum = min(minimum, float(d2.min()))
    return math.sqrt(minimum)


def safe_indices(values, count):
    result = [int(v) for v in values if v is not None]
    if any(i < 0 or i >= count for i in result):
        raise ValueError('Manifest atom index is outside solute atom range')
    return result


def helix_values(xyz, backbone):
    values = {}
    for i in range(1, len(backbone)-1):
        prev, cur, nxt = backbone[i-1:i+2]
        indices = [prev.get('C'), cur.get('N'), cur.get('CA'), cur.get('C'), nxt.get('N')]
        if any(v is None for v in indices):
            continue
        p, n, ca, c, nn = [xyz[int(v)] for v in indices]
        phi, psi = dihedral(p, n, ca, c), dihedral(n, ca, c, nn)
        if np.isfinite(phi) and np.isfinite(psi):
            values[str(cur['resnr'])] = int(-100 <= phi <= -30 and -80 <= psi <= -5)
    return values


def stereocenter_retention(xyz, centers):
    """Check alpha/side-chain centers with three CIP-ordered heavy neighbors.

    The manifest orders the three heavy neighbors by descending CIP priority;
    the omitted hydrogen is the lowest-priority substituent. For this convention
    a positive determinant is S and a negative determinant is R. This monitors
    geometric retention; it does not assign CIP priorities from coordinates.
    """
    failures, volumes = [], []
    for item in centers:
        center = int(item['center'])
        neighbors = [int(i) for i in item['neighbors']]
        if len(neighbors) != 3 or len(set(neighbors + [center])) != 4:
            raise ValueError('A stereocenter requires three distinct ordered neighbors')
        if min(neighbors + [center]) < 0 or max(neighbors + [center]) >= len(xyz):
            raise ValueError('Stereocenter index outside solute atom range')
        expected = item['expected']
        if expected not in ('S', 'R'):
            raise ValueError('Expected stereocenter descriptor must be S or R')
        volume = float(np.linalg.det(xyz[neighbors] - xyz[center]))
        volumes.append(abs(volume))
        observed = 'S' if volume > 1e-9 else 'R' if volume < -1e-9 else 'degenerate'
        if observed != expected:
            failures.append({'center': center, 'expected': expected, 'observed': observed,
                             'signed_volume_nm3': volume})
    return failures, min(volumes) if volumes else None


def xvg_summary(path):
    if not path.is_file():
        return {}
    legends, rows = {}, []
    for line in path.read_text(errors='replace').splitlines():
        match = re.match(r'@\s+s(\d+)\s+legend\s+"(.*?)"', line)
        if match:
            legends[int(match.group(1))] = match.group(2)
        if not line.strip() or line.lstrip().startswith(('#', '@')):
            continue
        try:
            rows.append([float(x) for x in line.split()])
        except ValueError:
            continue
    if not rows:
        return {}
    arr = np.asarray(rows)
    result = {}
    for i in range(1, arr.shape[1]):
        name = legends.get(i-1, 'column_%d' % i)
        values = arr[:, i]
        def finite_mean(v):
            good = v[np.isfinite(v)]
            return float(np.mean(good)) if len(good) else None
        good = values[np.isfinite(values)]
        result[name] = {'mean': finite_mean(values), 'std': float(np.std(good)) if len(good) else None,
                        'first_half_mean': finite_mean(values[:max(1, len(values)//2)]),
                        'second_half_mean': finite_mean(values[len(values)//2:]),
                        'finite': bool(np.isfinite(values).all()), 'n': len(values)}
    return result


def mean_stats(values):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if not len(arr):
        return None
    return {'mean': float(arr.mean()), 'min': float(arr.min()), 'max': float(arr.max()),
            'std': float(arr.std()), 'n': len(arr)}


def write_plot(path, rows):
    """Small dependency-free SVG; axes and units are explicit."""
    if len(rows) < 2:
        return
    panels = [('backbone_rmsd_nm', 'Backbone RMSD to first sampled frame (nm)'),
              ('rg_nm', 'Whole-solute mass-weighted Rg (nm)'),
              ('tail_peptide_contacts', 'C18 / peptide heavy-atom contacts < 0.45 nm'),
              ('alpha_region_fraction', 'Backbone alpha-region fraction (proxy, not DSSP)')]
    width, height = 900, 890
    content = ['<svg xmlns="http://www.w3.org/2000/svg" width="900" height="890" viewBox="0 0 900 890">',
               '<rect width="900" height="890" fill="white"/>',
               '<style>text{font-family:Arial,sans-serif;fill:#263445;font-size:13px}</style>',
               '<text x="55" y="27" style="font-size:19px;font-weight:bold">Single semaglutide: short-run diagnostics</text>',
               '<text x="55" y="48">These curves do not establish equilibrium or experimental force-field accuracy.</text>']
    xmin, xmax = rows[0]['time_ps']/1000., rows[-1]['time_ps']/1000.
    if xmax <= xmin:
        xmax = xmin + .001
    for panel, (key, label) in enumerate(panels):
        top, left, plotw, ploth = 88 + panel*195, 80, 770, 130
        values = [(r['time_ps']/1000., float(r[key])) for r in rows
                  if r[key] is not None and math.isfinite(float(r[key]))]
        content.append('<text x="55" y="%d" style="font-weight:bold">%s</text>' % (top-12, escape(label)))
        if not values:
            content.append('<text x="80" y="%d">Not available</text>' % (top+40))
            continue
        ymin, ymax = min(v[1] for v in values), max(v[1] for v in values)
        pad = max((ymax-ymin)*.1, .01)
        ymin, ymax = ymin-pad, ymax+pad
        content.append('<path d="M%d,%d V%d H%d" fill="none" stroke="#728092"/>' % (left,top,left*0+top+ploth,left+plotw))
        points = ' '.join('%.2f,%.2f' % (left+(t-xmin)/(xmax-xmin)*plotw,
                                       top+ploth-(v-ymin)/(ymax-ymin)*ploth) for t,v in values)
        content.append('<polyline points="%s" fill="none" stroke="#237a9c" stroke-width="1.8"/>' % points)
        for fraction in (0., .5, 1.):
            xx = left + fraction*plotw
            content.append('<text x="%.1f" y="%d" text-anchor="middle">%.3g</text>' % (xx,top+ploth+18,xmin+fraction*(xmax-xmin)))
        for fraction in (0., 1.):
            yy = top+ploth-fraction*ploth
            content.append('<text x="72" y="%.1f" text-anchor="end">%.3g</text>' % (yy+4,ymin+fraction*(ymax-ymin)))
        content.append('<text x="450" y="%d" text-anchor="middle">Time (ns)</text>' % (top+ploth+35))
    content.append('</svg>')
    path.write_text('\n'.join(content)+'\n')


def analyze(run, gmx=None):
    run = Path(run).resolve()
    analysis = run / 'analysis'
    analysis.mkdir(exist_ok=True)
    manifest = read_json(run / 'input/manifest.json', {})
    run_status = read_json(run / 'run_status.json', {'status': 'unknown'})
    messages = []
    stage = None
    if gmx:
        stage, messages = gromacs_export(run, analysis, gmx)
    if not stage:
        stage = run_status.get('stage', 'unknown')
    if stage in ('nvt','npt') and manifest.get('position_restraints'):
        messages.append('This trajectory has solute heavy-atom restraints. Its structural stability is not evidence of unrestrained molecular stability.')
    atoms = manifest.get('atoms', [])
    audit = manifest.get('chemical_audit', manifest.get('audit', {}))
    param = manifest.get('parameter_audit', {})
    if not isinstance(audit, dict):
        audit = {'detail': audit}
    backbone = manifest.get('peptide_backbone', [])
    expected_charge = manifest.get('expected_formal_charge', manifest.get('expectedformalcharge'))
    charges = [a.get('charge') for a in atoms]
    charge = sum(charges) if atoms and all(q is not None and math.isfinite(float(q)) for q in charges) else None
    charge_ok = (charge is not None and expected_charge is not None
                 and abs(charge-float(expected_charge)) < 1e-4)
    count = len(atoms)
    masses = np.array([float(a.get('mass', 0.)) for a in atoms])
    if count and (np.any(masses <= 0) or not np.isfinite(masses).all()):
        raise ValueError('Manifest contains non-positive/non-finite atom masses')
    heavy = np.array([i for i,a in enumerate(atoms) if a.get('element', '').upper() not in ('H','D')], dtype=int)
    bb = safe_indices([b.get(k) for b in backbone for k in ('N','CA','C')], count)
    tail = safe_indices(manifest.get('tail_atom_indices', []), count)
    peptide = safe_indices(manifest.get('peptide_atom_indices', []), count)
    heavy_set = set(heavy.tolist())
    tail = [i for i in tail if i in heavy_set]
    peptide = [i for i in peptide if i in heavy_set and i not in set(tail)]
    excluded = set()
    neighbors = {i:set() for i in range(count)}
    for a,b in manifest.get('bonds', []):
        neighbors[int(a)].add(int(b))
        neighbors[int(b)].add(int(a))
    for a in tail:
        one = neighbors[a]
        two = set().union(*(neighbors[i] for i in one)) if one else set()
        for b in one | two | {a}:
            excluded.add((a,b))
    pairmask = np.array([[(a,b) not in excluded for b in peptide] for a in tail], dtype=bool)
    rows, per_residue, reference, trajectory_error = [], {}, None, None
    stereocenters = manifest.get('stereocenters', [])
    stereo_failures, minimum_stereo_volume = [], None
    gro = analysis / 'frames.gro'
    if gro.is_file() and count:
        try:
            for time, xyz, box in gro_frames(gro):
                if len(xyz) != count:
                    raise ValueError('Exported trajectory has %d atoms; manifest expects %d' % (len(xyz), count))
                if reference is None:
                    reference = xyz.copy()
                contacts = None
                if tail and peptide:
                    dist2 = np.sum((xyz[tail,None,:] - xyz[None,peptide,:])**2, axis=2)
                    contacts = int(np.count_nonzero((dist2 < .45**2) & pairmask))
                helix = helix_values(xyz, backbone)
                for resnr, value in helix.items():
                    per_residue.setdefault(resnr, []).append(value)
                stereo, volume = stereocenter_retention(xyz, stereocenters)
                stereo_failures.extend(dict(failure, time_ps=time) for failure in stereo)
                if volume is not None:
                    minimum_stereo_volume = volume if minimum_stereo_volume is None else min(minimum_stereo_volume, volume)
                rows.append({'time_ps': time,
                             'added_NaCl_mM': (1000.*manifest.get('ions',{}).get('NaCl_pairs',0)/(.602214076*abs(float(np.linalg.det(box))))),
                             'backbone_rmsd_nm': kabsch_rmsd(xyz[bb], reference[bb]) if len(bb)>2 else None,
                             'rg_nm': radius_of_gyration(xyz, masses),
                             'tail_peptide_contacts': contacts,
                             'alpha_region_fraction': float(np.mean(list(helix.values()))) if helix else None,
                             'stereochemistry_mismatches': len(stereo) if stereocenters else None,
                             'min_periodic_copy_distance_nm': minimum_image_copy_distance(xyz[heavy], box) if len(heavy) else None})
        except (ValueError, np.linalg.LinAlgError) as error:
            trajectory_error = str(error)
            messages.append('Trajectory analysis stopped: '+str(error))
    cutoff = float(manifest.get('cutoff_nm', 1.0))
    images = [r['min_periodic_copy_distance_nm'] for r in rows if r['min_periodic_copy_distance_nm'] is not None]
    image_min = min(images) if images else None
    image_ok = image_min is not None and image_min > cutoff
    thermo = xvg_summary(analysis / 'thermodynamics.xvg')
    log_findings = []
    for s in ('em','nvt','npt','production'):
        path = run / (s+'.log')
        if not path.is_file():
            continue
        for line in path.read_text(errors='replace').splitlines():
            if re.search(r'LINCS WARNING|constraint.*(error|fail)|Fatal error|nan\b|not finite|could not be settled', line, re.I):
                log_findings.append({'stage':s, 'line':line.strip()[:240]})
    numeric_ok = bool(rows) and bool(thermo) and not trajectory_error and not log_findings and image_ok and not stereo_failures
    if thermo and not all(v['finite'] for v in thermo.values()):
        numeric_ok = False
    sanity_alerts = []
    temperature = thermo.get('Temperature',{}).get('second_half_mean')
    target_temperature = float(manifest.get('temperature_K',298.15))
    if temperature is not None and abs(temperature-target_temperature)>15.:
        sanity_alerts.append('Mean temperature differs from target by more than 15 K.')
    density = thermo.get('Density',{}).get('second_half_mean')
    if density is not None and not 800.<density<1200.:
        sanity_alerts.append('Mean density lies outside the broad 800-1200 kg/m3 water-box sanity range.')
    salt_values = [r['added_NaCl_mM'] for r in rows]
    if salt_values and manifest.get('ions',{}).get('NaCl_pairs'):
        target_salt = 1000.*manifest['ions'].get('requested_added_salt_M',.150)
        if abs(float(np.mean(salt_values))-target_salt)>15.:
            sanity_alerts.append('Sampled added NaCl differs from target by more than 15 mM; inspect box relaxation and ion count.')
    if sanity_alerts:
        numeric_ok = False
        messages.extend(sanity_alerts)
    requested_ns = manifest.get('production_ns', manifest.get('requested', {}).get('production_ns', 1.0))
    duration_ns = (rows[-1]['time_ps'] - rows[0]['time_ps'])/1000. if len(rows)>1 else 0.
    complete = run_status.get('status') == 'complete'
    goals = [
        {'id':1,'name':'Chemical identity/connectivity','status':'attention_required' if stereo_failures else audit.get('status','not_verified'),
         'evidence':audit,
         'stereochemistry_retention': {'centers_per_frame':len(stereocenters), 'frames_checked':len(rows) if stereocenters else 0,
                                      'mismatches':stereo_failures, 'minimum_abs_volume_nm3':minimum_stereo_volume,
                                      'note':'Geometry checked against manifest CIP ordering; absent centers are not implicitly validated.'}},
        {'id':2,'name':'Assigned protonation/solute charge','status':'checked' if charge_ok else 'not_verified',
         'expected_solute_charge_e':expected_charge,'actual_solute_charge_e':charge,
         'note':'This checks the assigned state, not experimental pKa values. See manifest for ion counts and system neutrality.'},
        {'id':3,'name':'Nonstandard force-field parameters','status':'coverage_checked_physics_unvalidated' if param.get('coverage_complete') else 'not_verified',
         'evidence':param,'note':'Short MD cannot validate torsional energetics, hydration free energies or binding kinetics.'},
        {'id':4,'name':'Numerical behavior / sampling','status':'short_run_checks_passed_sampling_unvalidated' if numeric_ok else 'attention_required',
         'run_complete':complete,'sampled_stage':stage,'min_periodic_copy_distance_nm':image_min,
         'nonbonded_cutoff_nm':cutoff,'log_findings':log_findings,
         'note':'The absence of direct periodic-image contact is not a finite-size convergence test. One short trajectory cannot establish an equilibrium ensemble.'},
        {'id':5,'name':'Solution-structure agreement with experiment','status':'not_validated',
         'note':'Metrics are reported, but no matched monomer experiment is fitted or declared passed. Alpha-region fraction is a torsion proxy, not DSSP or a CD spectrum.'},
    ]
    summary = {'schema_version':1, 'run_directory':str(run), 'run_status':run_status,
               'sampled_stage':stage,'requested_production_ns':requested_ns,
               'sampled_time_span_ns':duration_ns,'sampled_frames':len(rows),
               'trajectory_error':trajectory_error,'goals':goals,'thermodynamics':thermo,
               'sanity_alerts':sanity_alerts,
               'sanity_filter_note':'Broad implementation alerts, not experimental validation tolerances; pressure is reported without a short-trajectory pass/fail target.',
               'metrics':{key:mean_stats([r[key] for r in rows if r[key] is not None])
                          for key in ('backbone_rmsd_nm','rg_nm','tail_peptide_contacts','alpha_region_fraction','min_periodic_copy_distance_nm','added_NaCl_mM')},
               'alpha_region_definition':{'phi_deg':[-100,-30], 'psi_deg':[-80,-5], 'method':'dihedral-region proxy; not DSSP'},
               'residue_alpha_region_fraction':{k:float(np.mean(v)) for k,v in per_residue.items()},
               'provenance':manifest.get('provenance',{}), 'messages':messages,
               'interpretation':'A short single-molecule run is an implementation and local-structure screen, not complete force-field validation or a PLGA-release calculation.'}
    (analysis / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    if rows:
        with (analysis / 'timeseries.csv').open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)
        write_plot(analysis / 'diagnostics.svg', rows)
    text = ['SINGLE SEMAGLUTIDE: SHORT-RUN REPORT',
            'Run status: '+str(run_status.get('status','unknown')),
            'Sampled stage: '+str(stage),
            'Requested production: %s ns; sampled span in this stage: %.4f ns (%d frames).' % (requested_ns,duration_ns,len(rows)),
            '']
    for goal in goals:
        text.append('%d. %s: %s' % (goal['id'],goal['name'],goal['status']))
        if 'note' in goal:
            text.append('   '+goal['note'])
    text += ['', 'PRIMARY DIAGNOSTICS', 'Minimum periodic-copy distance: %s nm; nonbonded cutoff: %.3f nm.' % (image_min,cutoff),
             'Total solute charge: %s e; expected assigned charge: %s e.' % (charge,expected_charge),
             'Added NaCl concentration uses the instantaneous box volume; see added_NaCl_mM in timeseries.csv and summary.json.',
             'Stereocenters checked per frame: %d; mismatches across sampled frames: %d.' % (len(stereocenters),len(stereo_failures)),
             'Metrics: timeseries.csv / diagnostics.svg; full details: summary.json.',
             '', 'A stable 1-ns trajectory does not establish equilibrium or experimental force-field accuracy.']
    text.extend(messages)
    (analysis / 'report.txt').write_text('\n'.join(text)+'\n')
    if (analysis/'semaglutide.gro').is_file() and (analysis/'semaglutide.xtc').is_file():
        (run/'view_single_sema.vmd').write_text('''# VMD Tk Console: source /path/to/run/view_single_sema.vmd
set sema_dir [file dirname [file normalize [info script]]]
mol new [file join $sema_dir analysis semaglutide.gro] type gro waitfor all
mol addfile [file join $sema_dir analysis semaglutide.xtc] type xtc waitfor all
mol delrep 0 top
mol representation Licorice 0.2 12.0 12.0
mol selection {all}
mol color Name
mol addrep top
display projection Orthographic
axes location Off
''')
    print('\n'.join(text))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', required=True, type=Path)
    parser.add_argument('--gmx', help='GROMACS executable; omit to reuse existing analysis/frames.gro')
    args = parser.parse_args()
    analyze(args.run_dir, args.gmx)


if __name__ == '__main__':
    main()
