"""Build a seed-specific, solvent-explicit W1/O/W2 input on the compute node."""
import argparse
import json
import math
import tempfile
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist

from initial_preview import write_preview
from molecule_packing import pack_molecules
from pocket_geometry import place_pockets
from particle_series import initial_pocket_domains


# Composition, units and molecular mapping of the archived folded-C18 case.
# Huilly WO2025108020A1, Example 1: 70 mg drug / 3 mL W1,
# 2.4 g 50:50 PLGA / 10 mL DCM. Only these four species are represented.
D_REF = 39.2                     # Wet outer diameter (nm), reported precision
BOX_REF = 44.65246606818088      # One finite external bath, not the feed reservoir
RC = 0.70
PITCH = RC / np.cbrt(3.0)
W1_RADIUS_REF = 5.152           # Each of the ten reference water pockets
TYPE_NAMES = ('water', 'DCM', 'lactyl', 'glycolyl', 'polar peptide',
              'apolar peptide', 'C18 hydrocarbon')
PEPTIDE_MASS = 4109.616         # Supplied, ionized molecular topology in Da
WATER_MASS = 4 * 18.01528
MIN_OIL_GAP = 1.05                # 1.5 cutoffs; verify visibility on the actual analysis grid


def sphere_volume(r):
    return 4 * math.pi * r ** 3 / 3


def composition(diameter, count):
    if not math.isfinite(diameter) or diameter <= 0:
        raise ValueError('Wet droplet diameter must be positive and finite')
    if not 1 <= count <= 10:
        raise ValueError('W1 count must be an integer from 1 to 10')
    scale = (diameter / D_REF) ** 3
    wet = diameter / 2
    length = diameter + (BOX_REF - D_REF)
    bath_ref = BOX_REF ** 3 - sphere_volume(D_REF / 2)
    bath_volume = length ** 3 - sphere_volume(wet)
    def nearest(number):
        return math.floor(number + 0.5)
    amounts = dict(chains=nearest(140 * scale), drug=nearest(20 * scale),
                   dcm=nearest(182322 * scale), w1=nearest(48724 * scale),
                   w2=nearest(503130 * bath_volume / bath_ref))
    if min(amounts.values()) < 1:
        raise ValueError('At this diameter at least one component rounds to zero. '
                         'Increase the wet diameter to model a drug-loaded particle.')
    r = W1_RADIUS_REF * (diameter / D_REF) * (10 / count) ** (1 / 3)
    return wet, length, r, amounts


def solvent_lattice(length, wet, centers, radius, rng):
    # Slightly oversample a jittered cubic grid. Remove surplus sites at random
    # within each phase so rotating/repacking the pockets cannot change the recipe.
    # The target density applies after this thinning, not to the candidate grid.
    n = math.ceil(length / PITCH) + 1
    grid = np.indices((n, n, n), dtype=np.int32).reshape(3, -1).T
    xyz = (grid.astype(float) + 0.5 - n / 2) * (length / n)
    xyz += rng.uniform(-0.035, 0.035, size=xyz.shape) * (length / n)
    inside = np.einsum('ij,ij->i', xyz, xyz) < wet ** 2
    distance, nearest = cKDTree(centers).query(xyz[inside])
    pocket = np.zeros(len(xyz), dtype=np.int16)
    pocket[inside] = np.where(distance < radius, nearest + 1, 0)
    return xyz, inside, pocket, n


def pick_sites(rng, candidates, count, phase):
    if len(candidates) < count:
        raise ValueError(f'{phase} has {len(candidates)} available sites for {count} beads. '
                         'The requested geometry and molecular clearances cannot be filled.')
    return rng.choice(candidates, size=count, replace=False)


def solvent_sites(rng, lattice, inside, pocket, molecules, target):
    oil = inside & (pocket == 0)
    oil[molecules['polymer_site_indices'].ravel()] = False
    dcm = pick_sites(rng, np.flatnonzero(oil), target['dcm'], 'O')
    w2 = pick_sites(rng, np.flatnonzero(~inside), target['w2'], 'W2')
    w1_ids = np.flatnonzero(pocket)
    available = np.ones(len(w1_ids), dtype=bool)
    drug_xyz = molecules['xyz_nm'][molecules['types'] >= 5]
    for xyz, gap in ((drug_xyz, 0.38), (molecules['peptide_atoms_nm'].reshape(-1, 3), 0.22)):
        distances = cKDTree(xyz).query(lattice[w1_ids], distance_upper_bound=gap)[0]
        available &= distances >= gap
    w1 = pick_sites(rng, w1_ids[available], target['w1'], 'W1')
    return np.concatenate((dcm, w1, w2))


def assemble(diameter, count, seed, attempt_seed):
    wet, length, radius, target = composition(diameter, count)
    # Independent streams isolate pocket positions from molecular and solvent packing.
    streams = [np.random.default_rng(s) for s in attempt_seed.spawn(5)]
    centers = place_pockets(streams[0], count, wet, radius, MIN_OIL_GAP)
    rim_gaps = wet - radius - np.linalg.norm(centers, axis=1)
    pair_gaps = pdist(centers) - 2 * radius
    if np.any(rim_gaps < MIN_OIL_GAP - 1e-9) or np.any(pair_gaps < MIN_OIL_GAP - 1e-9):
        raise ValueError('W1 placement violated the oil film clearance')
    per_pocket, extra = divmod(target['drug'], count)
    occupancy = np.full(count, per_pocket, dtype=int)
    occupancy[streams[1].permutation(count)[:extra]] += 1

    grid, inside, pocket, n_axis = solvent_lattice(length, wet, centers, radius, streams[2])
    mol = pack_molecules(streams[3], centers, radius, wet, occupancy, target['chains'],
                         grid, inside & (pocket == 0), n_axis)
    sites = solvent_sites(streams[4], grid, inside, pocket, mol, target)
    # All archived geometry arrays and LAMMPS input positions are in the positive box.
    xyz = (np.vstack((mol['xyz_nm'], grid[sites])) + length / 2).astype(np.float32)
    types = np.r_[mol['types'], np.where(inside[sites] & (pocket[sites] == 0), 2, 1)].astype(np.uint8)
    group = np.r_[mol['mol'], np.zeros(len(sites), dtype=np.int32)].astype(np.int32)
    pockets = np.r_[mol['pocket'], pocket[sites]].astype(np.int16)
    counts = np.bincount(types, minlength=8)
    waters = np.bincount(pockets[types == 1], minlength=count + 1)
    expected = [target['w1'] + target['w2'], target['dcm'], target['chains'] * 154,
                target['chains'] * 154, target['drug'] * 27, target['drug'] * 13, target['drug'] * 3]
    if not np.array_equal(counts[1:], expected) or waters[0] != target['w2']:
        raise ValueError('Assembled geometry does not match component counts')
    bonds = np.asarray(mol['bonds'], dtype=np.int32)
    data = dict(xyz_nm=xyz, types=types, mol=group, pocket=pockets,
                bonds=bonds[:, 1:] - 1, bond_types=bonds[:, 0].astype(np.uint8),
                peptide_atoms_nm=(mol['peptide_atoms_nm'] + length / 2).astype(np.float32))
    for k in ('peptide_bead_names', 'peptide_bead_regions', 'peptide_atom_to_bead',
              'peptide_reference_charge', 'peptide_bond_rest_nm'):
        data[k] = mol[k]
    # The archived case has 42 chemical and 29 local-shape springs per peptide.
    data['bond_kinds'] = np.where(bonds[:, 0] == 1, 'polymer',
                                  np.where(bonds[:, 0] <= 43, 'peptide', 'elastic'))
    plga_mass = counts[3] * 72.063 + counts[4] * 58.036
    drug_mass = target['drug'] * PEPTIDE_MASS
    measured = dict(water_W1_g_per_g_PLGA=target['w1'] * WATER_MASS / plga_mass,
                    DCM_g_per_g_PLGA=target['dcm'] * 84.933 / plga_mass,
                    semaglutide_g_per_g_PLGA=drug_mass / plga_mass,
                    semaglutide_mg_per_mL_W1=1.66053906660 * drug_mass /
                    (count * sphere_volume(radius)))
    layout = dict(seed=seed, requested_wet_diameter_nm=diameter,
                  rc_nm=RC, water_molecules_per_bead=4, target_reduced_number_density=3.0,
                  site_pitch_nm=PITCH, candidate_lattice_pitch_nm=length / n_axis,
                  box_length_nm=length, oil_center_nm=[length / 2] * 3, oil_radius_nm=wet,
                  min_W2_thickness_nm=length / 2 - wet, pocket_radius_nm=radius,
                  pocket_centers_nm=(centers + length / 2).tolist(), required_oil_gap_nm=MIN_OIL_GAP,
                  pocket_center_radius_nm=np.linalg.norm(centers, axis=1).tolist(),
                  pocket_outer_oil_gap_nm=rim_gaps.tolist(),
                  min_outer_pocket_oil_cover_nm=float(rim_gaps.min()),
                  min_inter_pocket_oil_gap_nm=float(pair_gaps.min()) if len(pair_gaps) else None,
                  pockets=count, peptides_per_pocket=occupancy.tolist(),
                  water_beads_per_pocket=waters[1:].tolist(), W1_water_sites=target['w1'],
                  W2_water_sites=target['w2'], DCM_sites=target['dcm'],
                  peptide_count=target['drug'], peptide_beads_per_molecule=43,
                  PLGA_chain_count=target['chains'], repeats_per_PLGA_chain=308,
                  peptide_template_files=list(map(str, mol['peptide_template_files'])),
                  peptide_template_index_per_molecule=mol['peptide_template_indices'].tolist(),
                  type_names=dict(zip(range(1, 8), TYPE_NAMES)), n_sites=len(xyz),
                  site_counts={str(i): int(counts[i]) for i in range(1, 8)},
                  n_bonds_including_local_shape=len(bonds),
                  target_feed=dict(water_W1_g_per_g_PLGA=3 / 2.4, DCM_g_per_g_PLGA=13.26 / 2.4,
                                   semaglutide_g_per_g_PLGA=0.07 / 2.4,
                                   semaglutide_mg_per_mL_W1=70 / 3),
                  realized_composition=measured)
    return layout, data, mol['bond_coeffs']


def make_case(diameter, count, seed):
    if seed < 0:
        raise ValueError('Initial seed must be nonnegative')
    # The particle analysis smooths the water/organic interface. Geometrically
    # disjoint pockets must also be distinct on that analysis grid at step zero.
    max_attempts = 8
    for attempt, child in enumerate(np.random.SeedSequence(seed).spawn(max_attempts), 1):
        layout, data, coefficients = assemble(diameter, count, seed, child)
        audit = initial_pocket_domains(data['types'], data['xyz_nm'], layout['box_length_nm'],
                                       data['pocket'], count)
        if audit['all_distinct_and_closed']:
            layout['initial_aqueous_domains'] = audit
            layout['placement_attempt'] = attempt
            return layout, data, coefficients
        print(f"Placement {attempt}: resolved {audit['resolved_w1_count']}/{count} "
              "closed W1 pockets on the analysis grid.", flush=True)
    raise ValueError(f'Could not resolve all {count} separate W1 pockets at step zero '
                     f'after {max_attempts} seeded placements at {diameter:g} nm. '
                     'Increase the diameter, reduce the pocket count or change the initial seed.')


def write_data(path, layout, data, coefficients):
    with path.open('w', encoding='ascii') as f:
        f.write('Seeded semaglutide / PLGA / DCM W1/O/W2, reduced DPD units\n\n')
        f.write(f"{len(data['types'])} atoms\n{len(data['bonds'])} bonds\n\n")
        f.write(f'7 atom types\n{len(coefficients)} bond types\n\n')
        for axis in 'xyz':
            f.write(f"0 {layout['box_length_nm'] / RC:.10f} {axis}lo {axis}hi\n")
        f.write('\nMasses\n\n')
        for i, name in enumerate(TYPE_NAMES, 1):
            f.write(f'{i} 1.0 # {name}\n')
        f.write('\nBond Coeffs\n\n')
        for i, (k, r0) in enumerate(coefficients, 1):
            f.write(f'{i} {k:.5f} {r0:.8f}\n')
        f.write('\nAtoms # molecular\n\n')
        for i, (m, t, (x, y, z)) in enumerate(zip(data['mol'], data['types'],
                                                   data['xyz_nm'].astype(float) / RC), 1):
            f.write(f'{i} {m} {t} {x:.8f} {y:.8f} {z:.8f}\n')
        f.write('\nBonds\n\n')
        for i, (t, (a, b)) in enumerate(zip(data['bond_types'], data['bonds']), 1):
            f.write(f'{i} {t} {a + 1} {b + 1}\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--droplet-diameter-nm', type=float, required=True)
    p.add_argument('--w1-count', type=int, required=True)
    p.add_argument('--init-seed', type=int, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    if args.out.exists():
        p.error(f'Initial directory already exists: {args.out}')
    try:
        layout, data, coefficients = make_case(args.droplet_diameter_nm, args.w1_count, args.init_seed)
        args.out.parent.mkdir(exist_ok=True, parents=True)
        with tempfile.TemporaryDirectory(prefix='.wow-build-', dir=args.out.parent) as staging:
            out = Path(staging) / 'initial'
            out.mkdir()
            write_data(out / 'system.data', layout, data, coefficients)
            np.savez_compressed(out / 'geometry.npz', **data)
            (out / 'layout.json').write_text(json.dumps(layout, indent=2, allow_nan=False) + '\n')
            write_preview(out / 'initial_preview.svg',
                          np.asarray(layout['pocket_centers_nm']) - layout['box_length_nm'] / 2,
                          layout['pocket_radius_nm'], layout['oil_radius_nm'], layout['box_length_nm'],
                          layout['peptides_per_pocket'], args.init_seed)
            out.rename(args.out)
    except (ValueError, RuntimeError, OverflowError) as e:
        p.exit(2, f'Cannot build W/O/W input: {e}\n')
    keys = ('seed', 'requested_wet_diameter_nm', 'pockets', 'n_sites', 'PLGA_chain_count',
            'DCM_sites', 'W1_water_sites', 'W2_water_sites', 'peptides_per_pocket',
            'min_inter_pocket_oil_gap_nm', 'min_outer_pocket_oil_cover_nm',
            'initial_aqueous_domains', 'placement_attempt', 'target_feed', 'realized_composition')
    print(json.dumps({k: layout[k] for k in keys}, indent=2))
    print(f"Initial layout and preview: {args.out.resolve()}")


if __name__ == '__main__':
    main()
