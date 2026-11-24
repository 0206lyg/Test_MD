"""Pack rigid aqueous semaglutide conformers and connected 50:50 PLGA chains.

folded_templates.npz contains the five previously distributed project-topology
starts adapted from Agarwal & Haworth (doi:10.1007/s40203-025-00393-7).
It stores the source GRO hashes, local atom/bead positions in nm, and the
72 harmonic bond coefficients from the reference LAMMPS input. Type 1 is a
PLGA bond. Types 2 through 72 encode the peptide bonds and short backbone
shape restraints shared by all five conformers. No tail–backbone restraints
are added during placement.
"""
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist
from scipy.spatial.transform import Rotation


REPEATS = 308
WALL_GAP_NM = 0.55
HEAVY_GAP_NM = 0.65
ATOM_GAP_NM = 0.40


def _balanced_templates(rng, total, n_templates):
    """Sample without replacement per group, then shuffle across all pockets."""
    full, extra = divmod(total, n_templates)
    indices = np.tile(np.arange(n_templates), full)
    if extra:
        indices = np.concatenate((indices, rng.permutation(n_templates)[:extra]))
    return rng.permutation(indices)


def _place_peptides(rng, template, centers, radius, occupancies):
    base_atoms = template['atom_xyz_nm']
    base_beads = template['bead_xyz_nm']
    heavy = template['atom_heavy']
    count = int(sum(occupancies))
    assignment = _balanced_templates(rng, count, len(base_atoms))
    all_atoms = np.empty((count, base_atoms.shape[1], 3))
    all_beads = np.empty((count, base_beads.shape[1], 3))
    allowed_radius = radius - WALL_GAP_NM
    span = np.array([pdist(x).max() for x in base_atoms])

    offset = 0
    for pocket, (center, occupancy) in enumerate(zip(centers, occupancies), start=1):
        if occupancy == 0:
            continue
        indices = assignment[offset:offset+occupancy]
        if 2 * allowed_radius < span[indices].max():
            raise ValueError(f'W1 pocket {pocket} (radius {radius:.3f} nm) is too '
                             f'small for its folded peptide with {WALL_GAP_NM} nm clearance')
        # Redraw this pocket if early placements leave a cavity of the wrong shape.
        for _ in range(6):
            atom_tree = heavy_tree = None
            for molecule, conformer in enumerate(indices, start=offset):
                for _ in range(2000):
                    # Isotropic orientation and COM candidates uniform in accessible volume.
                    rotation = Rotation.from_quat(rng.normal(size=4)).as_matrix()
                    direction = rng.normal(size=3)
                    direction *= allowed_radius * rng.random()**(1/3) / np.linalg.norm(direction)
                    positions = base_atoms[conformer] @ rotation + direction
                    if (np.einsum('ij,ij->i', positions, positions) > allowed_radius**2).any():
                        continue
                    positions += center
                    if atom_tree is not None:
                        if (atom_tree.query(positions, distance_upper_bound=ATOM_GAP_NM)[0] < ATOM_GAP_NM).any():
                            continue
                        if (heavy_tree.query(positions[heavy], distance_upper_bound=HEAVY_GAP_NM)[0] < HEAVY_GAP_NM).any():
                            continue
                    all_atoms[molecule] = positions
                    all_beads[molecule] = base_beads[conformer] @ rotation + (direction + center)
                    break
                else:
                    break
                if molecule < offset + occupancy - 1:
                    packed = all_atoms[offset:molecule+1]
                    atom_tree = cKDTree(packed.reshape(-1, 3))
                    heavy_tree = cKDTree(packed[:, heavy].reshape(-1, 3))
            else:
                break
        else:
            raise ValueError(f'Could not pack {occupancy} folded peptides in W1 pocket '
                             f'{pocket} (radius {radius:.3f} nm) at the requested clearances')
        offset += occupancy
    return all_atoms, all_beads, assignment


def _polymer_walks(rng, oil, n_chains, n_axis):
    """Grow site-disjoint, nearest-neighbour walks in the organic region."""
    available = np.asarray(oil, dtype=bool).reshape(-1).copy()
    if len(available) != n_axis**3:
        raise ValueError('Oil mask must have one entry per site of the cubic grid')
    # The external W2 bath covers all box faces. This guarantees that adding
    # a flattened neighbour offset never accidentally wraps across a box face.
    cube = available.reshape((n_axis,) * 3)
    if any(np.any(face) for k in range(3) for face in
           (cube.take(0, axis=k), cube.take(-1, axis=k))):
        raise ValueError('The organic grid touches the box boundary')
    oil_ids = np.flatnonzero(available)
    if n_chains * REPEATS > len(oil_ids):
        raise ValueError('The organic grid has fewer sites than the requested PLGA repeats')
    offsets = np.array([1, -1, n_axis, -n_axis, n_axis**2, -n_axis**2])
    chains = np.empty((n_chains, REPEATS), dtype=np.intp)
    for chain_number, chain in enumerate(chains):
        for _ in range(150):
            starts = oil_ids[available[oil_ids]]
            chain[0] = rng.choice(starts)
            available[chain[0]] = False
            used = 1
            backtracks = 0
            while used < REPEATS and backtracks < 2 * REPEATS:
                candidates = chain[used-1] + offsets
                candidates = candidates[available[candidates]]
                if len(candidates):
                    chain[used] = candidates[rng.integers(len(candidates))]
                    available[chain[used]] = False
                    used += 1
                elif used > 12:
                    available[chain[used-12:used]] = True
                    used -= 12
                    backtracks += 1
                else:
                    break
            if used == REPEATS:
                break
            available[chain[:used]] = True
        else:
            raise ValueError(f'Could not fit connected PLGA chain {chain_number+1}/{n_chains} '
                             f'({REPEATS} repeats) in the available organic volume')
    return chains


def pack_molecules(rng, centers, radius, wet_radius, occupancies, n_chains,
                   lattice, oil_mask, n_axis):
    """Return molecule-major molecular beads, topology and occupied grid sites.

    Centers, lattice positions and returned positions use nm about the wet-
    particle origin. Solvent sites follow the molecular beads in the output
    writer. Each bond row is [LAMMPS bond type, 1-based ID, 1-based ID].
    """
    centers = np.asarray(centers, dtype=float)
    occupancies = np.asarray(occupancies)
    if centers.shape != (len(occupancies), 3) or np.any(occupancies < 0):
        raise ValueError('Expected one nonnegative peptide occupancy per pocket')
    if np.any(np.linalg.norm(centers, axis=1) + radius > wet_radius + 1.e-9):
        raise ValueError('W1 pocket extends outside the wet particle')

    path = Path(__file__).with_name('folded_templates.npz')
    with np.load(path, allow_pickle=False) as data:
        template = dict(data)
    atoms, peptides, assignment = _place_peptides(rng, template, centers, radius, occupancies)
    n_drug, n_sites, _ = peptides.shape
    chains = _polymer_walks(rng, oil_mask, n_chains, n_axis)
    # Each chain independently shuffles exactly equal numbers of LA and GA.
    repeats = np.repeat(np.array([3, 4], dtype=np.uint8), REPEATS // 2)
    chain_types = np.empty(chains.shape, dtype=np.uint8)
    for row in chain_types:
        row[:] = rng.permutation(repeats)
    grid_ids = chains.ravel()
    polymer_ids = np.arange(n_drug*n_sites+1, n_drug*n_sites+1+len(grid_ids))
    polymer_ids = polymer_ids.reshape(chains.shape)
    chain_bonds = np.column_stack((np.ones(n_chains * (REPEATS-1), dtype=int),
                                   polymer_ids[:, :-1].ravel(), polymer_ids[:, 1:].ravel()))
    local_bonds = template['bonds'].astype(int)
    endpoints = (local_bonds[None, :, :] + 1 +
                 n_sites * np.arange(n_drug)[:, None, None]).reshape(-1, 2)
    kinds = np.tile(np.arange(2, 2+len(local_bonds)), n_drug)
    peptide_bonds = np.column_stack((kinds, endpoints))
    bonds = np.concatenate((chain_bonds, peptide_bonds)).astype(np.int32)
    return {
        'xyz_nm': np.vstack((peptides.reshape(-1, 3), lattice[grid_ids])),
        'types': np.r_[np.tile(template['bead_types'], n_drug), chain_types.ravel()],
        'mol': np.r_[np.repeat(np.arange(1, n_drug+1), n_sites),
                     np.repeat(np.arange(n_drug+1, n_drug+n_chains+1), REPEATS)].astype(np.int32),
        'pocket': np.r_[np.repeat(np.arange(1, len(occupancies)+1), occupancies*n_sites),
                        np.zeros(len(grid_ids), dtype=np.int16)].astype(np.int16),
        'bonds': bonds,
        'bond_coeffs': template['bond_coeffs'],
        'polymer_site_indices': chains,
        'peptide_atoms_nm': atoms,
        'peptide_template_indices': assignment,
        'peptide_template_files': template['template_files'],
        'peptide_bead_names': template['bead_names'],
        'peptide_bead_regions': template['bead_regions'],
        'peptide_atom_to_bead': template['atom_to_bead'],
        'peptide_reference_charge': template['bead_reference_charge'],
        'peptide_bond_rest_nm': template['bond_rest_nm'],
    }
