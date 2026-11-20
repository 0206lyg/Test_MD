"""Identify DCM in the external aqueous reservoir from an instantaneous periodic frame."""
import math

import numpy as np
from scipy.ndimage import gaussian_filter, label


def eligible_ids(ids, types, xyz, box_length, *, grid_target=1.0, water_min=0.80, guard=1.0):
    """Return eligible IDs and diagnostics. All lengths are in DPD cutoffs.

    Center a cubic periodic box on its PLGA. A DCM bead must lie in water-rich
    fluid connected to the outside and beyond every polymer bead by a radial
    safety margin. Connected water channels inside the particle are protected.
    """
    if box_length <= 0 or grid_target <= 0 or guard < 0 or not 0 < water_min < 1:
        raise ValueError("Invalid box, grid, water fraction or guard")
    ids = np.asarray(ids)
    types = np.asarray(types)
    xyz = np.asarray(xyz)
    if types.ndim != 1 or ids.shape != types.shape or xyz.shape != (len(ids), 3):
        raise ValueError("Invalid coordinates or labels")
    if not np.isfinite(xyz).all():
        raise ValueError("Nonfinite coordinates")
    polymer = np.isin(types, [3, 4])
    if not np.any(polymer):
        raise ValueError("The droplet contains no PLGA")
    # PBC-safe centroid; a nearly uniform circular distribution has no
    # well-defined center and cannot support this single-droplet classification.
    moment = np.exp(2j * np.pi * xyz[polymer] / box_length).mean(axis=0)
    if min(abs(moment)) < 0.05:
        raise ValueError("Unable to center the droplet on its polymer distribution")
    center = (np.angle(moment) / (2 * np.pi) % 1.0) * box_length
    offset = (xyz - center + box_length / 2) % box_length - box_length / 2
    radii_sq = np.einsum('ij,ij->i', offset, offset)
    radius = float(np.sqrt(np.max(radii_sq[polymer])))
    candidates = (types == 2) & (radii_sq > (radius + guard)**2)
    info = {'polymer_center_x_rc': float(center[0]), 'polymer_center_y_rc': float(center[1]),
            'polymer_center_z_rc': float(center[2]), 'polymer_max_radius_rc': radius,
            'outside_radial_guard': int(np.count_nonzero(candidates))}
    if not np.any(candidates):
        info['eligible'] = 0
        return np.empty(0, dtype=ids.dtype), info

    nc = max(4, math.ceil(box_length / grid_target))
    ijk = np.floor((offset / box_length + 0.5) * nc).astype(np.intp).clip(0, nc - 1)
    voxel = np.ravel_multi_index(ijk.T, (nc, nc, nc))
    def density(indices):
        counts = np.bincount(indices, minlength=nc**3).reshape(nc, nc, nc).astype(float)
        # Standard deviation = one cutoff, not one fixed grid cell.
        return gaussian_filter(counts, sigma=nc / box_length, mode='wrap')
    water = density(voxel[types == 1])
    all_material = density(voxel)
    wet = (water >= water_min * all_material) & (all_material > 0)
    # The recentered box edge belongs to the remote aqueous bath. All six
    # faces seed connectivity, which also handles a percolating periodic bath.
    domains, n_domains = label(wet)
    external = np.zeros(n_domains + 1, dtype=bool)
    for face in (domains[0], domains[-1], domains[:, 0], domains[:, -1],
                 domains[:, :, 0], domains[:, :, -1]):
        external[face] = True
    external[0] = False
    selected = ids[candidates][external[domains.ravel()[voxel[candidates]]]]
    info['eligible'] = len(selected)
    return np.sort(selected), info


def exchange_ids(eligible, rate, elapsed, rng):
    """Discrete first-order external removal, including the closed control."""
    if not (math.isfinite(rate) and math.isfinite(elapsed) and rate >= 0 and elapsed >= 0):
        raise ValueError("Rate and elapsed time must be nonnegative and finite")
    p = -math.expm1(-rate * elapsed)
    eligible = np.asarray(eligible)
    return eligible[rng.random(len(eligible)) < p] if p else eligible[:0]
