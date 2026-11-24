"""Measure instantaneous particle water, solvent mass and aqueous connectivity in a DPD run."""
import argparse
import csv
import gzip
from itertools import islice
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import binary_closing, binary_fill_holes, gaussian_filter, label, map_coordinates
from scipy.spatial import ConvexHull, cKDTree

WATER, DCM, LA, GA = 1, 2, 3, 4
M_W, M_D, M_L, M_G = 18.01528, 84.933, 72.063, 58.036


def dump_frames(path, rc, stride):
    """Stream the run_dpd.py orthorhombic cubic dump, retaining and verifying atom IDs."""
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as f:
        for frame_index, header in enumerate(iter(f.readline, '')):
            if header.strip() != 'ITEM: TIMESTEP':
                raise ValueError('Expected ITEM: TIMESTEP')
            step = int(next(f))
            if next(f).strip() != 'ITEM: NUMBER OF ATOMS':
                raise ValueError('Missing atom count')
            n = int(next(f))
            if next(f).split() != 'ITEM: BOX BOUNDS pp pp pp'.split():
                raise ValueError('Requires orthogonal periodic boundaries')
            lo, hi = np.array([list(map(float, next(f).split())) for _ in range(3)]).T
            if not np.allclose(hi-lo, hi[0]-lo[0]) or hi[0] <= lo[0]:
                raise ValueError('Requires a cubic box')
            header = next(f).split()
            if header[:2] != ['ITEM:', 'ATOMS']:
                raise ValueError('Missing atom table')
            rows = islice(f, n)
            if frame_index % stride:
                if sum(1 for _ in rows) != n:
                    raise ValueError('Incomplete frame')
                continue
            cols = [header[2:].index(key) for key in ('id', 'type', 'x', 'y', 'z')]
            data = np.loadtxt(rows, usecols=cols, ndmin=2)
            if len(data) != n:
                raise ValueError('Incomplete frame')
            if not np.array_equal(data[:, 0], np.arange(1, n+1)):
                data = data[np.argsort(data[:, 0])]
                if not np.array_equal(data[:, 0], np.arange(1, n+1)):
                    raise ValueError('IDs do not match contiguous original geometry IDs')
            yield step, data[:, 1].astype(int), (data[:, 2:] - lo)*rc, (hi[0]-lo[0])*rc


def box_faces(grid):
    return (grid[0], grid[-1], grid[:, 0], grid[:, -1], grid[:, :, 0], grid[:, :, -1])


def surface_vertices(organic, contrast, cell):
    """Interpolate the zero-contrast contour on edges bordering the selected organic component."""
    pieces = []
    for axis in range(3):
        left = [slice(None)]*3
        right = left.copy()
        left[axis], right[axis] = slice(None, -1), slice(1, None)
        left, right = tuple(left), tuple(right)
        crosses = organic[left] != organic[right]
        a, b = contrast[left][crosses], contrast[right][crosses]
        p = np.argwhere(crosses).astype(float) + 0.5
        p[:, axis] += a/(a-b)
        pieces.append(p*cell)
    return np.vstack(pieces)


def hull_test(points, hull, depth):
    """Exact convex halfspace test, accelerated using a maximum dot-product search."""
    normals, intercepts = hull.equations[:, :3], hull.equations[:, 3]
    center = hull.points[hull.vertices].mean(axis=0)
    widths = -(normals @ center + intercepts) - depth
    if widths.min() <= 0:
        raise ValueError('Deep-core margin exceeds the usable particle radius')
    dual = normals/widths[:, None]
    lengths = np.einsum('ij,ij->i', dual, dual)
    # Pad to equal squared length in 4D. Nearest vector maximizes the dot product.
    dual4 = np.column_stack((dual, np.sqrt(np.maximum(0, lengths.max()-lengths))))
    tree = cKDTree(dual4)
    rel = points-center
    radii = np.linalg.norm(rel, axis=1)
    queries = np.zeros((len(points), 4))
    queries[:, :3] = rel * (np.sqrt(lengths.max())/np.maximum(radii, 1e-15))[:, None]
    nearest = tree.query(queries)[1]
    return np.einsum('ij,ij->i', rel, dual[nearest]) <= 1


def probe_envelope(component, points, pitch, radius):
    """Close sub-probe openings and fill disconnected voids. Keep broad concavities."""
    k = int(np.ceil(radius/pitch))
    x, y, z = np.ogrid[-k:k+1, -k:k+1, -k:k+1]
    ball = x*x + y*y + z*z <= (radius/pitch)**2
    # Allow the dilation to extend past a box face before the erosion.
    padded = np.pad(component, k+1)
    closed = binary_closing(padded, structure=ball)
    trimmed = closed[k+1:-(k+1), k+1:-(k+1), k+1:-(k+1)]
    filled = binary_fill_holes(trimmed)
    # Interpolate the voxel mask to avoid quantizing all beads to voxel centers.
    weight = map_coordinates(filled.astype(float), (points/pitch-0.5).T,
                             order=1, mode='nearest', prefilter=False)
    return weight >= 0.5, filled.sum()*pitch**3


def composition(inside, types, layout):
    w, d, la, ga = (np.count_nonzero(inside & (types == k)) for k in (WATER, DCM, LA, GA))
    p_mass = la*M_L + ga*M_G
    w_mass = w*layout['water_molecules_per_bead']*M_W
    d_mass = d*M_D
    if not p_mass:
        raise ValueError('No PLGA within the detected carrier envelope')
    # Additive bulk-volume convention used in the project; peptide is excluded.
    water_volume = w_mass/1.000
    return dict(water_inside=w, DCM_inside=d, DCM_outside=(types == DCM).sum()-d,
                LA_inside=la, GA_inside=ga, peptide_sites_inside=np.count_nonzero(inside & (types >= 5)),
                phi_water_all_sites=w/inside.sum(),
                phi_water_additive_carrier=water_volume/(water_volume + d_mass/1.326 + p_mass/1.340),
                dcm_wt_fraction_particle=d_mass/(d_mass+p_mass))


def carrier_grid(types, xyz, length, grid_nm, smooth_nm):
    """Recenter on PLGA and compute organic minus water density on a periodic grid."""
    polymer = (types == LA) | (types == GA)
    moment = np.exp((2j*np.pi/length)*xyz[polymer]).mean(axis=0)
    if np.min(abs(moment)) < 0.05:
        raise ValueError('No localized PLGA center for a single-particle analysis')
    center = np.angle(moment)*length/(2*np.pi)
    pos = (xyz-center+0.5*length) % length
    n = max(4, int(np.ceil(length/grid_nm)))
    shape, h = (n, n, n), length/n
    bins = np.ravel_multi_index(np.minimum((pos/h).astype(int), n-1).T, shape)

    def density(mask):
        grid = np.bincount(bins[mask], minlength=n**3).reshape(shape).astype(float)
        return gaussian_filter(grid, smooth_nm/h, mode='wrap')

    water = types == WATER
    dcm = types == DCM
    # Interface defined where water and organic carrier densities are equal. Peptide
    # is a tracer in this operation and can occupy either aqueous or organic domains.
    contrast = density(polymer | dcm)-density(water)
    return pos, bins, h, contrast


def aqueous_domains(contrast):
    """Label water-rich voxels and identify components connected to a box face."""
    wet_domains, n_wet = label(contrast < 0)
    external = np.zeros(n_wet+1, dtype=bool)
    for face in box_faces(wet_domains):
        external[face] = True
    external[0] = False
    return wet_domains, external


def initial_pocket_domains(types, xyz_nm, length, pocket, n_expected):
    """Map the water in each initial W1 pocket onto the trajectory analysis grid.

    A pocket is closed and distinct when its water-rich voxels all belong to
    one isolated component and no other pocket uses the same component. Beads
    in interfacial/organic voxels (label zero) do not define an aqueous path.
    """
    _, bins, _, contrast = carrier_grid(types, np.asarray(xyz_nm, dtype=float), length, 0.5, 0.7)
    domains, exterior = aqueous_domains(contrast)
    inner = (types == WATER) & (pocket > 0)
    origin = pocket[inner]
    labels = domains.ravel()[bins[inner]]
    domain_ids, bath_connected = [], []
    unsplit = True
    for index in range(1, n_expected + 1):
        ids, counts = np.unique(labels[(origin == index) & (labels != 0)], return_counts=True)
        dominant = int(ids[counts.argmax()]) if len(ids) else 0
        domain_ids.append(dominant)
        bath_connected.append(bool(exterior[ids].any()))
        unsplit = unsplit and len(ids) == 1
    closed = {d for d, bath in zip(domain_ids, bath_connected) if d and not bath}
    return dict(resolved_w1_count=len(closed), domain_ids=domain_ids,
                bath_connected=bath_connected,
                all_distinct_and_closed=unsplit and len(closed) == n_expected)


def measure(step, types, xyz, length, orig, pockets, layout, grid_nm, smooth_nm, deep_nm, probe_nm):
    if len(types) != len(orig) or not np.all((types == orig) | ((orig == DCM) & (types == WATER))):
        raise ValueError('Frame IDs or types incompatible with initial geometry and DCM-to-water sink')
    pos, bins, h, contrast = carrier_grid(types, xyz, length, grid_nm, smooth_nm)
    polymer = (types == LA) | (types == GA)
    water = types == WATER
    dcm = types == DCM
    domains, n_organic = label(contrast > 0)
    if n_organic == 0:
        raise ValueError('No organic-rich carrier')
    sizes = np.bincount(domains.ravel(), minlength=n_organic+1)[1:]
    largest = domains == (np.argmax(sizes)+1)
    if any(face.any() for face in box_faces(largest)):
        raise ValueError('Main organic phase reaches a recentered box face; no isolated envelope')
    hull = ConvexHull(surface_vertices(largest, contrast, h))
    inside = hull_test(pos, hull, 0)
    wet_domains, external = aqueous_domains(contrast)
    wet_in = inside & water
    wi = np.flatnonzero(wet_in)
    domain_at_water = wet_domains.ravel()[bins[wi]]
    bath_linked = external[domain_at_water]
    isolated = (domain_at_water != 0) & ~bath_linked
    deep = hull_test(pos[wi], hull, deep_nm)
    local = composition(inside, types, layout)
    group = np.where(orig[wi] == DCM, 2, np.where(pockets[wi] > 0, 0, 1))
    w1, w2, replacement = np.bincount(group, minlength=3)
    probe_inside, probe_volume = probe_envelope(largest, pos, h, probe_nm)
    probe_local = composition(probe_inside, types, layout)
    probe_wet = wet_domains.ravel()[bins[water & probe_inside]]
    probe_fields = dict(envelope_volume_nm3=probe_volume,
                        polymer_sites_inside_fraction=(probe_inside & polymer).sum()/polymer.sum(),
                        water_in_bath_linked_wet_voxels=external[probe_wet].sum(), **probe_local)
    return dict(step=step, envelope_volume_nm3=hull.volume, grid_pitch_nm=h, smooth_nm=smooth_nm,
                deep_nm=deep_nm, probe_nm=probe_nm, largest_organic_voxel_fraction=sizes.max()/sizes.sum(),
                polymer_sites_inside_fraction=(inside & polymer).sum()/polymer.sum(),
                initial_W1_inside=w1, initial_W2_inside=w2, replacement_water_inside=replacement,
                water_in_bath_linked_wet_voxels=bath_linked.sum(),
                water_in_isolated_wet_voxels=isolated.sum(),
                water_in_interfacial_or_organic_voxels=(domain_at_water == 0).sum(),
                isolated_occupied_wet_domains=len(np.unique(domain_at_water[isolated])),
                deep_water_in_bath_linked_wet_voxels=(bath_linked & deep).sum(),
                dcm_wt_fraction_box=dcm.sum()*M_D/(dcm.sum()*M_D+(types == LA).sum()*M_L+(types == GA).sum()*M_G),
                **local, **{'probe_'+key: value for key, value in probe_fields.items()})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--initial', type=Path, required=True, help='Original layout.json and geometry.npz directory')
    p.add_argument('--dump', type=Path, help='run_dpd.py id mol type x y z dump, optionally gzipped')
    p.add_argument('--out', type=Path, required=True, help='CSV path')
    p.add_argument('--stride', type=int, default=1, help='Evaluate every nth saved frame')
    p.add_argument('--grid-nm', type=float, default=0.5)
    p.add_argument('--smooth-nm', type=float, default=0.7)
    p.add_argument('--deep-nm', type=float, default=1.4, help='Inner depth for external-channel diagnostic')
    p.add_argument('--probe-nm', type=float, default=2.1, help='Radius for the concavity-sensitive envelope')
    args = p.parse_args()
    if args.stride < 1 or not all(np.isfinite(v) and v>0 for v in (args.grid_nm,args.smooth_nm,args.deep_nm,args.probe_nm)):
        p.error('Stride and spatial lengths must be positive')
    layout = json.loads((args.initial/'layout.json').read_text())
    with np.load(args.initial/'geometry.npz') as g:
        orig, pockets = g['types'], g['pocket']
        series = dump_frames(args.dump, layout['rc_nm'], args.stride) if args.dump else [
            (0, orig, g['xyz_nm'], layout['box_length_nm'])]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w', newline='') as f:
        writer = None
        for step, types, xyz, length in series:
            row = measure(step, types, xyz.astype(float), length, orig, pockets, layout,
                          args.grid_nm, args.smooth_nm, args.deep_nm, args.probe_nm)
            if writer is None:
                writer = csv.DictWriter(f, fieldnames=list(row))
                writer.writeheader()
            writer.writerow(row)
            f.flush()
        if writer is None:
            raise ValueError('No sampled frames')


if __name__ == '__main__':
    main()
