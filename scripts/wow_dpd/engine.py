"""Native LAMMPS input and coordinate access for the solvent-explicit W/O/W particle."""
import math

import numpy as np


def quote(path):
    path = str(path.resolve())
    if any(c in path for c in ('"', '\n', '$')):
        raise ValueError("LAMMPS file paths cannot contain a quote, newline or dollar sign")
    return f'"{path}"'


def static_commands(args, interactions):
    matrix = np.asarray(interactions['a_ij'], dtype=float)
    if matrix.shape != (7, 7) or not np.all(np.isfinite(matrix)) or not np.all(matrix == matrix.T):
        raise ValueError("a_ij must be a finite symmetric 7 by 7 matrix")
    gamma = float(interactions['gamma'])
    if gamma <= 0 or not math.isfinite(gamma):
        raise ValueError("gamma must be positive")
    setup = ['units lj', 'atom_style molecular', 'boundary p p p', 'newton on',
             'atom_modify map array', 'bond_style harmonic']
    if args.restart:
        setup += [f'read_restart {quote(args.restart)}']
    else:
        setup += ['special_bonds lj 1.0 1.0 1.0', f'read_data {quote(args.data)}']
    setup += ['comm_modify vel yes cutoff 4.5', 'neighbor 0.3 bin',
              'neigh_modify every 1 delay 0 check yes',
              f'pair_style dpd 1.0 1.0 {args.seed}']
    for i in range(len(matrix)):
        for j in range(i, len(matrix)):
            setup.append(f'pair_coeff {i+1} {j+1} {matrix[i,j]:g} {gamma:g} 1.0')
    setup += ['fix integrate all nve', f'timestep {args.dt:.12g}',
              'thermo_style custom step atoms temp pe ke press',
              'thermo_modify lost error flush yes']
    return setup


def frame(lmp):
    n = int(lmp.get_natoms())
    if int(lmp.extract_setting('nlocal')) != n:
        raise RuntimeError("Run one Python process and one MPI rank per simulation")
    ids = lmp.numpy.extract_atom('id')[:n]
    types = lmp.numpy.extract_atom('type')[:n]
    xyz = lmp.numpy.extract_atom('x')[:n, :3]
    lo, hi, xy, yz, xz, periodic, _ = lmp.extract_box()
    lengths = np.asarray(hi) - np.asarray(lo)
    if (not np.allclose(lengths, lengths[0]) or any(abs(v) > 1e-10 for v in (xy, yz, xz))
            or not all(periodic)):
        raise ValueError("The sink expects a cubic, triply periodic cell")
    return ids, types, xyz - np.asarray(lo), float(lengths[0])
