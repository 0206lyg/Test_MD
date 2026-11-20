#!/bin/bash
# Default layout for lee6456 on Negishi. Source from the other scripts.
export MD_ROOT="${MD_ROOT:-$HOME/MD}"
export MD_SCRATCH_ROOT="${MD_SCRATCH_ROOT:-/scratch/negishi/$(id -un)}"

# Fixed release pair; the installer verifies the matching official checksums.
export MD_GROMACS_VERSION=2026.4
export MD_CMAKE_VERSION=3.31.6
export MD_GCC_MODULE="${MD_GCC_MODULE:-gcc/12.2.0}"
export MD_CONDA_MODULE="${MD_CONDA_MODULE:-anaconda/2024.02-py311}"
export MD_AMBER_MODULE="${MD_AMBER_MODULE:-amber/24}"
export MD_AMBER_GCC_MODULE="${MD_AMBER_GCC_MODULE:-gcc/12.2.0}"
export MD_AMBER_MPI_MODULE="${MD_AMBER_MPI_MODULE:-openmpi/4.1.4}"

export MD_GMX_PREFIX="$MD_ROOT/software/gromacs-$MD_GROMACS_VERSION-cpu"
export MD_CMAKE_PREFIX="$MD_ROOT/envs/cmake-$MD_CMAKE_VERSION"
export MD_PREP_PREFIX="$MD_ROOT/envs/prep"
export MD_LAMMPS_VERSION=2025.7.22.4.0
export MD_DPD_PREFIX="$MD_ROOT/envs/lammps-dpd-$MD_LAMMPS_VERSION"
# Reuse the original installation in place; Python venvs cannot be moved.
if [[ ! -d $MD_DPD_PREFIX && $MD_LAMMPS_VERSION == 2025.7.22.4.0 &&
      -x $MD_ROOT/envs/wow-dpd-2025.7.22.4/bin/python && -r $MD_ROOT/envs/wow-dpd-2025.7.22.4/lib/libmpi.so.12 ]]; then
    export MD_DPD_PREFIX="$MD_ROOT/envs/wow-dpd-2025.7.22.4"
fi
export MD_DPD_PYTHON="$MD_DPD_PREFIX/bin/python"
export MD_BUILD_DIR="$MD_ROOT/build/gromacs-$MD_GROMACS_VERSION-cpu"
export MD_SOURCE_DIR="$MD_ROOT/src/gromacs-$MD_GROMACS_VERSION"
export MD_TEST_DIR="$MD_ROOT/src/regressiontests-$MD_GROMACS_VERSION"
