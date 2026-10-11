#!/bin/bash
# Usage: source /home/lee6456/MD/env.sh md|prep|amber|build
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    printf 'Use: source ~/MD/env.sh md|prep|amber|build\n' >&2
    exit 2
fi
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)/config.sh" || return
case "${1:-md}" in
    md|prep|amber|build) ;;
    *) printf 'Choose md, prep, amber, or build.\n' >&2; return 2 ;;
esac
if declare -F conda >/dev/null; then
    while (( ${CONDA_SHLVL:-0} > 0 )); do
        conda deactivate || return
    done
fi
module purge || return
case "${1:-md}" in
    md)
        module load "$MD_GCC_MODULE" || return
        if [[ ! -x "$MD_GMX_PREFIX/bin/gmx" || ! -f "$MD_GMX_PREFIX/.validated" ]]; then
            printf 'Validated GROMACS installation missing. Run sbatch ~/MD/setup_cpu.sbatch.\n' >&2
            return 1
        fi
        # Upstream GMXRC reads optional, possibly unset environment variables.
        # Keep the caller's nounset setting after sourcing the external hook.
        MD_ENV_NOUNSET=0
        if [[ $- == *u* ]]; then MD_ENV_NOUNSET=1; fi
        set +u
        MD_ENV_SOURCE_STATUS=0
        source "$MD_GMX_PREFIX/bin/GMXRC" || MD_ENV_SOURCE_STATUS=$?
        if (( MD_ENV_NOUNSET )); then set -u; fi
        if (( MD_ENV_SOURCE_STATUS )); then return "$MD_ENV_SOURCE_STATUS"; fi
        unset MD_ENV_NOUNSET MD_ENV_SOURCE_STATUS
        export GMX="$MD_GMX_PREFIX/bin/gmx"
        ;;
    build)
        module load "$MD_GCC_MODULE" || return
        [[ -x "$MD_CMAKE_PREFIX/bin/cmake" ]] || { printf 'Run setup_cpu.sbatch first.\n' >&2; return 1; }
        export PATH="$MD_CMAKE_PREFIX/bin:$PATH"
        ;;
    prep)
        module load "$MD_CONDA_MODULE" || return
        [[ -f "$MD_PREP_PREFIX/conda-meta/history" ]] || { printf 'Run setup_cpu.sbatch first.\n' >&2; return 1; }
        source "$(conda info --base)/etc/profile.d/conda.sh" || return
        conda activate "$MD_PREP_PREFIX" || return
        export PYTHONNOUSERSITE=1
        ;;
    amber)
        module load "$MD_AMBER_GCC_MODULE" "$MD_AMBER_MPI_MODULE" "$MD_AMBER_MODULE" || return
        ;;
esac
