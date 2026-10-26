#!/bin/bash
#SBATCH --job-name=md_prepare
#SBATCH --account=ardekani
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=md_prepare_%j.out

# Usage: sbatch prepare_poc.sbatch cases/poc_lawsone_dmso.json [OUTPUT_DIR]
# Preparation runs on a compute node; reusable input files remain under MD.
set -Eeuo pipefail
export MD_ROOT="${MD_ROOT:-$HOME/MD}"
[[ -n "${SLURM_JOB_ID:-}" ]] || { printf 'Use sbatch to run this script.\n' >&2; exit 1; }
[[ $# -ge 1 && $# -le 2 ]] || { printf 'Usage: sbatch prepare_poc.sbatch CASE.json [OUTPUT_DIR]\n' >&2; exit 2; }
[[ -f "$1" ]] || { printf 'Missing case configuration: %s\n' "$1" >&2; exit 1; }
MD_CASE_CONFIG="$(realpath -- "$1")"
MD_CASE_NAME="$(basename -- "${MD_CASE_CONFIG%.json}")"
MD_CASE_OUTPUT="$(realpath -m -- "${2:-$MD_ROOT/cases/${MD_CASE_NAME}_input}")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export TMPDIR="$MD_ROOT/tmp/prepare_${SLURM_JOB_ID}"
mkdir -p -- "$TMPDIR"
printf 'Configuration: %s\nPrepared inputs: %s\n' "$MD_CASE_CONFIG" "$MD_CASE_OUTPUT"

# Each subprocess starts with one coherent module/Conda environment. Amber's
# Python and the RDKit/ParmEd environment need not be loaded simultaneously.
MD_PREP_STAGE=structures
trap 'printf "Preparation failed during %s; see the Slurm log above.\n" "$MD_PREP_STAGE" >&2' ERR
(
    source "$MD_ROOT/env.sh" prep
    # Coordinate push-off uses SciPy. Older installations of this bundle did
    # not explicitly request it. Add only this missing preparation dependency
    # on the compute node, without rebuilding GROMACS or upgrading the env.
    if ! python3 -c 'import scipy.optimize, scipy.spatial' >/dev/null 2>&1; then
        printf 'Adding the missing SciPy coordinate-preparation dependency.\n'
        conda install --yes --prefix "$MD_PREP_PREFIX" --override-channels \
            --channel conda-forge --solver libmamba --freeze-installed scipy
        python3 -c 'import scipy.optimize, scipy.spatial'
        mkdir -p "$MD_ROOT/software-info"
        conda list --prefix "$MD_PREP_PREFIX" --explicit > "$MD_ROOT/software-info/prep-explicit.txt"
    fi
    python3 "$MD_ROOT/tools/prepare_case.py" structures --config "$MD_CASE_CONFIG" --output "$MD_CASE_OUTPUT"
)
MD_PREP_STAGE=parameters
(
    source "$MD_ROOT/env.sh" amber
    python3 "$MD_ROOT/tools/prepare_case.py" parameters --config "$MD_CASE_CONFIG" --output "$MD_CASE_OUTPUT"
)
MD_PREP_STAGE=assemble
(
    source "$MD_ROOT/env.sh" prep
    python3 "$MD_ROOT/tools/prepare_case.py" assemble --config "$MD_CASE_CONFIG" --output "$MD_CASE_OUTPUT"
)
rmdir -- "$TMPDIR" 2>/dev/null || true
printf '\nPreparation finished. Submit the CPU test with:\n'
printf 'sbatch %q %q\n' "$MD_ROOT/run_poc.sbatch" "$MD_CASE_OUTPUT"
