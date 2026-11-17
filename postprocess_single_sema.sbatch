#!/bin/bash
#SBATCH --job-name=sema_check
#SBATCH --account=ardekani
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=sema_check_%j.out
set -Eeuo pipefail
[[ $# == 1 ]] || { printf 'Usage: sbatch postprocess_single_sema.sbatch RUN_DIRECTORY\n' >&2; exit 2; }
export MD_ROOT="${MD_ROOT:-$HOME/MD}"
SEMA_RUN=$(realpath -- "$1")
[[ -f "$SEMA_RUN/input/manifest.json" ]] || { printf 'Missing single-semaglutide manifest.\n' >&2; exit 1; }
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
source "$MD_ROOT/env.sh" prep
SEMA_PYTHON="$MD_PREP_PREFIX/bin/python3"
source "$MD_ROOT/env.sh" md
"$SEMA_PYTHON" -E -s -u "$MD_ROOT/tools/analyze_single_sema.py" --run-dir "$SEMA_RUN" --gmx "$GMX"
