#!/bin/bash
#SBATCH --job-name=md_cpu_poc
#SBATCH --account=ardekani
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=md_cpu_poc_%j.out

# Usage: sbatch ~/MD/run_cpu.sbatch INPUT.tpr [additional mdrun options]
# Relative INPUT.tpr paths are resolved from the submission directory.
set -Eeuo pipefail
export MD_ROOT="${MD_ROOT:-$HOME/MD}"
[[ -n "${SLURM_JOB_ID:-}" ]] || { printf 'Use sbatch to run this script.\n' >&2; exit 1; }
[[ $# -ge 1 ]] || { printf 'Usage: sbatch run_cpu.sbatch INPUT.tpr [mdrun options]\n' >&2; exit 2; }
[[ -f "$1" ]] || { printf 'Missing TPR: %s\n' "$1" >&2; exit 1; }
MD_TPR="$(realpath -- "$1")"
shift
source "$MD_ROOT/env.sh" md
[[ "${SLURM_NTASKS:-1}" == 1 ]] || { printf 'Use --ntasks=1 for this thread-MPI build.\n' >&2; exit 1; }

MD_CPUS="${SLURM_CPUS_PER_TASK:-8}"
MD_NTMPI="${MD_NTMPI:-1}"
[[ "$MD_NTMPI" =~ ^[1-9][0-9]*$ ]] || { printf 'MD_NTMPI must be a positive integer.\n' >&2; exit 1; }
(( MD_CPUS % MD_NTMPI == 0 )) || { printf 'CPU count must be divisible by MD_NTMPI.\n' >&2; exit 1; }
export OMP_NUM_THREADS="$(( MD_CPUS / MD_NTMPI ))"

# All large simulation outputs are written in scratch. The Slurm output path
# above remains relative to the sbatch submission directory.
[[ -d "$MD_SCRATCH_ROOT" && -w "$MD_SCRATCH_ROOT" ]] || { printf 'Scratch directory missing or not writable: %s\n' "$MD_SCRATCH_ROOT" >&2; exit 1; }
MD_RUN_DIR="$MD_SCRATCH_ROOT/md_${SLURM_JOB_ID}_$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -- "$MD_RUN_DIR"
cp -- "$MD_TPR" "$MD_RUN_DIR/input.tpr"
cp -- "$MD_ROOT/config.sh" "$MD_RUN_DIR/config.sh"
cp -- "$MD_ROOT/run_cpu.sbatch" "$MD_RUN_DIR/run_cpu.sbatch"
printf '%s\n' "$MD_TPR" > "$MD_RUN_DIR/source_tpr.txt"
"$GMX" --version > "$MD_RUN_DIR/gromacs-version.txt" 2>&1
printf 'Results directory: %s\n' "$MD_RUN_DIR"
cd -- "$MD_RUN_DIR"
MD_COMMAND=("$GMX" mdrun -s input.tpr -deffnm md \
    -ntmpi "$MD_NTMPI" -ntomp "$OMP_NUM_THREADS" \
    -nb cpu -pme cpu -bonded cpu -update cpu \
    -cpt 10 -maxh "${MD_MAXH:-5.5}" "$@")
printf '%q ' srun --ntasks=1 --cpu-bind=cores "${MD_COMMAND[@]}" > command.sh
printf '\n' >> command.sh
exec srun --ntasks=1 --cpu-bind=cores "${MD_COMMAND[@]}"
