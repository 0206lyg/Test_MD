#!/bin/bash
#SBATCH --job-name=wow_sema_dpd
#SBATCH --account=ardekani
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --signal=USR1@300
#SBATCH --output=wow_sema_dpd_%j.out

set -Eeuo pipefail
usage() {
    cat <<'EOF'
Usage:
  sbatch run_wow_sema.sbatch
  sbatch run_wow_sema.sbatch --resume /scratch/negishi/USER/wow_sema_JOBID_UTC
  sbatch run_wow_sema.sbatch --check-input
A fresh submission uses the distributed 778,156-bead input, 1 MPI rank, and 8 CPU
threads by default. Results and a copy of the code and inputs are stored together.
Optional controls: --hold-steps, --max-steps, --seed, --sample-every, --chunk,
                   --dt, --target, --evap-rate (each takes one value).
On resume the original settings are reused. To extend a step cap, pass --max-steps N.
EOF
}
[[ ${1:-} != --help && ${1:-} != -h ]] || { usage; exit 0; }
[[ -n ${SLURM_JOB_ID:-} && ${SLURM_NTASKS:-1} == 1 ]] || {
    usage >&2; echo 'Submit with sbatch and one MPI task.' >&2; exit 2;
}
export MD_ROOT="${MD_ROOT:-$HOME/MD}"
source "$MD_ROOT/config.sh"
WOW_RESUME=
WOW_CHECK=0
WOW_REQUEST=()
while (( $# )); do
    case "$1" in
        --resume)
            [[ $# -ge 2 && -n $2 && -z $WOW_RESUME ]] || { usage >&2; exit 2; }
            WOW_RESUME="$2"; shift 2 ;;
        --check-input) WOW_CHECK=1; shift ;;
        --hold-steps|--max-steps|--seed|--sample-every|--chunk|--dt|--target|--evap-rate)
            (( $# >= 2 )) && [[ $2 != --* ]] || { usage >&2; exit 2; }
            WOW_REQUEST+=("$1" "$2"); shift 2 ;;
        *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done

WOW_CODE="$MD_ROOT/scripts/wow_dpd"
if [[ -n $WOW_RESUME ]]; then
    (( ! WOW_CHECK )) || { echo '--resume and --check-input are separate operations.' >&2; exit 2; }
    [[ $WOW_RESUME == /* ]] || { echo 'Use an absolute run path with --resume.' >&2; exit 2; }
    WOW_RUN="$(cd -- "$WOW_RESUME" && pwd -P)"
    [[ -f $WOW_RUN/WOW_SEMA_RUN && -f $WOW_RUN/requested_options.nul ]] || {
        echo "Missing W/O/W run marker or options in $WOW_RUN" >&2; exit 1;
    }
    WOW_CODE="$WOW_RUN/code/wow_dpd"
fi
source "$MD_ROOT/env.sh" dpd
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}" OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
WOW_LAUNCH=(env -u SLURM_TRES_PER_TASK srun --mpi=none --ntasks=1 --cpus-per-task="$OMP_NUM_THREADS" --cpu-bind=cores)
if (( WOW_CHECK )); then
    exec "${WOW_LAUNCH[@]}" "$MD_DPD_PYTHON" -E -s -u "$WOW_CODE/run_extraction.py" \
        --initial "$WOW_CODE/initial" --interactions "$WOW_CODE/interactions.json" \
        --threads "$OMP_NUM_THREADS" --check-input "${WOW_REQUEST[@]}"
fi
if [[ -z $WOW_RESUME ]]; then
    [[ -d $MD_SCRATCH_ROOT && -w $MD_SCRATCH_ROOT ]] || { echo "Scratch is unavailable: $MD_SCRATCH_ROOT" >&2; exit 1; }
    WOW_RUN="$MD_SCRATCH_ROOT/wow_sema_${SLURM_JOB_ID}_$(date -u +%Y%m%dT%H%M%SZ)"
    mkdir -- "$WOW_RUN"
    mkdir -- "$WOW_RUN/code" "$WOW_RUN/tmp"
    cp -R -- "$WOW_CODE" "$WOW_RUN/code/wow_dpd"
    cp -- "${BASH_SOURCE[0]}" "$WOW_RUN/code/run_wow_sema.sbatch"
    cp -- "$MD_ROOT/config.sh" "$MD_ROOT/env.sh" "$WOW_RUN/code/"
    WOW_CODE="$WOW_RUN/code/wow_dpd"
    : > "$WOW_RUN/requested_options.nul"
    if (( ${#WOW_REQUEST[@]} )); then printf '%s\0' "${WOW_REQUEST[@]}" > "$WOW_RUN/requested_options.nul"; fi
    date -u +%FT%TZ > "$WOW_RUN/WOW_SEMA_RUN"
else
    for ((i=0; i<${#WOW_REQUEST[@]}; i+=2)); do
        [[ ${WOW_REQUEST[i]} == --max-steps ]] || {
            echo 'A resumed run only accepts --max-steps; its original physics settings are restored.' >&2; exit 2;
        }
    done
    WOW_SAVED=()
    mapfile -d '' -t WOW_SAVED < "$WOW_RUN/requested_options.nul"
    WOW_REQUEST=("${WOW_SAVED[@]}" "${WOW_REQUEST[@]}")
fi
export TMPDIR="$WOW_RUN/tmp"
printf 'Results directory: %s\n' "$WOW_RUN"
cd -- "$WOW_RUN"
WOW_END="${SLURM_JOB_END_TIME:-}"
if [[ ! $WOW_END =~ ^[0-9]+$ ]]; then
    WOW_RECORD=$(scontrol show job -o "$SLURM_JOB_ID")
    [[ $WOW_RECORD =~ (^|[[:space:]])EndTime=([^[:space:]]+) ]] || { echo 'Cannot determine Slurm end time.' >&2; exit 1; }
    WOW_END=$(date -d "${BASH_REMATCH[2]}" +%s)
fi
WOW_BUDGET=$((WOW_END - $(date +%s) - 300))
(( WOW_BUDGET > 0 )) || { echo 'Less than five minutes remain; resubmit with a longer allocation.' >&2; exit 1; }
WOW_CMD=("${WOW_LAUNCH[@]}" "$MD_DPD_PYTHON" -E -s -u "$WOW_CODE/run_extraction.py"
    --initial "$WOW_CODE/initial" --interactions "$WOW_CODE/interactions.json"
    --out "$WOW_RUN/result" --threads "$OMP_NUM_THREADS" --time-budget "$WOW_BUDGET" "${WOW_REQUEST[@]}")
printf '%q ' "${WOW_CMD[@]}" > "$WOW_RUN/command_${SLURM_JOB_ID}.txt"
printf '\n' >> "$WOW_RUN/command_${SLURM_JOB_ID}.txt"
"${WOW_CMD[@]}"
printf '\nRun state: %s/result/state.json\n' "$WOW_RUN"
printf 'Resume a paused run: sbatch run_wow_sema.sbatch --resume %q\n' "$WOW_RUN"
