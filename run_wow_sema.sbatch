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

# EDIT THESE INPUTS, then submit with: sbatch run_wow_sema.sbatch
DROPLET_DIAMETER_NM=39.2   # Initial W1+O diameter at the O/W2 boundary
W1_COUNT=10                # Number of equal-volume inner water pockets (1-10)
INIT_SEED=20261120         # Initial coordinates and drug allocation
CHECK_INPUT=0              # 1 = generate, preview and check step 0 only

# DPD and DCM extraction controls. Seeds for coordinates and dynamics are independent.
DYNAMICS_SEED=48727
DT=0.01
HOLD_STEPS=1000
CHUNK=100                  # Steps per DCM exchange decision
SAMPLE_EVERY=1000          # Steps per trajectory/composition sample
MAX_STEPS=2000000          # Cumulative extraction steps, excluding the hold
TARGET=0.05                # In-particle DCM/(DCM + PLGA) mass fraction
EVAP_RATE=0.02             # Per reduced DPD time unit

# To resume, set an absolute existing run directory. Only the optional step cap applies.
RESUME_DIR=""
RESUME_MAX_STEPS=""

set -Eeuo pipefail
usage() {
    cat <<'EOF'
Edit the inputs at the top and submit: sbatch run_wow_sema.sbatch
Optional overrides for a fresh run:
  --droplet-diameter-nm D  --w1-count N  --init-seed S  --check-input
  --hold-steps N  --max-steps N  --seed S  --sample-every N
  --chunk N  --dt DT  --target F  --evap-rate K
Resume the saved geometry, code and settings:
  sbatch run_wow_sema.sbatch --resume /absolute/run/path [--max-steps N]
EOF
}
[[ ${1:-} != --help && ${1:-} != -h ]] || { usage; exit 0; }
[[ -n ${SLURM_JOB_ID:-} && ${SLURM_NTASKS:-1} == 1 ]] || {
    usage >&2; echo 'Submit with sbatch and one MPI task.' >&2; exit 2;
}
WOW_OVERRIDES=()
WOW_FRESH_OPTIONS=0
while (( $# )); do
    if [[ $1 != --check-input ]]; then
        [[ $# -ge 2 && -n $2 && $2 != --* ]] || { usage >&2; exit 2; }
    fi
    case "$1" in
        --resume) RESUME_DIR="$2" ;;
        --droplet-diameter-nm) DROPLET_DIAMETER_NM="$2"; WOW_FRESH_OPTIONS=1 ;;
        --w1-count) W1_COUNT="$2"; WOW_FRESH_OPTIONS=1 ;;
        --init-seed) INIT_SEED="$2"; WOW_FRESH_OPTIONS=1 ;;
        --check-input) CHECK_INPUT=1; WOW_FRESH_OPTIONS=1; shift; continue ;;
        --hold-steps|--seed|--sample-every|--chunk|--dt|--target|--evap-rate)
            WOW_OVERRIDES+=("$1" "$2"); WOW_FRESH_OPTIONS=1 ;;
        --max-steps) WOW_OVERRIDES+=("$1" "$2") ;;
        *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
    shift 2
done
[[ $CHECK_INPUT == 0 || $CHECK_INPUT == 1 ]] || { echo 'Set CHECK_INPUT to 0 or 1.' >&2; exit 2; }

export MD_ROOT="${MD_ROOT:-$HOME/MD}"
source "$MD_ROOT/config.sh"
source "$MD_ROOT/env.sh" dpd
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}" OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
WOW_LAUNCH=(env -u SLURM_TRES_PER_TASK srun --mpi=none --ntasks=1 --cpus-per-task="$OMP_NUM_THREADS" --cpu-bind=cores)
if [[ -z $RESUME_DIR ]]; then
    [[ -d $MD_SCRATCH_ROOT && -w $MD_SCRATCH_ROOT ]] || { echo "Scratch is unavailable: $MD_SCRATCH_ROOT" >&2; exit 1; }
    WOW_RUN="$MD_SCRATCH_ROOT/wow_sema_${SLURM_JOB_ID}_$(date -u +%Y%m%dT%H%M%SZ)"
    mkdir -- "$WOW_RUN"
    WOW_CODE="$WOW_RUN/code/wow_dpd"
    mkdir -p -- "$WOW_CODE" "$WOW_RUN/tmp"
    export TMPDIR="$WOW_RUN/tmp"
    printf 'Scratch run directory: %s\n' "$WOW_RUN"
    # A previous overlay may have left a large initial/ directory in the checkout.
    tar -C "$MD_ROOT/scripts/wow_dpd" --exclude='./initial' --exclude='__pycache__' -cf - . | tar -C "$WOW_CODE" -xf -
    cp -- "${BASH_SOURCE[0]}" "$WOW_RUN/code/run_wow_sema.sbatch"
    cp -- "$MD_ROOT/config.sh" "$MD_ROOT/env.sh" "$WOW_RUN/code/"
    WOW_REQUEST=(--hold-steps "$HOLD_STEPS" --max-steps "$MAX_STEPS" --seed "$DYNAMICS_SEED"
                 --sample-every "$SAMPLE_EVERY" --chunk "$CHUNK" --dt "$DT" --target "$TARGET"
                 --evap-rate "$EVAP_RATE" "${WOW_OVERRIDES[@]}")
    printf '%s\0' "${WOW_REQUEST[@]}" > "$WOW_RUN/requested_options.nul"
    "${WOW_LAUNCH[@]}" "$MD_DPD_PYTHON" -E -s -u "$WOW_CODE/build_initial.py" \
        --out "$WOW_CODE/initial" --droplet-diameter-nm "$DROPLET_DIAMETER_NM" \
        --w1-count "$W1_COUNT" --init-seed "$INIT_SEED"
    printf 'Initial coordinates, composition and preview: %s/initial/\n' "$WOW_CODE"
    if (( CHECK_INPUT )); then
        exec "${WOW_LAUNCH[@]}" "$MD_DPD_PYTHON" -E -s -u "$WOW_CODE/run_extraction.py" \
            --initial "$WOW_CODE/initial" --interactions "$WOW_CODE/interactions.json" \
            --threads "$OMP_NUM_THREADS" --check-input "${WOW_REQUEST[@]}"
    fi
    date -u +%FT%TZ > "$WOW_RUN/WOW_SEMA_RUN"
else
    (( ! WOW_FRESH_OPTIONS )) || { echo 'Only --max-steps can override settings on resume.' >&2; exit 2; }
    [[ $RESUME_DIR == /* ]] || { echo 'Use an absolute run path to resume.' >&2; exit 2; }
    WOW_RUN="$(cd -- "$RESUME_DIR" && pwd -P)"
    [[ -f $WOW_RUN/WOW_SEMA_RUN && -f $WOW_RUN/requested_options.nul ]] || {
        echo "Missing W/O/W run marker or options in $WOW_RUN" >&2; exit 1;
    }
    WOW_CODE="$WOW_RUN/code/wow_dpd"
    mapfile -d '' -t WOW_REQUEST < "$WOW_RUN/requested_options.nul"
    if [[ -n $RESUME_MAX_STEPS ]]; then WOW_REQUEST+=(--max-steps "$RESUME_MAX_STEPS"); fi
    WOW_REQUEST+=("${WOW_OVERRIDES[@]}")
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
