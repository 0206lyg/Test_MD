#!/bin/bash
#SBATCH --job-name=md_single_sema
#SBATCH --account=ardekani
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=md_single_sema_%j.out

# Each submission creates a NEW run. Existing Test_MD installations are reused.
# Default: one semaglutide, neutralizing ions + 0.150 M NaCl, 1 ns production.
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Usage: sbatch run_single_sema.sbatch [preparation options]
Defaults: --production-ns 1 --nvt-ps 50 --npt-ps 200
          --temperature 298.15 --seed 1729

Each submission starts a fresh calculation in MD_SCRATCH_ROOT.
The six-hour allocation is shared by preparation, all MD stages, and analysis.
An incomplete run is labelled PARTIAL and its available trajectory is analyzed.
This short run checks implementation and initial behavior, not equilibrium or
experimental force-field validation.
EOF
}
if [[ ${1:-} == --help || ${1:-} == -h ]]; then usage; exit 0; fi
[[ -n ${SLURM_JOB_ID:-} ]] || { usage >&2; printf 'Submit this script with sbatch.\n' >&2; exit 2; }
[[ ${SLURM_NTASKS:-1} == 1 ]] || { printf 'Use --ntasks=1 for the existing thread-MPI installation.\n' >&2; exit 2; }
for MD_ARG in "$@"; do
    case "$MD_ARG" in
        --output|--output=*|--assets|--assets=*|--gmx|--gmx=*|--continue|--continue=*)
            printf 'The launcher owns output/assets/gmx; each submission starts fresh.\n' >&2; exit 2 ;;
    esac
done

export MD_ROOT="${MD_ROOT:-$HOME/MD}"
[[ -f "$MD_ROOT/config.sh" && -f "$MD_ROOT/env.sh" ]] || {
    printf 'Missing config.sh or env.sh under %s. Set MD_ROOT to the existing Test_MD installation.\n' "$MD_ROOT" >&2; exit 1;
}
source "$MD_ROOT/config.sh"
MD_STAGE=startup
MD_RUN_DIR=not_created
MD_PYTHON=
trap 'printf "Single-semaglutide run failed during %s. Inspect %s and the Slurm output.\n" "$MD_STAGE" "$MD_RUN_DIR" >&2' ERR
export PYTHONUNBUFFERED=1 PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
MD_STARTED=$(date +%s)
MD_DEADLINE=${SLURM_JOB_END_TIME:-$((MD_STARTED + 21600))}
[[ $MD_DEADLINE =~ ^[0-9]+$ ]] && (( MD_DEADLINE > MD_STARTED )) || {
    printf 'SLURM_JOB_END_TIME is invalid or already expired.\n' >&2; exit 1;
}
MD_CPUS=${SLURM_CPUS_PER_TASK:-8}
[[ $MD_CPUS =~ ^[1-9][0-9]*$ ]] || { printf 'Invalid allocated CPU count.\n' >&2; exit 1; }
for MD_FILE in tools/prepare_single_sema.py tools/analyze_single_sema.py; do
    [[ -f "$MD_ROOT/$MD_FILE" ]] || { printf 'Missing overlay file: %s/%s\n' "$MD_ROOT" "$MD_FILE" >&2; exit 1; }
done
[[ -d "$MD_ROOT/cases/single_sema" ]] || { printf 'Missing cases/single_sema assets. Extract the overlay in %s.\n' "$MD_ROOT" >&2; exit 1; }
[[ -d "$MD_SCRATCH_ROOT" && -w "$MD_SCRATCH_ROOT" ]] || {
    printf 'Scratch directory is not writable: %s\n' "$MD_SCRATCH_ROOT" >&2; exit 1;
}

MD_RUN_DIR="$MD_SCRATCH_ROOT/single_sema_${SLURM_JOB_ID}_$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -- "$MD_RUN_DIR"
mkdir -- "$MD_RUN_DIR/code" "$MD_RUN_DIR/code/tools" "$MD_RUN_DIR/code/cases" "$MD_RUN_DIR/tmp"
cp -- "${BASH_SOURCE[0]}" "$MD_RUN_DIR/code/run_single_sema.sbatch"
cp -- "$MD_ROOT/config.sh" "$MD_ROOT/env.sh" "$MD_RUN_DIR/code/"
cp -- "$MD_ROOT/tools/prepare_single_sema.py" "$MD_ROOT/tools/analyze_single_sema.py" "$MD_RUN_DIR/code/tools/"
cp -R -- "$MD_ROOT/cases/single_sema" "$MD_RUN_DIR/code/cases/"
printf '%q ' "$@" > "$MD_RUN_DIR/requested_options.txt"
printf '\n' >> "$MD_RUN_DIR/requested_options.txt"
printf 'Single semaglutide screening run, created %s\n' "$(date -u +%FT%TZ)" > "$MD_RUN_DIR/SINGLE_SEMA_RUN"
cd -- "$MD_RUN_DIR"
export TMPDIR="$MD_RUN_DIR/tmp"
printf 'Results directory: %s\n' "$MD_RUN_DIR"

MD_STAGE=environment
# Keep preparation Python explicit: module changes must not select old Python.
source "$MD_ROOT/env.sh" prep
MD_PYTHON="$MD_PREP_PREFIX/bin/python3"
[[ -x "$MD_PYTHON" ]] || { printf 'Preparation Python missing: %s\n' "$MD_PYTHON" >&2; exit 1; }
"$MD_PYTHON" -E -s -c 'import sys, numpy
if sys.version_info < (3, 8):
    raise SystemExit("Python >= 3.8 is required")
print("Preparation Python:", sys.executable, sys.version.split()[0])'
source "$MD_ROOT/env.sh" md
export OMP_NUM_THREADS="$MD_CPUS"
# Site-adjusted allocations can leave a stale duplicate TRES request.
MD_LAUNCH=(env -u SLURM_TRES_PER_TASK srun --ntasks=1 --cpus-per-task="$MD_CPUS" --cpu-bind=cores)
"$GMX" --version > gromacs-version.txt 2>&1

MD_STAGE=preparation
"$MD_PYTHON" -E -s -u code/tools/prepare_single_sema.py \
    --output input --gmx "$GMX" --assets code/cases/single_sema "$@"
for MD_FILE in system.gro topol.top index.ndx manifest.json mdp/em.mdp mdp/nvt.mdp mdp/npt.mdp mdp/production.mdp; do
    [[ -s "input/$MD_FILE" ]] || { printf 'Preparation did not produce input/%s\n' "$MD_FILE" >&2; exit 1; }
done

MD_COMPLETED=()
write_status() {
    "$MD_PYTHON" -E -s - "$1" "$MD_STAGE" "$2" "${MD_COMPLETED[@]}" <<'PY'
import json, sys
from pathlib import Path
Path("run_status.json").write_text(json.dumps({
    "status": sys.argv[1], "stage": sys.argv[2], "reason": sys.argv[3],
    "completed_stages": sys.argv[4:],
    "scope": "short implementation/initial-behavior screen; not full force-field validation"
}, indent=2) + "\n")
PY
}
run_analysis() {
    printf '\nAnalyzing available results...\n'
    "$MD_PYTHON" -E -s -u code/tools/analyze_single_sema.py --run-dir "$MD_RUN_DIR" --gmx "$GMX"
}
finish_partial() {
    write_status partial "$1"
    printf '\nPARTIAL RUN: %s\n' "$1"
    run_analysis
    printf '\nPARTIAL results retained: %s\n' "$MD_RUN_DIR"
    printf 'See run_status.json and analysis/. This is not a completed production run.\n'
    exit 0
}

MD_PREVIOUS_GRO=input/system.gro
MD_PREVIOUS_CPT=
for MD_STAGE in em nvt npt production; do
    # Reserve five minutes for checkpoint output and analysis, for all stages combined.
    MD_REMAINING=$((MD_DEADLINE - $(date +%s) - 300))
    (( MD_REMAINING >= 120 )) || finish_partial "Wall-time budget reached before $MD_STAGE."
    MD_GROMPP=("$GMX" grompp -f "input/mdp/$MD_STAGE.mdp"
        -c "$MD_PREVIOUS_GRO" -r "$MD_PREVIOUS_GRO" -p input/topol.top
        -n input/index.ndx -o "$MD_STAGE.tpr" -po "$MD_STAGE.out.mdp")
    [[ -z $MD_PREVIOUS_CPT ]] || MD_GROMPP+=(-t "$MD_PREVIOUS_CPT")
    printf '\nPreparing stage: %s\n' "$MD_STAGE"
    if ! "${MD_GROMPP[@]}" > "$MD_STAGE.grompp.log" 2>&1; then
        cat "$MD_STAGE.grompp.log" >&2
        write_status failed "grompp failed in $MD_STAGE; inspect $MD_STAGE.grompp.log."
        exit 1
    fi
    MD_REMAINING=$((MD_DEADLINE - $(date +%s) - 300))
    (( MD_REMAINING >= 120 )) || finish_partial "Wall-time budget reached before $MD_STAGE mdrun."
    MD_MAXH=$(awk -v seconds="$MD_REMAINING" 'BEGIN {printf "%.6f", seconds/3600}')
    MD_COMMAND=("$GMX" mdrun -s "$MD_STAGE.tpr" -deffnm "$MD_STAGE"
        -ntmpi 1 -ntomp "$MD_CPUS" -nb cpu -pme cpu -bonded cpu -update cpu
        -cpt 5 -maxh "$MD_MAXH")
    printf '\nRunning stage: %s; CPU threads: %s; shared MD budget: %s seconds\n' "$MD_STAGE" "$MD_CPUS" "$MD_REMAINING"
    printf '%q ' "${MD_LAUNCH[@]}" "${MD_COMMAND[@]}" >> commands.sh
    printf '\n' >> commands.sh
    if ! "${MD_LAUNCH[@]}" "${MD_COMMAND[@]}"; then
        write_status failed "mdrun failed in $MD_STAGE; inspect $MD_STAGE.log."
        # Do not call a failed simulation successful just because files exist.
        run_analysis || true
        exit 1
    fi
    if [[ $MD_STAGE == em ]]; then
        if ! grep -q 'Steepest Descents converged to Fmax' em.log; then
            write_status failed 'Energy minimization did not reach its force tolerance; later MD was not started.'
            run_analysis || true
            printf 'Minimization did not reach its force tolerance. See em.log.\n' >&2
            exit 1
        fi
    else
        [[ -s $MD_STAGE.cpt ]] || { write_status failed "Missing $MD_STAGE.cpt"; exit 1; }
        MD_LAST_STEP=$("$GMX" dump -cp "$MD_STAGE.cpt" 2> "$MD_STAGE.checkpoint.log" | awk '$1=="step" && $2=="=" {print $3}')
        MD_TARGET_STEP=$(awk '$1=="nsteps" {n=$3} $1=="init-step" || $1=="init_step" {i=$3} END {if(n!="") print n+i}' "$MD_STAGE.out.mdp")
        [[ $MD_LAST_STEP =~ ^[0-9]+$ && $MD_TARGET_STEP =~ ^[0-9]+$ ]] || {
            write_status failed "Cannot verify final step for $MD_STAGE."; exit 1;
        }
        (( MD_LAST_STEP >= MD_TARGET_STEP )) || finish_partial "$MD_STAGE reached $MD_LAST_STEP / $MD_TARGET_STEP steps."
    fi
    [[ -s $MD_STAGE.gro ]] || { write_status failed "Missing $MD_STAGE.gro"; exit 1; }
    date -u +%FT%TZ > "$MD_STAGE.complete"
    MD_COMPLETED+=("$MD_STAGE")
    MD_PREVIOUS_GRO="$MD_STAGE.gro"
    [[ $MD_STAGE == em ]] || MD_PREVIOUS_CPT="$MD_STAGE.cpt"
done
write_status complete 'All requested MD steps completed; physical validation may remain incomplete.'
MD_STAGE=analysis
run_analysis
printf '\nREQUESTED SHORT RUN COMPLETE: %s\n' "$MD_RUN_DIR"
printf 'Read analysis/ for the five validation items; completion does not mean all five passed.\n'
grep 'Performance:' production.log || true
