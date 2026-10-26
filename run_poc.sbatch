#!/bin/bash
#SBATCH --job-name=md_poc
#SBATCH --account=ardekani
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=md_poc_%j.out

# Usage: sbatch run_poc.sbatch CASE_INPUT_DIR
#        sbatch run_poc.sbatch --restart EXISTING_SCRATCH_RUN_DIR
set -Eeuo pipefail
export MD_ROOT="${MD_ROOT:-$HOME/MD}"
[[ -n "${SLURM_JOB_ID:-}" ]] || { printf 'Use sbatch to run this script.\n' >&2; exit 1; }
[[ "${SLURM_NTASKS:-1}" == 1 ]] || { printf 'Use --ntasks=1 for this thread-MPI build.\n' >&2; exit 1; }
source "$MD_ROOT/env.sh" md
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
[[ "$OMP_NUM_THREADS" =~ ^[1-9][0-9]*$ ]] || { printf 'Invalid CPU count.\n' >&2; exit 1; }
# A site adjustment can leave SLURM_TRES_PER_TASK with an older CPU count.
# Remove that duplicate request only for this CPU-only job step, and pass the
# allocation's current CPU count explicitly to match mdrun's OpenMP threads.
MD_LAUNCH=(env -u SLURM_TRES_PER_TASK srun --ntasks=1
    --cpus-per-task="$OMP_NUM_THREADS" --cpu-bind=cores)

if [[ $# == 2 && "$1" == --restart ]]; then
    MD_RUN_DIR="$(realpath -- "$2")"
    [[ -f "$MD_RUN_DIR/input/system.gro" && -f "$MD_RUN_DIR/source_case.txt" ]] || { printf 'Not a PoC run directory: %s\n' "$MD_RUN_DIR" >&2; exit 1; }
elif [[ $# == 1 && "$1" != --* ]]; then
    MD_CASE_INPUT="$(realpath -- "$1")"
    for MD_REQUIRED in READY.json system.gro topol.top mdp/em.mdp mdp/nvt.mdp mdp/npt.mdp mdp/benchmark.mdp; do
        [[ -f "$MD_CASE_INPUT/$MD_REQUIRED" ]] || { printf 'Missing prepared input: %s\n' "$MD_CASE_INPUT/$MD_REQUIRED" >&2; exit 1; }
    done
    [[ -d "$MD_SCRATCH_ROOT" && -w "$MD_SCRATCH_ROOT" ]] || { printf 'Scratch directory missing or not writable: %s\n' "$MD_SCRATCH_ROOT" >&2; exit 1; }
    MD_RUN_DIR="$MD_SCRATCH_ROOT/md_${SLURM_JOB_ID}_$(date -u +%Y%m%dT%H%M%SZ)"
    mkdir -- "$MD_RUN_DIR"
    cp -a -- "$MD_CASE_INPUT" "$MD_RUN_DIR/input"
    cp -- "$MD_ROOT/run_poc.sbatch" "$MD_RUN_DIR/run_poc.sbatch"
    printf '%s\n' "$MD_CASE_INPUT" > "$MD_RUN_DIR/source_case.txt"
else
    printf 'Usage: sbatch run_poc.sbatch CASE_INPUT_DIR\n       sbatch run_poc.sbatch --restart RUN_DIR\n' >&2
    exit 2
fi
cd -- "$MD_RUN_DIR"
printf 'Results directory: %s\n' "$MD_RUN_DIR"
"$GMX" --version > "gromacs-version-${SLURM_JOB_ID}.txt" 2>&1
MD_STARTED="$(date +%s)"
# One shared wall-time budget for all stages, including grompp. Reserve ten
# minutes for checkpoint/output writing. Slurm supplies its actual end time.
MD_DEADLINE="${SLURM_JOB_END_TIME:-$((MD_STARTED + 21600))}"
[[ "$MD_DEADLINE" =~ ^[0-9]+$ ]] || { printf 'Invalid SLURM_JOB_END_TIME.\n' >&2; exit 1; }

pause_run() {
    printf '\n%s\nFiles and checkpoints retained in: %s\n' "$1" "$MD_RUN_DIR"
    printf 'Continue with: sbatch %q --restart %q\n' "$MD_ROOT/run_poc.sbatch" "$MD_RUN_DIR"
    exit 0
}
trap 'printf "PoC failed in %s; inspect %s and the Slurm output.\n" "${MD_STAGE:-startup}" "$MD_RUN_DIR" >&2' ERR

MD_PREVIOUS_GRO=input/system.gro
MD_PREVIOUS_CPT=
for MD_STAGE in em nvt npt benchmark; do
    if [[ ! -f "$MD_STAGE.complete" ]]; then
        MD_REMAINING=$((MD_DEADLINE - $(date +%s) - 600))
        (( MD_REMAINING >= 120 )) || pause_run "Wall-time budget reached before $MD_STAGE."
        if [[ ! -f "$MD_STAGE.tpr" ]]; then
            MD_GROMPP=("$GMX" grompp -f "input/mdp/$MD_STAGE.mdp" -c "$MD_PREVIOUS_GRO"
                -r input/system.gro -p input/topol.top -o "$MD_STAGE.tpr" -po "$MD_STAGE.out.mdp")
            if [[ -n "$MD_PREVIOUS_CPT" ]]; then MD_GROMPP+=(-t "$MD_PREVIOUS_CPT"); fi
            printf '\nPreparing stage: %s\n' "$MD_STAGE"
            "${MD_GROMPP[@]}" > "$MD_STAGE.grompp.log" 2>&1 || {
                cat "$MD_STAGE.grompp.log" >&2
                exit 1
            }
        fi
        MD_REMAINING=$((MD_DEADLINE - $(date +%s) - 600))
        (( MD_REMAINING >= 120 )) || pause_run "Wall-time budget reached before $MD_STAGE mdrun."
        MD_MAXH="$(awk -v seconds="$MD_REMAINING" 'BEGIN { printf "%.6f", seconds/3600 }')"
        MD_COMMAND=("$GMX" mdrun -s "$MD_STAGE.tpr" -deffnm "$MD_STAGE"
            -ntmpi 1 -ntomp "$OMP_NUM_THREADS" -nb cpu -pme cpu -bonded cpu -update cpu
            -cpt 5 -maxh "$MD_MAXH")
        if [[ -f "$MD_STAGE.cpt" ]]; then MD_COMMAND+=(-cpi "$MD_STAGE.cpt" -append); fi
        printf '\nRunning stage: %s (shared budget remaining: %s seconds)\n' "$MD_STAGE" "$MD_REMAINING"
        printf '%q ' "${MD_LAUNCH[@]}" "${MD_COMMAND[@]}" >> commands.sh
        printf '\n' >> commands.sh
        "${MD_LAUNCH[@]}" "${MD_COMMAND[@]}"

        if [[ "$MD_STAGE" == em ]]; then
            if ! grep -q 'Steepest Descents converged to Fmax' em.log; then
                printf 'Minimization has not converged; MD stages were not started. Inspect em.log.\n' >&2
                exit 1
            fi
        else
            # mdrun can return success after a wall-time stop. Verify the saved
            # step before advancing, so partial equilibration is never skipped.
            [[ -f "$MD_STAGE.cpt" ]] || { printf 'Missing checkpoint for %s.\n' "$MD_STAGE" >&2; exit 1; }
            MD_LAST_STEP="$("$GMX" dump -cp "$MD_STAGE.cpt" 2> "$MD_STAGE.checkpoint.log" | awk '$1 == "step" && $2 == "=" { print $3 }')"
            MD_TARGET_STEP="$(awk '$1 == "nsteps" { n=$3 } $1 == "init-step" || $1 == "init_step" { i=$3 } END { if (n != "") print n+i }' "$MD_STAGE.out.mdp")"
            [[ "$MD_LAST_STEP" =~ ^[0-9]+$ && "$MD_TARGET_STEP" =~ ^[0-9]+$ ]] || { printf 'Cannot verify completed step for %s.\n' "$MD_STAGE" >&2; exit 1; }
            (( MD_LAST_STEP >= MD_TARGET_STEP )) || pause_run "$MD_STAGE stopped at step $MD_LAST_STEP of $MD_TARGET_STEP."
        fi
        [[ -f "$MD_STAGE.gro" ]] || { printf 'Missing final coordinates for %s.\n' "$MD_STAGE" >&2; exit 1; }
        date -u +%Y-%m-%dT%H:%M:%SZ > "$MD_STAGE.complete"
    fi
    MD_PREVIOUS_GRO="$MD_STAGE.gro"
    if [[ "$MD_STAGE" != em ]]; then MD_PREVIOUS_CPT="$MD_STAGE.cpt"; fi
done
printf '\nPoC stages finished. Results: %s\n' "$MD_RUN_DIR"
grep 'Performance:' benchmark.log || true
printf 'This confirms completion of the configured short run, not equilibrium or a converged free energy.\n'
