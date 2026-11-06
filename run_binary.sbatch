#!/bin/bash
#SBATCH --job-name=md_binary
#SBATCH --account=ardekani
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=md_binary_%j.out

# Fresh calculation: sbatch run_binary.sbatch DMSO [preparation options]
# Explicit checkpoint continuation: sbatch run_binary.sbatch --continue RUN_DIR
set -Eeuo pipefail
export MD_ROOT="${MD_ROOT:-$HOME/MD}"
usage() {
    cat <<'EOF'
Usage: sbatch run_binary.sbatch DMSO|DMF|Acetone|Acetonitrile [options]
       sbatch run_binary.sbatch --continue RUN_DIR
Options: --box-nm 6 --organic-volume 1.5 --water-volume 10
         --temperature 298.15 --production-ns 20 --seed 1729
         --nvt-ps 100 --npt-ps 1000
Each solvent invocation creates a NEW calculation with a uniformly packed mixture.
EOF
}
if [[ ${1:-} == --help || ${1:-} == -h ]]; then usage; exit 0; fi
[[ -n ${SLURM_JOB_ID:-} ]] || { usage >&2; printf 'Submit with sbatch.\n' >&2; exit 2; }
[[ ${SLURM_NTASKS:-1} == 1 ]] || { printf 'Use --ntasks=1 for thread-MPI.\n' >&2; exit 2; }
[[ -f "$MD_ROOT/env.sh" ]] || { printf 'Missing %s/env.sh\n' "$MD_ROOT" >&2; exit 1; }
source "$MD_ROOT/config.sh"
MD_STAGE=startup
MD_RUN_DIR=not_created
trap 'printf "Binary run failed during %s. Inspect %s and the Slurm output.\n" "$MD_STAGE" "$MD_RUN_DIR" >&2' ERR
export PYTHONUNBUFFERED=1 PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
MD_STARTED=$(date +%s)
MD_DEADLINE=${SLURM_JOB_END_TIME:-$((MD_STARTED + 86400))}
[[ $MD_DEADLINE =~ ^[0-9]+$ ]] || { printf 'Invalid SLURM_JOB_END_TIME.\n' >&2; exit 1; }
MD_PREPARE=1
if [[ $# == 2 && $1 == --continue ]]; then
    MD_RUN_DIR=$(realpath -- "$2")
    [[ -f "$MD_RUN_DIR/BINARY_RUN" && -f "$MD_RUN_DIR/input/READY.json" ]] || {
        printf 'Continuation requires a prepared binary run directory.\n' >&2; exit 2;
    }
    MD_PREPARE=0
elif [[ $# -ge 1 && $1 != --* ]]; then
    MD_SOLVENT=$1
    shift
    # Restrict names before using them in a directory name. Python canonicalizes.
    case ${MD_SOLVENT,,} in dmso|dmf|acetone|acn|acetonitrile) ;;
        *) printf 'Choose DMSO, DMF, Acetone, or Acetonitrile.\n' >&2; exit 2 ;;
    esac
    # The output path and solvent belong to this launcher, not passthrough flags.
    for MD_ARG in "$@"; do
        case $MD_ARG in --output|--output=*|--solvent|--solvent=*)
            printf 'Pass the solvent as the first argument; output is assigned automatically.\n' >&2; exit 2 ;;
        esac
    done
    [[ -d "$MD_SCRATCH_ROOT" && -w "$MD_SCRATCH_ROOT" ]] || {
        printf 'Scratch directory is not writable: %s\n' "$MD_SCRATCH_ROOT" >&2; exit 1;
    }
    MD_RUN_DIR="$MD_SCRATCH_ROOT/binary_${MD_SOLVENT}_${SLURM_JOB_ID}_$(date -u +%Y%m%dT%H%M%SZ)"
    mkdir -- "$MD_RUN_DIR"
    mkdir -- "$MD_RUN_DIR/code" "$MD_RUN_DIR/code/tools"
    cp -- "${BASH_SOURCE[0]}" "$MD_RUN_DIR/code/run_binary.sbatch"
    cp -- "$MD_ROOT/env.sh" "$MD_ROOT/config.sh" "$MD_RUN_DIR/code/"
    for MD_FILE in prepare_binary.py analyze_binary.py prepare_case.py mdprep_forcefield.py; do
        cp -- "$MD_ROOT/tools/$MD_FILE" "$MD_RUN_DIR/code/tools/"
    done
    printf 'Binary solvent/water run; created %s\n' "$(date -u +%FT%TZ)" > "$MD_RUN_DIR/BINARY_RUN"
    printf '%q ' "$MD_SOLVENT" "$@" > "$MD_RUN_DIR/requested_options.txt"
    printf '\n' >> "$MD_RUN_DIR/requested_options.txt"
else
    usage >&2; exit 2
fi
cd -- "$MD_RUN_DIR"
printf 'Results directory: %s\n' "$MD_RUN_DIR"
export TMPDIR="$MD_RUN_DIR/tmp"
mkdir -p -- "$TMPDIR"
if (( MD_PREPARE )); then
    MD_STAGE=structures
    (
        export OMP_NUM_THREADS=1
        source "$MD_ROOT/env.sh" prep
        python3 code/tools/prepare_binary.py structures --solvent "$MD_SOLVENT" --output input "$@"
    )
    MD_STAGE=parameters
    (
        export OMP_NUM_THREADS=1
        source "$MD_ROOT/env.sh" amber
        python3 code/tools/prepare_binary.py parameters --output input
    )
    MD_STAGE=assemble
    (
        export OMP_NUM_THREADS=1
        source "$MD_ROOT/env.sh" prep
        python3 code/tools/prepare_binary.py assemble --output input
    )
fi
for MD_REQUIRED in READY.json manifest.json system.gro topol.top index.ndx mdp/em.mdp mdp/nvt.mdp mdp/npt.mdp mdp/production.mdp; do
    [[ -f input/$MD_REQUIRED ]] || { printf 'Missing input/%s\n' "$MD_REQUIRED" >&2; exit 1; }
done
source "$MD_ROOT/env.sh" md
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
[[ $OMP_NUM_THREADS =~ ^[1-9][0-9]*$ ]] || { printf 'Invalid CPU count.\n' >&2; exit 1; }
# Match the actual Slurm allocation, including a site-adjusted CPU count.
MD_LAUNCH=(env -u SLURM_TRES_PER_TASK srun --ntasks=1 --cpus-per-task="$OMP_NUM_THREADS" --cpu-bind=cores)
"$GMX" --version > "gromacs-version-${SLURM_JOB_ID}.txt" 2>&1
pause_run() {
    printf '\nNOT COMPLETE: %s\nCheckpoints retained: %s\n' "$1" "$MD_RUN_DIR"
    printf 'To continue this same calculation explicitly:\nsbatch %q --continue %q\n' "$MD_ROOT/run_binary.sbatch" "$MD_RUN_DIR"
    exit 0
}
MD_PREVIOUS_GRO=input/system.gro
MD_PREVIOUS_CPT=
for MD_STAGE in em nvt npt production; do
    if [[ ! -f $MD_STAGE.complete ]]; then
        MD_REMAINING=$((MD_DEADLINE - $(date +%s) - 600))
        (( MD_REMAINING >= 120 )) || pause_run "Wall-time budget reached before $MD_STAGE."
        if [[ ! -f $MD_STAGE.tpr ]]; then
            MD_GROMPP=("$GMX" grompp -f "input/mdp/$MD_STAGE.mdp" -c "$MD_PREVIOUS_GRO"
                -p input/topol.top -n input/index.ndx -o "$MD_STAGE.tpr" -po "$MD_STAGE.out.mdp")
            [[ -z $MD_PREVIOUS_CPT ]] || MD_GROMPP+=(-t "$MD_PREVIOUS_CPT")
            printf '\nPreparing stage: %s\n' "$MD_STAGE"
            "${MD_GROMPP[@]}" > "$MD_STAGE.grompp.log" 2>&1 || { cat "$MD_STAGE.grompp.log" >&2; exit 1; }
        fi
        MD_REMAINING=$((MD_DEADLINE - $(date +%s) - 600))
        (( MD_REMAINING >= 120 )) || pause_run "Wall-time budget reached before $MD_STAGE mdrun."
        MD_MAXH=$(awk -v seconds="$MD_REMAINING" 'BEGIN { printf "%.6f", seconds/3600 }')
        MD_COMMAND=("$GMX" mdrun -s "$MD_STAGE.tpr" -deffnm "$MD_STAGE"
            -ntmpi 1 -ntomp "$OMP_NUM_THREADS" -nb cpu -pme cpu -bonded cpu -update cpu
            -cpt 5 -maxh "$MD_MAXH")
        [[ ! -f $MD_STAGE.cpt ]] || MD_COMMAND+=(-cpi "$MD_STAGE.cpt" -append)
        printf '\nRunning stage: %s; allocated CPU threads: %s\n' "$MD_STAGE" "$OMP_NUM_THREADS"
        printf '%q ' "${MD_LAUNCH[@]}" "${MD_COMMAND[@]}" >> commands.sh
        printf '\n' >> commands.sh
        "${MD_LAUNCH[@]}" "${MD_COMMAND[@]}"
        if [[ $MD_STAGE == em ]]; then
            if ! grep -q 'Steepest Descents converged to Fmax' em.log; then
                printf 'Minimization did not reach its force tolerance; MD not advanced. See em.log.\n' >&2
                exit 1
            fi
        else
            [[ -f $MD_STAGE.cpt ]] || { printf 'Missing checkpoint: %s\n' "$MD_STAGE" >&2; exit 1; }
            MD_LAST_STEP=$("$GMX" dump -cp "$MD_STAGE.cpt" 2> "$MD_STAGE.checkpoint.log" | awk '$1 == "step" && $2 == "=" { print $3 }')
            MD_TARGET_STEP=$(awk '$1 == "nsteps" {n=$3} $1 == "init-step" || $1 == "init_step" {i=$3} END {if(n != "") print n+i}' "$MD_STAGE.out.mdp")
            [[ $MD_LAST_STEP =~ ^[0-9]+$ && $MD_TARGET_STEP =~ ^[0-9]+$ ]] || {
                printf 'Cannot verify final step for %s.\n' "$MD_STAGE" >&2; exit 1;
            }
            (( MD_LAST_STEP >= MD_TARGET_STEP )) || pause_run "$MD_STAGE reached $MD_LAST_STEP / $MD_TARGET_STEP steps."
        fi
        [[ -f $MD_STAGE.gro ]] || { printf 'Missing final coordinates: %s\n' "$MD_STAGE" >&2; exit 1; }
        date -u +%FT%TZ > "$MD_STAGE.complete"
    fi
    MD_PREVIOUS_GRO="$MD_STAGE.gro"
    [[ $MD_STAGE == em ]] || MD_PREVIOUS_CPT="$MD_STAGE.cpt"
done
MD_STAGE=analysis
mkdir -p analysis
# A water molecule contributes once, just like an organic molecule. The Python
# analysis reconstructs molecular centers; it never substitutes atom fractions.
for MD_TRAJECTORY in nvt npt production; do
    MD_ANALYSIS_DIR="analysis/$MD_TRAJECTORY"
    mkdir -p "$MD_ANALYSIS_DIR"
    printf 'System\n' | "$GMX" trjconv -f "$MD_TRAJECTORY.xtc" -s "$MD_TRAJECTORY.tpr" \
        -n input/index.ndx -o "$MD_ANALYSIS_DIR/frames.gro" -dt 100 -pbc mol \
        > "$MD_ANALYSIS_DIR/trjconv.log" 2>&1
    (
        export OMP_NUM_THREADS=1
        source "$MD_ROOT/env.sh" prep
        python3 code/tools/analyze_binary.py --trajectory "$MD_ANALYSIS_DIR/frames.gro" \
            --manifest input/manifest.json --output "$MD_ANALYSIS_DIR"
    )
    # The XTC is the original trajectory; sampled GRO frames are reproducible.
    rm -- "$MD_ANALYSIS_DIR/frames.gro"
done
printf 'Temperature\nPressure\nDensity\nPotential\n0\n' | "$GMX" energy \
    -f production.edr -o analysis/production_thermodynamics.xvg > analysis/production_energy.log 2>&1
printf 'System\n' | "$GMX" trjconv -f production.xtc -s production.tpr -n input/index.ndx \
    -dump 0 -pbc mol -o visualization.gro > analysis/visualization.log 2>&1
cat > view_binary.vmd <<'VMD'
# Open VMD's Tk Console and run: source /path/to/run/view_binary.vmd
set binary_dir [file dirname [file normalize [info script]]]
mol new [file join $binary_dir visualization.gro] type gro waitfor all
set binary_mol [molinfo top]
mol addfile [file join $binary_dir production.xtc] type xtc waitfor all molid $binary_mol
mol delrep 0 $binary_mol
mol representation VDW 0.65 12.0
mol selection {not water}
mol color ColorID 3
mol material Opaque
mol addrep $binary_mol
mol representation Lines 1.0
mol selection {water}
mol color ColorID 0
mol material Transparent
mol addrep $binary_mol
display projection Orthographic
axes location Off
VMD
date -u +%FT%TZ > analysis.complete
printf '\nBINARY RUN COMPLETE. Results: %s\n' "$MD_RUN_DIR"
printf 'View: visualization.gro + production.xtc (or source view_binary.vmd in VMD).\n'
printf 'Composition: analysis/production/; early stages: analysis/nvt/ and analysis/npt/.\n'
grep 'Performance:' production.log || true
printf 'This is a finite-time binary-mixture screening run, not an encapsulation-efficiency measurement.\n'
