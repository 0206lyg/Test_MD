#!/bin/bash
#SBATCH --job-name=binary_analysis
#SBATCH --account=ardekani
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=binary_analysis_%j.out

# Postprocess an existing binary run without starting or continuing MD.
set -Eeuo pipefail
export MD_ROOT="${MD_ROOT:-$HOME/MD}"
usage() {
    cat <<'EOF'
Usage: sbatch postprocess_binary.sbatch RUN_DIR [options]
  --stage production|npt|nvt   default: production
  --stride-ps 100             sample interval in ps
  --grid 4                    cells per box direction (default 1.5 nm at L=6 nm)
  --begin-ps TIME             optional absolute trajectory start time
  --end-ps TIME               optional absolute trajectory end time
Outputs: a new analysis_spatial_STAGE_TIMESTAMP_JOBID directory inside RUN_DIR.
No simulation is started. The original trajectory and topology are only read.
EOF
}
if [[ ${1:-} == --help || ${1:-} == -h ]]; then usage; exit 0; fi
[[ -n ${SLURM_JOB_ID:-} ]] || { usage >&2; printf 'Submit this analysis with sbatch.\n' >&2; exit 2; }
[[ $# -ge 1 && $1 != --* ]] || { usage >&2; exit 2; }
BIN_RUN=$(realpath -- "$1")
shift
BIN_STAGE=production
BIN_STRIDE=100
BIN_GRID=4
BIN_BEGIN=
BIN_END=
while (( $# )); do
    [[ $# -ge 2 ]] || { printf 'Missing value for %s\n' "$1" >&2; exit 2; }
    case $1 in
        --stage) BIN_STAGE=$2 ;;
        --stride-ps) BIN_STRIDE=$2 ;;
        --grid) BIN_GRID=$2 ;;
        --begin-ps) BIN_BEGIN=$2 ;;
        --end-ps) BIN_END=$2 ;;
        *) printf 'Unknown option: %s\n' "$1" >&2; exit 2 ;;
    esac
    shift 2
done
case $BIN_STAGE in production|npt|nvt) ;; *) printf 'Choose production, npt or nvt.\n' >&2; exit 2 ;; esac
[[ $BIN_GRID =~ ^[0-9]+$ ]] && (( 10#$BIN_GRID >= 2 && 10#$BIN_GRID <= 20 )) || { printf 'Grid must be 2..20.\n' >&2; exit 2; }
# AWK validates without depending on the preparation Python before it is loaded.
for BIN_VALUE in "$BIN_STRIDE" "${BIN_BEGIN:-0}" "${BIN_END:-0}"; do
    [[ $BIN_VALUE =~ ^[0-9]+([.][0-9]+)?$ ]] || { printf 'Times must be finite nonnegative decimal ps.\n' >&2; exit 2; }
done
awk -v x="$BIN_STRIDE" 'BEGIN {exit !(x>0)}' || { printf 'Stride must be positive.\n' >&2; exit 2; }
if [[ -n $BIN_BEGIN && -n $BIN_END ]]; then
    awk -v b="$BIN_BEGIN" -v e="$BIN_END" 'BEGIN {exit !(e>b)}' || { printf 'End time must exceed begin time.\n' >&2; exit 2; }
fi
for BIN_FILE in "$BIN_STAGE.tpr" "$BIN_STAGE.xtc" input/manifest.json input/index.ndx; do
    [[ -f "$BIN_RUN/$BIN_FILE" ]] || { printf 'Missing %s/%s\n' "$BIN_RUN" "$BIN_FILE" >&2; exit 1; }
done
[[ -f "$MD_ROOT/tools/analyze_binary.py" ]] || { printf 'Install the analysis overlay under %s first.\n' "$MD_ROOT" >&2; exit 1; }
BIN_OUT="$BIN_RUN/analysis_spatial_${BIN_STAGE}_$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_JOB_ID}"
mkdir -- "$BIN_OUT"
cp -- "$MD_ROOT/tools/analyze_binary.py" "$BIN_OUT/analyze_binary.py"
cp -- "$BIN_RUN/input/manifest.json" "$BIN_OUT/manifest.json"
printf 'Postprocessing results: %s\n' "$BIN_OUT"
trap 'printf "Binary postprocessing failed. Inspect %s and its extraction.log.\n" "$BIN_OUT" >&2' ERR
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
(
    source "$MD_ROOT/env.sh" md
    BIN_EXTRACT=("$GMX" trjconv -f "$BIN_RUN/$BIN_STAGE.xtc" -s "$BIN_RUN/$BIN_STAGE.tpr"
        -n "$BIN_RUN/input/index.ndx" -o "$BIN_OUT/frames.gro" -dt "$BIN_STRIDE" -pbc mol)
    [[ -z $BIN_BEGIN ]] || BIN_EXTRACT+=(-b "$BIN_BEGIN")
    [[ -z $BIN_END ]] || BIN_EXTRACT+=(-e "$BIN_END")
    printf '%q ' "${BIN_EXTRACT[@]}" > "$BIN_OUT/extraction_command.txt"
    printf '\n' >> "$BIN_OUT/extraction_command.txt"
    printf 'System\n' | "${BIN_EXTRACT[@]}" > "$BIN_OUT/extraction.log" 2>&1
)
(
    source "$MD_ROOT/env.sh" prep
    python3 "$BIN_OUT/analyze_binary.py" --trajectory "$BIN_OUT/frames.gro" \
        --manifest "$BIN_OUT/manifest.json" --output "$BIN_OUT" --grid "$BIN_GRID"
)
rm -- "$BIN_OUT/frames.gro"
date -u +%FT%TZ > "$BIN_OUT/analysis.complete"
printf '\nBINARY POSTPROCESSING COMPLETE: %s\n' "$BIN_OUT"
printf '2D: projection_last.svg, projection_time_average.svg\n'
printf 'Dispersion: concentration_timeseries.csv, dispersion_timeseries.svg\n'
printf '3D ParaView: concentration_3d.pvd and organic_molecules_3d.pvd (keep vtk/ alongside).\n'
