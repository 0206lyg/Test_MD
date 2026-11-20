#!/bin/bash
#SBATCH --job-name=setup_md_cpu
#SBATCH --account=ardekani
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=setup_md_cpu_%j.out

# Common installer for GROMACS, molecular preparation and LAMMPS DPD.
set -Eeuo pipefail
export MD_ROOT="${MD_ROOT:-$HOME/MD}"
[[ -f "$MD_ROOT/tools/install_cpu.sh" ]] || { printf 'Extract the ZIP into your home directory to create ~/MD.\n' >&2; exit 1; }
exec bash "$MD_ROOT/tools/install_cpu.sh"
