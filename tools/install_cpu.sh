#!/bin/bash
# Invoked by setup_cpu.sbatch on a CPU compute node.
set -Eeuo pipefail
trap 'printf "Installation failed at line %s. See the setup job log.\n" "$LINENO" >&2' ERR
export MD_ROOT="${MD_ROOT:-$HOME/MD}"
source "$MD_ROOT/config.sh"
[[ -n "${SLURM_JOB_ID:-}" ]] || { printf 'Submit with sbatch ~/MD/setup_cpu.sbatch.\n' >&2; exit 1; }
MD_BUILD_JOBS="${SLURM_CPUS_PER_TASK:-4}"
# Keep regression decomposition independent of compilation parallelism.
# Upstream thread-MPI CI also defaults to two ranks. The full check still runs.
MD_TEST_NTMPI="${MD_TEST_NTMPI:-2}"
if [[ ! "$MD_TEST_NTMPI" =~ ^[1-9][0-9]*$ ]] || (( MD_TEST_NTMPI > MD_BUILD_JOBS )); then
    printf 'MD_TEST_NTMPI must be a positive integer no larger than SLURM_CPUS_PER_TASK.\n' >&2
    exit 1
fi
printf 'Build jobs: %s; regression tests: %s thread-MPI ranks x 1 OpenMP thread\n' \
    "$MD_BUILD_JOBS" "$MD_TEST_NTMPI"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONUNBUFFERED=1 PYTHONNOUSERSITE=1
mkdir -p "$MD_ROOT"/{software,src,build,downloads,envs,cache/conda-pkgs,tmp,software-info,cases}
export TMPDIR="$MD_ROOT/tmp"
export CONDA_PKGS_DIRS="$MD_ROOT/cache/conda-pkgs"
export CONDA_ENVS_PATH="$MD_ROOT/envs"

# Preserve the failed force comparison before a retry can overwrite it.
save_check_diagnostics() {
    local phase="$1" relative
    local name="check-${SLURM_JOB_ID}-${phase}"
    local destination="$MD_ROOT/software-info/$name"
    mkdir -p "$destination"
    for relative in \
        CMakeCache.txt Testing/Temporary/LastTest.log \
        Testing/Temporary/LastTestsFailed.log tests/CTestTestfile.cmake \
        tests/complex/acetonitrilRF/{checkforce.out,checkforce.err,checkpot.out,mdrun.out,md.log,grompp.out,grompp.err,topol.tpr,traj.trr}; do
        [[ -f "$MD_BUILD_DIR/$relative" ]] || continue
        mkdir -p "$destination/$(dirname "$relative")"
        cp -- "$MD_BUILD_DIR/$relative" "$destination/$relative"
    done
    if [[ -f "$MD_ROOT/software-info/check.log" ]]; then
        cp -- "$MD_ROOT/software-info/check.log" "$destination/check.log"
    fi
    tar -czf "$destination.tar.gz" -C "$MD_ROOT/software-info" "$name"
    printf 'Saved check diagnostics: %s.tar.gz\n' "$destination"
}
if [[ -f "$MD_BUILD_DIR/Testing/Temporary/LastTestsFailed.log" ]]; then
    save_check_diagnostics before-retry
fi

# Check existing AmberTools without coupling its MPI module to the new engine.
(
    source "$MD_ROOT/env.sh" amber
    module list 2>&1
    for tool in antechamber parmchk2 tleap; do
        command -v "$tool" || { printf 'Missing AMBER utility: %s\n' "$tool" >&2; exit 1; }
    done
) 2>&1 | tee "$MD_ROOT/software-info/amber.txt"

# CMake is a build-only Conda environment; preparation tools are separate.
(
    if declare -F conda >/dev/null; then
        while (( ${CONDA_SHLVL:-0} > 0 )); do conda deactivate; done
    fi
    module purge
    module load "$MD_CONDA_MODULE"
    if [[ ! -f "$MD_CMAKE_PREFIX/conda-meta/history" ]]; then
        conda create --yes --prefix "$MD_CMAKE_PREFIX" \
            --override-channels --channel conda-forge --strict-channel-priority \
            --solver libmamba "cmake=$MD_CMAKE_VERSION" python=3.11
    fi
    if [[ ! -f "$MD_PREP_PREFIX/conda-meta/history" ]]; then
        conda create --yes --prefix "$MD_PREP_PREFIX" \
            --override-channels --channel conda-forge --strict-channel-priority \
            --solver libmamba python=3.11 packmol rdkit parmed
    fi
    conda list --prefix "$MD_CMAKE_PREFIX" --explicit > "$MD_ROOT/software-info/cmake-explicit.txt"
    conda list --prefix "$MD_PREP_PREFIX" --explicit > "$MD_ROOT/software-info/prep-explicit.txt"
    conda env export --prefix "$MD_PREP_PREFIX" --from-history > "$MD_ROOT/software-info/prep-history.yml"
    source "$MD_ROOT/env.sh" prep
    command -v packmol
    python - <<'PY'
from importlib.metadata import version
from rdkit import Chem
import parmed
assert Chem.MolFromSmiles('CS(=O)C') is not None
for package in ('rdkit', 'ParmEd'):
    print(f'{package}: {version(package)}')
print('Preparation imports OK')
PY
)

# Configure with the native GCC compiler and without an activated Conda env.
source "$MD_ROOT/env.sh" build
for tool in gcc g++ make perl curl tar md5sum sha256sum; do command -v "$tool" >/dev/null; done
MD_CMAKE="$MD_CMAKE_PREFIX/bin/cmake"
"$MD_CMAKE" --version | tee "$MD_ROOT/software-info/cmake.txt"
[[ "$("$MD_CMAKE" --version)" == "cmake version $MD_CMAKE_VERSION"* ]] || { printf 'Unexpected CMake version.\n' >&2; exit 1; }
gcc --version > "$MD_ROOT/software-info/compiler.txt"
module list > "$MD_ROOT/software-info/build-modules.txt" 2>&1

download_checked() {
    local url="$1" filename="$2" expected_md5="$3"
    local target="$MD_ROOT/downloads/$filename"
    if [[ -f "$target" ]] && printf '%s  %s\n' "$expected_md5" "$target" | md5sum --check --status; then
        printf 'Using verified download: %s\n' "$filename"
        return
    fi
    curl --fail --location --retry 3 --connect-timeout 30 --output "$target.part" "$url"
    printf '%s  %s\n' "$expected_md5" "$target.part" | md5sum --check
    mv -- "$target.part" "$target"
}
download_checked "https://ftp.gromacs.org/gromacs/gromacs-2026.4.tar.gz" \
    gromacs-2026.4.tar.gz fbe7d5071dee17e3f448502b1fd30d55
download_checked "https://ftp.gromacs.org/regressiontests/regressiontests-2026.4.tar.gz" \
    regressiontests-2026.4.tar.gz 7098c9e520c3ba3654a7bb5cc4e8adf5
sha256sum "$MD_ROOT/downloads/gromacs-2026.4.tar.gz" \
          "$MD_ROOT/downloads/regressiontests-2026.4.tar.gz" > "$MD_ROOT/software-info/downloads.sha256"

if [[ ! -f "$MD_SOURCE_DIR/.unpacked" ]]; then
    tar -xzf "$MD_ROOT/downloads/gromacs-2026.4.tar.gz" -C "$MD_ROOT/src"
    touch "$MD_SOURCE_DIR/.unpacked"
fi
if [[ ! -f "$MD_TEST_DIR/.unpacked" ]]; then
    tar -xzf "$MD_ROOT/downloads/regressiontests-2026.4.tar.gz" -C "$MD_ROOT/src"
    touch "$MD_TEST_DIR/.unpacked"
fi

# Zen3 CPU, mixed precision, thread-MPI + OpenMP, one node. No GPU dependency.
# Unix Makefiles is required for this bundled FFTW build workflow.
"$MD_CMAKE" -S "$MD_SOURCE_DIR" -B "$MD_BUILD_DIR" -G "Unix Makefiles" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="$(command -v gcc)" \
    -DCMAKE_CXX_COMPILER="$(command -v g++)" \
    -DPython3_EXECUTABLE="$MD_CMAKE_PREFIX/bin/python" \
    -DCMAKE_INSTALL_PREFIX="$MD_GMX_PREFIX" \
    -DGMX_GPU=OFF -DGMX_MPI=OFF -DGMX_THREAD_MPI=ON -DGMX_OPENMP=ON \
    -DGMX_SIMD=AVX2_256 -DGMX_DOUBLE=OFF \
    -DGMX_BUILD_OWN_FFTW=ON -DGMX_FFT_LIBRARY=fftw3 \
    -DREGRESSIONTEST_DOWNLOAD=OFF -DREGRESSIONTEST_PATH="$MD_TEST_DIR" \
    "-DREGRESSIONTEST_EXTRA_ARGS=-nt;$MD_TEST_NTMPI;-ntomp;1" \
    -DBUILD_TESTING=ON \
    2>&1 | tee "$MD_ROOT/software-info/configure.log"
"$MD_CMAKE" --build "$MD_BUILD_DIR" --parallel "$MD_BUILD_JOBS" \
    2>&1 | tee "$MD_ROOT/software-info/build.log"

# Do not report a successful installation unless the official checks pass.
# Tests choose their own -ntomp (including 2); inherited OMP_NUM_THREADS=1
# conflicts with those settings. Only the test subprocess needs it unset.
if (
    unset OMP_NUM_THREADS
    # XALT 3.0.2 sampling initializes the application's global C RNG with
    # srand(time/pid/clock). NbnxmTests uses std::rand without resetting it.
    # Keep XALT tracking, but prevent sampling from changing test inputs.
    # Scope this setting to checks, after all module loads have completed.
    printf 'XALT sampling before checks: primary=%s; scalar=%s; scalar/spsr=%s\n' \
        "${XALT_SAMPLING-<unset>}" "${XALT_SCALAR_SAMPLING-<unset>}" \
        "${XALT_SCALAR_AND_SPSR_SAMPLING-<unset>}"
    export XALT_SAMPLING=no
    export CTEST_PARALLEL_LEVEL=1 CTEST_OUTPUT_ON_FAILURE=1
    "$MD_CMAKE" --build "$MD_BUILD_DIR" --target check --parallel "$MD_BUILD_JOBS"
) 2>&1 | tee "$MD_ROOT/software-info/check.log"; then
    printf 'Official GROMACS checks passed.\n'
else
    save_check_diagnostics failed
    for diagnostic in "$MD_BUILD_DIR/tests/complex/acetonitrilRF/"checkforce.{out,err}; do
        if [[ -s "$diagnostic" ]]; then
            printf '\nForce-comparison diagnostic: %s\n' "$diagnostic"
            sed -n '1,160p' "$diagnostic"
        fi
    done
    printf '\nOfficial checks failed; installation was not marked complete.\n' >&2
    exit 1
fi
# A failed or partial reinstall must not leave an old validation marker.
if [[ -f "$MD_GMX_PREFIX/.validated" ]]; then
    mv -- "$MD_GMX_PREFIX/.validated" "$MD_GMX_PREFIX/.validated.previous"
fi
"$MD_CMAKE" --install "$MD_BUILD_DIR" \
    2>&1 | tee "$MD_ROOT/software-info/install.log"
"$MD_GMX_PREFIX/bin/gmx" --version 2>&1 | tee "$MD_ROOT/software-info/gromacs.txt"
cp "$MD_BUILD_DIR/CMakeCache.txt" "$MD_ROOT/software-info/CMakeCache.txt"
date -u +%Y-%m-%dT%H:%M:%SZ > "$MD_GMX_PREFIX/.validated"

printf '\nMD CPU SETUP COMPLETE\n'
printf 'GROMACS: source "%s/env.sh" md\n' "$MD_ROOT"
printf 'Preparation: source "%s/env.sh" prep\n' "$MD_ROOT"
printf 'AmberTools: source "%s/env.sh" amber\n' "$MD_ROOT"
printf 'Results will be created under: %s\n' "$MD_SCRATCH_ROOT"
