# Negishi MD: GROMACS 2026.4 CPU 연구 환경

이 ZIP을 `/home/lee6456`에서 풀면 `/home/lee6456/MD`가 생깁니다. 소스, 빌드, 설치된 소프트웨어, 분자 준비 도구와 입력 파일은 MD 아래에서 관리합니다. MD 계산 결과는 `/scratch/negishi/lee6456/md_JOBID_UTC시간/`에 자동 저장하며, Slurm 로그는 `sbatch`를 실행한 디렉터리에 남습니다.

## 1. 설치

다운로드한 `negishi_md_cpu.zip`을 Negishi의 `/home/lee6456`에 업로드하고 실행합니다.

```bash
cd /home/lee6456
unzip negishi_md_cpu.zip
mv negishi_md_cpu.zip MD/
cd /home/lee6456/MD
sbatch setup_cpu.sbatch
```

ZIP 안에 MD 폴더가 이미 있으므로, `MD` 안에서 압축을 풀어 `MD/MD`를 만들지 않습니다. 압축 해제 후 ZIP도 MD로 옮깁니다. 기존 `OpenLB`와 `ondemand`는 수정하지 않습니다. 이전 ZIP의 `negishi_md_cpu` 폴더는 이번 설치에 사용하지 않습니다.

설치 job은 `ardekani` 계정의 `cpu` 파티션에서 1 node, 1 task, 4 CPU, 16 GB, 최대 6시간을 요청합니다. 이 시간은 설치 job의 상한이며 설치 소요시간의 예측은 아닙니다. 로그인 노드에서는 위 `sbatch`만 제출하고, 실제 다운로드·빌드·검사는 계산 노드에서 실행합니다.

설치 스크립트가 순서대로 수행하는 일:

1. 기존 `amber/24` 모듈의 `antechamber`, `parmchk2`, `tleap` 확인.
2. CMake 3.31.6 전용 Conda 환경과 Packmol·RDKit·ParmEd 준비 환경 생성.
3. GROMACS 2026.4 소스와 동일 버전 regression tests 다운로드 및 공식 MD5 확인.
4. GCC 12.2.0으로 GROMACS 및 필요한 FFTW 빌드.
5. 공식 `check` target을 실행하고 통과하면 사용자 디렉터리에 설치.
6. 실제 버전, CMake 설정, 패키지 목록, 검사 로그 기록.

환경 및 패키지의 재현을 위해 GROMACS와 CMake 버전은 고정했습니다. 매번 자동으로 최신 버전을 선택하지 않습니다. CPU 엔진의 빌드 옵션은 `Release`, `AVX2_256`, mixed precision, thread-MPI 및 OpenMP이며 GPU 지원은 꺼져 있습니다. 단일 노드 내 CPU 병렬 계산을 지원합니다. 외부 OpenMPI는 AMBER 모듈을 불러올 때만 사용합니다.

네트워크 접근이 필요한 곳은 conda-forge, GROMACS 다운로드 서버와 FFTW 다운로드 서버입니다. 관리자 권한이나 `sudo`는 사용하지 않습니다. CMake 3.31.6은 GROMACS 2026.4의 요구사항인 CMake 3.28 이상을 충족합니다.

## 2. 설치 완료 확인

```bash
squeue -u lee6456
```

완료된 `setup_md_cpu_JOBID.out`의 마지막에 `MD CPU SETUP COMPLETE`가 있어야 합니다. 이어서:

```bash
cd /home/lee6456/MD
source ./env.sh md
gmx --version
```

출력에서 버전 `2026.4`, MPI library `thread_mpi`, OpenMP 지원, SIMD `AVX2_256`, GPU support 비활성 상태를 확인합니다. 실제 실행 파일은 `/home/lee6456/MD/software/gromacs-2026.4-cpu/bin/gmx`입니다.

`make check`에 해당하는 검사가 실패하면 설치 완료를 출력하지 않습니다. 오류는 Slurm 로그와 `software-info/check.log`에서 확인합니다. 설치가 중단됐을 때 같은 `sbatch setup_cpu.sbatch`를 재제출하면 기존 다운로드와 빌드 파일을 재사용합니다. Conda 환경도 재사용하며, 기존 환경의 자동 업그레이드는 하지 않습니다.

### 기존 설치의 검사 실패 후 재시도

최신 `setup_md_cpu_44122858.out`에서는 `acetonitrilRF`를 포함한 나머지 검사가 통과하고, `NbnxmTests`의 `SingleEnergyGroupSimd2xMM` 비교 하나가 실패했습니다. 이 공식 검사는 전역 `std::rand()`로 입력을 만들고, float 순차 합산을 기준으로 SIMD 합산과 비교합니다. 이번 차이는 5 ULP이며 검사 허용값은 4 ULP입니다. 별도 환경에서 공식 소스와 검사를 수정하지 않은 동일 실행 파일은 정상 환경에서 연속 3회 통과했고, 외부 초기화 코드가 난수 상태를 바꾸면 같은 검사 실패가 재현됐습니다. 진단용으로 같은 float 입력의 double 합을 계산했을 때 SIMD 결과가 더 정확해도 실패하는 사례가 확인됐습니다.

Negishi 로그에는 `module purge` 뒤에도 `xalt/3.0.2`가 남아 있습니다. XALT 3.0.2의 sampling 초기화 코드는 전역 `srand()`를 호출하므로 이 검사 입력을 바꿀 수 있습니다. 해당 job의 sampling 설정은 로그에 없어 현장 원인으로 확정하지는 않았습니다. 수정본은 모든 모듈 로드가 끝난 뒤 검사 subprocess에만 `XALT_SAMPLING=no`를 설정해 이 간섭을 방지하고 이전 sampling 설정을 기록합니다. XALT tracking은 유지하며 실제 MD 실행 환경은 바꾸지 않습니다. GROMACS 소스·기준값·허용오차·전체 검사 목록은 그대로이며 전체 공식 `check`가 통과해야 설치를 완료합니다. 수동으로 `.validated`를 만들지 않습니다.

`setup_md_cpu_44121885.out`에서는 앞서 수정한 OpenMP 환경변수 충돌이 해결돼 106개 검사 중 105개가 통과했습니다. 남은 `acetonitrilRF` 실패는 기준 trajectory와 이번 trajectory의 미세한 좌표 차이가 1.0 nm cutoff를 가로지르면서 발생한 힘 비교 차이입니다. 실패 자료에서 두 원자의 최종 거리는 기준값 0.999972808 nm, 이번 값 1.000001565 nm입니다. 해당 테스트의 RF/LJ 상호작용은 cutoff에서 힘이 불연속이며, 계산되는 힘 변화 약 10.09 kJ/mol/nm가 관찰된 차이 약 10.1과 일치합니다.

이전 스크립트는 빌드에 할당한 CPU 수를 regression test의 thread-MPI rank 수에도 그대로 사용했습니다. 첨부된 CMake 설정과 실행 로그에서는 9 ranks였습니다. 수정본은 컴파일 병렬도와 검사 병렬도를 분리하고, regression tests에 기본 2 ranks × 1 OpenMP thread를 사용합니다. upstream CI의 thread-MPI 기본 rank 수도 2입니다. 공식 reference, 허용오차, 테스트 목록은 그대로 유지하며 전체 `check`가 통과해야 설치를 완료합니다. 컴파일에는 계속 Slurm이 할당한 CPU를 사용하고, 실제 MD 실행은 `run_cpu.sbatch`의 설정을 따릅니다.

수정된 ZIP을 `/home/lee6456/MD/negishi_md_cpu.zip`에 업로드해 기존 ZIP을 교체한 다음, 변경된 두 파일만 추출합니다. 현재 설정과 설치 job의 CPU 수 등 사용자가 변경한 내용은 유지됩니다.

```bash
cd /home/lee6456
unzip -o MD/negishi_md_cpu.zip MD/tools/install_cpu.sh MD/README.md
cd /home/lee6456/MD
sbatch setup_cpu.sbatch
```

`build`, `envs`, `downloads`를 지우지 않습니다. 기존 컴파일 결과를 재사용하고, 검사 설정을 갱신한 뒤 전체 검사를 다시 실행합니다. 이전 실패 자료는 `software-info/check-새JOBID-before-retry.tar.gz`에 보존합니다. 다시 실패하면 새 `setup_md_cpu_새JOBID.out`과 `software-info/check-새JOBID-failed.tar.gz`를 함께 확인합니다. 진단 압축 파일에는 CTest 로그, CMake 설정, `acetonitrilRF`의 비교 결과·실행 로그 및 존재하는 TPR/TRR 파일이 포함됩니다.

## 3. 디렉터리 역할

| 경로 | 내용 |
|---|---|
| `/home/lee6456/MD/config.sh` | 버전·모듈·설치 위치·scratch 기본 위치 |
| `/home/lee6456/MD/env.sh` | 작업 환경 선택 |
| `/home/lee6456/MD/setup_cpu.sbatch` | 설치 job 제출 스크립트 |
| `/home/lee6456/MD/tools/install_cpu.sh` | 다운로드·빌드·검사 작업 |
| `/home/lee6456/MD/run_cpu.sbatch` | 준비된 TPR의 CPU 계산 job |
| `/home/lee6456/MD/software/gromacs-2026.4-cpu/` | 설치된 GROMACS와 자체 빌드 의존성 |
| `/home/lee6456/MD/envs/cmake-3.31.6/` | 빌드용 CMake 환경 |
| `/home/lee6456/MD/envs/prep/` | Packmol·RDKit·ParmEd 환경 |
| `/home/lee6456/MD/downloads/` | 다운로드한 GROMACS·test archive |
| `/home/lee6456/MD/src/` | 풀어 놓은 소스 및 regression tests |
| `/home/lee6456/MD/build/` | 컴파일 파일, FFTW 빌드 파일 |
| `/home/lee6456/MD/cache/` | Conda 패키지 cache |
| `/home/lee6456/MD/tmp/` | 설치 중 임시 파일 |
| `/home/lee6456/MD/software-info/` | 패키지 목록, 버전, 설정, 설치 검사 기록 |
| `/home/lee6456/MD/cases/` | 이후 작성할 분자 구조·force field·MDP·TPR |
| `/scratch/negishi/lee6456/md_JOBID_UTC시간/` | 각 계산의 trajectory, energy, checkpoint 등 |

기존 시스템 모듈은 공유 설치본을 그대로 참조합니다. Conda가 자체 환경 등록 정보를 홈에 남길 수 있지만, 새 소프트웨어와 package cache의 설치 위치는 위 표대로 MD 아래입니다. Conda 환경을 만든 뒤 MD 디렉터리를 다른 곳으로 옮기지 않습니다.

## 4. 작업 환경 선택

각 셸에서 필요한 환경만 선택합니다. `env.sh`는 기존 모듈 및 활성화된 Conda 환경을 정리한 뒤 지정한 환경을 불러옵니다. `.bashrc`는 수정하지 않습니다.

```bash
# MD 엔진
source /home/lee6456/MD/env.sh md
gmx --version

# 구조 생성, 분자 배치, topology 변환
source /home/lee6456/MD/env.sh prep
python -c 'from rdkit import Chem; import parmed; print("OK")'

# AmberTools
source /home/lee6456/MD/env.sh amber
command -v antechamber parmchk2 tleap

# CMake 확인이 필요한 경우
source /home/lee6456/MD/env.sh build
cmake --version
```

## 5. 분자 모델이 준비된 뒤의 계산

이 ZIP은 설치 환경과 실행 스크립트입니다. PLGA·약물·용매의 좌표, 파라미터, MDP, 실제 `model.tpr`는 아직 포함하지 않습니다. 소프트웨어 설치가 force field의 적합성을 검증하는 것은 아닙니다.

예를 들어 `cases/poc/model.tpr`를 만든 뒤:

```bash
cd /home/lee6456/MD
sbatch run_cpu.sbatch cases/poc/model.tpr
```

기본 요청은 1 node, 1 task, 8 CPU, 16 GB, 6시간입니다. thread-MPI 1 rank × OpenMP 8 threads로 시작하며, 모든 계산은 CPU에서 수행합니다. `-maxh 5.5`로 Slurm의 6시간 제한 전에 checkpoint를 기록하고 종료하도록 합니다. TPR에 지정된 계산 길이가 더 짧으면 먼저 종료됩니다.

job log에는 실제 결과 경로가 출력됩니다. 예를 들어:

`/scratch/negishi/lee6456/md_12345678_20261011T180000Z/`

입력 TPR도 `input.tpr`로 해당 폴더에 복사하고, 실제 실행 명령과 GROMACS 버전을 함께 기록합니다. `md.xtc`/`md.trr`, `md.edr`, `md.log`, `md.cpt` 등은 MDP 설정과 계산 진행 상태에 따라 이 폴더에 생성됩니다. `md_cpu_poc_12345678.out`은 제출한 `/home/lee6456/MD`에 남습니다.

짧은 성능 측정:

```bash
MD_MAXH=0.15 sbatch --time=00:15:00 run_cpu.sbatch cases/poc/model.tpr -nsteps 10000
```

최대 약 9분 또는 10,000 step 중 먼저 도달한 시점까지 실행합니다. `md.log`의 `Performance:`에 나오는 ns/day로 필요한 wall time을 추정합니다. 짧은 실행의 처리속도는 근사이며, 6시간 제한은 mechanism의 수렴 보장이 아닙니다.

CPU 병렬 배치를 비교해야 한다면:

```bash
MD_NTMPI=2 sbatch run_cpu.sbatch cases/poc/model.tpr
```

이는 같은 8 CPU에서 2 thread-MPI ranks × 4 OpenMP threads를 사용합니다. thread-MPI 빌드이므로 Slurm의 `--ntasks`는 항상 1로 유지합니다. scratch 상위 경로는 `config.sh` 또는 `MD_SCRATCH_ROOT` 환경변수로 지정할 수 있고, 기본 경로가 없거나 쓸 수 없으면 홈으로 우회하지 않고 오류를 출력합니다. 현재 스크립트는 새 계산 폴더를 생성하는 용도입니다.

## 검증 범위와 출처

Bash/Python 문법과 경로·실행 인자 전달을 로컬에서 검사했습니다. 이번 수정에서는 별도 Linux 환경에서 GROMACS 2026.4를 실제 빌드해, 공식 `acetonitrilRF` 검사와 첨부된 TPR을 각각 1·2·4·9 ranks로 실행했습니다. 9 ranks에서는 동일한 힘 비교 오류가 재현됐고, 1·2·4 ranks에서는 통과했습니다. 기본값인 2 ranks × 1 OpenMP thread로 공식 complex 42개와 essentialdynamics 7개도 모두 통과했습니다.

이 재현 환경은 Intel Xeon, GCC 13.3, AVX2_256, mixed precision, thread-MPI, FFTPACK입니다. 사용자 설치본의 GCC 12.2·FFTW 설정은 변경하지 않았습니다. `acetonitrilRF`는 RF 방식을 사용하므로 해당 힘 계산에 FFT를 사용하지 않습니다. Negishi에서 최종 검증은 재제출한 설치 job의 전체 공식 `check`가 수행합니다. PLGA·약물의 본 분자계를 계산한 결과는 아닙니다.

- GROMACS 2026.4와 공식 source/test MD5: https://manual.gromacs.org/2026.4/download.html
- GCC, CMake, FFTW, SIMD, thread-MPI 및 검사: https://manual.gromacs.org/2026.4/install-guide/index.html
- Compiler 관련 알려진 문제: https://manual.gromacs.org/2026.4/user-guide/known-issues.html
- XALT 3.0.2 sampling 초기화: https://github.com/xalt/xalt/blob/xalt-3.0.2/src/libxalt/xalt_initialize.c#L725-L742
- Negishi AMBER 모듈: https://docs.rcac.purdue.edu/software/apps_md/amber/
- Packmol: https://m3g.github.io/packmol/download.shtml
- RDKit: https://rdkit.org/docs/Install.html
- ParmEd: https://parmed.github.io/ParmEd/html/readwrite.html
