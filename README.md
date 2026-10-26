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
2. CMake 3.31.6 전용 Conda 환경과 Packmol·RDKit·ParmEd·SciPy 준비 환경 생성.
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
| `/home/lee6456/MD/prepare_poc.sbatch` | 조성 JSON에서 분자 좌표·force field·GROMACS 입력 생성 |
| `/home/lee6456/MD/run_poc.sbatch` | 최소화·초기 MD·성능 측정 및 checkpoint 재시작 |
| `/home/lee6456/MD/tools/prepare_case.py` | 입력 생성의 세 단계 실행 |
| `/home/lee6456/MD/software/gromacs-2026.4-cpu/` | 설치된 GROMACS와 자체 빌드 의존성 |
| `/home/lee6456/MD/envs/cmake-3.31.6/` | 빌드용 CMake 환경 |
| `/home/lee6456/MD/envs/prep/` | Packmol·RDKit·ParmEd 환경 |
| `/home/lee6456/MD/downloads/` | 다운로드한 GROMACS·test archive |
| `/home/lee6456/MD/src/` | 풀어 놓은 소스 및 regression tests |
| `/home/lee6456/MD/build/` | 컴파일 파일, FFTW 빌드 파일 |
| `/home/lee6456/MD/cache/` | Conda 패키지, 재사용하는 PLGA 좌표·전하 계산 cache |
| `/home/lee6456/MD/tmp/` | 설치 중 임시 파일 |
| `/home/lee6456/MD/software-info/` | 패키지 목록, 버전, 설정, 설치 검사 기록 |
| `/home/lee6456/MD/cases/` | 조성 JSON 및 생성한 분자 구조·force field·MDP |
| `/home/lee6456/MD/validation/` | 로컬에서 실제 수행한 입력 생성 검증의 요약 |
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

## 5. PLGA / 약물 / 혼합용매 입력 생성

기존 설치가 완료됐다면 GROMACS를 재설치하지 않습니다. 새 ZIP을 `/home/lee6456/MD/negishi_md_cpu.zip`에 업로드한 뒤 입력 생성 관련 파일만 추가·갱신합니다. 아래 명령은 기존 `config.sh`, `env.sh`와 설치 파일을 덮어쓰지 않습니다.

```bash
cd /home/lee6456
unzip -o MD/negishi_md_cpu.zip 'MD/tools/mdprep_*.py' MD/tools/prepare_case.py MD/prepare_poc.sbatch MD/run_poc.sbatch 'MD/cases/*' 'MD/validation/*' MD/README.md
cd /home/lee6456/MD
```

ZIP에는 기본 DMSO 계의 **이미 생성한 `cases/poc_lawsone_dmso_input/`**도 포함돼 있습니다. 첫 CPU 실행은 준비 작업을 반복하지 않고 바로 제출할 수 있습니다.

```bash
sbatch run_poc.sbatch cases/poc_lawsone_dmso_input
```

아래의 `prepare_poc.sbatch`는 입력을 재생성하거나 다른 조건을 만드는 방법입니다. 용매별 비교 계산을 시작할 때는 기준 DMSO 계를 포함한 모든 조건을 Negishi의 같은 준비 환경에서 생성합니다.

기본 설정은 `cases/poc_lawsone_dmso.json`입니다. **PLGA 30사슬 × 90 반복단위**, LA:GA 50:50, 무작위 D/L 배열, hydroxyl/carboxyl 말단, lawsone 1분자, DMSO/물, 전체 박스 8 × 8 × 14 nm를 사용합니다. 기존 6 × 450과 총 반복단위 수 2,700개를 유지하면서, PoC에서 사슬당 약 5.87 kDa로 줄인 모델입니다. 말단 수와 사슬 이완 시간은 실험 분자량의 사슬과 달라지며, 분자량 의존성을 제거했다고 가정하지 않습니다. 실험 재료 Resomer RG 503 H는 50:50 D,L-PLGA 산 말단 고분자입니다. 반복단위의 실제 LA/GA 배열과 입체화학은 알려져 있지 않아 고정 seed의 무작위 배열로 모델링합니다.

초기 고분자 배치 영역은 **8 × 8 × 8 nm**입니다. 처음부터 4 nm 두께의 조밀한 고분자를 강제로 만들지 않고, 중앙의 느슨한 고분자 영역과 양쪽 약 3 nm의 외부 용매 영역에서 시작합니다. `slab_nm`은 초기 배치 영역이며 MD에서 고분자를 가두는 벽이나 위치 restraint가 아닙니다. 주변 용매는 고분자 사이 빈 공간에도 배치됩니다. 평형 고분자 두께·밀도를 이 초기 설정으로 가정하지 않으며, 이완 이후 고분자 내부·외부의 조성 영역이 실제로 유지되는지 확인해야 합니다.

첫 조성의 `water_mole_fraction = 0.9`는 물+유기용매 분자 수 기준 90 mol% 물이며, 실험의 최종 부피비를 그대로 재현한 값이 아닙니다. 기작 탐색을 위한 시작 조성입니다. 분자 수는 순수액체 밀도를 이용한 부피 가산 근사로 정하고 manifest에 기록합니다. 동일 몰분율에서 유기용매 종류를 바꾸면 분자 부피가 달라지므로 물·유기용매 분자 수도 재계산합니다. 초기 밀도는 NPT에서 이완합니다.

기본 DMSO 계는 PLGA 20,340원자, lawsone 1분자, DMSO 1,628분자, 물 14,655분자로 **총 80,604원자**입니다.

입력 생성은 계산 노드에서 수행합니다. 기존 RDKit, Packmol, ParmEd와 AmberTools를 사용합니다. 사슬 간 초기 겹침을 해소하는 좌표 최적화에 SciPy가 필요하며, 기존 준비 환경에 없으면 `prepare_poc.sbatch`가 계산 노드에서 이 의존성만 추가합니다. 기존 패키지 버전을 동결한 상태로 설치하며 GROMACS를 재빌드하지 않습니다.

```bash
cd /home/lee6456/MD
sbatch prepare_poc.sbatch cases/poc_lawsone_dmso.json
```

완료되면 다음 파일들이 `cases/poc_lawsone_dmso_input/`에 생깁니다.

| 파일 | 내용 |
|---|---|
| `case.json` | 사용한 설정의 사본 |
| `polymer.json` | 사슬별 LA/GA 배열, D/L 입체화학, 원자 대응표, 생성 기록 |
| `molecules/` | 분자별 좌표, 결합, 전하, Amber topology 및 계산 로그 |
| `charge_fragments/` | 짧은 올리고머 전하 계산과 긴 사슬로의 전하 이전 근거 |
| `manifest.json` | 원자·분자 수, 실제 조성, 설정 식별값 |
| `packmol.inp`, `packmol.log`, `packed.pdb` | 혼합계 배치 입력·로그·좌표 |
| `system.gro`, `topol.top` | GROMACS 초기 좌표와 topology |
| `mdp/` | 최소화·초기 이완·짧은 성능 측정 설정 |
| `READY.json` | 입력 조립 완료 기록; 평형화 완료를 의미하지 않음 |

좌표와 파라미터를 직접 손으로 입력하지 않습니다. PLGA는 결합 길이·각도·입체화학을 가진 반복단위를 연결해 사슬을 각각 생성합니다. 이어서 결합·각도에 해당하는 거리와 입체화학을 유지하는 기하학적 최적화로 사슬 간 겹침을 해소합니다. 이것은 초기 좌표 준비이며 GAFF2에서 평형화된 구조가 아닙니다. 비결합 heavy atom 1.4 Å, 모든 원자 1.2 Å의 초기 접근거리 검사를 통과한 뒤 실제 GAFF2 에너지 최소화가 수행돼야 합니다. MD의 force field나 최소화 기준을 바꾸는 soft potential은 사용하지 않습니다.

긴 사슬 전체에 ETKDG나 양자화학 계산을 적용하지 않습니다. 약물과 유기용매는 SMILES에서 좌표를 생성하고, Packmol은 준비된 사슬의 내부 결합을 변경하지 않은 채 나머지 분자를 배치합니다. **Packmol에서 고정하는 것은 초기 배치 중의 좌표뿐이며, MD에서 PLGA는 움직입니다.**

### 5.1 Force field와 화학 상태

약물·유기용매·PLGA에는 GAFF2를 사용하고 물에는 TIP3P를 사용합니다. 모두 원자 수준의 모델입니다. 작은 분자는 AmberTools AM1-BCC로 전하를 계산합니다. PLGA는 5개 반복단위로 된 짧은 참조 사슬 10종에서 전하를 계산합니다. 8종은 가운데 반복단위와 양옆 LA/GA 조합을 모두 포함하고, 추가 2종은 LA 말단 조합을 보완합니다. 내부 전하는 각 참조 사슬의 중앙 반복단위에서, 말단 전하는 해당 말단에서 가져옵니다. 인접 LA/GA 종류·말단·원자 역할을 대응시켜 긴 사슬로 이전하며, D/L 중심에는 같은 전하를 사용합니다. 동등한 수소 전하는 평균화하고 각 중성 반복단위의 전하 합을 0으로 맞춘 보정량을 별도 기록합니다. 이 방식은 국소 환경의 전하 이전 근사를 사용하는 **PoC 파라미터화**이며, 실험 loading 또는 분배 자유에너지를 재현한다고 검증된 PLGA force field는 아닙니다.

고분자의 추가 GAFF2 결합 항은 짧은 참조 사슬에서 `parmchk2`로 추출하고, 같은 항의 계수가 서로 다른 경우 중단합니다. 긴 사슬마다 같은 국소 항을 반복 검색하지 않으며, 전체 사슬의 모든 필요한 결합·각도·비틀림 항은 `tleap`에서 실제로 확인합니다. 각 조건은 동일한 준비 도구·force-field 버전으로 생성해 비교합니다.

현재 모든 약물과 PLGA 말단은 **중성 상태**로 설정합니다. 논문에 명시되지 않은 실제 pH·이온화 비율을 모델이 자동으로 결정하지 않습니다. 특히 lawsone과 산 말단의 전하 상태는 본 기작 비교에서 별도로 검토해야 합니다. poloxamer, 염, 용매 제거·세척 과정은 이 국소 분배 모델에 포함하지 않습니다.

분자 생성·양자화학 출력·전하 이전·중성화 기록은 보존합니다. 누락된 force-field 항, 원자 수 불일치, Packmol 미수렴, 심한 분자 간 겹침 또는 GROMACS 입력 오류가 있으면 중단하며, 테스트 허용오차를 완화하거나 `grompp -maxwarn`으로 통과시키지 않습니다.

### 5.2 첫 CPU 계산

입력 준비가 완료된 뒤 제출합니다.

```bash
cd /home/lee6456/MD
sbatch run_poc.sbatch cases/poc_lawsone_dmso_input
```

기본 요청은 8 CPU, 16 GB, 최대 6시간입니다. 하나의 job에서 에너지 최소화 → NVT 100 ps → NPT 200 ps → 성능 측정용 MD 100 ps를 순서대로 실행합니다. NVT는 1 fs, 이후 MD는 수소 결합 길이를 구속한 2 fs를 사용합니다. 고분자·약물의 위치는 구속하지 않습니다. NPT는 slab의 x/y 면적을 유지하고 z 방향 박스 길이를 이완합니다. 이 시간들은 초기 실행 범위이며, 평형화·분배 자유에너지 수렴을 보장하는 값이 아닙니다.

전체 단계가 **하나의 6시간 예산**을 공유하고 종료 전에 checkpoint를 저장합니다. 미완료 단계를 완료로 표시하지 않습니다. 결과는 `/scratch/negishi/lee6456/md_JOBID_UTC시간/`에, Slurm 로그는 제출한 MD 폴더에 남습니다. 이어서 실행할 때는 로그에 출력된 실제 경로를 사용합니다.

```bash
sbatch run_poc.sbatch --restart /scratch/negishi/lee6456/md_JOBID_UTC시간
```

재시작은 동일한 결과 폴더에서 완료된 단계를 건너뛰고 미완료 단계의 checkpoint를 이어갑니다. 정상 완료된 마지막 `benchmark.log`의 `Performance:`가 해당 계의 실제 ns/day입니다.

### 5.3 용매·약물을 바꾸는 방법

DMSO 기준 설정을 복사해 새 이름으로 저장합니다.

```bash
cp cases/poc_lawsone_dmso.json cases/poc_lawsone_acn.json
```

복사한 JSON의 `mixture.solvent`를 `"ACN"`으로 바꾼 뒤 제출합니다.

```bash
sbatch prepare_poc.sbatch cases/poc_lawsone_acn.json
```

지원하는 용매 값은 `DMSO`, `DMF`, `acetone`, `ACN`이며 약물 `drug.name`은 `lawsone`, `menadione`, `plumbagin`입니다. 새로운 조건은 새로운 출력 폴더를 사용해 기존 입력을 덮어쓰지 않습니다.

**용매만 바꾸고 polymer 설정·box·seed를 유지하면 같은 초기 PLGA 좌표를 재사용합니다.** 반복단위 배열·사슬 입체화학·고분자 topology도 유지되며, 유기용매의 구조·파라미터·분자 수·배치와 혼합계 topology를 새로 만듭니다. PLGA 좌표와 전하 계산은 `MD/cache/`에 저장해 재사용합니다. cache가 없거나 RDKit 버전이 달라도, `cases/*_input/`에 같은 고분자 설정의 완성된 입력이 있으면 그 PLGA PDB 좌표를 그대로 가져옵니다. 새 cache는 좌표를 double precision으로 보존합니다. 용매마다 달라져야 하는 평형 PLGA 구조는 각각의 MD 이완에서 얻습니다. 동일한 초기 구조를 쓰는 것이 평형 구조까지 같게 고정하는 것은 아닙니다.

서로 다른 초기 구조를 비교할 때는 `polymer.seed`를 바꿉니다. 용매 배치만 바꾸려면 `mixture.packing_seed`를 바꿉니다. 약물을 처음부터 고분자 영역에 놓으려면 `drug.initial_region`을 `"polymer"`로 바꿀 수 있습니다. 시작 위치가 다른 두 trajectory의 짧은 체류 시간만으로 평형 분배 자유에너지를 결정하지 않습니다.

## 6. 이미 준비된 TPR의 단독 실행

기존 `run_cpu.sbatch`도 유지됩니다. 이는 준비된 TPR 하나를 실행하는 스크립트이며, 위 여러 단계의 준비·재시작 처리는 `run_poc.sbatch`가 담당합니다.

```bash
sbatch run_cpu.sbatch cases/my_case/model.tpr
```

## 검증 범위와 출처

입력 생성 도구는 별도 Linux 환경의 실제 AmberTools 24.8, RDKit, Packmol, ParmEd와 GROMACS 2026.4로 확인했습니다.

- 지원하는 약물 3종과 유기용매 4종 모두 AM1-BCC/GAFF2 파라미터화를 완료했습니다. 원자 종류·순서와 모든 공유결합의 원자 쌍이 원래 분자 그래프와 일치하고, 누락된 force-field 항이 없으며 중성 topology를 생성했습니다.
- PLGA 12-mer + lawsone + DMSO + TIP3P의 8,878원자 혼합계를 실제 조립했습니다. 최소화는 459 step에서 Fmax = 879.966 kJ/mol/nm로 지정 기준 1000 아래에 도달했습니다. NVT 0.1 ps, NPT 0.2 ps, 후속 MD 0.2 ps가 완료됐고, 기록된 에너지는 유한하며 LINCS 경고가 없었습니다. 이 짧은 실행은 파일 연결과 초기 안정성 검사입니다.
- CPU job의 중간 종료·checkpoint 재시작은 작은 물 계에서 실제 GROMACS로 확인했습니다. 미완료 단계에서 다음 단계로 넘어가지 않고, 같은 결과 폴더에서 이어집니다.
- 기본 30 × DP90 고분자 구조의 모든 1,350개 입체중심과 비결합 접근거리를 확인했습니다. DMSO → ACN 구조 생성에서는 기존 binary cache가 없는 경우에도 완성된 DMSO 입력을 참조해 30개 PLGA PDB 파일을 바이트 단위까지 동일하게 재사용했습니다. ACN 계의 원자 수는 78,877개로 재계산됐습니다.
- **기본 80,604원자 전체 계**는 실제 GROMACS 최소화가 1,954 step에서 Fmax = 989.41772 kJ/mol/nm로 수렴했습니다. 이어진 NVT 100 step × 1 fs가 완료됐고 에너지가 유한하며 LINCS 경고가 없었습니다. 최소화 및 NVT 이후에도 1,350개 입체중심의 부호가 모두 유지됐습니다. 검사에서만 NVT를 0.1 ps로 줄였으며, 배포 MDP는 100 ps입니다. 이는 실행 가능한 초기 계의 확인이며 온도·밀도 평형화 또는 분배 기작의 검증이 아닙니다.

SQM이 출력하는 개별 Mulliken 전하는 소수 셋째 자리여서, 일부 중성 분자의 AM1-BCC 출력 전하 합에 약 0.001–0.002 e 잔차가 있었습니다. 코드에서는 SQM의 총전하가 중성임을 확인한 뒤 출력 정밀도 범위 안의 잔차만 보정하고, 원본·보정량을 함께 보존합니다. PLGA의 반복단위 전하 이전에 사용하는 중성화 보정은 별도로 기록합니다. QM 수렴 기준, GROMACS 허용오차 또는 경고 처리 기준은 완화하지 않습니다.

Bash/Python 문법과 경로·실행 인자 전달도 검사했습니다. 이전 설치 문제 진단에서는 별도 Linux 환경에서 GROMACS 2026.4를 실제 빌드해, 공식 `acetonitrilRF` 검사와 첨부된 TPR을 각각 1·2·4·9 ranks로 실행했습니다. 9 ranks에서는 동일한 힘 비교 오류가 재현됐고, 1·2·4 ranks에서는 통과했습니다. 기본값인 2 ranks × 1 OpenMP thread로 공식 complex 42개와 essentialdynamics 7개도 모두 통과했습니다.

이 재현 환경은 Intel Xeon, GCC 13.3, AVX2_256, mixed precision, thread-MPI, FFTPACK입니다. 사용자 설치본의 GCC 12.2·FFTW 설정은 변경하지 않았습니다. `acetonitrilRF`는 RF 방식을 사용하므로 해당 힘 계산에 FFT를 사용하지 않습니다. Negishi에서 설치본의 최종 검증은 설치 job의 전체 공식 `check`가 수행합니다. 실제 PLGA 분배 자유에너지나 실험 loading을 검증한 결과는 아닙니다.

- GROMACS 2026.4와 공식 source/test MD5: https://manual.gromacs.org/2026.4/download.html
- GCC, CMake, FFTW, SIMD, thread-MPI 및 검사: https://manual.gromacs.org/2026.4/install-guide/index.html
- Compiler 관련 알려진 문제: https://manual.gromacs.org/2026.4/user-guide/known-issues.html
- XALT 3.0.2 sampling 초기화: https://github.com/xalt/xalt/blob/xalt-3.0.2/src/libxalt/xalt_initialize.c#L725-L742
- Negishi AMBER 모듈: https://docs.rcac.purdue.edu/software/apps_md/amber/
- Packmol: https://m3g.github.io/packmol/download.shtml
- RDKit: https://rdkit.org/docs/Install.html
- ParmEd: https://parmed.github.io/ParmEd/html/readwrite.html
- Resomer RG 503 H 재료 사양: https://www.evonik.com/en/products/hc/pr_52021485.html
- 논문: https://doi.org/10.1002/mame.202500269
- 고분자 전하 이전 방법론: https://doi.org/10.1021/acs.jcim.3c01691 (이 코드의 GAFF2/Amber AM1-BCC 파라미터 자체를 제공하는 논문은 아님)
- SQM 개별 Mulliken 전하 출력 정밀도: https://github.com/Amber-MD/AmberClassic/blob/main/src/sqm/qm2_print_charges.F90
