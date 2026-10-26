# Negishi MD: GROMACS 2026.4 CPU Research Environment

Extracting the distribution ZIP in `/home/lee6456` creates `/home/lee6456/MD`. Source code, build files, installed software, molecular preparation tools, and simulation inputs are managed under `MD`. MD results are automatically written to `/scratch/negishi/lee6456/md_JOBID_UTC_TIMESTAMP/`, while Slurm logs remain in the directory from which `sbatch` was submitted.

## 1. Installation

Upload `negishi_md_cpu.zip` to `/home/lee6456` on Negishi, then run:

```bash
cd /home/lee6456
unzip negishi_md_cpu.zip
mv negishi_md_cpu.zip MD/
cd /home/lee6456/MD
sbatch setup_cpu.sbatch
```

The ZIP already contains the `MD` directory. Extract it from `/home/lee6456` to avoid creating a nested `MD/MD` directory, then move the ZIP into `MD`. The existing `OpenLB` and `ondemand` directories are left unchanged. The `negishi_md_cpu` directory from an earlier distribution is not used for this installation.

The setup job requests one node, one task, four CPUs, 16 GB of memory, and a maximum wall time of six hours on the `cpu` partition under the `ardekani` account. This is a job time limit, not an estimate of installation time. Submit the `sbatch` command from the login node; downloading, compilation, and testing run on a compute node.

The installation script performs the following steps:

1. Check that `antechamber`, `parmchk2`, and `tleap` are available through the existing `amber/24` module.
2. Create a dedicated Conda environment for CMake 3.31.6 and a preparation environment containing Packmol, RDKit, ParmEd, and SciPy.
3. Download the GROMACS 2026.4 source and matching regression tests, and verify their official MD5 checksums.
4. Build GROMACS and the required FFTW dependency with GCC 12.2.0.
5. Run the official `check` target and install into the user directory only after it passes.
6. Record the actual software versions, CMake configuration, package lists, and test logs.

GROMACS and CMake versions are pinned for reproducibility. The installer does not automatically select the latest release. The CPU engine uses a `Release` build, `AVX2_256` SIMD, mixed precision, thread-MPI, and OpenMP, with GPU support disabled. It supports CPU parallelism within a single node. External OpenMPI is loaded only as a dependency of the AMBER module.

Network access is required for conda-forge and the GROMACS and FFTW download servers. Administrator privileges and `sudo` are not required. CMake 3.31.6 satisfies the GROMACS 2026.4 requirement of CMake 3.28 or later.

## 2. Verify the Installation

Check the job queue:

```bash
squeue -u lee6456
```

A successful `setup_md_cpu_JOBID.out` ends with `MD CPU SETUP COMPLETE`. Then run:

```bash
cd /home/lee6456/MD
source ./env.sh md
gmx --version
```

Verify version `2026.4`, MPI library `thread_mpi`, OpenMP support, SIMD `AVX2_256`, and disabled GPU support. The installed executable is `/home/lee6456/MD/software/gromacs-2026.4-cpu/bin/gmx`.

If the tests corresponding to `make check` fail, the installer does not report successful completion. Inspect the Slurm log and `software-info/check.log` for the error. Resubmitting `sbatch setup_cpu.sbatch` after an interrupted installation reuses existing downloads and build files. Existing Conda environments are also reused without automatic upgrades.

### Retry After a Failed Installation Test

In `setup_md_cpu_44122858.out`, the remaining tests, including `acetonitrilRF`, passed, but one `SingleEnergyGroupSimd2xMM` comparison in `NbnxmTests` failed. This official test generates its input using the global `std::rand()` state and compares SIMD summation against sequential float summation. The observed difference was 5 ULP, while the test tolerance was 4 ULP. In a separate environment, the same executable, built without modifying the official source or tests, passed three consecutive runs under normal conditions. Changing the random-number state through external initialization code reproduced the failure. A diagnostic calculation using a double-precision sum of the same float inputs also identified cases in which the test failed even though the SIMD result was more accurate.

The Negishi logs show that `xalt/3.0.2` remains loaded after `module purge`. The sampling initialization code in XALT 3.0.2 calls the global `srand()` function and can therefore change the test input. The sampling configuration of the failed job was not recorded, so this was not established as the cause of that specific failure. The revised installer sets `XALT_SAMPLING=no` only for the test subprocess, after all modules have been loaded, and records the previous sampling setting. XALT tracking remains enabled, and the MD runtime environment is unchanged. The GROMACS source, reference values, tolerances, and full test list are preserved. Installation completes only after the entire official `check` target passes. Do not create `.validated` manually.

In the earlier `setup_md_cpu_44121885.out`, the OpenMP environment-variable conflict had been resolved and 105 of 106 tests passed. The remaining `acetonitrilRF` failure resulted from a small coordinate difference between the reference and computed trajectories crossing the 1.0 nm cutoff. In the failure data, the final distance between the affected atoms was 0.999972808 nm in the reference and 1.000001565 nm in the computed trajectory. The RF/LJ interactions used by this test have a force discontinuity at the cutoff. The calculated force change of approximately 10.09 kJ/mol/nm agrees with the observed difference of approximately 10.1 kJ/mol/nm.

The earlier script used the CPU allocation for compilation as the number of thread-MPI ranks for regression tests. The supplied CMake configuration and execution logs showed nine ranks. The revised installer separates compilation parallelism from test parallelism and runs regression tests with two ranks and one OpenMP thread per rank by default. The upstream CI default for thread-MPI is also two ranks. Official reference data, tolerances, and the test list are unchanged, and the full `check` target must pass before installation completes. Compilation continues to use the CPUs allocated by Slurm, while production MD follows the settings in `run_cpu.sbatch`.

Upload the revised ZIP as `/home/lee6456/MD/negishi_md_cpu.zip`, replacing the previous ZIP, then extract only the two updated files. Existing configuration and user changes, such as the CPU allocation in the setup job, are preserved.

```bash
cd /home/lee6456
unzip -o MD/negishi_md_cpu.zip MD/tools/install_cpu.sh MD/README.md
cd /home/lee6456/MD
sbatch setup_cpu.sbatch
```

Keep the existing `build`, `envs`, and `downloads` directories. The installer reuses compiled files, updates the test configuration, and reruns the full test suite. Previous failure data are preserved in `software-info/check-NEW_JOBID-before-retry.tar.gz`. If the retry fails, inspect both `setup_md_cpu_NEW_JOBID.out` and `software-info/check-NEW_JOBID-failed.tar.gz`. The diagnostic archive includes CTest logs, CMake configuration, `acetonitrilRF` comparison results and execution logs, and any available TPR/TRR files.

## 3. Directory Layout

| Path | Contents |
|---|---|
| `/home/lee6456/MD/config.sh` | Versions, modules, installation paths, and the default scratch location |
| `/home/lee6456/MD/env.sh` | Environment selection |
| `/home/lee6456/MD/setup_cpu.sbatch` | Installation job submission script |
| `/home/lee6456/MD/tools/install_cpu.sh` | Download, build, and test workflow |
| `/home/lee6456/MD/run_cpu.sbatch` | CPU job for an existing TPR |
| `/home/lee6456/MD/prepare_poc.sbatch` | Generate molecular coordinates, force-field parameters, and GROMACS inputs from a composition JSON |
| `/home/lee6456/MD/run_poc.sbatch` | Minimization, initial MD, benchmarking, and checkpoint restart |
| `/home/lee6456/MD/tools/prepare_case.py` | Three-stage input generation workflow |
| `/home/lee6456/MD/software/gromacs-2026.4-cpu/` | Installed GROMACS and locally built dependencies |
| `/home/lee6456/MD/envs/cmake-3.31.6/` | CMake build environment |
| `/home/lee6456/MD/envs/prep/` | Packmol, RDKit, and ParmEd environment |
| `/home/lee6456/MD/downloads/` | Downloaded GROMACS and test archives |
| `/home/lee6456/MD/src/` | Extracted source code and regression tests |
| `/home/lee6456/MD/build/` | Compilation files and FFTW build files |
| `/home/lee6456/MD/cache/` | Conda packages and reusable PLGA coordinates and charge calculations |
| `/home/lee6456/MD/tmp/` | Temporary installation files |
| `/home/lee6456/MD/software-info/` | Package lists, versions, configuration, and installation test records |
| `/home/lee6456/MD/cases/` | Composition JSON files and generated molecular structures, force fields, and MDP files |
| `/home/lee6456/MD/validation/` | Summaries of input generation checks performed locally |
| `/scratch/negishi/lee6456/md_JOBID_UTC_TIMESTAMP/` | Trajectories, energies, checkpoints, and other outputs for each run |

Existing system modules continue to reference their shared installations. Conda may store environment registration metadata in the home directory, but new software and package caches are installed under `MD` as listed above. Keep the `MD` directory in its original location after creating the Conda environments.

## 4. Select a Working Environment

Select the required environment in each shell. `env.sh` clears existing modules and the active Conda environment before loading the requested environment. It does not modify `.bashrc`.

```bash
# MD engine
source /home/lee6456/MD/env.sh md
gmx --version

# Structure generation, molecular packing, and topology conversion
source /home/lee6456/MD/env.sh prep
python -c 'from rdkit import Chem; import parmed; print("OK")'

# AmberTools
source /home/lee6456/MD/env.sh amber
command -v antechamber parmchk2 tleap

# CMake, when needed
source /home/lee6456/MD/env.sh build
cmake --version
```

## 5. Prepare PLGA / Drug / Mixed-Solvent Inputs

An existing, completed GROMACS installation can be reused. Upload the new ZIP as `/home/lee6456/MD/negishi_md_cpu.zip`, then add or update only the input preparation files. The following command preserves the existing `config.sh`, `env.sh`, and installed software.

```bash
cd /home/lee6456
unzip -o MD/negishi_md_cpu.zip 'MD/tools/mdprep_*.py' MD/tools/prepare_case.py MD/prepare_poc.sbatch MD/run_poc.sbatch 'MD/cases/*' 'MD/validation/*' MD/README.md
cd /home/lee6456/MD
```

The distribution ZIP also includes **pregenerated inputs for the default DMSO system in `cases/poc_lawsone_dmso_input/`**. The first CPU run can therefore be submitted directly:

```bash
sbatch run_poc.sbatch cases/poc_lawsone_dmso_input
```

Use `prepare_poc.sbatch` as described below to regenerate inputs or prepare another condition. For solvent comparisons, generate every condition, including the reference DMSO system, with the same preparation environment on Negishi.

The default configuration is `cases/poc_lawsone_dmso.json`. It uses **30 PLGA chains with 90 repeat units per chain**, an LA:GA ratio of 50:50, a random D/L sequence, hydroxyl/carboxyl end groups, one lawsone molecule, and DMSO/water in an 8 × 8 × 14 nm box. This proof of concept (PoC) retains the 2,700 total repeat units of the earlier 6 × 450 configuration while shortening each chain to approximately 5.87 kDa. The number of chain ends and the chain relaxation time differ from those of the experimental polymer; molecular-weight effects are not assumed to be eliminated. The experimental material, Resomer RG 503 H, is an acid-terminated 50:50 D,L-PLGA. Its exact LA/GA sequence and stereochemistry are unknown and are represented by random sequences generated with a fixed seed.

The initial polymer placement region is **8 × 8 × 8 nm**. The system starts with a loosely packed central polymer region and approximately 3 nm of external solvent on each side, rather than a dense polymer slab constrained to an initial thickness of 4 nm. `slab_nm` defines the initial placement region; it does not impose walls or positional restraints during MD. Solvent molecules are also placed in voids between the polymer chains. This initial geometry does not prescribe the equilibrium polymer thickness or density. After relaxation, verify whether distinct compositions inside and outside the polymer region are maintained.

The initial `water_mole_fraction = 0.9` specifies 90 mol% water among the water and organic-solvent molecules. It is a starting composition for exploring the mechanism, rather than a direct reproduction of the experimental final volume ratio. Molecular counts are estimated using additive volumes based on pure-liquid densities and recorded in the manifest. Changing the organic solvent at the same mole fraction changes the molecular volumes, so both water and organic-solvent counts are recalculated. The initial density relaxes during NPT.

The default DMSO system contains 20,340 PLGA atoms, one lawsone molecule, 1,628 DMSO molecules, and 14,655 water molecules, for a **total of 80,604 atoms**.

Input preparation runs on a compute node using the existing RDKit, Packmol, ParmEd, and AmberTools installations. SciPy is required for the coordinate optimization that removes initial interchain overlaps. If it is missing from the preparation environment, `prepare_poc.sbatch` adds this dependency on a compute node while keeping existing package versions frozen. GROMACS is not rebuilt.

```bash
cd /home/lee6456/MD
sbatch prepare_poc.sbatch cases/poc_lawsone_dmso.json
```

On completion, `cases/poc_lawsone_dmso_input/` contains:

| File | Contents |
|---|---|
| `case.json` | Copy of the configuration used |
| `polymer.json` | LA/GA sequence, D/L stereochemistry, atom mapping, and generation records for each chain |
| `molecules/` | Molecular coordinates, bonds, charges, Amber topologies, and calculation logs |
| `charge_fragments/` | Short-oligomer charge calculations and records supporting charge transfer to longer chains |
| `manifest.json` | Atom and molecule counts, actual composition, and configuration identifiers |
| `packmol.inp`, `packmol.log`, `packed.pdb` | Packing input, log, and coordinates for the mixture |
| `system.gro`, `topol.top` | Initial GROMACS coordinates and topology |
| `mdp/` | Minimization, initial relaxation, and short benchmark settings |
| `READY.json` | Record of completed input assembly; this does not indicate completed equilibration |

Coordinates and parameters are generated programmatically. Each PLGA chain is built by connecting repeat units with defined bond lengths, angles, and stereochemistry. Geometric optimization then removes interchain overlaps while preserving distances associated with bonds and angles, as well as stereochemistry. This prepares initial coordinates; it does not produce a structure equilibrated with GAFF2. After passing initial nonbonded distance checks of 1.4 Å for heavy atoms and 1.2 Å for all atoms, the system must undergo actual GAFF2 energy minimization. No soft potential is introduced to alter the MD force field or minimization criterion.

ETKDG and quantum-chemical calculations are not applied to entire long chains. Drug and organic-solvent coordinates are generated from SMILES. Packmol places the remaining molecules without changing the internal bonds of the prepared chains. **The polymer coordinates are fixed only during the initial Packmol placement; PLGA is free to move during MD.**

### 5.1 Force Field and Chemical States

GAFF2 is used for the drug, organic solvent, and PLGA; TIP3P is used for water. All components are represented at the atomic level. Small-molecule charges are calculated with AmberTools AM1-BCC. PLGA charges are calculated from ten short reference chains, each containing five repeat units. Eight reference chains cover all combinations of the central repeat unit and its LA/GA neighbors; two additional chains cover LA terminal combinations. Internal charges are taken from the central repeat unit of each reference chain, and terminal charges from the corresponding chain ends. Charges are transferred to the long chains by matching neighboring LA/GA identities, terminal environments, and atom roles. The same charges are used for D/L centers. Equivalent hydrogen charges are averaged, and corrections that set the total charge of each neutral repeat unit to zero are recorded separately. This is a **PoC parameterization** based on approximate charge transfer between local environments. It has not been validated as a PLGA force field that reproduces experimental drug loading or partition free energies.

Additional GAFF2 bonded parameters for the polymer are obtained from the short reference chains using `parmchk2`. Preparation stops if different coefficients are found for the same term. The same local terms are not repeatedly searched for each long chain; `tleap` checks that all required bond, angle, and torsion terms are present for every full chain. Generate comparison cases with the same preparation tools and force-field versions.

All drugs and PLGA end groups are currently modeled in their **neutral states**. The model does not infer the experimental pH or ionization fractions when these are unspecified in the paper. In particular, the charge states of lawsone and the acid end groups require separate assessment for mechanistic comparisons. Poloxamer, salts, solvent removal, and washing are outside the scope of this local partitioning model.

Molecular generation records, quantum-chemical output, charge-transfer records, and neutralization corrections are retained. Preparation stops on missing force-field terms, atom-count mismatches, Packmol nonconvergence, severe intermolecular overlaps, or GROMACS input errors. Test tolerances are not relaxed, and `grompp -maxwarn` is not used to bypass errors.

### 5.2 First CPU Run

After input preparation is complete, submit:

```bash
cd /home/lee6456/MD
sbatch run_poc.sbatch cases/poc_lawsone_dmso_input
```

The default request is eight CPUs, 16 GB of memory, and up to six hours. A single job runs energy minimization, 100 ps of NVT, 200 ps of NPT, and a 100 ps MD benchmark in sequence. NVT uses a 1 fs timestep; subsequent MD uses a 2 fs timestep with bond lengths involving hydrogen constrained. Polymer and drug positions are not restrained. NPT maintains the slab area in the x/y plane while allowing the box length in z to relax. These durations define an initial run and do not guarantee equilibration or convergence of partition free energies.

All stages share **one six-hour time budget**, with checkpoints written before the job exits. Incomplete stages are not marked as completed. Results are written to `/scratch/negishi/lee6456/md_JOBID_UTC_TIMESTAMP/`, and Slurm logs remain in the submission directory. To resume a run, use the actual results path printed in its log:

```bash
sbatch run_poc.sbatch --restart /scratch/negishi/lee6456/md_JOBID_UTC_TIMESTAMP
```

Restarting uses the same results directory, skips completed stages, and resumes an incomplete stage from its checkpoint. The `Performance:` entry in the successfully completed `benchmark.log` reports the measured simulation throughput in ns/day.

### 5.3 Change the Solvent or Drug

Copy the reference DMSO configuration under a new name:

```bash
cp cases/poc_lawsone_dmso.json cases/poc_lawsone_acn.json
```

Set `mixture.solvent` to `"ACN"` in the copied JSON, then submit:

```bash
sbatch prepare_poc.sbatch cases/poc_lawsone_acn.json
```

Supported solvent values are `DMSO`, `DMF`, `acetone`, and `ACN`. Supported `drug.name` values are `lawsone`, `menadione`, and `plumbagin`. Each new condition uses a separate output directory to preserve existing inputs.

**Changing only the solvent while keeping the polymer settings, box, and seeds unchanged reuses the same initial PLGA coordinates.** The repeat-unit sequences, chain stereochemistry, and polymer topology are also preserved. The organic-solvent structures, parameters, molecular counts, packing, and mixture topology are regenerated. PLGA coordinates and charge calculations are cached in `MD/cache/`. If the cache is missing or the RDKit version differs, completed inputs in `cases/*_input/` with matching polymer settings can still supply the same PLGA PDB coordinates. New caches preserve coordinates in double precision. Solvent-dependent equilibrium PLGA structures are obtained through MD relaxation of each system; sharing initial coordinates does not constrain the systems to the same equilibrium structure.

Change `polymer.seed` to compare different initial polymer structures. Change `mixture.packing_seed` to vary only the solvent placement. Set `drug.initial_region` to `"polymer"` to place the drug initially in the polymer region. Short residence times in trajectories started from different locations are insufficient to determine an equilibrium partition free energy.

## 6. Run an Existing TPR Independently

`run_cpu.sbatch` runs a single prepared TPR. The multistage workflow and its restart handling are provided by `run_poc.sbatch`.

```bash
sbatch run_cpu.sbatch cases/my_case/model.tpr
```

## Validation Scope and References

The input generation tools were checked in a separate Linux environment using actual installations of AmberTools 24.8, RDKit, Packmol, ParmEd, and GROMACS 2026.4.

- AM1-BCC/GAFF2 parameterization completed for all three supported drugs and four organic solvents. Atom identities and ordering, and the atom pairs in all covalent bonds, matched the original molecular graphs. The generated topologies were neutral and contained no missing force-field terms.
- An 8,878-atom mixture containing a PLGA 12-mer, lawsone, DMSO, and TIP3P water was assembled. Minimization reached Fmax = 879.966 kJ/mol/nm in 459 steps, below the specified threshold of 1000 kJ/mol/nm. NVT for 0.1 ps, NPT for 0.2 ps, and subsequent MD for 0.2 ps completed with finite recorded energies and no LINCS warnings. This short run checked input integration and initial stability.
- Interruption and checkpoint restart of CPU jobs were tested with GROMACS on a small water system. The workflow resumed in the same results directory and did not advance past an incomplete stage.
- All 1,350 stereocenters and the nonbonded distance criteria were checked for the default 30 × DP90 polymer structure. During DMSO-to-ACN input generation, all 30 PLGA PDB files were reused byte for byte from the completed DMSO inputs, even without the previous binary cache. The ACN system size was recalculated as 78,877 atoms.
- The **full default system of 80,604 atoms** reached Fmax = 989.41772 kJ/mol/nm in 1,954 GROMACS minimization steps. A subsequent NVT run of 100 steps at 1 fs completed with finite energies and no LINCS warnings. All 1,350 stereocenters retained their signs after minimization and NVT. NVT was shortened to 0.1 ps only for this check; the distributed MDP specifies 100 ps. This validates an executable initial system, not temperature or density equilibration or the partitioning mechanism.

SQM prints individual Mulliken charges to three decimal places. Consequently, the AM1-BCC output charge sums of some neutral molecules had residuals of approximately 0.001–0.002 e. After confirming that the SQM total charge is neutral, the code corrects only residuals within the output precision and preserves both the original values and the corrections. Neutralization corrections used for PLGA repeat-unit charge transfer are recorded separately. QM convergence criteria, GROMACS tolerances, and warning-handling criteria are unchanged.

Bash and Python syntax, paths, and command-line argument handling were also checked. During investigation of earlier installation failures, GROMACS 2026.4 was built in a separate Linux environment, and the official `acetonitrilRF` test and the supplied TPR were run with 1, 2, 4, and 9 ranks. The same force-comparison failure was reproduced with 9 ranks, while runs with 1, 2, and 4 ranks passed. All 42 official complex tests and 7 essentialdynamics tests also passed with the default configuration of 2 ranks and 1 OpenMP thread per rank.

That reproduction environment used an Intel Xeon CPU, GCC 13.3, AVX2_256, mixed precision, thread-MPI, and FFTPACK. The user's GCC 12.2/FFTW installation settings were unchanged. The `acetonitrilRF` test uses a reaction-field method and therefore does not use FFTs for the relevant force calculation. On Negishi, final installation validation is performed by the full official `check` target in the setup job. These checks do not validate PLGA partition free energies or experimental drug loading.

- [GROMACS 2026.4 downloads and official source/test MD5 checksums](https://manual.gromacs.org/2026.4/download.html)
- [GCC, CMake, FFTW, SIMD, thread-MPI, and testing](https://manual.gromacs.org/2026.4/install-guide/index.html)
- [Known compiler-related issues](https://manual.gromacs.org/2026.4/user-guide/known-issues.html)
- [XALT 3.0.2 sampling initialization](https://github.com/xalt/xalt/blob/xalt-3.0.2/src/libxalt/xalt_initialize.c#L725-L742)
- [AMBER module on Negishi](https://docs.rcac.purdue.edu/software/apps_md/amber/)
- [Packmol](https://m3g.github.io/packmol/download.shtml)
- [RDKit](https://rdkit.org/docs/Install.html)
- [ParmEd](https://parmed.github.io/ParmEd/html/readwrite.html)
- [Resomer RG 503 H material specifications](https://www.evonik.com/en/products/hc/pr_52021485.html)
- [Experimental paper](https://doi.org/10.1002/mame.202500269)
- [Polymer charge-transfer methodology](https://doi.org/10.1021/acs.jcim.3c01691) — methodological reference; this paper does not provide the GAFF2/Amber AM1-BCC parameters implemented here.
- [SQM output precision for individual Mulliken charges](https://github.com/Amber-MD/AmberClassic/blob/main/src/sqm/qm2_print_charges.F90)
