# Single semaglutide in water: short first run

This overlay adds a **one-molecule, fixed-protonation MD screen** to the existing
`/home/lee6456/MD` installation. It reuses `config.sh`, `env.sh`, the preparation
Python environment, and the installed CPU GROMACS. No software reinstall or
online parameter generation is performed when the job starts.

## Submit

Place the overlay archive directly in `/home/lee6456/MD`, then run:

```bash
cd /home/lee6456/MD
unzip -o single_sema_overlay.zip
sbatch run_single_sema.sbatch
```

Each submission creates a **fresh** directory under
`/scratch/negishi/lee6456/single_sema_JOBID_TIMESTAMP`. The Slurm output prints
the exact directory. Existing binary-solvent and PLGA scripts are not replaced.

| Setting | Default |
|---|---|
| Solute | One complete semaglutide: 31 residues and 291 heavy atoms |
| Water and ions | OPC water; neutralizing Na+ plus approximately 150 mM NaCl |
| Protonation | Fixed state representing pH 7.4; neutral HIE histidine; solute charge −4 e |
| Temperature / pressure | 298.15 K / 1 bar in NPT and production |
| Box | Prebuilt truncated octahedron; exact dimensions, water/ion counts and concentration are in `input/manifest.json` |
| Stages | Energy minimization → 50 ps NVT → 200 ps NPT → **1 ns production** |
| Resources | 8 CPU cores, 16 GB, **6-hour shared job limit** |

The bundled box has 5,795 water molecules, 20 Na+ and 16 Cl- ions: 16 salt
pairs plus four neutralizing Na+. It contains **17,999 physical atoms** and
23,794 integration sites including OPC virtual sites. Initial box vectors are
6.443 nm long with 109.471-degree angles (205.87 nm³).

The initial packing volume contracts during NPT. Salt count is therefore
calibrated against a short pilot's NPT volume, 179.44 nm³: 16 pairs correspond
to 148.1 mM there. The initial packing concentration is 129.1 mM, not 150 mM;
production volume determines the actual concentration, which analysis reports
as `added_NaCl_mM`. Neutralizing Na+ is counted separately. Integer ion counts
and NPT volume fluctuations prevent an exactly constant 150.000 mM value.

NVT/NPT restrain solute heavy atoms while the solvent relaxes. Production has no
positional restraints. The wall-time limit is shared across all stages; it does
not guarantee that 1 ns will finish on a particular node. If time runs short,
available results are analyzed and explicitly marked **PARTIAL**. Checkpoints
are retained; this launcher does not silently resume an earlier calculation.

Optional fresh velocity-seed replicate:

```bash
sbatch run_single_sema.sbatch --seed 1730
```

A different seed changes initial velocities, not the starting conformation.
This alone does not test convergence from distinct conformational states.
`--production-ns`, `--nvt-ps`, `--npt-ps`, and `--temperature` are also accepted.

## What was reused from the paper

The coordinate template comes from Supplementary Data 1 of Colston, Faivre and
Schneebeli, *Nature Communications* **16**, 9512 (2025):
<https://doi.org/10.1038/s41467-025-64891-0> (CC BY 4.0).

**Only the coordinate template is reused. This is not a reproduction of that
paper's CHARMM36m/CGenFF constant-pH force field.** The supplied source contains
all 31 residues and 291 heavy atoms. Source proton-bearing/dummy hydrogen
positions are not transferred to the fixed-state model. The linker γ-Glu
coordinate center was corrected from R to S to match the official L-γ-Glu
chemical structure; the modification is recorded in the parameter audit.

The bundled topology is a **provisional AMBER-based model**: ff14SB for canonical
residues; custom AMBER atom types and **single-point AM1-BCC charges on
converged MMFF94s ACE–X–NME geometries** for Aib and lipidated Lys; OPC water; and Li–Merz 12–6 OPC ion
parameters. Parameter
coverage and chemical identity checks do **not** establish the accuracy of
custom torsions, solvation energetics or binding kinetics. Source coordinates,
parameter provenance and checksums are included with the assets.

The scripts under `tools/build_sema_parameters.py` and
`tools/solvate_sema_parameters.py` document asset construction. They are not
called during normal job submission and require additional developer tools.
The single-point calculation uses converged AM1 electronic SCF (`scfconv=1e-10`)
and `maxcyc=0`; it does not claim AM1 geometry optimization. This documented
Antechamber option and its geometry dependence are recorded in the audit.
The parameter source directory contains the SQM inputs/outputs, equivalent-atom
charge checks, and every charge adjustment made when removing the caps.

## Read the results

Start with `analysis/report.txt` and `run_status.json`. Then inspect:

| File | Contents |
|---|---|
| `analysis/summary.json` | Five validation items, model provenance, numerical diagnostics and completion status |
| `analysis/diagnostics.svg` | Backbone RMSD, Rg, C18–peptide contacts and a backbone helical-region proxy |
| `analysis/timeseries.csv` | Numerical values, including the minimum distance to periodic copies |
| `analysis/thermodynamics.xvg` | Temperature, pressure, potential energy and density where available |
| `analysis/semaglutide.gro` + `analysis/semaglutide.xtc` | Whole, solute-only structure and trajectory; water/ions are excluded |

For VMD, keep the run directory layout and enter the following in **Tk Console**:

```tcl
source /path/to/run/view_single_sema.vmd
```

Analysis uses the latest available stage, including an incomplete production
trajectory. The report identifies the sampled stage. To regenerate analysis:

```bash
sbatch postprocess_single_sema.sbatch /scratch/negishi/lee6456/single_sema_JOBID_TIMESTAMP
```

The five goals are **not** automatically marked PASS after a successful run:
chemical identity and assigned charge can be checked; parameter coverage can be
audited; numerical problems and initial structural behavior can be screened.
**A 1 ns trajectory cannot validate equilibrium sampling or experimental
force-field accuracy.** The helical-region metric is a φ/ψ classification, not
DSSP and not a calculated CD spectrum. The analysis does not fit or claim
agreement with a matched monomer experiment, and this run does not measure
PLGA release or ECM binding.

## Execution checked before packaging

The supplied solute/OPC/ion system completed **EM + 10 ps NVT + 20 ps NPT +
20 ps production** and all analysis exports in local GROMACS 2024.5, using
four CPU threads. The final salt setting gave a sampled mean of **148.4 mM**;
mean temperature was **297.95 K**; minimum solute-to-periodic-copy distance
was **2.67 nm**; all 30 stereocenters retained their assignments in the three
analyzed production frames. No LINCS warning or numerical sanity alert was
found. Evidence is in `validation/single_sema/local_execution_check.json`.

This verifies short execution, not the full default 1 ns trajectory, a remote
Negishi/Slurm run, or physical accuracy of the provisional force field.
