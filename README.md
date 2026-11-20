# Test_MD on Purdue Negishi

## Shared CPU setup

From the repository root (default `~/MD`), run `sbatch setup_cpu.sbatch` once after installing or updating the dependencies. If the checkout is elsewhere, export `MD_ROOT="$PWD"` first. Wait for `MD CPU SETUP COMPLETE` in the setup log before submitting a case.

`tools/install_cpu.sh` provisions GROMACS 2026.4, the molecular preparation environment and the LAMMPS DPD runtime. Existing validated GROMACS and compatible DPD installations are reused. `config.sh` defines paths and module names. `env.sh` selects `md`, `prep`, `amber`, `build` or `dpd`. Jobs inherit `MD_SCRATCH_ROOT`, which defaults to `/scratch/negishi/$USER`. Slurm logs stay in the submission directory.

## Submit a case

| Command | Calculation |
|---|---|
| `sbatch run_wow_sema.sbatch` | DCM extraction from a PLGA W₁/O/W₂ particle with 10 inner water droplets and 20 semaglutide molecules. |
| `sbatch run_single_sema.sbatch` | One semaglutide molecule in explicit water with 0.150 M added NaCl. Default production is 1 ns. Coordinates and parameter provenance are in `cases/single_sema/NOTICE.txt`. |
| `sbatch run_binary.sbatch DMSO` | Water mixed with DMSO, DMF, acetone or acetonitrile. The job prepares the mixture and measures its spatial composition. |
| `sbatch prepare_poc.sbatch cases/poc_lawsone_dmso.json` | Generate the PLGA/lawsone inputs from the case specification. Run once per input composition. |
| `sbatch run_poc.sbatch cases/poc_lawsone_dmso_input` | Run the prepared PLGA/lawsone case after preparation succeeds. |
| `sbatch run_cpu.sbatch model.tpr` | Run an existing GROMACS TPR on CPU. |

## W/O/W extraction and results

The default job requests 8 CPU threads and up to 48 h for 778,156 beads. After a 1,000-step closed relaxation, the code attempts DCM exchange in the external water every 100 steps. Every 1,000 extraction steps it saves the particle composition and coordinates. The first sampled particle with **DCM/(DCM + PLGA) ≤ 5 wt% by mass** is the endpoint.

| File under the scratch run directory | What to read |
|---|---|
| `result/particle.csv` | Local `dcm_wt_fraction_particle`, internal water fraction, water origins and connectivity to the bath. |
| `result/sink.csv` | Eligible and exchanged DCM counts at every exchange decision. |
| `result/trajectory.lammpstrj` | Full bead coordinates including the endpoint frame. Multiply coordinates by 0.70 for nm. |
| `result/state.json` | `target_reached`, `paused` or `max_steps`, plus the committed step and restart file. |

Copy the **Results directory** path printed in the job log to resume, for example `sbatch run_wow_sema.sbatch --resume /absolute/path/to/wow_sema_run`. For a run at its cumulative step cap, append `--max-steps 3000000` to extend it. Use `--check-input` for a zero-step input check. The source, inputs and water-analysis code are in `scripts/wow_dpd/`.

For binary mixtures, use `sbatch run_binary.sbatch --continue RUN_DIR` to resume or `sbatch postprocess_binary.sbatch RUN_DIR` to analyze saved frames. Single-peptide results are in `analysis/summary.json`; `sbatch postprocess_single_sema.sbatch RUN_DIR` repeats that analysis. Prepared PLGA/drug runs resume with `sbatch run_poc.sbatch --restart RUN_DIR`.
