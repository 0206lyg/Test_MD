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

Edit the input block at the top of `run_wow_sema.sbatch`, then submit it without arguments. `DROPLET_DIAMETER_NM=39.2` is the wet diameter at the O/W₂ boundary, `W1_COUNT=10` selects 1–10 equal-volume pockets, and `INIT_SEED=20261120` controls the initial configuration. At fixed diameter, pocket count and seed preserve the amounts of PLGA, DCM, W₁ water and drug. The generator places pockets randomly within the available volume with a 1.05 nm minimum oil film. It checks the actual bead density for separate W₁ regions disconnected from W₂ at step zero. Each fresh job builds its coordinates on scratch from compact templates. Set `CHECK_INPUT=1` to generate the same files and check step zero before production. Command-line overrides remain available (`bash run_wow_sema.sbatch --help`).

The default job requests 8 CPU threads and up to 48 h for 778,156 beads. After a 1,000-step closed relaxation, the code attempts DCM exchange in the external water every 100 steps. Every 1,000 extraction steps it saves the particle composition and coordinates. The first sampled particle with **DCM/(DCM + PLGA) ≤ 5 wt% by mass** is the endpoint.

| File under the scratch run directory | What to read |
|---|---|
| `code/wow_dpd/initial/layout.json`, `initial_preview.svg` in the same folder | Generated composition, pocket geometry and preview. |
| `result/particle.csv` | Local `dcm_wt_fraction_particle`, internal water fraction, water origins and connectivity to the bath. |
| `result/sink.csv` | Eligible and exchanged DCM counts at every exchange decision. |
| `result/trajectory.lammpstrj` | Full bead coordinates including the endpoint frame. Multiply coordinates by 0.70 for nm. |
| `result/state.json` | `target_reached`, `paused` or `max_steps`, plus the committed step and restart file. |

Set `RESUME_DIR` to the **Results directory** printed in the job log to continue a production run, including a run created with the older fixed input. The job reuses its saved coordinates, code, and options. Set `RESUME_MAX_STEPS` to extend the cumulative cap if needed.

For binary mixtures, use `sbatch run_binary.sbatch --continue RUN_DIR` to resume or `sbatch postprocess_binary.sbatch RUN_DIR` to analyze saved frames. Single-peptide results are in `analysis/summary.json`; `sbatch postprocess_single_sema.sbatch RUN_DIR` repeats that analysis. Prepared PLGA/drug runs resume with `sbatch run_poc.sbatch --restart RUN_DIR`.
