# Binary organic-solvent / water screening on Negishi

This v2 overlay adds `run_binary.sbatch`, `postprocess_binary.sbatch`,
`tools/prepare_binary.py`, `tools/analyze_binary.py`, and this guide to the existing Test_MD installation.
It uses the installed environments configured by `env.sh` and `config.sh` and
the existing `tools/prepare_case.py` and `tools/mdprep_forcefield.py` helpers.
The inspected upstream commit was `205139762df71a58a03b2ec66a1dffd567a3d0ef`.

## Purpose and expected result

At 298.15 K and 1 bar, pure binary mixtures of water with DMSO, DMF, acetone,
or acetonitrile are expected to remain a single liquid phase. Local clustering
and composition fluctuations are compatible with complete macroscopic miscibility.
This calculation checks the binary reference system and the behavior of its
force field. It does not measure drug encapsulation efficiency.

In the nanoprecipitation experiment, a water-miscible solvent transports PLGA
and drug into water. The deteriorating solvent quality produces polymer-rich
particles while solvent and water mix. This is a different process from
liquid-liquid demixing of the two solvents alone. The original PoC (proof of
concept) checked the preparation and MD pipeline; a preplaced polymer slab
and one drug molecule do not reproduce nanoparticle formation or quantify EE.

## Install and start a new calculation

Upload `binary_solvent_overlay_v2.zip` into `/home/lee6456/`, then run:

```bash
unzip -o /home/lee6456/binary_solvent_overlay_v2.zip -d /home/lee6456/MD
sbatch --chdir=/home/lee6456/MD /home/lee6456/MD/run_binary.sbatch DMSO
```

The second command can use `DMF`, `Acetone`, or `Acetonitrile` instead.
Names are case-insensitive; `ACN` is also accepted. Each invocation creates a
fresh, uniformly packed calculation, with no PLGA, drug, surfactant, ions, or
pre-existing liquid interface. It never selects an earlier run automatically.
Preparation runs inside the same compute job, not on the login node.

To submit all four independent solvents:

```bash
for solvent in DMSO DMF Acetone Acetonitrile; do
    sbatch --chdir=/home/lee6456/MD /home/lee6456/MD/run_binary.sbatch "$solvent"
done
```

Preview molecule counts without preparing molecules or submitting a job:

```bash
python3 /home/lee6456/MD/tools/prepare_binary.py plan --solvent DMSO
```

Change a setting explicitly, for example an independent packing/velocity seed:

```bash
sbatch --chdir=/home/lee6456/MD /home/lee6456/MD/run_binary.sbatch DMSO --seed 1730
sbatch --chdir=/home/lee6456/MD /home/lee6456/MD/run_binary.sbatch DMSO --box-nm 8 --production-ns 50
```

The default request is one CPU task, 8 CPU threads, 16 GB and 24 hours.
Slurm options such as `--time=2-00:00:00` must precede the script path and must
fit the account/partition limits. These are requested resources, not a measured
runtime prediction. The launcher uses the actual allocated CPU count and clears
the duplicated `SLURM_TRES_PER_TASK` request for each `srun` step, avoiding the
earlier 9-versus-8 CPU conflict.

## Composition and box

Sabura et al., Experimental Section 4.3, uses **1.5 mL organic solution into
10 mL aqueous solution**. The aqueous phase contains 0.5% poloxamer 407 in the
experiment. This binary reference deliberately omits that third component.
The ratio represents the nominal composition after injection, before solvent
evaporation or purification. It does not reproduce the injection process or
the changing local composition during nanoprecipitation.

The default pure-component volume ratio is therefore **1.5:10**, giving a
nominal organic volume fraction of 0.130435. It is not an organic mole fraction
of 0.130435, and it is not 90 mol% water. For an initial cube of volume V:

    N_i = round(602.214076 * density_i * V_nm3 * volume_fraction_i / molar_mass_i)

The conversion assumes additive initial pure-liquid volumes. NPT then relaxes
the mixture volume and density without changing the molecule counts. Densities
are room-temperature initialization estimates, not simulated mixture densities.

| Solvent | Water molecules | Organic molecules | Organic mole fraction | Total atoms |
| --- | ---: | ---: | ---: | ---: |
| DMSO | 6260 | 238 | 0.03663 | 21160 |
| DMF | 6260 | 219 | 0.03380 | 21408 |
| Acetone | 6260 | 229 | 0.03529 | 21070 |
| Acetonitrile | 6260 | 325 | 0.04935 | 20730 |

The **6 x 6 x 6 nm periodic cube** has no preferred interfacial direction and
is a practical first screen of a homogeneous bulk liquid. Its length comfortably
exceeds twice the 1 nm real-space cutoff; PME treats long-range electrostatics.
Isotropic pressure coupling preserves the cubic shape while allowing density
relaxation. Approximately 21,000 atoms keep the first comparison manageable on
the existing CPU installation. A cube is not a thermodynamic-limit calculation.
If a persistent domain appears, repeat with independent seeds and an 8 nm cube
(about 2.37 times as many atoms) before interpreting size-dependent segregation.
An elongated box is useful for a deliberately prepared planar-interface study,
which is not this homogeneous-mixture test.

## Molecular model and stages

Organic molecules use the existing GAFF2 / AM1-BCC workflow; water is TIP3P.
Charges, topology assignment, atom ordering, and severe periodic contacts are
checked. No positional restraints are applied. There is one thermostat group
for the complete system. The model is a consistent screening starting point,
not a mixture-specific thermodynamic validation.

Amber is the parameter-generation software used here, not the complete name of
the physical model. `antechamber -at gaff2 -c bcc` assigns organic parameters and
charges; `tleap` loads `leaprc.gaff2` and `leaprc.water.tip3p`; ParmEd converts
the topology, and GROMACS integrates the equations of motion.

The original PLGA PoC uses the same small-molecule and water models. Its polymer
charges are a custom transfer from AM1-BCC pentamers, with neutrality corrections.
This is **not** the published PLGA-specific GAFF/RESP parameterization. The
binary workflow contains no polymer and does not call that charge-transfer path.

For a long-term PLGA production model, the published GAFF/RESP PLGA models are
a better-supported candidate than the custom pentamer-BCC transfer: Andrews
et al. compared RESP and BCC for PLGA(50:50), and RESP better reproduced the
caloric curve/glass transition. This evidence concerns PLGA properties, not
the accuracy of all drug/solvent cross interactions. Reproducing that model
requires its full parameter set, not simply attaching its charges to GAFF2.
The relevant sources are https://doi.org/10.1016/j.polymer.2020.122903 and
https://doi.org/10.1021/acs.jpcb.9b06681.

For binary reference calculations, Honti et al.'s optimized CS2-based DMF with
TIP4P/2005 directly targets mixture free energies; the acetone study below
supports CGenFF acetone with TIP4P/epsilon over the other tested combinations.
These are **paired-model references**, not an interchangeable collection for
encapsulation. Using different water models for the four solvents would confound
the solvent comparison. Selecting one common production framework requires
checking solvent-water mixing thermodynamics and drug transfer between mixed
solvent and PLGA, rather than relying only on neat-liquid density or visual mixing.
This v2 overlay upgrades postprocessing; its simulation model remains the
explicitly provisional GAFF2/AM1-BCC + TIP3P baseline.

| Stage | Default duration | Integration / coupling |
| --- | ---: | --- |
| Energy minimization | up to 50,000 steps | steepest descent; must meet force tolerance |
| NVT | 100 ps | 1 fs, 298.15 K, velocity-rescale thermostat |
| NPT | 1 ns | 2 fs, isotropic C-rescale, 1 bar |
| Production | 20 ns | 2 fs, same NPT conditions; coordinates every 10 ps |

All MD stages are retained, including early relaxation. Initial velocities and
packing depend on `--seed`. A 20 ns production window is a starting observation
period, not proof of equilibrium or sufficient sampling of a phase transition.
Inspect density and composition over time before judging whether the mixture
has stabilized. A mixed initial condition alone cannot exclude slow nucleation;
an intentionally segregated starting condition would be a separate convergence
test if necessary. No evaporation is modeled in this closed periodic system.

Model dependence matters: Honti et al. compared 15 DMF/water model combinations
and only four recovered the known complete miscibility. That paper did not
establish the accuracy of this GAFF2/AM1-BCC + TIP3P combination. Unexpected
demixing is therefore a reason to validate the force field and sampling, not
evidence that experimental ambient DMF and water are immiscible.

## Where to look

The Slurm output prints the exact results directory:

    /scratch/negishi/lee6456/binary_SOLVENT_JOBID_UTC_TIMESTAMP

| Output | Use |
| --- | --- |
| `input/binary_config.json`, `input/manifest.json` | actual settings, counts and composition |
| `input/system.gro` | initial uniformly packed configuration |
| `nvt.xtc`, `npt.xtc`, `production.xtc` | trajectories for all MD stages |
| `visualization.gro`, `production.xtc` | topology/coordinates and production movie for VMD |
| `view_binary.vmd` | ready-made organic/water representations in VMD |
| `analysis/production/concentration_timeseries.csv` | long-wavelength composition fluctuations over time |
| `analysis/production/composition_profiles.csv` | per-frame molecular composition along x, y and z |
| `analysis/production/concentration_diagnostics.svg` | ready-to-open trend graph and final directional profiles |
| `analysis/production/dispersion_timeseries.svg` | random-mixing-corrected dispersion trend |
| `analysis/production/projection_last.svg` | last-frame XY/XZ/YZ projected composition maps |
| `analysis/production/projection_time_average.svg` | time-averaged projected composition maps |
| `analysis/production/concentration_3d.pvd` and `vtk/` | time-resolved 3D coarse-grained composition for ParaView |
| `analysis/production/organic_molecules_3d.pvd` and `vtk/` | time-resolved organic molecule centres for ParaView |
| `analysis/production/summary.json` | numerical summary and metric definitions |
| `analysis/nvt/`, `analysis/npt/` | the same composition analysis for early relaxation |
| `analysis/production_thermodynamics.xvg` | temperature, pressure, density and potential energy |
| `analysis/production_energy.log` | thermodynamic summary reported by GROMACS |
| `production.log` | MD performance and numerical diagnostics |
| `code/` | exact preparation, analysis and launcher scripts used for this run |

Analysis samples one frame every 100 ps and uses molecular centers, so an organic
molecule is counted once just like a water molecule. Atom fractions would bias
this comparison because the organic molecules contain more atoms than water.
Look for growing, persistent composition domains together with increasing
long-wavelength fluctuations. Transient small clusters or one uneven snapshot
do not establish phase coexistence. The analyzer intentionally does not produce
an automatic binary "phase separated" verdict.

For VMD, download `visualization.gro`, `production.xtc`, and `view_binary.vmd`
into the same local folder. In VMD's Tk Console, source the `.vmd` file using
its local path, or load the GRO and then add the XTC to that same molecule.
Organic molecules are shown in orange and water in blue. Inspecting only the
organic molecules can make local clusters easier to see, but compare with the
quantitative composition outputs before interpreting them.

## Spatial analysis of an existing trajectory

For runs created with this v2 overlay, spatial outputs are generated automatically
for NVT, NPT and production. For a run already started or completed with the
previous overlay, use this separate analysis job. Its input is the **existing
binary result directory**, not a PLGA PoC directory:

```bash
sbatch /home/lee6456/MD/postprocess_binary.sbatch /absolute/path/to/binary_run
```

The script prints a new `analysis_spatial_production_TIMESTAMP_JOBID` directory
inside that run. It reads the existing TPR/XTC; it neither reruns MD nor uses the
old analyzer saved in the run's code snapshot. To restrict the time window or
change the sampling/grid, append options, for example:

```bash
sbatch /home/lee6456/MD/postprocess_binary.sbatch /absolute/path/to/binary_run --begin-ps 11100 --end-ps 21100 --stride-ps 100 --grid 4
```

Times are absolute GROMACS trajectory times: with the default stage offsets,
production spans 1100 to 21100 ps. An omitted time window analyzes all available
frames, including an unfinished production trajectory if it is readable.

### Dispersion metric

For each of the `grid^3` cells, let n_c be its molecular count and x_c its organic
mole fraction; x is the global organic mole fraction and N the total molecules.

    V = sum_c n_c (x_c - x)^2 / N
    V_random = x(1-x) (occupied_cells - 1) / (N - 1)
    segregation_excess = (V - V_random) / (x(1-x) - V_random)
    mixing_index_random_corrected = 1 - segregation_excess

The random reference holds total molecule counts fixed. Random spatial mixing
has mean segregation_excess near 0 and mixing_index near 1; component-pure
occupied cells give segregation_excess 1 and mixing_index 0 when the denominator
is nonzero. Negative segregation values (mixing index above 1) are retained:
they indicate fluctuations below the random-label reference, not a physical
mixing fraction above 100%. If every occupied cell contains a single molecule,
the normalization cannot distinguish random from segregated labels and the index
is undefined. No arbitrary phase-separation threshold is imposed.

The default 4^3 grid is intentionally coarse (1.5 nm cells for an initial 6 nm
cube). With only about 220-325 organic molecules, a very fine grid mostly measures
empty cells and counting noise. Compare time evolution at the same grid and
check nearby resolutions when drawing a physical conclusion. The Fourier
concentration metrics remain useful because they do not change under a common
periodic translation of a moving domain.

### Two-dimensional maps

Open `projection_last.svg` or `projection_time_average.svg` in a browser.
They show XY, XZ and YZ **projections through the full box**, not midplane slices.
Each projected mole fraction is the organic count divided by the total molecular
count in that column. Colour is **local organic mole fraction / global organic
mole fraction**: 1 is the bulk composition, 2 is twice the bulk composition.
It is not a volume fraction. The common colour scale is 0-3; values above 3
remain in numerical fields and are shown with the saturated upper colour.
Time averages are at fractional box coordinates and use the time-average box lengths
for the displayed physical axes. Inspect individual frames too: averaging in
the fixed laboratory frame can conceal moving domains.

### Three-dimensional fields and molecules in ParaView

Copy each `.pvd` **together with its sibling `vtk/` directory** to your desktop.
The PVD files reference the individual time frames using relative paths.

1. Open `concentration_3d.pvd` in ParaView and click Apply.
2. Select the organic mole fraction or relative enrichment cell-data array.
   Use Slice for a cross-section, or Threshold to inspect enriched cells.
3. Open `organic_molecules_3d.pvd` as well; use Points/Point Gaussian or a Sphere
   Glyph to see organic molecular centres. These are not atomic van der Waals
   surfaces. Play the animation to inspect persistence and domain movement.

Empty cells carry a validity mask and undefined composition, rather than being
silently labelled pure water. Position units are nm; PVD time units are ns.
All solvent species are counted by molecules, so differences in atom count do
not bias the concentration comparison. No additional Python packages beyond
the existing NumPy preparation environment are required on Negishi.

Completion is printed as `BINARY RUN COMPLETE` only after every MD stage has
reached its target and analysis has finished. A successful `mdrun` return alone
is not considered completion: its checkpoint step is checked. If wall time is
exhausted, the script prints `NOT COMPLETE` and retains the checkpoint. Only
then, to continue that same calculation, use the exact command printed in the log:

```bash
sbatch --chdir=/home/lee6456/MD /home/lee6456/MD/run_binary.sbatch --continue /absolute/path/to/binary_run
```

## Sources

- Sabura et al., *Impact of Formulation Ingredient Compatibility on Encapsulation
  of Hydrophobic Substances in Polymeric Nanocarriers Prepared by Nanoprecipitation*,
  sections 1, 2.1 and 4.3. https://doi.org/10.1002/mame.202500269
- Idrissi et al. (2015), water-DMSO mixing thermodynamics and local structure.
  https://doi.org/10.1039/C4CP04839C
- Honti, Idrissi and Jedlovszky (2021), DMF-water free energy of mixing and model
  sensitivity. https://doi.org/10.1021/acs.jpcb.1c01749
- de Jesus-Gonzalez et al. (2018), acetone-water miscibility and solvent model
  dependence. https://doi.org/10.1021/acs.jced.7b00573
- Taha and Lee (2013), ambient complete miscibility of acetone and acetonitrile
  with water, and demixing induced by adding a third component (TES buffer).
  https://doi.org/10.1063/1.4809995
- GROMACS 2026.4 MD parameters:
  https://manual.gromacs.org/documentation/2026.4/user-guide/mdp-options.html

## Verification scope

The delivery is checked locally for Python and shell syntax, molecule-count and
configuration logic, actual RDKit structure generation for all four solvents,
actual full-size periodic Packmol packing and contact checks for all four,
trajectory-analysis behavior on synthetic reference data, and launcher control
flow with simulated command responses. These checks do not substitute for
AmberTools parameterization or a GROMACS simulation on Negishi. The first actual
simulation remains to be submitted there.
