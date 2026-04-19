# Repository for the Computational Medicine Project: In Silico Analysis of Sex Differences in Cardiomyocyte Repolarisation Response to Groundwater Arsenic Exposure

This repository extends the **Sex-Specific Human Electromechanics** codebase (ToR-ORd-Land male/female models) by Holmes et al. (2025) to study how arsenic-induced IKr / IKs block produces sex-differentiated APD prolongation, and to compare the model's predictions against the clinical QTc data of Mumford et al. (2007) and Chen et al. (2013).

The base repository and its male/female ToR-ORd-Land models are described in:

> Holmes M, Wang ZJ, Doste R, et al. *Sex-specific human electromechanical multiscale in-silico models for virtual therapy evaluation.* J Mol Cell Cardiol Plus 2025; 13: 100479. https://doi.org/10.1016/j.jmccpl.2025.100479

All files listed below are additions specific to this project. The original model files ([model_ToRORd_Land_Male.m](model_ToRORd_Land_Male.m), [model_ToRORd_Land_Female.m](model_ToRORd_Land_Female.m), [modelRunner.m](modelRunner.m), etc.) are retained unchanged from the base repo.

## New scripts

### Baseline verification

- **[scriptVerification.m](scriptVerification.m)** — Runs the unmodified male and female ToR-ORd-Land models for endo / epi / mid cell types at BCL = 1000 ms and compares baseline APD90 against Holmes et al. (2025). Sanity check that the sex-specific models reproduce the published values before any drug perturbation is applied.

### Concentration sweeps (produce `.mat` data consumed by the validation / ablation scripts)

- **[scriptSweep_SingleCell_Parallel.m](scriptSweep_SingleCell_Parallel.m)** — Single male and single female ToR-ORd-Land cell per cell type (endo / epi / mid), simulated across 20 log-spaced arsenic concentrations (IKr and IKs blocked via Hill equation, IC50 from Drolet et al. 2004). Writes `singlecell_sweep.mat`.
- **[scriptSweep_Population_Parallel.m](scriptSweep_Population_Parallel.m)** — Population-of-models version: 25 male + 25 female virtual cells per cell type with variability in 11 ionic conductances, same concentration grid. Writes `population_sweep.mat`.
- **[scriptSensitivity_Parallel.m](scriptSensitivity_Parallel.m)** — One-at-a-time perturbation of the Hill coefficient (h ∈ {0.7, 1.0, 1.5}) and the IKr / IKs IC50s (±1 SD around the Drolet values) to check whether the qualitative sex difference survives calibration uncertainty.

All three sweep scripts are parallelised with `parfor` over the flattened (scenario, concentration, cell-type) grid. The `.mat` outputs are not committed — regenerate them locally by running the sweep scripts.

### Validation / analysis (consume the sweep `.mat` files; produce the paper figures)

- **[scriptValidationSingleCell.m](scriptValidationSingleCell.m)** — Loads `singlecell_sweep.mat` and plots ΔAPD90 vs [ATO] and AP-trace panels (male vs female, per cell type) for the "Single-cell validation" section.
- **[scriptValidationPopulation.m](scriptValidationPopulation.m)** — Loads `population_sweep.mat` and produces the population-level figures: mean ± SD ΔAPD90 curves, boxplots at low / mid / high concentration, prevalence of ΔAPD90 > T ms for T ∈ {30, 40, 50}, and AP-trace panels at Bangladesh-realistic calibrated concentrations.
- **[scriptConcentrationMapping.m](scriptConcentrationMapping.m)** — Derives the water → cardiomyocyte free-arsenite mapping by matching model output against two clinical anchors (Chen 2013 female:male slope ratio of 4.3; Mumford 2007 prevalence of 3.9 / 11.1 / 20.6 % across exposure groups) and fits a linear-through-origin bridge `C_cell = k · C_water` per cell type.
- **[scriptVirtualAblation.m](scriptVirtualAblation.m)** — At three calibrated arsenic concentrations (Mumford low / mid / high), swaps each sex-difference conductance scaling factor from male to female one at a time and measures how much of the female-minus-male ΔAPD90 gap remains. Identifies which ion-channel sex difference drives the greater female vulnerability.

## Reproducing the figures

1. (Optional) Run `scriptVerification.m` to confirm the baseline male / female APD90 values match Holmes et al. (2025).
2. Run `scriptSweep_SingleCell_Parallel.m` and `scriptSweep_Population_Parallel.m` (requires the Parallel Computing Toolbox) to produce `singlecell_sweep.mat` and `population_sweep.mat`.
3. Run `scriptValidationSingleCell.m` and `scriptValidationPopulation.m` to produce the qualitative validation figures.
4. Run `scriptConcentrationMapping.m` to calibrate the water → cell concentration bridge used by the ablation script.
5. Run `scriptVirtualAblation.m` for the mechanistic decomposition.
6. Run `scriptSensitivity_Parallel.m` to verify robustness of the sex-difference prediction to Hill / IC50 uncertainty.
