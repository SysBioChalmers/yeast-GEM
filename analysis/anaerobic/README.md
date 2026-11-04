# Anaerobic revision analyses

Scripts, input data and result tables for the simulations that address the reviewer comments on the anaerobic study. The results in `results/` were computed with the model before the open curation pull requests (#423, #426) and have to be regenerated with the final model.

## Setup

| Variable | Content |
| --- | --- |
| `RAVEN_ROOT` | RAVEN 3.0.0-rc1 or later |
| `GUROBI_HOME` | Gurobi installation (`matlab` folder is added to the path) |
| `YEAST_GEM_ROOT` | root of this repository (default: derived from the script location) |
| `POLYTOPESAMPLER_ROOT` | [PolytopeSamplerMatlab](https://github.com/ConstrainedSamplingGroup/PolytopeSamplerMatlab), for the RHMC sampling |
| `MATLAB_PREFDIR` | separate MATLAB preference folder, so that `savepath` does not change the user's path |

`matlab/setupStudy.m` sets the paths. RAVEN 3.0.0-rc1 `sampleChebyshevCenter.m` needs a sparse constraint matrix and a Gurobi call to run the sampling on the full model; this patch is not included.

## Reviewer comments, analyses and outputs

| Comment | Analysis | Scripts | Output in `results/` |
| --- | --- | --- | --- |
| Predictions are degenerate (R1 major 5) | flux variability analysis of 13 scenarios on the optimal face and in a 1% band; RHMC and CHRR flux sampling | `fvaRun.m`, `runSampling.m`, `rhmcSample.m`, `analyzeSampling.m`, `compareSensitivity.m`, `validateRHMC.m`, `validateCHRR.m` | `fva`, `analysis_*`, `rhmc_validation`, `chrr_validation` |
| Validation with deletion strains | anaerobic batch predictions of deletion strains against `data/knockout_data*.tsv`; essentiality screen | `knockoutPanel.m`, `adh3Gpd2.m`, `geneDeletionScreen.m`, `gpdViability.m`, `frdKnockouts.m` | `ko_panel.tsv`, `adh3_gpd2.tsv`, `ko_screen`, `gpd_viability*.tsv` |
| Engineered strains | yields of strains with altered redox or energy metabolism | `engineeredPanel.m` | `engineered_panel.tsv` |
| Growth and 13C flux benchmarks | predictions against the Sjöberg and four 13C data sets | `growthNewSim.m`, `growthOldSim.m`, `fluxPredSim.m`, `sjobergSim.m`, `reproStep1.m`, `reproVariants.m` | `repro_*` |
| Effect of the individual curations (R1 major 3) | (1) model version (9.0.2, 9.1.0) by anaerobic function (original, updated); (2) each component of the anaerobic function undone on 9.1.0; (3) the 9.0.2 to 9.1.0 curation replayed section by section, with leave-one-out | `variantModel.m`, `curateAnaerobic.m`, `curationAblation.m` | `function_ablation`, `curation_ablation*.tsv`, `supplementary` |

`py/` contains the plotting and table scripts (`make_supplementary.py`, `fig_attribution.py`, biomass and chain diagnostics). `data/` holds the deletion-strain and engineered-strain measurements with their sources.
