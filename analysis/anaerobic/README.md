# Anaerobic revision analyses

Simulations that address the reviewer comments on the anaerobic study, on yeast-GEM 9.1.2 (develop with #423 and #426 merged). Each folder answers one comment of the letter to the reviewers and holds its scripts, a README and the summary tables in `results/`; `paper_figures/` reruns everything and builds the figures and Supplementary Tables. The raw samples are not stored: the working folder `ANAEROBIC_RESULTS` holds them and `paper_figures/collect_results.sh` copies the summary tables into the comment folders.

## Setup

| Variable | Content |
| --- | --- |
| `RAVEN_ROOT` | RAVEN 3.0.0-rc1 or later |
| `GUROBI_HOME` | Gurobi installation (its `matlab` folder is added to the path) |
| `YEAST_GEM_ROOT` | yeast-GEM checkout with the model and `code/` to use (default: this repository) |
| `POLYTOPESAMPLER_ROOT` | [PolytopeSamplerMatlab](https://github.com/ConstrainedSamplingGroup/PolytopeSamplerMatlab), for the RHMC sampling |
| `ANAEROBIC_RESULTS` | working folder for all result files, including the raw samples (default: `results/` here) |
| `MATLAB_PREFDIR` | separate MATLAB preference folder, so that RAVEN's `setpref` and `savepath` leave the user's settings alone; one folder for each job when jobs run at the same time |

Start MATLAB with `addpath('analysis/anaerobic/common'); setupStudy('3.0.0-rc1')`, which puts all folders on the path.

A model is named `<release>_<constraints>` (`common/variantModel.m`), for example `vdev_new`: the yeast-GEM release tag (`dev`: the model of `YEAST_GEM_ROOT`, here 9.1.2) and the anaerobic constraints, `new` for `applyEnvironment(model,'anaerobic')` and `old` for the constraints before 9.1.0 (for 9.1.2 plus thiamine uptake, which the model needs without oxygen). A `+<modifier>` suffix undoes or adds one constraint.

## Folders

| Folder | Comment | Content |
| --- | --- | --- |
| `r1_major1_data_use` | same data for curation and validation | no simulation; points to the results below |
| `r1_major2_constraints` | state of the previous anaerobic function | Supplementary Table SX; glycerol with each constraint undone comes from `r1_major3_attribution` |
| `r1_major3_attribution` | which change improves the predictions | `reproVariants`: 2×2 of model version and constraints, each constraint undone (`function_ablation`), two undone (`function_pairs`), constraints opened (`sectionB`). `curationAblation`: the 9.0.2 to 9.1.0 curation replayed section by section and applied on 9.0.2 together with the later changes; one part, one of the four categories of the comment, or two parts left out (`categories`, `pairs`). `biomassDR`, `biomassTermVariants` |
| `r1_major4_gam_protein` | GAM, NGAM and protein content | no simulation; the `oldGAM` row of `function_ablation` |
| `r1_major5_degeneracy` | pFBA degeneracy | `fvaRun` (13 scenarios, optimal face and 1% band), `runSampling` (RHMC, four chains), `analyzeSampling`, `compareSensitivity`, `chain_diagnostics.py`, `validateRHMC` |
| `r1_major6_gcy1` | Gcy1 | `gcy1_thermodynamics.py` |
| `r1_major7_reproducibility` | reproducibility and backward compatibility | no simulation; the parity test in `code/python/tests/reference` |
| `r2_major1_deletions` | deletion strains | `knockoutPanel`, `adh3Gpd2`, `engineeredPanel`, `gpdViability`, `gpdSinkChain`, `frdKnockouts`, `geneDeletionScreen` (single-gene essentiality), `deletion_comparison.py` (Supplementary Table S11 and Figure R2: measured against predicted growth and yields of all deletion strains with data) |
| `r2_major2_fig3` | what Figure 3C-D compares | `flux_comparison.py` |
| `r2_major3_protons` | proton diffusion and the plasma-membrane ATPase | `protonBudget` (the FVA of the proton fluxes is in `r1_major5_degeneracy`) |
| `paper_figures` | | `run_all.m`, `sample_array.sh`, `makeFigures.m`, `make_supplementary.py`, `combine_ablation.py`, `collect_results.sh` |
| `common` | | path set-up, model variants, and the v9.1.0 benchmark code (`growthNewSim`, `growthOldSim`, `fluxPredSim`, `sjobergSim`) |
| `data` | | deletion-strain and engineered-strain measurements, with their sources |

## Run order

1. `sbatch paper_figures/sample_array.sh` (148 chains: 4 models × 13 scenarios, then 3 more chains of 8 scenarios), and `sbatch --array=<tasks> paper_figures/sample_array.sh 0.99` for the 1% band (the Sjöberg, Bakker wild type, Jouhten 2008 and Wasylenko 2014 scenarios of each model).
2. `run_all()` in MATLAB. The curation ablation is split over jobs: `curationAblation('curation_ablation','categories',k)` for k = 1 to 38 and `curationAblation('curation_pairs','pairs',k:k+4)` for k = 1, 6, ..., 196; `python paper_figures/combine_ablation.py` joins them.
3. `python r1_major5_degeneracy/chain_diagnostics.py`, `python r2_major2_fig3/flux_comparison.py`, `python r1_major6_gcy1/gcy1_thermodynamics.py`, `python paper_figures/make_supplementary.py`.
4. `bash paper_figures/collect_results.sh`.

## Not yet in this folder

- Supplementary Table S12 (CRPS by data set).
