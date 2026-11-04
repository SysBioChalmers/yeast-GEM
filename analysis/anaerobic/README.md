# Anaerobic revision analyses

Simulations that address the reviewer comments on the anaerobic study. Each folder answers one comment (see the letter to the reviewers); `paper_figures/` reruns everything on one model and writes the tables and figures of the revision. No results are stored in the repository: `results/` is created by the scripts and is ignored by git.

## Setup

| Variable | Content |
| --- | --- |
| `RAVEN_ROOT` | RAVEN 3.0.0-rc1 or later |
| `GUROBI_HOME` | Gurobi installation (its `matlab` folder is added to the path) |
| `YEAST_GEM_ROOT` | yeast-GEM checkout with the model and `code/` to use (default: this repository) |
| `POLYTOPESAMPLER_ROOT` | [PolytopeSamplerMatlab](https://github.com/ConstrainedSamplingGroup/PolytopeSamplerMatlab), for the RHMC sampling |
| `MATLAB_PREFDIR` | separate MATLAB preference folder, so that RAVEN's `setpref` and `savepath` leave the user's settings alone |

Start MATLAB with `addpath('analysis/anaerobic/common')`; every script then calls `setupStudy`, which sets up the remaining paths. RAVEN 3.0.0-rc1 `sampleChebyshevCenter.m` needs a sparse constraint matrix and Gurobi to sample the full model; this patch is not included.

A model is named `<release>_<constraints>` (`common/variantModel.m`), for example `v910_new`: the yeast-GEM release tag and the anaerobic constraints, `new` for `applyEnvironment(model,'anaerobic')` and `old` for the constraints before 9.1.0. A `+<modifier>` suffix undoes or adds one constraint.

## Folders

| Folder | Comment | Content |
| --- | --- | --- |
| `r1_major3_attribution` | which change improves the predictions | `reproVariants` (2×2 of model version and constraints; each component of the constraints undone; constraints opened), `curationAblation` (9.0.2 to 9.1.0 curation replayed section by section, with leave-one-out), `biomassDR`, `biomassTermVariants` (biomass reducing-equivalent term) |
| `r1_major5_degeneracy` | pFBA degeneracy | `fvaRun` (13 scenarios, optimal face and 1% band), `runSampling` (RHMC), `analyzeSampling`, `compareSensitivity`, `chain_diagnostics.py`, `validateRHMC` |
| `r2_major1_deletions` | validation with deletion strains | `knockoutPanel`, `adh3Gpd2`, `engineeredPanel`, `gpdViability`, `gpdSinkChain`, `frdKnockouts`, `geneDeletionScreen` (single-gene essentiality) |
| `r2_major3_protons` | proton diffusion and the plasma-membrane ATPase | `protonBudget` (the FVA of the proton fluxes is in `r1_major5_degeneracy`) |
| `paper_figures` | | `run_all.m` (all analyses above), `sample_array.sh` (Slurm array for the sampling), `makeFigures.m`, `make_supplementary.py` (Supplementary Tables SX, SY, SZ1 to SZ3 and S10) |
| `common` | | path set-up, model variants, and the v9.1.0 benchmark code (`growthNewSim`, `growthOldSim`, `fluxPredSim`, `sjobergSim`) |
| `data` | | deletion-strain and engineered-strain measurements, with their sources |

## Run order

1. `sbatch paper_figures/sample_array.sh` (52 jobs: 4 models × 13 scenarios).
2. `run_all()` in MATLAB.
3. `python paper_figures/make_supplementary.py`.

## Not yet in this folder

- The statistics of the flux comparison in Figure 3C-D (R2 major 2) and the Supplementary Tables S11 (deletion strains) and S12 (CRPS by data set).
- The Gcy1 thermodynamics (R1 major 6).
- `make_supplementary.py` contains the glycerol values of the single-constraint ablation (`effect`) as constants from an earlier model; they have to be read from the rerun instead.
