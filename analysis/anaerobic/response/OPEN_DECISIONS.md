# Open decisions and unchecked statements

## Decisions for the co-authors

| # | Decision | Options and numbers (9.1.2, updated constraints) |
| --- | --- | --- |
| 1 | NADPH term in the biomass reaction (degree of reduction 4.37 per C-mol) | Keep it, with a justification, or remove it (degree of reduction 4.22, experimental 4.2). Removing it changes glycerol from 4.35 to 4.37 and the 13C error from 0.31 to 0.28; the letter reports the effect and does not say the term was removed. `MANUSCRIPT_EDITS.md`, l.223-234, options A and B |
| 2 | FADH2 production in the biomass reaction | 0.019 mmol/gDW (upper limit from the cysteine content) or 0.08 (fitted to the succinate production). Without the term glycerol is 4.36 and the 13C error 0.32. `MANUSCRIPT_EDITS.md`, l.196-201, options A and B |
| 3 | Restriction of Mdh1 to the oxidative direction | The letter says it was removed for lack of evidence. Restoring it changes glycerol by +0.13 and the 13C error from 0.30 to 0.37 |
| 4 | Definition of the original constraints on 9.1.2 | The original constraints give no growth on 9.1.2, because #423 blocks the oxygen-dependent HMP-P synthase and only the updated constraints provide thiamine. The 2×2 therefore uses the original constraints plus thiamine uptake. Alternative: use the released 9.1.0 for that cell |
| 5 | gpd1 gpd2, gpd2, osm1 and ald6 ald4 ald5 | No curation is proposed: each route the model uses is a real reaction, and no solid evidence supports a reaction change. The letter reports them as limitations and states that the explanation (regulation, enzyme capacity) is not tested |
| 6 | Release naming | The model of the analyses is called 9.1.2 (develop with #423 and #426 merged); the letter and manuscript edits use this name |

## Statements that are not yet checked

| Statement | Where | Status |
| --- | --- | --- |
| Origin of NGAM 0.7 and release 8.3.1 for the GAM fit | R1 major 4 | The model files have GAM 55.3 and NGAM 0.7 in every release from v8.5.0 to 9.1.1 (earlier ones not read); `fitGAM.m` fits GAM to the data of van Hoek 1998 and not NGAM |
| Anaerobic function of 9.0.2 | R1 major 4 | The override (GAM 30.49, NGAM 0, protein 0.461) was read in v8.7.1 and v9.0.0; 9.0.2 is assumed to be the same |
| Origin of the original values GAM 30.49 and NGAM 0 | R1 major 4 | The code cites Nissen et al. 1997 and Jouhten et al. 2012; the papers were not read for this |
| Glycolysis and alcoholic fermentation are proton-neutral; the proton of ammonium is released on assimilation | R2 major 3 | Carried over from an earlier analysis; the proton budget of 9.1.2 (`../r2_major3_protons/results/protonBudget.txt`) was not compared with the sentence |
| Diffusion reactions: 0.9 mmol/gDW/h in the extreme of the 1% band | R2 major 3 | FVA values of 9.1.2 (`../r1_major5_degeneracy/results/fva`); the sampled values at 99% growth were not recomputed |
| Largest relative deviations of the small fluxes (acetate, acetyl-CoA synthetase, NAD-dependent isocitrate dehydrogenase, threonine aldolase) | R2 major 2 | Carried over; the statistics of `../r2_major2_fig3` give the mean errors and the three largest absolute deviations |
| Biomass yields of engineered strains within 5-10% | R2 major 1 | Not compared with `../r2_major1_deletions/results/engineered_panel.tsv` |
| Largest projected internal flux 43 mmol/gDW/h; same quantiles with internal bounds ±100 and ±400 | Methods, sampling | Measured in an earlier run; not repeated |
| Software versions in the Methods (RAVEN 2.11.1, MATLAB 2026a, Gurobi 13.0.2) | Methods | These analyses used RAVEN 3.0.0-rc1, MATLAB 2024b and Gurobi 13.0.1 |

## Not in the material

- Supplementary Table S12 (CRPS by data set).
- The R1 minor and R2 minor corrections are listed in `MANUSCRIPT_EDITS.md` but not applied to a manuscript file.
