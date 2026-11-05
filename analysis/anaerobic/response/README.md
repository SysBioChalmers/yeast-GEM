# Response to the reviewers: material for the co-authors

All simulations are for yeast-GEM 9.1.2 (develop with #423 and #426 merged); the folders of `analysis/anaerobic` hold the scripts and the summary tables.

| File | Content |
| --- | --- |
| `RESPONSE_LETTER.md` | draft reply to each comment of both reviewers, with the numbers of the rerun |
| `MANUSCRIPT_EDITS.md` | the proposed text of the passages of the manuscript that the replies promise to change |
| `OPEN_DECISIONS.md` | decisions for the co-authors and statements that are not yet checked |

## Comment, reply and evidence

| Comment | Reply in the letter | Evidence |
| --- | --- | --- |
| R1 major 1 (same data for curation and validation) | distinguishes Figure 1 (data used in curation) from Figure 3 and the deletion strains | `../r1_major1_data_use` |
| R1 major 2 (state of the old function) | Supplementary Table SX | `../r1_major2_constraints`, `../paper_figures/results/supplementary/SX_anaerobic_constraints.tsv`, SY |
| R1 major 3 (which change improves the predictions) | 2×2 of model and constraints, each constraint and each curation part left out, four categories, pairs | `../r1_major3_attribution` (`results/`), SZ1 to SZ5, Supplementary Figures S_b and S_c |
| R1 major 4 (GAM, NGAM, protein) | values are those of yeast-GEM, not fitted here | `../r1_major4_gam_protein`, SZ2 |
| R1 major 5 (pFBA degeneracy) | RHMC sampling in four chains, FVA, 99% band | `../r1_major5_degeneracy` (`results/`), Supplementary Figure S_a |
| R1 major 6 (Gcy1) | thermodynamic bound on the oxidative direction | `../r1_major6_gcy1` |
| R1 major 7 (reproducibility) | one data file for the constraints, MATLAB and Python give the same model | `../r1_major7_reproducibility` |
| R2 major 1 (deletion strains) | growth and yields of 56 strain comparisons, adh3 and mae1 against Grigaitis 2023 | `../r2_major1_deletions` (`results/`), Table S11, Supplementary Figures S_d and S_e |
| R2 major 2 (Figure 3C-D) | statistics of the compared fluxes | `../r2_major2_fig3`, Table S10 |
| R2 major 3 (protons, H+-ATPase) | proton diffusion from FVA, ATPase not determined | `../r2_major3_protons`, `../r1_major5_degeneracy/results/fva` |
| R1 minor, R2 minor | corrected in the text | `MANUSCRIPT_EDITS.md`, section "Minor" |

## Supplementary material

| Item | File |
| --- | --- |
| Table SX: constraints of the previous and updated function, with the glycerol when each is undone | `../paper_figures/results/supplementary/SX_anaerobic_constraints.tsv` |
| Table SY: opening glycine cleavage, glutamate synthase and Gcy1, against the 13C data | `.../SY_opening_constraints_13C.tsv` |
| Table SZ1: model version × anaerobic constraints | `.../SZ1_model_x_constraints.tsv` |
| Table SZ2: each constraint undone | `.../SZ2_constraints_undone.tsv` |
| Table SZ3: each curation part and category left out or added | `.../SZ3_curation_parts.tsv` |
| Table SZ4: two curation parts left out | `.../SZ4_curation_pairs.tsv` |
| Table SZ5: two constraints undone | `.../SZ5_constraints_pairs.tsv` |
| Table S10: every compared 13C flux, pFBA and sampled interval | `.../S10_fluxomics_points.tsv` |
| Table S11: deletion strains, measured against predicted | `../r2_major1_deletions/results/deletion_comparison.tsv` |
| Table S12: CRPS by data set | not built; the per-data-set values are in `../r1_major5_degeneracy/results/chains/crps_opt0.99999.tsv` |
| Figure S_a: sampled intervals and CRPS | `../paper_figures/results/figures/figR1_sampling.pdf` |
| Figure S_b: effect of leaving out each change | `.../figAttribution.pdf` |
| Figure S_c: categories, without and alone | `.../figCategories.pdf` |
| Figure S_d: deletion strains | `.../figDeletions.pdf` |
| Figure S_e: Bakker 2000 wild type and adh3 | `.../figR2_knockouts.pdf` |

## Main results

- The updated constraints on 9.1.2 predict glycerol within the measurement error (4.35 against 4.5 ± 0.4 mmol/gDW/h), growth R² 0.90 and a 13C error of 30%; the other three combinations of model version and constraints are worse on each (SZ1).
- The growth improvement comes from the general GAM, NGAM and protein content, not from the curation (SZ2).
- Glycerol depends on the mitochondrial and nitrogen changes and on the reaction directions together; charge and element balance, the biomass equation and the changes after 9.1.0 do not change the predictions (SZ3, SZ4, Supplementary Figure S_c).
- The fluxes compared with the 13C data are nearly determined in 9.1.2 (95% intervals below 3.2% of the glucose uptake) and not in 9.0.2 (about 90%); the ranking of the models is the same in four chains and in the 99% band.
- The growth class of deletion strains is predicted for 10 of 18 strain and study combinations; the model does not predict the lack of growth of the gpd1 gpd2 mutant.
