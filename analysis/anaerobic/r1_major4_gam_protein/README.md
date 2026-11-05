# R1 major 4: GAM, NGAM and protein content

GAM, NGAM and protein content are not fitted in this study. Their effect on the anaerobic predictions is the `oldGAM` row of `function_ablation` (`../r1_major3_attribution/reproVariants.m`: the general values replaced by GAM 30.49, NGAM 0 and protein 0.461) and the chemostat growth errors of `repro`. The protein content of the nitrogen-limited conditions is set in `../common/growthNewSim.m` (environment variable `NLIM_PROTEIN`, default 0.28 g/gDW).
