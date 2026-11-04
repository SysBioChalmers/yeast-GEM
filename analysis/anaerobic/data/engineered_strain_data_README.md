# Anaerobic physiology of engineered *S. cerevisiae* strains — companion notes

Companion to `engineered_strain_data.tsv`. That file holds one row per
strain x condition x quantity, with values transcribed exactly as printed
(`<`, `N.D.`, `>99` kept verbatim). **Only anaerobic cultures are included**;
aerobic data from the same papers are described here but deliberately left out of
the TSV.

Sources (local PDFs and text extractions in `../papers/`):

| key | paper | DOI |
| --- | --- | --- |
| Kozak2014 | Kozak et al. 2014, Metab Eng 21:46–59 | 10.1016/j.ymben.2013.11.005 |
| Nissen2000 | Nissen et al. 2000, Metab Eng 2:69–77 | 10.1006/mben.1999.0140 |
| Bro2006 | Bro et al. 2006, Metab Eng 8:102–111 | 10.1016/j.ymben.2005.09.007 |
| deKok2011 | de Kok et al. 2011, Metab Eng 13:518–526 | 10.1016/j.ymben.2011.06.001 |
| Basso2011 | Basso et al. 2011, Metab Eng 13:694–703 | 10.1016/j.ymben.2011.09.005 |

## Extraction caveats that apply to all five papers

* `pdftoppm`/`pdftotext` and every Python PDF library are unavailable on this
  node, so the PDFs could **not** be re-rendered to resolve table layout. All
  values were read from the pre-extracted `.txt` files. Table columns were
  matched to strains by the header order printed in the text; where a table was
  garbled enough to be ambiguous this is flagged below and in the row `notes`.
* The extraction loses Greek letters and superscripts: `μ` becomes `m`,
  `Δ` becomes `1`/`2`/`D`, `±` becomes `7`, `−1` becomes `/C0 1`, `×` becomes
  `/C2`. Values were reconstructed accordingly (e.g. de Kok's enzyme activities
  printed as "mmol mg protein-1 min-1" are µmol mg protein⁻¹ min⁻¹). Where the
  original could be either mmol or µmol the choice was checked against the
  measured flux (see Basso invertase note).
* No value in the TSV is a figure read; the `source_type` column is `table` or
  `text` throughout. Where the text quotes a slope taken from a figure
  (Nissen Fig. 1 growth rates) the printed number in the text was used, and the
  location says so.
* Deletion genotypes use systematic ORF names. All were checked; the only
  unresolved ones are `MAL21`/`MAL22` (see deKok2011 below).

---

## 1. Kozak et al. 2014 — *acs1Δ acs2Δ* with heterologous PFL or A-ALD

### Cultivation

* Strain background: CEN.PK. Reference strain **CEN.PK113-7D**
  (*MATa MAL2-8C SUC2*, prototrophic). Engineered strains derive from
  **CEN.PK102-12A** (*MATa MAL2-8C SUC2 ura3-53 leu2-3,112 his3-Δ1*) and carry
  2 µm URA3 plasmids, so they are *ura3* complemented by the plasmid but remain
  *leu2*/*his3* complemented by the deletion markers.
* Medium: synthetic medium of Verduyn et al. 1992, ammonium sulfate as nitrogen
  source, 30 °C, pH 5.0 (2 M KOH), Antifoam C 0.2 g L⁻¹.
* Anaerobiosis: 2-L Applikon bioreactors, 1 L working volume, 800 rpm, sparged
  with 500 mL min⁻¹ N₂ (<10 ppm O₂), Nonprene tubing, Viton O-rings, medium
  reservoirs continuously N₂-flushed. Anaerobic growth factors:
  **ergosterol 10 mg L⁻¹ + Tween 80 420 mg L⁻¹**. Anaerobic shake flasks were
  run in a Bactron chamber (5% H₂, 6% CO₂, 89% N₂).
* Carbon source: anaerobic **batch** bioreactors 25 g L⁻¹ glucose; anaerobic
  **chemostat** D = 0.10 h⁻¹. **Ambiguity:** Methods says anaerobic chemostat
  medium contained 20 g L⁻¹ glucose, while the Table 4 caption says 25 g L⁻¹.
  Either way the culture is glucose-limited, so this only matters if the
  feed concentration is simulated explicitly.
* Chemostat steady state: ≥5 volume changes, CO₂ rate stable within 4% over
  2 volume changes, samples between 10 and 15 volume changes.
* Ethanol corrected for evaporation; formate, glycerol and residual glucose by
  enzymatic kits after rapid quenching.

### Engineered strains in the TSV (anaerobic only)

**IMZ383 — *acs1Δ acs2Δ* + *E. coli* PFL** (anaerobic chemostat + batch)
**IMZ384 — *acs1Δ acs2Δ* + *L. plantarum* PFL** (anaerobic shake flask only)

Model changes needed:

* **Delete**: `ACS1` = YAL054C and `ACS2` = YLR153C, i.e. remove all
  acetyl-CoA synthetase activity — both the cytosolic (Acs1/Acs2) reaction
  `acetate[c] + ATP[c] + CoA[c] -> acetyl-CoA[c] + AMP[c] + PPi[c]` and, in
  models that also place ACS in the mitochondrion, the corresponding
  mitochondrial copy if it is assigned to the same genes. Note the paper's ATP
  accounting: with subsequent PPi hydrolysis this is equivalent to 2 ATP → 2 ADP
  + 2 Pi per acetate activated.
* **Add (IMZ383/IMZ384)**: pyruvate-formate lyase, EC 2.3.1.54, **cytosol**:
  `pyruvate[c] + CoA[c] <=> acetyl-CoA[c] + formate[c]`
  — no ATP, no redox cofactor, mass- and charge-balanced as written.
  Authors treat it as reversible in principle but in vivo it runs towards
  acetyl-CoA + formate. Genes: *pflB* (PFL) + *pflA* (PFL-activating enzyme,
  EC 1.97.1.4) from *Escherichia coli* (IMZ383, plasmid pUDE204) or
  *Lactobacillus plantarum* (IMZ384, pUDE214), codon-optimised, TPI1 and FBA1
  promoters. PFL-AE carries an O₂-sensitive [4Fe–4S] cluster; it has no
  stoichiometric role but explains why these strains are viable only
  anaerobically (add no aerobic capability). A GPR for the new reaction should
  be the heterologous genes, not a yeast ORF.
* **Transporters**: none changed. Formate is excreted; the model needs a
  formate export/exchange reaction. Native FDH1/FDH2 (formate dehydrogenase,
  `formate + NAD+ -> CO2 + NADH`) were transcriptionally up >25-fold but were
  *not* able to reoxidise all formate in vivo, so a simulation that lets FDH
  run freely will under-predict formate excretion.
* **Authors' own stoichiometric statements**: acetyl-CoA production via PFL is
  stoichiometrically coupled to formate production, so the measured
  2.5 ± 0.1 mmol formate (g DW)⁻¹ is a direct in vivo measure of the cytosolic
  acetyl-CoA demand for biosynthesis. Their bottom-up estimate of that demand is
  **1.16 mmol (g DW)⁻¹** = 1.04 (lipids + lysine, Flikweert et al. 1999) + 0.05
  (Met + Cys in protein) + 0.05 (glutathione) + 0.02 (N-acetylglucosamine for
  chitin) + ~2.5 × 10⁻³ (protein lysine acetylation, assuming 39% protein and
  56 kDa average protein). The measured formate exceeds this, i.e. PFL flux is
  more than sufficient for the known sinks. **This 1.16 mmol (g DW)⁻¹ number is
  a useful check on the cytosolic acetyl-CoA drain of the yeast-GEM biomass
  equation.**
* Despite an ATP-neutral pathway, the anaerobic biomass yield **fell** 15%
  (0.096 → 0.082 g g⁻¹). The authors attribute this to weak-acid uncoupling by
  the excreted formate and acetate, not to pathway stoichiometry — a model will
  therefore *over*-predict IMZ383's yield unless extra maintenance ATP or an
  uncoupling term is added.

### Deliberately excluded

* **IMZ305** (*ald2-6Δ acs1Δ acs2Δ* + *E. coli eutE*, acetylating acetaldehyde
  dehydrogenase) was characterised only in **aerobic** glucose-limited
  chemostats (Table 6: Y_biomass 0.501 reference vs 0.429 IMZ305 g g⁻¹) and in
  aerobic shake flasks. Its A-ALD reaction, if needed later, is
  `acetaldehyde[c] + NAD+[c] + CoA[c] <=> acetyl-CoA[c] + NADH[c] + H+[c]`
  (EC 1.2.1.10, *E. coli eutE*), combined with deletion of `ALD2` (YMR170C),
  `ALD3` (YMR169C), `ALD4` (YOR374W), `ALD5` (YER073W), `ALD6` (YPL061W),
  `ACS1`, `ACS2`.
* A-ALD screening growth rates and enzyme activities (Table 3) and all
  intracellular metabolite data (Table 7) are aerobic.
* The anaerobic shake-flask growth rate of the Acs⁺ reference **IME140** is not
  printed; only the ratios are (IMZ383 = 73%, IMZ384 = 53% of IME140). No
  value was derived for the TSV.

---

## 2. Nissen et al. 2000 — *gdh1Δ* with *GLT1*/*GLN1* or *GDH2* over-expression

### Cultivation

* Strain background: haploid derivative of the diploid industrial model strain
  **CBS8066**. TN1 = *MATα ho-2Δ* wild type; TN2 = TN1 with *ura3* replaced by
  *SUC2*; all engineered strains descend from TN2. (Genotype strings are
  reproduced with the extraction's garbled Δ symbols flagged; the `ura3-2Δ::SUC2`
  and `gdh1-2Δ::URA3` readings are my reconstruction of `ura3-220::SUC2` /
  `gdh1-21::URA3` in the text — **uncertain in detail, unambiguous in
  substance**: *URA3* is disrupted, *GDH1* is replaced by *URA3*.)
* Medium: minimal medium of Verduyn et al. 1990, vitamins sterile-filtered.
  Initial **glucose 25 g L⁻¹**, **(NH₄)₂SO₄ 3.75 g L⁻¹** — ammonium is the sole
  nitrogen source, which is the whole point of the strategy. Ergosterol
  4.2 mg (g DW)⁻¹ and Tween 80 175 mg (g DW)⁻¹ (note: per gram dry weight, not
  per litre, as printed), dissolved in hot 96% ethanol. Antifoam 75 µL L⁻¹.
* Anaerobic **batch** only (no chemostats in this paper), 4.5 L working volume,
  30 °C, 600 rpm, pH 5.00 (2 M KOH), sparged with N₂ containing <5 ppm O₂
  (technical N₂ passed over copper at 400 °C) at 0.50 L N₂ min⁻¹ L⁻¹, Norprene
  tubing, off-gas condensers at 2 °C. Inoculated to 1 mg L⁻¹. Each cultivation
  performed three times "with identical results"; no SDs are printed except a
  stated 3–4% SD on the ethanol yield.

### Engineered strains and the model change each needs

Wild type ammonium assimilation (Gdh1, NADPH-dependent, `GDH1` = **YOR375C**):
`2-oxoglutarate + NH4+ + NADPH -> glutamate + NADP+`.

* **TN9 (*gdh1Δ*)** — delete `GDH1` (YOR375C) only. Flux must then route via
  the GS/GOGAT couple already in the model (`GLN1` = YPR035W, `GLT1` = YDL171C)
  or via Gdh2 (`GDH2` = YDL215C). Do **not** also delete `GDH3` (YAL062W):
  the authors state Gdh3p "does not function in *gdh1* strains", so in a model
  Gdh3 should be constrained to zero when reproducing these data, otherwise the
  NADPH route survives the knockout and the phenotype disappears. **This is the
  single most important modelling caveat for this paper.**
* **TN15 (*gdh1Δ* PGK1p-GLN1)** — over-express `GLN1`
  (`glutamate + NH4+ + ATP -> glutamine + ADP + Pi`, cytosol).
* **TN17 (*gdh1Δ* PGK1p-GLT1)** — over-express `GLT1`
  (`2-oxoglutarate + glutamine + NADH -> 2 glutamate + NAD+`, cytosol).
* **TN19 (*gdh1Δ* PGK1p-GLT1 PGK1p-GLN1)** — the target strain. Net reaction
  replacing Gdh1:
  `2-oxoglutarate + NH4+ + NADH + ATP -> glutamate + NAD+ + ADP + Pi`.
  So relative to wild type: −1 NADPH, +1 NADH oxidised, +1 ATP hydrolysed per
  glutamate. Nothing is added to the network — all three genes are native; only
  the promoters change (chromosomal `PGK1p` insertions, not plasmids), so in a
  model this is a *gdh1Δ* deletion plus removal of any flux bound on GS/GOGAT.
* **TN22 (*gdh1Δ* PGK1p-GDH2)** — over-express `GDH2`
  (`2-oxoglutarate + NH4+ + NADH -> glutamate + NAD+`, no ATP). In the wild type
  Gdh2 runs catabolically (glutamate → 2-oxoglutarate + NH4+ + NADH); to
  reproduce TN22 its direction must be freed / reversed.
* **Transporters**: unchanged.

### Authors' quantitative predictions vs measurement

* Prediction (from an earlier flux analysis at D = 0.3 h⁻¹, Nissen et al. 1997):
  the Gdh1 flux is 8.0 C-mmol (g biomass)⁻¹ h⁻¹ although the net
  2-oxoglutarate → glutamate flux is only 2.0. Replacing Gdh1 by Glt1+Gln1
  should cut glycerol formation by **49%**, from 9.8 to 5.0 C-mmol
  (g biomass)⁻¹ h⁻¹, and raise ethanol formation by **6%**, from 54.9 to
  58.2 C-mmol (g biomass)⁻¹ h⁻¹.
* Measured (anaerobic batch, C-mol C-mol⁻¹ glucose): glycerol yield fell 38%
  (0.097 → 0.060) and ethanol yield rose 8–10% (0.480 → 0.520 in TN9,
  0.531 in TN19). So the measured glycerol reduction was smaller and the
  ethanol gain larger than predicted.
* µmax: 0.41 h⁻¹ (TN1) → 0.21 (TN9) → 0.24 (TN15) → 0.31 (TN17) → 0.37 (TN19)
  → 0.39 (TN22).

### Internal inconsistencies flagged in the TSV `notes`

* TN9 µmax is 0.21 h⁻¹ in Results/Fig. 1 but "0.22" in the Discussion.
* Gdh2p activity in TN22 is "30-fold higher than in the wild type" in Results
  and "10-fold higher" in the Discussion; the tabulated ratio is
  0.625/0.020 ≈ 31.
* The abstract credits TN19 with "10% higher ethanol yield and 38% lower
  glycerol yield"; the Results credit TN9 with +8%/−38%. TN17 and TN19 are
  indistinguishable from each other.

---

## 3. Bro et al. 2006 — GAPN expression (and the *gpd* comparison)

### Cultivation

* **Glucose batch strains**: background **M4054** (*MATa ura3 gap1*; *GAP1* =
  YKR039W, general amino acid permease). Reference = M4054 + empty pYX212;
  gapN strain = **CPB1001** = M4054 + pYX212-*gapN*.
* **Glucose/xylose chemostat strains**: background **JSECPB1**
  (*MATa SUC2 MAL2-8 pADH-XYL1 pPGK-XYL2 pPGK-XKS1 ura3⁻*, CEN.PK background).
  Xyl-reference = CPBCB5 (+ empty pYX212); Xyl-gapN = CPBCB4 (+ pCPBGAPN1).
  Note the two experiments use **different strain backgrounds**.
* Medium (per litre): 40 g glucose (batch) or 20 g glucose + 50 g xylose
  (chemostat); 5.0 g (NH₄)₂SO₄ (ammonium as N source); 3.0 g KH₂PO₄; 0.5 g
  MgSO₄·7H₂O; Verduyn 1992 trace metals and vitamins; **Tween 80 420 mg L⁻¹ +
  ergosterol 10 mg L⁻¹** dissolved in boiling ethanol; antifoam.
* Anaerobiosis: 2 L working volume (batch) / 1.0 L (chemostat), 30 °C,
  600 rpm (batch) / 500 rpm (chemostat), N₂ at 400 mL min⁻¹ (batch) /
  0.2 L min⁻¹ (chemostat), Norprene tubing, dissolved O₂ below detection limit,
  pH 5.0 (4 M KOH batch / 2 M NaOH chemostat).
* Chemostat **D = 0.05 h⁻¹**; steady state after ≥5 residence times, biomass,
  metabolites and CO₂ constant within ±2%.
* Batch yields are slopes of metabolite vs glucose concentration during
  exponential growth (so they are differential yields, not end-point yields).

### Model change

* **Add** *Streptococcus mutans* **gapN** — non-phosphorylating,
  NADP⁺-dependent glyceraldehyde-3-phosphate dehydrogenase, EC 1.2.1.9,
  **cytosol**, treated as **irreversible**:
  `glyceraldehyde-3-phosphate[c] + NADP+[c] + H2O[c] -> 3-phosphoglycerate[c] + NADPH[c] + 2 H+[c]`
  (the paper writes it without water as
  `D-glyceraldehyde 3-phosphate + NADP+ -> 3-phospho-D-glycerate + NADPH`;
  water must be added to balance it in a model). Expressed from the TPI1
  promoter on the 2 µm URA3 plasmid pYX212. GPR = heterologous *gapN*.
* **Nothing is deleted.** The native route
  (Tdh1/2/3 `GAP + NAD+ + Pi -> 1,3-bisphosphoglycerate + NADH` followed by
  Pgk1 `+ ADP -> 3PG + ATP`) stays intact and carries most of the flux; GAPN
  bypasses it, giving NADPH instead of NADH + ATP. **Per unit of bypassed flux:
  −1 ATP, −1 NADH, +1 NADPH.**
* **Transporters**: unchanged. For the xylose strains the model additionally
  needs xylose reductase (both `xylose + NADPH -> xylitol + NADP+` and
  `xylose + NADH -> xylitol + NAD+`, *P. stipitis* XYL1), xylitol dehydrogenase
  (`xylitol + NAD+ -> xylulose + NADH`, XYL2), xylitol export, and
  over-expressed native XKS1.
* Nitrogen source is ammonium, and the authors stress that the GAPN strategy —
  unlike the Nissen ammonium-assimilation strategy — does **not** depend on the
  nitrogen source or on growth.

### Authors' predictions vs measurement (all anaerobic, glucose)

* Genome-scale model (Förster et al. 2003) at µ = 0.27 h⁻¹, wild type predicted:
  0.40 g ethanol g⁻¹ glucose, 0.09 g glycerol g⁻¹, 0.10 g DW g⁻¹ — vs measured
  0.41 / 0.081 / 0.098.
* GAPN insertion predicted: ethanol **+10.4%**, glycerol **−100%** (complete
  elimination), biomass **+7.6%**.
* Measured: glycerol −40% (0.081 → 0.046 g g⁻¹), ethanol +3%
  (0.41 → 0.42 g g⁻¹), biomass 0.098 → 0.101 g g⁻¹, µmax unchanged
  (0.27 → 0.28 h⁻¹).
* When the model is re-constrained to the *measured* glycerol production rate it
  predicts +0.02 g ethanol g⁻¹ and +0.005 g DW g⁻¹ — matching the experiment.
  The model says ~10% of the GAP → 3PG flux must go through GAPN for the full
  effect; in vitro GAPN was ~10% of total GAPDH activity (0.75 vs 6.8 U mg⁻¹)
  yet the in vivo flux share was evidently lower, which the authors ascribe to
  low intracellular NADP⁺/NADPH.
* Other strategies scored in Table 2 (predictions only, not in the TSV, all
  anaerobic glucose at µ = 0.27 h⁻¹, % change vs wild-type simulation —
  ethanol / glycerol / biomass): *gdh1Δ* + GLT1/GLN1 +5.2 / −49.3 / +5.2
  (experimental +8 / −38 / −1); *gdh1Δ* + GDH2 +4.2 / −46.8 / +9.7
  (experimental +2 / −31 / +14); NAD⁺-dependent glycine dehydrogenase
  +4.2 / −46.8 / +9.7; NADP⁺-dependent orotate reductase, transhydrogenase,
  NADH kinase, NADP⁺-dependent glycerol dehydrogenase, NADP⁺-dependent
  glycerol-3-phosphate dehydrogenase, NADP⁺-dependent (phosphorylating) GAPDH
  and GAPN all 9.0–10.4 / −100 / +7.6–16.5. The authors note the
  transhydrogenase prediction is experimentally wrong (runs the other way in
  vivo, Nissen et al. 2001) — a useful negative control for model validation.
* Xylose/glucose chemostat predictions (Table 4, in silico columns not in the
  TSV) span Y_biomass 0.071–0.100, Y_ethanol 0.302–0.454, Y_glycerol 0.000,
  Y_xylitol 0.232–0.303 g g⁻¹, D 0.054–0.075 h⁻¹, depending on which fluxes are
  constrained; measured D was fixed at 0.05 h⁻¹, so the in silico and in vivo
  columns are **not at the same growth rate** and must not be compared directly.

### On *gpd* deletions

Despite the task framing, **this paper contains no *gpd1Δ*/*gpd2Δ* strain of its
own**. It only cites the earlier work (Björkqvist et al. 1997; Nissen et al.
2000a; Valadi et al. 1998) to argue that deleting `GPD1` (**YDL022W**) and
`GPD2` (**YOL059W**) is not a workable strategy: the double mutant cannot grow
anaerobically because NADH cannot be reoxidised, and single deletions neither
reduce glycerol much nor raise ethanol. If a *gpd* phenotype is needed for
validation it has to come from those papers, not from Bro et al. 2006.

---

## 4. de Kok et al. 2011 — maltose phosphorylase + β-phosphoglucomutase

### Cultivation

* Strain background: CEN.PK, all strains carry the constitutive MAL activator
  **MAL2-8C**. Reference **CEN.PK113-7D** (*MAL1x MAL2x MAL3x MAL4x LEU2 URA3*);
  knockout platform **IMK289** from the doubly auxotrophic **CEN.PK102-3A**
  (*leu2-112 ura3-52*), with plasmids restoring LEU2/URA3.
* Medium: synthetic medium of Verduyn et al. 1990. Chemostats: **25 g L⁻¹
  maltose** as sole carbon and energy source, supplemented with **ergosterol
  10 mg L⁻¹ + Tween 80 420 mg L⁻¹** (dissolved in ethanol), Antifoam C
  0.2 g L⁻¹, pH 5.0 (2 M KOH), 30 °C. Shake flasks used a different medium:
  20 g L⁻¹ maltose, 6.6 g K₂SO₄, 3 g KH₂PO₄, **2.3 g urea as nitrogen source**,
  0.5 g MgSO₄·7H₂O, pH 5.0 (those cultures are aerobic and are not in the TSV).
* Anaerobiosis: 2 L bioreactors, 1 L working volume, 800 rpm, 500 mL min⁻¹ N₂
  (<10 ppm O₂), Norprene tubing, Viton O-rings, N₂-flushed medium vessels.
* Anaerobic growth rates were obtained from the slope of log(CO₂ production
  rate) in batch fermentations at pH 5.0 with 25 g L⁻¹ maltose.
* Chemostat **D = 0.05 h⁻¹** (chosen because IMZ226's anaerobic µmax is
  0.10 h⁻¹); steady state ≥5 volume changes with <4% change in dry weight and
  CO₂ rate over 2 volume changes, sampled between 10 and 15 volume changes to
  avoid evolutionary adaptation. Viability by CFDA/PI flow cytometry.

### Model changes

**Deletions (strain IMK289 platform, shared by IMZ199, IMZ203, IMZ226):** all
native maltose transporters and maltases —
`MAL11`/`AGT1` = **YGR289C**, `MAL12` = **YGR292W**,
`MAL31` = **YBR298C**, `MAL32` = **YBR299W**,
`MPH2` = **YDL247W**, `MPH3` = **YJR160C**, plus **MAL21 and MAL22** of the
MAL2 locus on chromosome III. **Uncertain:** the MAL2 locus is absent from the
S288C reference genome, so *MAL21*/*MAL22* have no systematic ORF names and will
usually have no counterpart in yeast-GEM; in practice deleting the YGR/YBR/YDL/
YJR genes above removes all maltose transport and maltase GPRs in the model.
Southern blotting confirmed CEN.PK102-3A carries three MAL loci (MAL1x, MAL2x,
MAL3x on chromosomes VII, III, II).

* **IMZ199 (hydrolytic reference)** = IMK289 + `MAL11` (TDH3 promoter,
  integrated) + `MAL12` (TDH3 promoter, 2 µm plasmid). Network-wise this is
  *wild-type maltose metabolism*: **maltose-H⁺ symport** plus intracellular
  hydrolysis.
  * Transport: `maltose[e] + H+[e] -> maltose[c] + H+[c]`; the imported proton
    is expelled by Pma1 at **1 H⁺ per ATP**, i.e. `ATP + H2O -> ADP + Pi + H+[e]`.
  * Maltase (EC 3.2.1.20): `maltose[c] + H2O -> 2 glucose[c]`
    (ΔG°′ = −15.5 kJ mol⁻¹), then 2 × hexokinase (2 ATP).
  * **Net: 3 ATP per maltose = 1.5 ATP per hexose.**
* **IMZ203** = IMK289 + `MAL11` + *L. sanfranciscensis* **mapA** only. Grows
  very poorly anaerobically (µ = 0.03 ± 0.01 h⁻¹, and even that without
  consistent exponential growth) because β-glucose-1-phosphate cannot be
  isomerised efficiently.
* **IMZ226 (phosphorolytic strain)** = IMK289 + `MAL11` + **mapA** + **pgmB**.
  * **Add** maltose phosphorylase (EC 2.4.1.8, *Lactobacillus sanfranciscensis*
    *mapA*, PGK1 promoter, 2 µm plasmid pUDE82), **cytosol**, reversible in
    principle (ΔG°′ = +5.5 kJ mol⁻¹) but used in the cleavage direction:
    `maltose[c] + phosphate[c] <=> D-glucose[c] + beta-D-glucose-1-phosphate[c]`
  * **Add** β-phosphoglucomutase (EC 5.4.2.6, *Lactococcus lactis* **pgmB**,
    codon-optimised, TDH3 promoter), **cytosol**, reversible:
    `beta-D-glucose-1-phosphate[c] <=> D-glucose-6-phosphate[c]`
    This is **essential**: native Pgm1 (YKL127W) and Pgm2 (YMR105C) act on
    α-glucose-1-phosphate only, and β-phosphoglucomutase activity was below
    detection (<0.005 µmol mg protein⁻¹ min⁻¹) in both CEN.PK113-7D and IMZ203.
    If a model already lumps glucose-1-phosphate without anomeric distinction,
    the pgmB requirement disappears and IMZ203's poor growth cannot be
    reproduced — worth stating explicitly when using these data.
  * Only **one** hexokinase reaction is now needed per maltose (for the free
    glucose), so **net: 4 ATP per maltose = 2 ATP per hexose**, the same per
    hexose as glucose.
  * Transport is **still maltose-H⁺ symport** (MAL11), so the 1 ATP proton cost
    remains.
* **Transporter alternatives discussed but not realised**: a maltose
  **uniporter** would give 5 ATP per maltose (2.5 per hexose). The MAL11 point
  mutants E120A, D123G and E167A reported to abolish proton coupling did **not**
  increase the biomass yield in IMZ226 (data not shown), so they apparently
  still symport — do not model them as uniporters.

### Authors' ATP-yield table (Table 6, anaerobic, ATP per substrate / per hexose)

glucose uniport 2 / 2; **maltose H⁺-symport + hydrolysis 3 / 1.5 (native)**;
maltose uniport + hydrolysis 4 / 2; **maltose H⁺-symport + phosphorolysis
4 / 2 (IMZ226)**; maltose uniport + phosphorolysis 5 / 2.5; maltotriose
H⁺-symport + hydrolysis 5 / 1.67; maltotriose uniport + hydrolysis 6 / 2;
maltotriose H⁺-symport + phosphorolysis 7 / 2.33; maltotriose uniport +
phosphorolysis 8 / 2.67. *P. furiosus* PF1535 maltodextrin phosphorylase was
expressed but had no detectable activity at 30 °C.

### Predicted vs measured biomass yield

* Theoretical (Weusthuis et al. 1993 model, maintenance 1 mmol ATP (g DW)⁻¹ h⁻¹,
  µ = 0.05 h⁻¹): 0.066 (3 ATP) vs 0.088 (4 ATP) g (g glucose equivalent)⁻¹, i.e.
  a **+33%** increase.
* Measured: 0.067 ± 0.000 (IMZ199) vs 0.085 ± 0.000 (IMZ226) = **+26%**.
  Validation: CEN.PK113-7D on glucose under the same anaerobic conditions gave
  0.088 ± 0.001, exactly the calculated 4-ATP value.
* The 26% vs 33% gap is attributed to a residual maltase activity
  (0.03 µmol mg protein⁻¹ min⁻¹, <1% of the reference's 3.6) still hydrolysing
  part of the maltose; an in vivo maltase flux of 0.35 out of
  1.58 mmol (g biomass)⁻¹ h⁻¹ would explain it exactly. Candidate sources named:
  Ima1–5, Suc2, Gdb1, Nth1, Ath1, Cwh41.

### Not in the TSV

IMZ133 (*malΔ mphΔ*), IMZ135 (+MAL12) and IMZ229 (+MAL11) have no anaerobic
values (`N.D.` in Table 4); their aerobic rates are <0.001, <0.001 and
0.007 ± 0.001 h⁻¹ respectively. All aerobic growth rates and all enzyme
activities in Table 4 were measured in aerobic shake flasks and are therefore
excluded; the key ones are quoted above.

---

## 5. Basso et al. 2011 — iSUC2, intracellular sucrose hydrolysis

### Cultivation

* Strain background: CEN.PK. Reference **CEN.PK113-7D** (*MATa URA3 TRP1 SUC2*).
  **IMI056** = *MATa URA3 trp1-289 TRP1-P_ADH1::iSUC2*; **IMM007** = IMI056 after
  90 generations of anaerobic sucrose-limited chemostat evolution (single-cell
  isolate); **IMM008/IMM009** = IMM007 with one/both *AGT1* copies deleted.
  Uracil prototrophy was restored in IMI056 because *ura3* interferes with
  chemostat characterisation.
* Medium: synthetic medium of Verduyn et al. 1992 with **25 g L⁻¹ sucrose**
  (filter-sterilised as 50% w/v), ammonium as nitrogen source, **ergosterol
  10 mg L⁻¹ + Tween 80 420 mg L⁻¹** in ethanol, Antifoam C 0.2 g L⁻¹, pH 5.0
  (2 M KOH), 30 °C.
* Anaerobiosis: 2 L bioreactors, 1.0 L working volume held by level sensor,
  800 rpm, 500 mL min⁻¹ N₂ (<10 ppm O₂), Norprene tubing, Viton O-rings,
  N₂-flushed medium vessels.
* **D = 0.10 h⁻¹**, preceded by a batch phase in the same vessel; steady state
  after ≥5 volume changes with <2% variation in dry weight and specific CO₂
  production over 2 volume changes. Ethanol corrected for evaporation. Sucrose,
  glucose and fructose by enzymatic assay; residual sugars rapidly quenched.
* Yields in Table 3 are per mol sucrose; the text also gives the reference
  strain as 0.094 g biomass and 0.38 g ethanol per g **sucrose monohydrate**.

### Model changes — this is a *topology* change, not a new reaction

**Reference (SUC2, extracellular hydrolysis), 4 ATP per sucrose = 2 per hexose:**

* Invertase (EC 3.2.1.26) in the **extracellular space / cell wall**:
  `sucrose[e] + H2O[e] -> D-glucose[e] + D-fructose[e]`
* Glucose and fructose enter by **facilitated diffusion** via Hxt transporters:
  `glucose[e] <=> glucose[c]`, `fructose[e] <=> fructose[c]` — no ATP, no proton.
* 2 × hexokinase (glucose) / hexokinase + Pfk route (fructose) consume 2 ATP,
  glycolysis yields 4, **net 4 ATP per sucrose**.

**iSUC2 strains (intracellular hydrolysis), 3 ATP per sucrose = 1.5 per hexose:**

* **Delete / relocate**: remove the extracellular invertase reaction. Genetically
  the first 60 nucleotides of `SUC2` = **YIL162W** (the 20-amino-acid N-terminal
  secretion signal) were deleted and the native inducible promoter replaced by
  the constitutive **ADH1** promoter, so the same gene now produces only the
  short, non-glycosylated cytosolic isoform. In a model: move the invertase
  reaction from the extracellular compartment to the **cytosol**:
  `sucrose[c] + H2O[c] -> D-glucose[c] + D-fructose[c]`, keeping YIL162W as GPR.
* **Add / enable**: **sucrose-H⁺ symport** across the plasma membrane,
  `sucrose[e] + H+[e] -> sucrose[c] + H+[c]`, encoded by **AGT1 = MAL11 =
  YGR289C** (α-glucoside-H⁺ symporter; Km for sucrose 6 ± 1 mM, unchanged by
  evolution). The imported proton must be expelled by Pma1 at **1 H⁺ per ATP**
  (`ATP + H2O -> ADP + Pi + H+[e]`), which is the entire source of the yield
  difference. Hexose transport by facilitated diffusion is no longer used for
  the sucrose carbon.
* **IMM007** additionally has *AGT1* **duplicated** and ~6-fold over-expressed
  (transcriptome), with general up-regulation of MAL genes (MPH2, MPH3, FSP2,
  IMA5, MAL31, MAL12, MAL32, MAL13); sequencing found no mutation in *AGT1* or
  800 bp upstream, so the cause is presumed to be a second-site regulator
  mutation (e.g. *MALx3*). For a stoichiometric model the only change vs IMI056
  is the **capacity** of the symporter (Vmax 1.74 → 20.8 mmol (g DW)⁻¹ h⁻¹),
  which matters because in IMI056 the symport capacity was *below* the sucrose
  consumption rate.
* **IMM009** = both *AGT1* copies deleted → symport capacity <0.3, yields back to
  the *SUC2* reference. The authors conclude a residual **extracellular**
  sucrose-hydrolysing activity of unknown identity exists (7 ± 32 and
  5 ± 3 mmol (g DW)⁻¹ h⁻¹ measured extracellularly in IMI056/IMM007) and
  competes with the intracellular route. **Ambiguity for modelling:** that
  unidentified activity has no gene assignment, so a model with invertase
  strictly cytosolic will over-predict the iSUC2 effect in IMI056 and IMM009.
* An alternative route mentioned but published elsewhere (Dario et al.,
  unpublished): intracellular maltase (Malx2p) plus `SUC2` deletion.

### Authors' prediction vs measurement

* Because intracellular hydrolysis gives 1.5 instead of 2 ATP per hexose,
  **one extra sucrose must be fermented for every three taken up**, predicting a
  **25% lower biomass yield**, a **25% lower glycerol yield** (glycerol is
  strictly growth-coupled anaerobically), and a higher ethanol yield given by
  `Y_intra = 0.75·Y_extra + 0.25·Y_max` with Y_max = 4 mol ethanol per mol
  sucrose, i.e. `R = 0.75 + 1/Y_extra`. With the measured
  Y_extra = 2.96 mol mol⁻¹ this gives **+8.7% ethanol**, essentially the 8%
  maltose-vs-glucose difference of Weusthuis et al. 1993.
* Measured: IMI056 (unevolved) only −6 ± 2% biomass and +4% ethanol, with a high
  residual sucrose of 1.8 g L⁻¹ — the transport step, not the topology, limited
  it. IMM007 (evolved) −30 ± 2% biomass, **+11 ± 2% ethanol**, −20% glycerol,
  residual sucrose <0.1 g L⁻¹ — in good agreement with the prediction. IMM009
  (both *AGT1* deleted) reverted to 1.01 / 1.01 / 0.96 of the reference.
* Extrapolation given by the authors: at µ = 0.025 h⁻¹ (industrial cell-recycle
  conditions) the predicted ethanol gain shrinks to **3%**, because a larger
  share of the sucrose is already dissimilated for maintenance
  (Y_extra = 3.56 mol mol⁻¹ derived from Boender et al. 2009).
* Anaerobic **batch** cultures of the evolved iSUC2 strain gave ethanol yields
  "consistently higher (up to 20%, which exceeded theoretical predictions)"
  than the reference, but carbon recoveries exceeded 100% in several of those
  cultures and the data are not shown — **not transcribed into the TSV**; use
  the chemostat data only.

### Unit check

Table 3 prints invertase activity and sucrose:H⁺ transport capacity as
"mmol g DW⁻¹ h⁻¹". This is consistent with the surrounding text (the reference
strain's symport capacity is stated to be 7-fold lower than its sucrose
consumption rate of 3.11 mmol g DW⁻¹ h⁻¹, and 3.11/7 = 0.44 ≈ 0.43), so mmol —
not µmol — was kept.

---

## Systematic ORF names used (all checked; uncertainties flagged)

| gene | ORF | note |
| --- | --- | --- |
| ACS1 | YAL054C | |
| ACS2 | YLR153C | |
| ALD2 / ALD3 / ALD4 / ALD5 / ALD6 | YMR170C / YMR169C / YOR374W / YER073W / YPL061W | Kozak aerobic strains only |
| GDH1 | YOR375C | |
| GDH2 | YDL215C | |
| GDH3 | YAL062W | not deleted, but must be constrained (see Nissen) |
| GLT1 | YDL171C | |
| GLN1 | YPR035W | |
| GPD1 / GPD2 | YDL022W / YOL059W | cited only, no strain in these five papers |
| MAL11 = AGT1 | YGR289C | sucrose/maltose-H⁺ symporter |
| MAL12 | YGR292W | |
| MAL13 | YGR288W | regulator, up-regulated in Basso IMM007 |
| MAL31 / MAL32 / MAL33 | YBR298C / YBR299W / YBR297W | |
| MAL21 / MAL22 | **none** | MAL2 locus, chromosome III; absent from the S288C reference genome — **uncertain**, no ORF id exists |
| MPH2 / MPH3 | YDL247W / YJR160C | |
| SUC2 | YIL162W | |
| GAP1 | YKR039W | background mutation in Bro's M4054 |
| PGM1 / PGM2 | YKL127W / YMR105C | α-phosphoglucomutases; cannot act on β-G1P |
| URA3 / TRP1 / LEU2 / HIS3 | YEL021W / YDR007W / YCL018W / YOR202W | auxotrophic markers only |
