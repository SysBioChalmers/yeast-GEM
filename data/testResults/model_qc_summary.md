## Model quality report

:warning: **9 pre-existing finding(s), no regressions vs `develop`.** Non-blocking.

**Model size:** 4105 reactions (0), 2748 metabolites (0), 1142 genes (-1) vs `develop`

_Each check name links to its explanation in the [testResults README](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md)._

### Model checks
_Growth, duplicate keys and reactions without metabolites are gates and block the merge; every other row is a non-blocking report._

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| [Growth (biomass producible)](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#growth-biomass-producible) | 0.0811 | 0 | :white_check_mark: |
| [Duplicate keys in the model file](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#duplicate-keys-in-the-model-file) | 0 | 0 | :white_check_mark: |
| [Reactions without metabolites](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#reactions-without-metabolites) | 0 | 0 | :white_check_mark: |
| [Exact-duplicate reaction groups](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#exact-duplicate-reaction-groups) | 0 | 0 | :white_check_mark: |
| [Unused metabolites](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#unused-metabolites) | 0 | 0 | :white_check_mark: |
| [Unused genes](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#unused-genes) | 0 | 0 | :white_check_mark: |
| [Metabolites missing formula](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#metabolites-missing-formula) | [5](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/qc_findings.md#metabolites-missing-formula) | 0 | :warning: |
| [Metabolites missing charge](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#metabolites-missing-charge) | 0 | 0 | :white_check_mark: |
| [Reaction bound and gene rule issues](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#reaction-bound-and-gene-rule-issues) | 0 | 0 | :white_check_mark: |
| [Reactions split across compartments](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#reactions-split-across-compartments) | [20](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/qc_findings.md#reactions-split-across-compartments) | 0 | :warning: |

### Mass/charge balance and MACAW

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| [Mass-imbalanced reactions](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#mass-imbalanced-reactions) | [8](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/qc_findings.md#mass-imbalanced-reactions) | 0 | :warning: |
| [Charge-imbalanced reactions](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#charge-imbalanced-reactions) | [10](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/qc_findings.md#charge-imbalanced-reactions) | 0 | :warning: |
| [Dead-end metabolites](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#dead-end-metabolites) | [943](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/qc_findings.md#dead-end-metabolites) | 0 | :warning: |
| [Reactions that can only carry flux one way](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#reactions-that-can-only-carry-flux-one-way) | [570](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/qc_findings.md#reactions-that-can-only-carry-flux-one-way) | 0 | :warning: |
| [Reactions flagged as MACAW duplicates](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#reactions-flagged-as-macaw-duplicates) | [18](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/qc_findings.md#reactions-flagged-as-macaw-duplicates) | 0 | :warning: |

### Annotations
_Reaction/metabolite name agreement with the tsvs, and model/annotation-table consistency, are also gates; every other row is a non-blocking report._

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| [Malformed cross-references](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#malformed-cross-references) | [34](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/qc_findings.md#malformed-cross-references) | 0 | :warning: |
| [Cross-refs inconsistent across compartments](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#cross-refs-inconsistent-across-compartments) | 0 | 0 | :white_check_mark: |
| [Reaction/metabolite names disagreeing with the tsvs](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#reactionmetabolite-names-disagreeing-with-the-tsvs) | 0 | 0 | :white_check_mark: |
| [Model/annotation-table consistency](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#modelannotation-table-consistency) | 0 | 0 | :white_check_mark: |
| [Removed identifiers not added to the deprecated lists](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#removed-identifiers-not-added-to-the-deprecated-lists) | 0 | 0 | :white_check_mark: |
| [Metabolite structure (SMILES) vs. formula/charge](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#metabolite-structure-smiles-vs-formulacharge) | [916](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/qc_findings.md#metabolite-structure-smiles-vs-formulacharge) | 0 | :warning: |
| [deltaG tsvs out of sync with the model](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#deltag-tsvs-out-of-sync-with-the-model) | 0 | 0 | :white_check_mark: |

### Model files
_All three are gates._

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| [YAML round-trip](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#yaml-round-trip) | pass |  | :white_check_mark: |
| [SBML round-trip](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#sbml-round-trip) | pass |  | :white_check_mark: |
| [YAML lint](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#yaml-lint) | pass |  | :white_check_mark: |

### Validation metrics
_Predictions against measured data, to 4 significant digits. A change within tolerance still gets a checkmark; direction matters, so a rising R2 is a gain and a rising error is not._

| Metric | This branch | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| [Growth prediction R2](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#growth-prediction-r2) | [0.9019](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/validation_findings.md#growth-prediction-r2) | 0 | :white_check_mark: |
| [Anaerobic flux median fold error](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#anaerobic-flux-median-fold-error) | [1.062](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/validation_findings.md#anaerobic-flux-median-fold-error) | 0 | :white_check_mark: |
| [Anaerobic exchange mean relative error](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#anaerobic-exchange-mean-relative-error) | [0.06455](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/validation_findings.md#anaerobic-exchange-mean-relative-error) | 0 | :white_check_mark: |
| [Ammonium per ATPase](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#ammonium-per-atpase) | [1.072](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/validation_findings.md#ammonium-per-atpase) | 0 | :white_check_mark: |
| [Gene essentiality accuracy](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#gene-essentiality-accuracy) | [0.9014](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/validation_findings.md#gene-essentiality-accuracy) | -8.903e-05 | :white_check_mark: |
| [True non-essential genes](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#true-non-essential-genes) | [932](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/validation_findings.md#gene-essentiality-accuracy) | -1 | :white_check_mark: |
| [True essential genes](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#true-essential-genes) | [65](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/validation_findings.md#gene-essentiality-accuracy) | 0 | :white_check_mark: |
| [False non-essential genes](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#false-non-essential-genes) | [94](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/validation_findings.md#gene-essentiality-accuracy) | 0 | :warning: |
| [False essential genes](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/README.md#false-essential-genes) | [15](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/validation_findings.md#gene-essentiality-accuracy) | 0 | :warning: |

### Gene essentiality
_Single-gene deletions against the Stanford deletion collection. A false non-essential gene is essential in the collection and viable in the model; a false essential gene is viable in the collection and not in the model. The per-gene calls are in [essentialGenes.tsv](essentialGenes.tsv)._

#### Changed since `develop`

| Gene | Name | Before | Now | Effect |
| --- | --- | --- | --- | --- |
| `YBR132C` | AGP2 | true non-essential | not in the model | removed |

#### False non-essential genes (94)

| Gene | Name |
| --- | --- |
| `YAL038W` | CDC19 |
| `YBL030C` | PET9 |
| `YBR002C` | RER2 |
| `YBR004C` | GPI18 |
| `YBR011C` | IPP1 |
| `YBR038W` | CHS2 |
| `YBR070C` | ALG14 |
| `YBR110W` | ALG1 |
| `YBR153W` | RIB7 |
| `YBR192W` | RIM2 |
| `YBR196C` | PGI1 |
| `YBR243C` | ALG7 |
| `YBR252W` | DUT1 |
| `YBR256C` | RIB5 |
| `YBR265W` | TSC10 |
| `YCL052C` | PBN1 |
| `YCR012W` | PGK1 |
| `YDL004W` | ATP16 |
| `YDL015C` | TSC13 |
| `YDL120W` | YFH1 |
| `YDL166C` | FAP7 |
| `YDR023W` | SES1 |
| `YDR050C` | TPI1 |
| `YDR062W` | LCB2 |
| `YDR208W` | MSS4 |
| `YDR236C` | FMN1 |
| `YDR302W` | GPI11 |
| `YDR353W` | TRR1 |
| `YDR367W` | KEI1 |
| `YDR373W` | FRQ1 |
| `YDR376W` | ARH1 |
| `YDR437W` | GPI19 |
| `YDR487C` | RIB3 |
| `YER023W` | PRO3 |
| `YGL001C` | ERG26 |
| `YGL008C` | PMA1 |
| `YGL047W` | ALG13 |
| `YGL065C` | ALG2 |
| `YGL142C` | GPI10 |
| `YGL169W` | SUA5 |
| `YGL225W` | VRG4 |
| `YGR046W` | TAM41 |
| `YGR060W` | ERG25 |
| `YGR065C` | VHT1 |
| `YGR175C` | ERG1 |
| `YGR191W` | HIP1 |
| `YGR216C` | GPI1 |
| `YGR277C` | CAB4 |
| `YHR007C` | ERG11 |
| `YHR042W` | NCP1 |
| `YHR068W` | DYS1 |
| `YHR072W` | ERG7 |
| `YHR128W` | FUR1 |
| `YHR190W` | ERG9 |
| `YIL048W` | NEO1 |
| `YJL005W` | CYR1 |
| `YJL026W` | RNR2 |
| `YJL091C` | GWT1 |
| `YJL097W` | PHS1 |
| `YJR013W` | GPI14 |
| `YJR016C` | ILV3 |
| `YJR057W` | CDC8 |
| `YKL004W` | AUR1 |
| `YKL019W` | RAM2 |
| `YKL035W` | UGP1 |
| `YKL060C` | FBA1 |
| `YKL141W` | SDH3 |
| `YKL152C` | GPM1 |
| `YKL165C` | MCD4 |
| `YKL182W` | FAS1 |
| `YLL031C` | GPI13 |
| `YLR100W` | ERG27 |
| `YLR153C` | ACS2 |
| `YLR305C` | STT4 |
| `YLR355C` | ILV5 |
| `YLR359W` | ADE13 |
| `YML125C` | PGA3 |
| `YMR013C` | SEC59 |
| `YMR108W` | ILV2 |
| `YMR281W` | GPI12 |
| `YMR296C` | LCB1 |
| `YMR298W` | LIP1 |
| `YMR301C` | ATM1 |
| `YNL036W` | NCE103 |
| `YNL038W` | GPI15 |
| `YNL267W` | PIK1 |
| `YOL066C` | RIB2 |
| `YOR143C` | THI80 |
| `YOR149C` | SMP3 |
| `YPL076W` | GPI2 |
| `YPL175W` | SPT14 |
| `YPL231W` | FAS2 |
| `YPL252C` | YAH1 |
| `YPR035W` | GLN1 |

#### False essential genes (15)

| Gene | Name |
| --- | --- |
| `YBR126C` | TPS1 |
| `YDR074W` | TPS2 |
| `YER014W` | HEM14 |
| `YER026C` | CHO1 |
| `YGR157W` | CHO2 |
| `YJR073C` | OPI3 |
| `YKL067W` | YNK1 |
| `YKR039W` | GAP1 |
| `YKR067W` | GPT2 |
| `YKR072C` | SIS2 |
| `YMR217W` | GUA1 |
| `YNL275W` | BOR1 |
| `YOL122C` | SMF1 |
| `YOR054C` | VHS3 |
| `YPL172C` | COX10 |

### MEMOTE score
_Fast subset -- flux-variability and matrix-rank tests are skipped, so this is not the same score `memote report snapshot` would give. Compared against `develop`'s own last-committed score, not a fresh run on it. Comment `/run memote` on the pull request for the full suite._

| Section | Score | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| **Total** | 66.2% | 0 | :white_check_mark: |
| consistency | 42.0% | 0 | :white_check_mark: |
| annotation_met | 72.6% | 0 | :white_check_mark: |
| annotation_rxn | 71.0% | 0 | :white_check_mark: |
| annotation_gene | 56.9% | 0 | :white_check_mark: |
| annotation_sbo | 88.5% | 0 | :white_check_mark: |

<details><summary>Per-test scores</summary>

| Section | Test | Score |
| --- | --- | ---: |
| Consistency | Stoichiometric Consistency | 100.0% |
| Consistency | Mass Balance | 5.3% |
| Consistency | Charge Balance | 0.4% |
| Consistency | Metabolite Connectivity | 0.0% |
| Consistency | Unbounded Flux In Default Medium | 100.0% |
| Annotation - Metabolites | Presence of Metabolite Annotation | 0.0% |
| Annotation - Metabolites | Metabolite Annotations Per Database | 73.2% |
| Annotation - Metabolites | Metabolite Annotation Conformity Per Database | 36.4% |
| Annotation - Metabolites | Uniform Metabolite Identifier Namespace | 0.0% |
| Annotation - Reactions | Presence of Reaction Annotation | 0.0% |
| Annotation - Reactions | Reaction Annotations Per Database | 82.6% |
| Annotation - Reactions | Reaction Annotation Conformity Per Database | 33.4% |
| Annotation - Reactions | Uniform Reaction Identifier Namespace | 0.0% |
| Annotation - Genes | Presence of Gene Annotation | 0.0% |
| Annotation - Genes | Gene Annotations Per Database | 89.3% |
| Annotation - Genes | Gene Annotation Conformity Per Database | 40.0% |
| Annotation - SBO Terms | Metabolite General SBO Presence | 0.0% |
| Annotation - SBO Terms | Metabolite SBO:0000247 Presence | 1.9% |
| Annotation - SBO Terms | Reaction General SBO Presence | 0.0% |
| Annotation - SBO Terms | Metabolic Reaction SBO:0000176 Presence | 0.9% |
| Annotation - SBO Terms | Transport Reaction SBO:0000185 Presence | 23.4% |
| Annotation - SBO Terms | Exchange Reaction SBO:0000627 Presence | 0.0% |
| Annotation - SBO Terms | Demand Reaction SBO:0000628 Presence | 100.0% |
| Annotation - SBO Terms | Sink Reactions SBO:0000632 Presence | 0.0% |
| Annotation - SBO Terms | Gene General SBO Presence | 0.0% |
| Annotation - SBO Terms | Gene SBO:0000243 Presence | 0.0% |
| Annotation - SBO Terms | Biomass Reactions SBO:0000629 Presence | 0.0% |

</details>

### Annotation verification (MetaNetX)
_Metabolite/reaction identifiers cross-checked against MetaNetX's reference data -- see the full report below for detail._

| Finding | Count | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| `wrong`, isolated | 30 | 0 | :white_check_mark: |
| `wrong`, recurring | 13 | 0 | :white_check_mark: |
| `drift` (safe to `--fix`) | 0 | 0 | :white_check_mark: |
| `missing` (safe to `--fix`) | 0 | 0 | :white_check_mark: |

[Full report](https://github.com/SysBioChalmers/yeast-GEM/blob/curation/gene-rules-tyr1-frd-carnitine/data/testResults/annotation_report.md)

:x: = a count rose vs the target branch (regression) &middot; :warning: = a non-zero or worse-than-baseline finding (non-blocking) &middot; :hourglass_flowing_sand: = still running &middot; :grey_question: = the check did not run, or has no baseline yet.

_Most deltas are computed in one run -- the target branch is checked out and measured the same way -- so they reflect this pull request and nothing else. MEMOTE and the MetaNetX annotation check are the exception: compared against the target branch's own last-committed result instead._
