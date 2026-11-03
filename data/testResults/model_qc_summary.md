## Model quality report

:warning: **9 pre-existing finding(s), no regressions vs `develop`.** Non-blocking.

**Model size:** 4105 reactions (0), 2748 metabolites (0), 1143 genes (0) vs `develop`

_Each check name links to its explanation in the [testResults README](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md)._

### Model checks
_Growth, duplicate keys and reactions without metabolites are gates and block the merge; every other row is a non-blocking report._

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| [Growth (biomass producible)](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#growth-biomass-producible) | 0.0811 | 0 | :white_check_mark: |
| [Duplicate keys in the model file](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#duplicate-keys-in-the-model-file) | 0 | 0 | :white_check_mark: |
| [Reactions without metabolites](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#reactions-without-metabolites) | 0 | 0 | :white_check_mark: |
| [Exact-duplicate reaction groups](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#exact-duplicate-reaction-groups) | 0 | 0 | :white_check_mark: |
| [Unused metabolites](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#unused-metabolites) | 0 | 0 | :white_check_mark: |
| [Unused genes](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#unused-genes) | 0 | 0 | :white_check_mark: |
| [Metabolites missing formula](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#metabolites-missing-formula) | [16](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/qc_findings.md#metabolites-missing-formula) | 0 | :warning: |
| [Metabolites missing charge](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#metabolites-missing-charge) | 0 | 0 | :white_check_mark: |
| [Reaction bound and gene rule issues](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#reaction-bound-and-gene-rule-issues) | 0 | 0 | :white_check_mark: |
| [Reactions split across compartments](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#reactions-split-across-compartments) | [20](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/qc_findings.md#reactions-split-across-compartments) | 0 | :warning: |

### Mass/charge balance and MACAW

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| [Mass-imbalanced reactions](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#mass-imbalanced-reactions) | [8](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/qc_findings.md#mass-imbalanced-reactions) | 0 | :warning: |
| [Charge-imbalanced reactions](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#charge-imbalanced-reactions) | [10](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/qc_findings.md#charge-imbalanced-reactions) | 0 | :warning: |
| [Dead-end metabolites](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#dead-end-metabolites) | [943](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/qc_findings.md#dead-end-metabolites) | 0 | :warning: |
| [Reactions that can only carry flux one way](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#reactions-that-can-only-carry-flux-one-way) | [570](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/qc_findings.md#reactions-that-can-only-carry-flux-one-way) | 0 | :warning: |
| [Reactions flagged as MACAW duplicates](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#reactions-flagged-as-macaw-duplicates) | [18](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/qc_findings.md#reactions-flagged-as-macaw-duplicates) | 0 | :warning: |

### Annotations
_Reaction/metabolite name agreement with the tsvs, and model/annotation-table consistency, are also gates; every other row is a non-blocking report._

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| [Malformed cross-references](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#malformed-cross-references) | [34](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/qc_findings.md#malformed-cross-references) | 0 | :warning: |
| [Cross-refs inconsistent across compartments](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#cross-refs-inconsistent-across-compartments) | 0 | 0 | :white_check_mark: |
| [Reaction/metabolite names disagreeing with the tsvs](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#reactionmetabolite-names-disagreeing-with-the-tsvs) | 0 | 0 | :white_check_mark: |
| [Model/annotation-table consistency](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#modelannotation-table-consistency) | 0 | 0 | :white_check_mark: |
| [Removed identifiers not added to the deprecated lists](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#removed-identifiers-not-added-to-the-deprecated-lists) | 0 | 0 | :white_check_mark: |
| [Metabolite structure (SMILES) vs. formula/charge](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#metabolite-structure-smiles-vs-formulacharge) | [916](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/qc_findings.md#metabolite-structure-smiles-vs-formulacharge) | 0 | :warning: |
| [deltaG tsvs out of sync with the model](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#deltag-tsvs-out-of-sync-with-the-model) | 0 | 0 | :white_check_mark: |

### Model files
_All three are gates._

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| [YAML round-trip](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#yaml-round-trip) | pass |  | :white_check_mark: |
| [SBML round-trip](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#sbml-round-trip) | pass |  | :white_check_mark: |
| [YAML lint](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#yaml-lint) | pass |  | :white_check_mark: |

### Validation metrics
_Predictions against measured data, to 4 significant digits. A change within tolerance still gets a checkmark; direction matters, so a rising R2 is a gain and a rising error is not._

| Metric | This branch | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| [Growth prediction R2](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#growth-prediction-r2) | [0.9019](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/validation_findings.md#growth-prediction-r2) | 0 | :white_check_mark: |
| [Anaerobic flux median fold error](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#anaerobic-flux-median-fold-error) | [1.062](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/validation_findings.md#anaerobic-flux-median-fold-error) | 0 | :white_check_mark: |
| [Anaerobic exchange mean relative error](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#anaerobic-exchange-mean-relative-error) | [0.06455](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/validation_findings.md#anaerobic-exchange-mean-relative-error) | 0 | :white_check_mark: |
| [Ammonium per ATPase](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#ammonium-per-atpase) | [1.072](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/validation_findings.md#ammonium-per-atpase) | 0 | :white_check_mark: |
| [Gene essentiality accuracy](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#gene-essentiality-accuracy) | [0.9015](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/validation_findings.md#gene-essentiality-accuracy) | 0 | :white_check_mark: |
| [True non-essential genes](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#true-non-essential-genes) | [933](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/validation_findings.md#gene-essentiality-accuracy) | 0 | :white_check_mark: |
| [True essential genes](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#true-essential-genes) | [65](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/validation_findings.md#gene-essentiality-accuracy) | 0 | :white_check_mark: |
| [False non-essential genes](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#false-non-essential-genes) | [94](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/validation_findings.md#gene-essentiality-accuracy) | 0 | :warning: |
| [False essential genes](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/README.md#false-essential-genes) | [15](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/validation_findings.md#gene-essentiality-accuracy) | 0 | :warning: |

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

[Full report](https://github.com/SysBioChalmers/yeast-GEM/blob/ci/model-qc-report/data/testResults/annotation_report.md)

:x: = a count rose vs the target branch (regression) &middot; :warning: = a non-zero or worse-than-baseline finding (non-blocking) &middot; :hourglass_flowing_sand: = still running &middot; :grey_question: = the check did not run, or has no baseline yet.

_Most deltas are computed in one run -- the target branch is checked out and measured the same way -- so they reflect this pull request and nothing else. MEMOTE and the MetaNetX annotation check are the exception: compared against the target branch's own last-committed result instead._
