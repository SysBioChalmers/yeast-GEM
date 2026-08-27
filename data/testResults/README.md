# Test results

Everything yeast-GEM checks about itself, what each check means, and the files
in this folder that hold the detail.

## 1. Where the results come from

Two places, and they answer different questions.

**Every pull request** that touches the model, the physiology data or the
Python package runs [`model-qc.yml`](../../.github/workflows/model-qc.yml),
which reports into a single comment on the pull request. It measures the
branch **and the branch it is targeting**, in the same run, and reports the
difference.

That comparison is the point. A model this size carries findings that predate
any one pull request, so an absolute count says little: 1753 reactions flagged
by MACAW is not something a curation branch caused or can be asked to fix. A
count that *rose* is. Because both sides are measured in the same job — same
code, same reference data, same solver, same dependency versions — the only
difference between the two columns is the model itself.

The icons read:

| | Meaning |
|---|---|
| :white_check_mark: | zero, or the model grows, or a metric improved |
| :warning: | a non-zero finding this pull request did not make worse — a report, not a blocker |
| :x: | a count that rose, or a metric that moved the wrong way, against the target branch |
| :grey_question: | the check did not run: its result is unknown, not clean |
| :hourglass_flowing_sand: | still being computed |

**Only growth blocks the merge.** Everything else reports, so that work is not
blocked on findings that were already there. A :x: asks for review; it does not
by itself mean the branch is wrong, since a curation can legitimately trade one
finding for another.

**At release**, `increaseVersion.m` regenerates the summary table below along
with `growth.md`, `growth.png` and `essentialGenes.tsv`, recording where the
model stood for that version.

## 2. Latest release results

- Model version: develop (post-9.1.0 curations, commit 58be7b3)
- Software: MATLAB 24.2.0.2923080 (R2024b) Update 6, RAVEN 2.11.2

| Metric | Value |
|---|---|
| Growth prediction R2 | 0.9019066691398621 |
| Anaerobic flux median fold error | 1.053461050692928 |
| Anaerobic fluxes with no comparable ratio | 11 |
| Anaerobic exchange mean relative error | 0.06454293960336684 |
| Anaerobic exchange within error | 0.75 |
| Ammonium per ATPase | 1.071951239704032 |
| Gene essentiality accuracy | 0.9024390243902439 |
| True non-essential genes | 934 |
| True essential genes | 65 |
| False non-essential genes | 94 |
| False essential genes | 14 |

Per measurement, at a glucose uptake rate of 23 mmol/gDW/h:

| Exchange | Measured | Predicted | Within error |
|---|---|---|---|
| glycerol | 4.5 +/- 0.4 | 4.3527 | yes |
| ethanol | 31 +/- 2 | 35.7072 | no |
| carbon dioxide | 38 +/- 10 | 38.7741 | yes |
| biomass | 0.36 +/- 0.02 | 0.3792 | yes |

Ethanol is predicted above the measured rate; the other three products fall
within their experimental error.

## 3. What each check means

The headings below match the check names in the pull-request comment exactly,
so a name there links straight to its explanation here — and to the section of
the same name in [`qc_findings.md`](qc_findings.md), which lists the entries
behind the count.

### Model checks

#### Growth (biomass producible)
**Gate.** The objective value under the model's own constraints, from
`slim_optimize`. A model that cannot produce biomass blocks the merge: nothing
else in the report means much if the model does not work.

#### Exact-duplicate reaction groups
Reactions sharing identical stoichiometry and direction. Keyed on the
metabolite/coefficient set together with reversibility, so a forward-only and a
reversible copy of the same conversion are counted as different reactions
rather than conflated. Reported as one row per reaction in each group larger
than one.

#### Unused metabolites
Metabolites taking part in no reaction. Usually the residue of a removed
reaction whose metabolites were left behind.

#### Unused genes
Genes associated with no reaction, normally left over after a gene rule was
edited.

#### Metabolites missing formula
Metabolites with no chemical formula. These also make every reaction they
appear in impossible to mass-balance, so a single missing formula shows up
twice: once here, and once in the mass-balance count.

#### Metabolites missing charge
Metabolites with no charge assigned, which likewise makes their reactions
impossible to charge-balance.

### Mass/charge balance and MACAW

#### Mass-imbalanced reactions
Elemental imbalance, from cobrapy's `Reaction.check_mass_balance()`.

Three classes are excluded because they are not expected to balance:

- **boundary reactions** — an exchange, demand or sink has one side by
  construction;
- **pseudoreactions** — biomass and its components are lumped compositions;
- **SLIME reactions** — yeast-GEM splits lipids into measurable entities using
  fractional, averaged coefficients.

The last exclusion is what makes this count usable. There are 186 SLIME
reactions, against 2 genuinely unexpected mass imbalances in the rest of the
model; including them would report 194 findings that nobody can act on, which
is the same as reporting nothing.

#### Charge-imbalanced reactions
Charge imbalance, from the same call and with the same exclusions.

#### Reactions flagged by MACAW dead-end test
Reactions that cannot carry flux, or can only carry it in one direction,
because a metabolite in them can only be produced or only consumed. From
[MACAW](https://github.com/Devlin-Moyer/macaw).

#### Reactions flagged as MACAW duplicates
Reactions MACAW considers redundant. It runs four sub-tests — exact duplicates,
and pairs differing only in direction, only in coefficients, or only in redox
partner — and writes one column per sub-test rather than a single verdict. The
count is the union: a reaction flagged by any sub-test is counted once, since
the sub-tests overlap.

A sub-test that does not apply to a reaction is reported by MACAW as the
literal string `N/A`, which is neither a null nor `ok`. Reading MACAW's output
back through pandas turns it into a null, so a file and the values a count is
taken from can disagree — which is why this count is taken from the frame in
memory and never from a written file.

### Annotations

#### Malformed cross-references
Cross-references that do not match the identifier pattern their namespace
declares at [identifiers.org](https://identifiers.org). Only namespaces with a
pattern recorded in `model_qc.py` are checked; a namespace without one is left
alone rather than guessed at, because a wrong pattern manufactures findings
that look exactly like real ones.

#### Cross-refs inconsistent across compartments
The same metabolite name appearing in several compartments with disagreeing
cross-references.

Read these before acting on them. A different ChEBI identifier across
compartments is sometimes correct, where the compartments genuinely hold
different protonation states or a zwitterion. The finding is that the name and
the identifier disagree about whether they are the same species: either the
names should differ, or the annotations should match.

### Validation metrics

Predictions against measured data. These are continuous rather than counts, so
each declares which direction is an improvement and a tolerance below which a
change is reported as unchanged — otherwise solver noise reads as a result.

#### Growth prediction R2
`growth.m` compares predicted against measured growth rates over the chemostat
data of [Österlund _et al._ (2013)](https://doi.org/10.1186/1752-0509-7-36),
across aerobic and anaerobic, carbon- and nitrogen-limited conditions. Reported
as the coefficient of determination about the line of identity. Higher is
better. The per-condition plot is in [growth.md](growth.md).

#### Anaerobic flux median fold error
`anaerobic_flux_predictions.m` fixes the glucose uptake rate to the measured
value for each dataset and compares predicted intracellular fluxes, scaled to
100 · v<sub>i</sub> / v<sub>glucose</sub>, against the measurements of
[Jouhten _et al._ (2008)](https://doi.org/10.1186/1752-0509-2-60).

Reported as a fold error, `10^median(|log10(predicted/measured)|)`. A value of
1.05 means the typical prediction is 5% out in either direction. Lower is
better, and 1.0 is exact.

Fold error rather than a difference, because these fluxes span three orders of
magnitude — 0.1 to 163 in the scaled units — and being 2× out matters as much
at a flux of 1 as at a flux of 100. A difference-based score answers "are the
big fluxes right", which is nearly a restatement of the glucose constraint:
the coefficient of determination reads 0.997 on a model whose smaller half of
measurements is out by about 50%, which is why it is not reported.

This is the only summary of how far out the comparable predictions are. A mean
fold error, a fraction within 2-fold and an R2 were all reported at one point,
but they are restatements of the same quantity and moved together, so the
table said the same thing four times.

#### Anaerobic fluxes with no comparable ratio
Measurements for which no fold error exists: a flux measured as non-zero but
predicted as zero, or predicted in the opposite direction.

Counted rather than dropped, and reported separately from the fold error
because it is a different kind of failure. Both are worse than any finite fold
error, and excluding them from the denominator would make the fold error
improve as the model got worse. Lower is better.

#### Anaerobic exchange mean relative error
Mean absolute relative deviation between predicted and measured exchange rates
for the main fermentation products, measured by [Sjöberg _et al._
(2024)](https://doi.org/10.1016/j.ymben.2024.01.007). Lower is better.

No single R<sup>2</sup> is reported for these: the product rates are in
mmol/gDW/h and growth is in 1/h, so pooling them into one coefficient of
determination would not mean anything, whereas the relative deviation per
measurement is comparable across units.

#### Anaerobic exchange within error
The fraction of those predictions falling inside the experimental error of the
measurement. Higher is better.

#### Ammonium per ATPase
Ammonium uptake divided by plasma-membrane ATPase flux. Ammonium uptake runs
against the proton gradient the ATPase maintains, so the two are coupled and
the measured ratio is close to 1. Neither direction is an improvement, so a
move is reported but never counted as a regression.

#### Gene essentiality accuracy
`essentialGenes.m` deletes each gene in turn on Kennedy synthetic complete
medium and compares predicted viability against the Stanford yeast deletion
collection. Higher is better. The per-gene lists are in
[essentialGenes.tsv](essentialGenes.tsv).

Note that the deletion collection was screened on complex medium supplemented
for genetic markers, so it is not ideal ground truth; it is most useful for
comparing model versions against each other.

#### True non-essential genes
Genes correctly predicted to be viable when deleted. Higher is better.

#### True essential genes
Genes correctly predicted to be lethal when deleted. Higher is better.

#### False non-essential genes
Genes predicted viable that are lethal in the collection. Lower is better.

#### False essential genes
Genes predicted lethal that are viable in the collection. Lower is better.

## 4. Files in this folder

Alphabetically. Regenerated by CI on each pull request and committed, so that a
change in a finding shows in the pull request's own diff rather than only in a
comment that later scrolls away:

| File | Holds |
|---|---|
| `qc_findings.md` | every entry behind every count, one section per check |
| `qc_metrics.tsv` | every structural, balance and annotation count, as `key<TAB>value` |
| `validation_metrics.tsv` | every validation metric, same format |

`qc_findings.md` keeps all its section headings whether or not the check found
anything: an empty section reads `_None._`. One file with a fixed set of
sections diffs far better than a CSV per check that appears and disappears —
a check going clean shows as rows vanishing under a heading that stays put,
rather than as a deleted file.

Regenerated at release by `increaseVersion.m`:

| File | Holds |
|---|---|
| `essentialGenes.tsv` | one row per gene, sorted, with `TP`/`TN`/`FP`/`FN` |
| `growth.md` | measured against predicted growth, per condition |
| `growth.png` | the same, as a figure |
| `README.md` | the summary table in section 2 |

`essentialGenes.tsv` is one sorted table rather than four lists of gene names,
so a gene changing category is a one-line diff. In the list form the gene moved
between sections, which showed as two edits far apart in the file and made a
single reclassification easy to miss.

## 5. Not yet implemented

Human-GEM additionally checks model/annotation-table consistency, that removed
identifiers were moved to a deprecated list, and structure (SMILES/InChI)
against formula and charge. All three need the YAML-based curation layout with
separate identifier TSVs, which yeast-GEM adopts in 9.2.0 (see
[#379](https://github.com/SysBioChalmers/yeast-GEM/issues/379)).
