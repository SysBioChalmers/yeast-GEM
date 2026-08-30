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

**At release**, `increase_version.py` regenerates the summary table below along
with `growth.md`, `growth.png` and `essentialGenes.tsv`, recording where the
model stood for that version.

## 2. Latest release results

- Model version: 9.1.1
- Software: Python 3.12.14, cobrapy 0.32.1

| Metric | Value |
|---|---|
| Growth prediction R2 | 0.9019063673672041 |
| Anaerobic flux median fold error | 1.061500419010032 |
| Anaerobic exchange mean relative error | 0.0645475428321038 |
| Ammonium per ATPase | 1.07195123970411 |
| Gene essentiality accuracy | 0.9015356820234869 |
| True non-essential genes | 933 |
| True essential genes | 65 |
| False non-essential genes | 94 |
| False essential genes | 15 |

Per measurement, at a glucose uptake rate of 23 mmol/gDW/h:

| Exchange | Measured | Predicted | Within error |
|---|---|---|---|
| glycerol | 4.5 +/- 0.4 | 4.3526 | yes |
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

### Model size
Reaction, metabolite and gene counts, with the change against the target
branch. Shown at the top of the comment, but not a check: unlike every count
below, a curation pull request is expected to change these, so it carries no
icon and is never a regression. It exists so a size change is visible, not so
it can be judged.

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

#### Dead-end metabolites
Metabolites that every reaction they appear in can only produce, or only
consume, so no steady-state flux through them is possible. From
[MACAW](https://github.com/Devlin-Moyer/macaw)'s dead-end test. The commonest
cause is a missing reaction; sometimes it is a metabolite annotated as two
different species that should be one.

Counted per metabolite rather than per blocked reaction. MACAW flags a reaction
for each dead-end metabolite in it, so one gap is counted many times: on the
committed model 943 dead-end metabolites account for 1753 flagged reactions,
and a single metabolite accounts for 49 of them. The metabolite is what gets
fixed, and the blocked reactions follow. `qc_findings.md` lists each metabolite
with its name and how many reactions it blocks, so the ones worth attention
sort to the top.

#### Reactions that can only carry flux one way
Reversible reactions that are not themselves dead ends, but that have a
reactant or product which every *other* reaction it appears in can only consume
or only produce. The reaction can therefore only ever run in one direction, and
MACAW reports which: `only when going forwards` or `only when going backwards`.

Either the reversibility is wrong and should be corrected, or something is
missing that would let the metabolite move the other way. This is reported
separately from the dead-end count because it is a different fix.

#### Reactions flagged as MACAW duplicates
Reactions MACAW considers possible duplicates of each other. It runs four
sub-tests and writes one column per sub-test rather than a single verdict; the
count is the union, since a reaction can be flagged by more than one.
`qc_findings.md` names the reaction it duplicates and why, in these terms:

| Reported as | Means |
|---|---|
| same metabolites and coefficients | identical stoichiometry, possibly different genes |
| same metabolites, opposite direction or reversibility | one runs the other way, or one is reversible and the other is not |
| same metabolites, different coefficients | same participants, different stoichiometry |
| same conversion using a different electron carrier | the same redox chemistry with, say, NAD(H) in one and NADP(H) in the other |

Not every one is a mistake: separate irreversible importer and exporter
reactions for the same metabolite are legitimate when different genes encode
them, as is an enzyme that genuinely uses either NAD(H) or NADP(H).

The redox sub-test needs a list of oxidised/reduced carrier pairs, which is not
currently configured, so MACAW writes `N/A` for it and it is reported as not
run rather than as zero findings.

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

This is the only summary of how far out the predictions are. A mean fold
error, a fraction within 2-fold and an R2 were all reported at one point, but
they are restatements of the same quantity and moved together, so the table
said the same thing four times.

One caveat, since it is no longer reported as its own number: a ratio needs
both values non-zero and pointing the same way, so measurements predicted as
zero or in the opposite direction have no fold error and are left out of the
median. There are 11 such measurements out of 90 on the committed model.
Because they are excluded rather than penalised, a change that turns a
badly-predicted flux into an unpredicted one improves this number. Read it
together with the flux plot rather than on its own.

#### Anaerobic exchange mean relative error
Mean absolute relative deviation between predicted and measured exchange rates
for the main fermentation products, measured by [Sjöberg _et al._
(2024)](https://doi.org/10.1016/j.ymben.2024.01.007). Lower is better.

No single R<sup>2</sup> is reported for these: the product rates are in
mmol/gDW/h and growth is in 1/h, so pooling them into one coefficient of
determination would not mean anything, whereas the relative deviation per
measurement is comparable across units.

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

### MEMOTE score

[`memote_snapshot.py`](../../code/python/tests/ci/memote_snapshot.py) runs the
[MEMOTE](https://memote.readthedocs.io) test suite in-process and scores this
branch, skipping the flux-variability and matrix-rank tests that dominate its
runtime on a genome-scale model — so this is a faster, lower number than
`memote report snapshot` (used at release, see below) would give. Reported as
a total score, one score per MEMOTE section (consistency, annotation
completeness by entity type), and the individual test scores behind those
sections in a collapsed table.

Unlike the other checks, it is not measured on both branches in this run:
even the fast subset roughly doubles the job's cost if run twice, for a score
that mostly reflects annotation completeness rather than something a typical
pull request moves. Instead it is compared against the target branch's own
last-committed score.

**At release**, `release.yml` additionally runs `memote report snapshot` (the
full suite, no tests skipped) and publishes it as an HTML report — see
[the latest snapshot](https://sysbiochalmers.github.io/yeast-GEM/release_report.html)
— which is the authoritative score for a given version; the pull-request
score above is a faster proxy for catching a regression early.

## 4. Files in this folder

Alphabetically. Regenerated by CI on each pull request and committed, so that a
change in a finding shows in the pull request's own diff rather than only in a
comment that later scrolls away:

| File | Holds |
|---|---|
| `qc_findings.md` | every entry behind every count, one section per check |
| `qc_metrics.tsv` | every structural, balance and annotation count, as `key<TAB>value` |
| `validation_findings.md` | the measurements behind each validation metric |
| `validation_metrics.tsv` | every validation metric, same format |
| `memote_score.md` | the [MEMOTE](https://memote.readthedocs.io) fast-subset score, total and per section |

Both findings files keep all their section headings whether or not the check found
anything: an empty section reads `_None._`. One file with a fixed set of
sections diffs far better than a CSV per check that appears and disappears —
a check going clean shows as rows vanishing under a heading that stays put,
rather than as a deleted file.

Regenerated at release by `increase_version.py`:

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
