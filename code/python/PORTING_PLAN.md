# yeast-GEM Python porting plan

Goal: provide Python equivalents of the MATLAB functions in `code/`, built on
[cobrapy](https://github.com/opencobra/cobrapy) and reusing
[ravengem](https://github.com/SysBioChalmers/ravengem) wherever it already
re-implements RAVEN functionality, so that the model can be loaded, curated,
saved and validated entirely from Python.

## Status

| Phase | Status | Notes |
|---|---|---|
| 1. Scaffold + comparator + reference fixture | **done** | `yeastgem` package importable, `read/write_yeast_model` ported, level-1 comparator + 15-test pytest suite passing, CI workflow in place (`matlab-reference-compare` job parked behind `if: false` until the reference bundle is seeded), reference-bundle scaffold + MATLAB regeneration stub. `code/io.py` is now a deprecated forwarding shim. |
| 2. Config-as-code refactor (both languages) | **done** | Data files (`data/yeastgem/ids.yml`, `data/conditions/{minimal_Y6,anaerobic,glycine_nitrogen,nitrogen_limitation}.yml`) created. MATLAB: `applyCondition.m`, `applyIDs.m`, `readYAML.m`; four legacy functions converted to one-line shims. Python: `yeastgem.config.load_ids()`, `yeastgem.conditions.apply()` cover prelude, cofactor-pseudoreaction edits, biomass-stoichiometry deltas, and bounds (33 tests passing). The `amino_acid_ratio` step in `anaerobic` is deferred to phase 4 (Tier 2) — calling `conditions.apply(model, 'anaerobic')` raises `NotImplementedError` with a clear pointer. **Verified end-to-end on MATLAB R2024b + RAVEN:** pre-refactor (`feat/anaerobic`) vs post-refactor (`feat/python-port`) is byte-identical on `rxns`, `mets`, `lb`, `ub`, `S` for all four conditions; Python `apply` matches MATLAB `lb`/`ub` for the three supported conditions; the two feasible SBML round-trips pass `yeastgem.compare`. Recipe + verification scripts in [tests/reference/README.md](tests/reference/README.md). |
| 3. Tier 1 — load/save parity (`commit_yeast_model`) | **done** | MATLAB: `saveYeastModel.m` → `commitYeastModel.m` (with deprecation shim); the cd-dance inside the pipeline replaced by `applyCondition('minimal_Y6')` / `applyCondition('anaerobic')`. Python: `yeastgem.io.commit_yeast_model` ports the release pipeline (apply minimal_Y6 → add SBO terms → SBML-validity gate → aerobic growth check → write SBML → ΔG CSVs → README update); `write_yeast_model` is a deprecated forwarding shim. Companion ports: `yeastgem.missing_fields.add_sbo_terms`, `load_delta_g`, `save_delta_g`. **20 new tests, 53 total passing.** `addSBOterms` faithfully replicates the legacy MATLAB pseudoreaction-loop bug (`for i=numel(model.rxns)` iterating only the last reaction); a fix is tracked as a future behaviour-change PR. Anaerobic growth check is deferred to phase 4 — emits a warning by default, raises `NotImplementedError` when `allow_no_growth=False`. Multi-format export (`.yml`/`.txt`/`.xlsx`/`.mat`) stays MATLAB-only for now; `.xml` is the contract. **Verified end-to-end:** pre-rename `saveYeastModel` vs post-rename `commitYeastModel` semantically equal; MATLAB `commitYeastModel` vs Python `commit_yeast_model` semantically equal. Recipe + verification driver in [tests/reference/README.md](tests/reference/README.md) and [`runPhase3.m`](tests/reference/runPhase3.m). |
| 3.5. Upstream restructure (raven-python + RAVEN) | **done** | Decision #1 reversed: generic helpers move upstream rather than living locally. Moved to raven-python (on `feat/yeast-gem-shared`): `raven_python.comparison.diff_models` + `DiffReport` (renamed from the local `compare_models`/`ComparisonReport`), `raven_python.annotation.{add_sbo_terms, load_delta_g_csv, save_delta_g_csv}`, `raven_python.conditions.{apply_condition, load_condition, set_reaction_bounds}`. Moved to RAVEN (on `feat/yeast-gem-shared`): `io/readYAML.m`, `core/applyCondition.m`. yeast-GEM now: depends on `raven-python` (git URL pinned to the feature branch), `yeastgem.compare`/`yeastgem.missing_fields`/`yeastgem.conditions` become thin wrappers that configure upstream defaults with yeast-specific data; MATLAB `code/readYAML.m` deleted, `code/applyCondition.m` → `code/applyYeastCondition.m` (handles the `amino_acid_ratio` pre-step then delegates to RAVEN). yeast-GEM uses the legacy `only_last_reaction_for_pseudo=True` flag on the upstream `add_sbo_terms` to stay byte-equivalent during the transition. **46 new raven-python tests + 46 yeast-GEM tests passing.** Verified: all 4 conditions byte-equivalent pre vs post restructure on MATLAB; Python `commit_yeast_model` (now through raven-python) semantically equal to MATLAB `commitYeastModel`. |
| 4. Tier 2 — biomass + conditions in Python | **done (core)** | Biomass subsystem moved upstream as `raven_python.biomass` (`BiomassConfig` + `BiomassComponent` + `sum_biomass` / `scale_biomass` / `rescale_pseudoreaction` / `set_gam`; 19 new tests on synthetic models). yeast-GEM ids.yml gained a `biomass_components` section; `yeastgem.biomass` exposes `sum_biomass`, `scale_biomass`, `rescale_pseudoreaction` (with the yeast `lipid` → backbone+chain aggregation), `set_gam` (auto-locates the NGAM reaction by name), and `change_amino_acid_ratio` (reads `data/physiology/aminoAcid_Bjorkeroth2020.tsv`). `yeastgem.conditions.apply` now handles `amino_acid_ratio` before delegating to upstream; `yeastgem.io.commit_yeast_model` runs the anaerobic growth check on a copy. **Verified** end-to-end on the real model: Python `conditions.apply('anaerobic')` produces SBML semantically equal to MATLAB `applyYeastCondition('anaerobic')`; Python `commit_yeast_model` (with anaerobic check active) produces SBML semantically equal to MATLAB `commitYeastModel`. 54 yeast-GEM tests + 38 new raven-python tests passing. **Deferred:** chemostat sweep + `fit_gam` (analysis/calibration, not part of the commit pipeline; tracked in UPSTREAM_CANDIDATES.md). |
| 5. Tier 3 — test suite | **done** | Ported the four ``code/modelTests/`` routines to ``yeastgem.model_tests``: ``growth`` (Tobias 2013 chemostat R² across 4 conditions), ``essential_genes`` (cobrapy ``single_gene_deletion`` + Stanford KO collection, returns ``EssentialGeneResult`` dataclass with accuracy / sensitivity / specificity / MCC), ``anaerobic_flux_predictions`` (Jouhten 2008 + Frick & Wittmann flux R² + mean relative error), ``plot_anaerobic`` (fermentation-product bar plot), ``find_duplicated_rxns`` (wrapper over the new ``raven_python.manipulation.find_duplicate_reactions``). Stanford ORF lists extracted from ``essentialGenes.m`` to ``data/essentialGenes/{inviable,verified}_orfs.txt`` so both languages read the same source. 7 new yeast-GEM tests + 6 new raven-python tests; full Python suite 61/61 passing. Verified vs MATLAB on the real model (`runPhase5Metrics.m`): growth R² matches at 1e-7; anaerobic flux R² and essential-gene accuracy/MCC match within 5e-3; single 1-gene difference in the essential-gene confusion matrix is a Gurobi/HiGHS solver-tolerance borderline at the 1e-6 ratio threshold. |
| 6. Tier 4 — curation framework | **done** | Generic `curateModelFromTables` engine moved to RAVEN (with `metPrefix` / `rxnPrefix` parameters defaulted to BiGG `M_`/`R_`); equivalent `raven_python.curation.{batch_curate, batch_curate_from_tsv}` in raven-python with the same schema (DataFrames + a `from_tsv` convenience). yeast-GEM keeps the user-facing `curateMetsRxnsGenes` MATLAB function as a 50-line shim that pins yeast's `s_`/`r_` prefixes and forwards upstream; the historical v8_*/v9_* curation scripts and `TEMPLATEcuration` keep working without change. New `yeastgem.curation.curate_mets_rxns_genes` Python entry point with the same prefix pinning. "Everything after the listed core columns is MIRIAM" — yeast-GEM's existing TSVs (12+10+9 MIRIAM columns) work unchanged. 13 new raven-python tests + 4 new yeast-GEM tests; full Python suite 65/65 passing. **MATLAB shim verified** to forward correctly (no-op call leaves the model unchanged). Direct MATLAB-vs-Python end-to-end parity check is blocked by pre-existing flakiness in the legacy `curateMetsRxnsGenes` (errors on the v8_6_3 VolPolyP schema and the v8_7_0 DBnewRxns pack); the Python implementation is more permissive than the legacy MATLAB on these edge cases. |
| 7. Docs + CI | **done** | Top-level README updated: "Contribution via Python is supported" section explaining the `yeastgem` + `raven-python` split + `saveYeastModel` → `commitYeastModel` rename. `code/python/README.md` rewritten with a getting-started block and an API map across the seven modules. CI workflow has three required jobs: `test` (matrix Python 3.10/3.11/3.12 + ruff + pytest), `parity-level-1-round-trip` (Python SBML read+write must round-trip the committed model semantically equal — `tests/ci/check_round_trip.py`), and `parity-level-2-metrics` (Python validation metrics must match the committed MATLAB reference within tolerance — `tests/ci/check_metrics.py` against `tests/reference/metrics.json`). Reference tolerances absorb the known Gurobi-vs-HiGHS solver drift (1 gene on the essential-gene confusion matrix, ≤ 5e-3 on R² metrics). Both parity scripts pass locally. **Prerequisite for CI to pass**: `raven-python`'s `feat/yeast-gem-shared` branch must be pushed to GitHub so the `pip install` URL dep resolves. |

## Design principles

- **Canonical object is `cobra.Model`** (ravengem's convention). No parallel
  RAVEN-style struct. RAVEN-only fields (`metDeltaG`, `rxnConfidenceScores`,
  SBO terms, MIRIAM) live in cobra `annotation`/`notes`.
- **Depend on the upstream toolboxes; do not duplicate.** *Revised after
  phase 3.* In MATLAB, yeast-GEM builds on RAVEN (`importModel`,
  `exportModel`, `solveLP`, the new `readYAML` / `applyCondition`, …). In
  Python, yeast-GEM builds on `raven-python` (which itself builds on
  cobrapy) — `diff_models`, `add_sbo_terms`, `apply_condition`, ΔG CSV
  helpers all live upstream. yeastgem keeps only the *yeast-specific
  configuration* of those generics: the data files under `data/`, the
  `applyYeastCondition` wrapper that handles the yeast-only
  `amino_acid_ratio` step, the legacy-bug-compat flag on
  `add_sbo_terms`, and the repo orchestration in `commit_yeast_model`
  (paths, README rewrite).
- **Package layout:** a proper importable package under `code/python/`
  (working name `yeastgem`), with flat submodules. The existing `code/io.py`
  is folded into `yeastgem.io`.
- **Parity over rewrite.** Each port must reproduce the MATLAB numeric result
  on the current model (validated against a MATLAB-produced reference) before
  it is considered done.

## Decisions taken

These four choices shape every section below:

1. **Upstream stance — generic helpers live upstream (revised after phase 3).**
   yeast-GEM depends on RAVEN (MATLAB) and `raven-python` (Python), and
   contributes the organism-agnostic helpers there rather than keeping
   them in-tree. Phase 3.5 moved the first batch (`diff_models`,
   `add_sbo_terms`, ΔG CSV persistence, `apply_condition`, `readYAML`).
   yeastgem keeps only the yeast-specific configuration / wrappers.
   See [UPSTREAM_CANDIDATES.md](UPSTREAM_CANDIDATES.md) for the
   remaining tracked items.

   *(Earlier, pre-phase-3 stance:)* keep everything in yeast-GEM for now;
   no new
   dependencies on RAVEN (MATLAB) or ravengem (Python) beyond what yeast-GEM
   already uses. Existing RAVEN usage in MATLAB stays; Python remains on a
   plain cobrapy baseline (no ravengem dependency added). Generic-looking
   functions that would *eventually* be useful upstream are tracked in
   [UPSTREAM_CANDIDATES.md](UPSTREAM_CANDIDATES.md) with their proposed API
   and rationale, but are implemented locally inside yeastgem. This decouples
   yeast-GEM's release schedule from those toolboxes' maturity (ravengem is
   pre-alpha), avoids designing upstream APIs around a single-organism use
   case, and lets the yeastgem implementations settle before they become
   anyone's public API.
2. **MATLAB direction — lock-step parity.** Every behavior change is made in
   MATLAB and Python in the same PR. CI verifies both toolchains produce the
   same model and the same analysis metrics. Idiomatic differences are allowed;
   observable outputs are not.
3. **Config-as-code — refactor in both languages now.** `minimal_Y6`,
   `anaerobicModel`, `glycineNitrogenSource`, `nitrogenLimitation`, the
   biomass/GAM yeast IDs, and the amino-acid ratios are demoted to data files
   (YAML/TSV) under `data/conditions/` and `data/yeastgem/`, consumed by thin
   loaders in both MATLAB and Python. Single source of truth.
4. **Validation contract — semantic + metric.** Two CI gates: (a) semantic
   equality on the committed model (rxns, mets, genes, S, bounds, GPRs, key
   annotations) and (b) metric parity within tolerance on the analyses
   (growth R², essential-gene accuracy, anaerobic flux R²). Detail in
   *Validation strategy* below.

## Proposed package layout

```
code/python/
  yeastgem/
    __init__.py
    io.py            # commit_yeast_model (release pipeline) + read helper;
                     #   replaces code/io.py. save_yeast_model kept as
                     #   deprecated shim for one release cycle.
    config.py        # loads data/yeastgem/ids.yml (biomass rxn, H+, GAM mets, …)
    conditions.py    # thin loader: apply_condition(model, name) reads
                     #   data/conditions/<name>.yml; covers media + presets
    missing_fields.py# loadDeltaG, saveDeltaG, addSBOterms, addConfidenceScores
    biomass.py       # sumBioMass, scaleBioMass, rescalePseudoReaction, GAM, AA-ratio
                     #   (tracked as upstream candidate; stays in yeastgem now)
    tests/           # ported model tests (growth, essentialGenes, …)
    curation.py      # curateMetsRxnsGenes + QC checks
  pyproject.toml     # deps: cobra, pandas, pyyaml, matplotlib (NO ravengem)
  UPSTREAM_CANDIDATES.md  # what should eventually move to ravengem / RAVEN
```

New data files (consumed by **both** MATLAB and Python loaders):

```
data/
  yeastgem/
    ids.yml                       # canonical yeast IDs (biomass rxn id,
                                  #   H+ met id, GAM cofactor met ids, …)
  conditions/
    minimal_Y6.yml                # exchange bounds for minimal media
    anaerobic.yml                 # bound + biomass changes for anaerobic
    glycine_nitrogen.yml          # glycine-as-N-source preset
    nitrogen_limitation.yml       # N-limitation preset
  physiology/
    aminoacid_Bjorkeroth2020.tsv  # (existing) AA ratios, aerobic + anaerobic
```

## Function triage

`.deprecated/` and the historical `v8_x_x` / `v9_x_x` curation scripts are
**out of scope** (frozen change-logs, not reusable code).

### Tier 1 — Core infrastructure (load/save path)

| MATLAB | Python target | Reuse (baseline only) |
|---|---|---|
| `loadYeastModel` | `io.read_yeast_model` (extend existing) | cobrapy SBML; local YAML reader if needed; `loadDeltaG` |
| `saveYeastModel` → `commit_yeast_model` | `io.commit_yeast_model` | cobrapy SBML/validator, local multi-format writer, local growth check |
| `loadDeltaG` / `saveDeltaG` | `missing_fields` | pandas; ΔG stored in cobra `annotation` |
| `addSBOterms` | `missing_fields.add_sbo_terms` | cobrapy reaction inspection (no ravengem) |
| `addConfidenceScores` | `missing_fields.add_confidence_scores` | pure logic |
| `minimal_Y6` | `conditions.apply('minimal_Y6')` (post-refactor) | pure bound-setting |

`saveYeastModel` also: SBML validity check (cobrapy validator), README
size/date update, `e-005`→`e-05` normalization — all reproducible in Python.

### Tier 2 — Biomass & condition modifiers (mostly pure math, low RAVEN coupling)

`sumBioMass`, `scaleBioMass`, `rescalePseudoReaction`, `changeGAM`, `fitGAM`,
`changeAminoAcidRatio`, `anaerobicModel`, `glycineNitrogenSource`,
`nitrogenLimitation`. These are stoichiometry/bound manipulations; they port
directly to cobrapy. `fitGAM` additionally needs the chemostat simulation
helper (see Tier 3).

### Tier 3 — Tests / analysis (the `increaseVersion` CI gate)

`growth`, `essentialGenes`, `findDuplicatedRxns`, `anaerobic_flux_predictions`,
`plotAnaerobic`, `anaerobiosis`.

- `essentialGenes` → cobrapy `single_gene_deletion` + Stanford KO comparison
  (cobrapy is baseline; this is the only "drop, just call the toolbox" case).
- `findDuplicatedRxns` → small local port (~20 lines: stoichiometry-signature
  detection over reactions). Logged as an upstream candidate.
- `growth`/`anaerobic_flux_predictions`/`plotAnaerobic` need a shared
  `simulate_chemostat` helper (FBA loop over dilution rates) + matplotlib;
  local now, upstream candidate.

### Tier 4 — Curation tooling (port the framework, not the one-offs)

- `curateMetsRxnsGenes` → `curation.curate` (TSV-driven), built on cobrapy
  `model.add_reactions` / `model.add_metabolites` / GPR setters + a local
  equation parser. Logged as an upstream candidate.
- QC checks `CheckBalanceforSce`, `CheckEnergyProduction`, `checkMetBalance`
  → local mass-balance helper (cobrapy's `Reaction.check_mass_balance`)
  + local FBA leak test for energy cycles. Logged as upstream candidates.
- `TEMPLATEcuration` → a Python notebook template.

### Deferred / likely skip

- `getEarlierModelVersion`, `increaseVersion` — git/release glue; port last or
  keep in MATLAB.
- `GetMNXID`, `mapIDsViaMNXref`, `addYMDBconcentrations` — depend on RAVEN
  MNXref tables / web downloads and overlap ravengem's reconstruction layer;
  low priority.

## Critical assessment: purpose, reuse, duplication

Reading the actual bodies (not just signatures) changes the picture. Most of
these functions fall into three nature-classes, and the right home depends on
the class:

- **Real generic algorithm** — organism-agnostic logic. Implemented locally in
  **yeastgem** for now; the abstraction + proposed upstream API is logged in
  [UPSTREAM_CANDIDATES.md](UPSTREAM_CANDIDATES.md) so the work isn't lost when
  ravengem (or RAVEN) is ready to absorb it.
- **Config-as-code** — a "function" whose entire body is hardcoded yeast
  reaction/metabolite IDs (`r_0501`, `s_0794`, exchange-ID lists). These are
  *data*, not algorithms; they belong in **yeastgem**, and several should be
  demoted to a data/config file rather than a function.
- **Thin wrapper / orchestration** — repo glue (paths, README edits, run order)
  or a one-liner over a toolbox call. Keep in **yeastgem**; for the one-liner
  case, call the underlying cobrapy primitive directly.

A recurring trap: functions that *look* generic (`sumBioMass`, `changeGAM`,
`rescalePseudoReaction`) are in fact laced with yeast-specific identifiers
(pseudoreaction names, the H⁺ metabolite `s_0794`, `r_4047`). The **algorithm**
is reusable; the **identifiers** are not. The local-yeastgem implementation
keeps this split internally: a generic helper inside the module, parameterised
by a small ID config, so the same code is upstream-ready when the time comes.

### Per-function verdict

| Function | Nature | Use freq | Verdict / home |
|---|---|---|---|
| `loadYeastModel` | thin shim — default path + ΔG fix-up for legacy formats | routine | **drop** (or ultra-thin shim); call `read_yaml_model`/`readYAMLmodel` on the default path directly. See *Load vs save asymmetry* below. |
| `saveYeastModel` → **rename to `commitYeastModel` / `commit_yeast_model`** | release pipeline — canonical state, validation gates, multi-format export, README metadata | routine | **yeastgem** — reframe as the commit function (run before `git commit`), not a wrapper. Keep `saveYeastModel` as a deprecated shim for one release cycle. See *Load vs save asymmetry* below. |
| `loadDeltaG` / `saveDeltaG` | annotation⇄CSV persistence | occasional | **yeastgem** (mechanism logged as upstream candidate) |
| `addSBOterms` | mostly generic rule-based annotation | routine | **yeastgem** (generic skeleton + yeast pseudoreaction tweak; upstream candidate) |
| `addConfidenceScores` | generic 0–3 scheme + yeast naming heuristics | occasional | **yeastgem** (heuristics too yeast-flavoured to upstream cleanly) |
| `minimal_Y6` | **config-as-code** (hardcoded exchange IDs) | routine | **yeastgem** — demote to a media data file |
| `sumBioMass` | real MW-weighted algorithm + yeast IDs | occasional (core to curation) | **yeastgem** (parameterised by yeast ID config internally; upstream candidate) |
| `scaleBioMass` / `rescalePseudoReaction` | real S-matrix rescale + yeast `s_0794` | occasional | **yeastgem** (same as `sumBioMass`; upstream candidate) |
| `changeGAM` | real but hardcoded met-name set | occasional | **yeastgem** (parameterised by cofactor met set; upstream candidate) |
| `fitGAM` | real fitting loop + yeast chemostat data + plot | rare | **yeastgem** (parameterised; upstream candidate with documented chemostat-data schema) |
| `changeAminoAcidRatio` | config-as-code (yeast data file, `r_4047`) | rare | **yeastgem** |
| `anaerobicModel` | config-as-code (large hardcoded yeast change set) | occasional | **yeastgem** — demote to condition data file |
| `glycineNitrogenSource` / `nitrogenLimitation` | **config-as-code**, hardcoded IDs | **rare** | **yeastgem** — demote to condition presets |
| `growth` | real chemostat validation + yeast data | routine (CI) | **yeastgem** (`simulate_chemostat` helper logged as upstream candidate) |
| `essentialGenes` | thin over cobrapy `single_gene_deletion` + yeast benchmark | routine (CI) | **yeastgem** glue over cobrapy (the only "use the toolbox directly" case) |
| `findDuplicatedRxns` | small generic algorithm | rare | **yeastgem** (~20-line local port; upstream candidate) |
| `anaerobic_flux_predictions` / `plotAnaerobic` / `anaerobiosis` | yeast validation + plots | occasional | **yeastgem** |
| `curateMetsRxnsGenes` | generic TSV→add/change engine | occasional | **yeastgem** (over cobrapy primitives; upstream candidate) |
| `CheckBalanceforSce` | thin wrapper over elemental balance | occasional | **yeastgem** — local helper using cobrapy mass-balance; thin results table |
| `CheckEnergyProduction` | real energy-cycle leak test (generic technique) | occasional | **yeastgem** (upstream candidate) |
| `checkMetBalance` | display helper (print rxns touching a met) | rare | **yeastgem** util, or a cobrapy idiom inline |
| `getEarlierModelVersion` / `increaseVersion` | git/release glue | rare | **yeastgem** (repo-specific) |
| `GetMNXID` / `mapIDsViaMNXref` / `addYMDBconcentrations` | ID mapping / downloads | rare | defer (overlap with future upstream reconstruction work) |

### Load vs save asymmetry

`loadYeastModel` and `saveYeastModel` look like a matched pair but they are
not:

- **`loadYeastModel`** decomposes into (i) a default path constant, (ii) a
  format dispatch already handled by the toolbox (`.yml` →
  `readYAMLmodel`/`read_yaml_model`, else SBML), and (iii) a conditional
  `loadDeltaG` fix-up **only on the non-YAML branch**. On the canonical YAML
  path it does literally nothing beyond calling the toolbox loader. → It is
  not a meaningful abstraction; **drop it** (or keep as a 3-line shim and
  deprecate) and document the one-liner alternative in the README. The
  legacy-format ΔG fix-up is a separate, explicitly-called step
  (`loadDeltaG`), not something to hide behind a "load" name.
- **`saveYeastModel`** is a release pipeline, not a wrapper: enforce
  canonical state (`minimal_Y6`, `addSBOterms`) → validation gates (SBML
  valid, aerobic + anaerobic growth) → multi-format export → repo metadata
  (README size/date) → MATLAB `e-005`→`e-05` text patch. Most of this is
  yeast/repo-specific policy that `exportModel`/`write_yaml_model` does not
  and should not do. **Keep it, but rename to `commitYeastModel` /
  `commit_yeast_model`** — the current name "save" implies a casual write,
  while the function is the heavy ceremony you run before opening a curation
  PR. The new name signals the workflow: *run this, then `git commit`*. The
  docstring is explicit that the function does not perform the git commit
  itself. In Python, two pieces drop out: the RAVEN-format conversion and
  the `e-005` patch (Python's SBML writers don't produce that string).

  **Deprecation:** keep `saveYeastModel` as a 3-line shim for one release
  cycle that emits a deprecation warning and forwards to `commitYeastModel`,
  in both MATLAB and Python. This prevents breaking external curation
  scripts that call the current name. Remove the shim at the next minor
  version bump after the rename ships.

Open question to verify early: does `ravengem.io.read_yaml_model` accept
yeast-GEM's RAVEN-style YAML (with top-level `metDeltaG`/`rxnDeltaG` arrays),
or does it expect cobrapy's YAML convention? If the latter, either teach
ravengem to read RAVEN YAML or migrate ΔG into `annotation` in the committed
file — this is the prerequisite for "drop loadYeastModel" to be clean.

### Duplication to eliminate (do not port)

Because we're not adding a ravengem dependency, the only "drop, use the
toolbox" case is cobrapy:

- `essentialGenes` deletion loop ≈ cobrapy `single_gene_deletion` — keep only
  the Stanford-KO comparison around the cobrapy call.
- `scaleBioMass` / `rescalePseudoReaction` / `sumBioMass` share a single
  biomass subsystem — port them together as one module, not three independent
  copies.

Functions that are *also* implemented in ravengem (e.g. `findDuplicatedRxns`,
`CheckBalanceforSce`'s elemental balance) are still re-implemented locally
here rather than reused; their existence in ravengem just makes them stronger
upstream candidates. See [UPSTREAM_CANDIDATES.md](UPSTREAM_CANDIDATES.md).

### Net effect on scope

Roughly half of the "functions to port" are not really algorithms to
re-implement: they are either yeast **data** (`minimal_Y6`, `glycine…`,
`nitrogenLimitation`, `anaerobicModel`, `changeAminoAcidRatio`) best expressed
as config, or one cobrapy call away. The genuinely new Python code yeastgem
owns is: the biomass subsystem, the chemostat helper, the `io` orchestration
(`commit_yeast_model`), and yeast-specific validation. With the no-new-deps
policy, yeastgem also owns local re-implementations of the generic algorithms
(duplicate detection, batch curation, energy-cycle test, mass-balance helper)
that would otherwise come from ravengem — these are kept as small, well-scoped
helpers so they can graduate to ravengem later without rewrites.

## MATLAB strategy and dual-language maintenance

The two ecosystems are **not symmetric**, but the chosen policy treats them
the same way: don't deepen the toolbox dependency.

- **Python:** everyone uses cobrapy. ravengem (built on cobrapy) is the
  natural eventual upstream, but it is pre-alpha and adding it as a
  dependency now would couple yeast-GEM's release schedule to ravengem's
  churn. Stay on cobrapy alone; track ravengem candidates in
  [UPSTREAM_CANDIDATES.md](UPSTREAM_CANDIDATES.md).
- **MATLAB:** yeast-GEM already *hard-depends* on RAVEN (`importModel`,
  `exportModel`, `ravenCobraWrapper`, `solveLP`, `getElementalBalance`,
  `constructEquations`, `getExchangeRxns`, `getTransportRxns`,
  `findGeneDeletions`). That existing baseline stays. We do not add *new*
  RAVEN-only requirements (i.e. don't move yeast-GEM helpers into RAVEN
  even when they look generic — the upstream candidates list covers both
  ecosystems).

### Should MATLAB functions move into RAVEN?

**Not now.** The generic primitives yeast-GEM needs already exist in RAVEN
(`getElementalBalance`, `constructEquations`, `findGeneDeletions`, …) and
yeast-GEM uses them via its existing baseline; we keep that as-is. The
organism-agnostic helpers yeast-GEM itself owns (`changeGAM`, biomass math,
energy-cycle test, `findDuplicatedRxns`, batch curation) are kept in
yeast-GEM for now and tracked in
[UPSTREAM_CANDIDATES.md](UPSTREAM_CANDIDATES.md). Moving them into RAVEN is
deferred — the rationale, proposed API, and what would trigger the move all
live in the candidates document.

### If we refactor in Python, must MATLAB change too?

Distinguish two kinds of change:

- **Behavioral / numeric** (changes the model that gets written, or a reported
  metric like growth R² or essential-gene accuracy): **must stay in lock-step**
  across MATLAB and Python. The committed `yeast-GEM.yml`/`.xml` must be
  reproducible and identical from either toolchain; a Python-only behavior
  change would silently fork the model. These changes are made (and reviewed) in
  both languages in the same PR, with a cross-check that both produce the same
  model/metrics.
- **Structural / idiomatic** (renaming, splitting a module, replacing a loop
  with a vectorised call): **need not be mirrored.** Each ecosystem may use
  its own idioms; only the observable result is contracted.

### Maintenance policy: lock-step parity

The chosen policy. Concretely:

- **Every PR that changes behavior touches both languages.** A Python-only or
  MATLAB-only behavioral change is a CI failure, not a code-review note.
- **CI gate per PR** runs MATLAB save → reference artifact; Python save →
  candidate artifact; a comparator enforces *semantic equality* on the model
  and *metric parity within tolerance* on the analyses (see *Validation
  strategy*). Detail there.
- **Idiomatic differences are free.** Splitting a module, vectorising a loop,
  renaming a private helper — none of these need a MATLAB counterpart, as
  long as outputs are unchanged.
- **One language "owns" producing the committed artifact** at a time to avoid
  races (initially MATLAB, until Python is at full parity). Both languages
  must be able to produce it; only one writes the committed file per release.

The hard invariant: **the model artifact is the contract.** Both toolchains
must read and write the same yeast-GEM files and agree on the validation
metrics.

## Upstream candidates (tracked separately)

Per Decision #1, no function is upstreamed in this porting effort. Generic
algorithms are implemented locally in yeastgem (and stay where they are in
MATLAB), and the candidates for future upstreaming — with their proposed
APIs, rationale, and triggers — live in
[UPSTREAM_CANDIDATES.md](UPSTREAM_CANDIDATES.md).

A summary, for orientation only:

- **High-confidence upstream candidates** (clean abstraction, broadly useful):
  biomass subsystem (`sumBioMass`/`scaleBioMass`/`rescalePseudoReaction`),
  `changeGAM`, `fitGAM` + chemostat sweep helper, `CheckEnergyProduction`,
  TSV-driven batch curation, duplicate-reaction detector,
  cross-language model comparator (used by CI).
- **Boundary cases** (probably worth upstreaming, may need API refinement):
  `addSBOterms`, ΔG annotation⇄CSV persistence.
- **Unlikely to upstream**: `addConfidenceScores` (heuristics too
  yeast-flavoured), all explicit yeast benchmarks and physiology, repo
  orchestration (`commitYeastModel`, `loadYeastModel`).
- **Drop, not upstream**: cobrapy already provides
  `single_gene_deletion`/mass balance — yeastgem just calls these directly.

A function being on the "tracked" list does not gate its yeastgem
implementation; it just means we keep an eye on the API shape so it can move
cleanly later.

## Config-as-code refactor (both languages)

Per the chosen stance, the yeast condition presets and the yeast ID set become
data files, with thin loaders in both MATLAB and Python. This is a real
workstream that also touches the existing MATLAB code — not a Python-only
change.

### Files

- `data/yeastgem/ids.yml` — canonical yeast IDs that are *parameters* to
  organism-agnostic algorithms: biomass reaction id, H⁺ metabolite id (`s_0794`),
  GAM cofactor mets, protein/carb/lipid/etc. pseudoreaction ids, protein
  reaction id (`r_4047`).
- `data/conditions/minimal_Y6.yml` — exchange bounds for minimal media (the
  body of `minimal_Y6.m`).
- `data/conditions/anaerobic.yml` — bound changes + heme-a removal + AA-ratio
  switch + FAD recycling, all expressed as a structured diff.
- `data/conditions/glycine_nitrogen.yml`, `data/conditions/nitrogen_limitation.yml` —
  the existing 3–5-line bound flips, as data.
- `data/physiology/aminoacid_Bjorkeroth2020.tsv` — already exists; the AA-ratio
  function becomes a thin loader.

### Loader API (mirrored in both languages)

- MATLAB: `model = applyCondition(model, 'anaerobic')` reads
  `data/conditions/anaerobic.yml` and applies the diff. `applyMedia`,
  `applyIDs` similarly. The current functions (`minimal_Y6`, `anaerobicModel`,
  `glycineNitrogenSource`, `nitrogenLimitation`) become 3-line shims that call
  `applyCondition` with a fixed name — kept for backwards compatibility, with
  a deprecation note.
- Python: `yeastgem.conditions.apply(model, 'anaerobic')` does the same. Same
  YAML files, same semantics.

### Format sketch

```yaml
# data/conditions/anaerobic.yml
name: anaerobic
description: Convert aerobic yeast-GEM to anaerobic.
biomass:
  amino_acid_ratio: anaerobic   # column selector for aminoacid_Bjorkeroth2020.tsv
  remove_cofactors:
    - heme a
bounds:
  - rxn: r_1992      # O2 exchange
    lb: 0
  # ergosterol/fatty acid uptake, MDH2/IDP2 block, FAD recycling …
```

The committed model artifact stays the contract: after the refactor, MATLAB
and Python loaders applied to the same condition file must produce the same
model (verified by the CI gate).

## Phased roadmap

1. **Scaffold + comparator + reference fixture.** Create `code/python/` with
   the `yeastgem` package and `pyproject.toml` (deps: cobra, pandas, pyyaml,
   matplotlib — no ravengem); fold existing `io.py` in. Build the model
   comparator (level-1) and the reference-bundle generator. CI infrastructure
   for lock-step parity is in place before any function is ported. Create
   [UPSTREAM_CANDIDATES.md](UPSTREAM_CANDIDATES.md) seeded with the candidates
   identified during this planning.
2. **Config-as-code refactor (both languages).** Create
   `data/yeastgem/ids.yml` and the four `data/conditions/*.yml` files; add
   `applyCondition` / `applyIDs` loaders in MATLAB; existing functions
   (`minimal_Y6.m`, `anaerobicModel.m`, `glycineNitrogenSource.m`,
   `nitrogenLimitation.m`) become 3-line shims. CI verifies the model is
   unchanged by the refactor. **This is a MATLAB PR before the Python work
   really begins**, so the Python side has a clean target to load.
3. **Tier 1 — load/save parity.** Python `read_yeast_model`/`commit_yeast_model`
   handle YAML + ΔG + SBO + confidence scores + the condition loaders +
   README/version updates. Both level-1 and level-2 CI gates active. Includes
   the `saveYeastModel` → `commitYeastModel` rename (both languages, with
   deprecation shim).
4. **Tier 2 — biomass + conditions in Python.** Local `yeastgem.biomass`
   (sum/scale/rescale/GAM/AA-ratio); `fitGAM` orchestrates a local
   `simulate_chemostat` helper. Validated against the metric tolerances.
5. **Tier 3 — test suite.** `growth`, `essentialGenes`,
   `anaerobic_flux_predictions`, `plotAnaerobic` in Python; the
   `increaseVersion` PR gate now runs in either language.
6. **Tier 4 — curation framework.** Local `yeastgem.curation` over cobrapy
   primitives; `TEMPLATEcuration` as a Python notebook.
7. **Docs + CI.** Update README's "Python contribution not yet functional"
   note; add the Python CI job. Cross-language parity job becomes required.

Out-of-band, not gating any phase above: the upstream candidates document is
revisited whenever a yeastgem helper is touched, to keep its API shape
aligned with the proposed upstream signature. Actual upstreaming happens in
ravengem / RAVEN repos when those projects are ready to receive it, on their
schedule, not this one.

## Validation strategy

Two-level CI gate per PR, matching the lock-step parity policy.

### Level 1 — semantic model equality

A comparator (in `code/python/yeastgem/tests/compare_models.py`, mirrored by a
MATLAB equivalent) checks, on the same model loaded by each toolchain:

- set of reaction ids, metabolite ids, gene ids (exact)
- S matrix (within numerical tolerance, e.g. `1e-9`)
- lower/upper bounds, objective coefficients (exact)
- GPR rules (parsed and normalised; insensitive to whitespace / operator
  spelling)
- key annotation fields: SBO terms, MIRIAM, ΔG, confidence scores

Formatting differences (key ordering, whitespace, float repr) are explicitly
**not** failures. The comparator is single-source — used both for PR CI and
for the release "models are identical" check.

### Level 2 — metric parity within tolerance

The analyses in `modelTests` are run in both languages on the same condition
preset; metrics must agree within ε:

| Metric | Tolerance (proposal) |
|---|---|
| Aerobic growth (optimal flux through `r_4041`) | `1e-6` |
| Chemostat fit R² (`growth` over 4 conditions) | `1e-4` |
| Essential-gene confusion matrix (TP/TN/FP/FN counts) | exact |
| Anaerobic flux R² (`anaerobic_flux_predictions`) | `1e-4` |
| Biomass component fractions (`sumBioMass` X, P, C, R, D, L, I, F) | `1e-6` |

Tolerances are revisited once we see real solver-vs-solver drift between
MATLAB (Gurobi via RAVEN) and Python (Gurobi/HiGHS via cobrapy).

### Reference fixtures

A MATLAB job (run at the start of each release cycle, not per PR) regenerates
a reference bundle in `code/python/tests/reference/`: the committed model
loaded + saved by MATLAB, plus the metric values above. PR CI compares the
Python toolchain against this bundle; any drift forces a deliberate refresh
of the reference, reviewed in the PR.

## Open items to confirm later

- **Package name: `yeastgem`** (settled). **PyPI: deferred** — yeastgem is
  tightly coupled to the in-repo model and condition files, so a PyPI package
  would either pin a model version per release (awkward) or be unusable
  without a clone. The repo is the distribution, mirroring the MATLAB side.
  `pyproject.toml` still ships so `pip install -e code/python/` works from a
  clone. Revisit if external demand appears.
- **Plot library: matplotlib** (settled). The plots are static publication
  artifacts saved into `data/testResults/`; matplotlib gives PDF/PNG output,
  the lightest dependency, deterministic rendering for CI image-diffs, and
  visual continuity with the current MATLAB figures. Plotly may be added
  later as an optional `yeastgem[notebooks]` extra if a curation notebook
  wants interactivity — not as a core dep.
- Exact tolerances in the metric-parity table (start with the values above,
  loosen only with evidence of solver drift).
- Future upstreaming of any tracked candidate — handled in
  [UPSTREAM_CANDIDATES.md](UPSTREAM_CANDIDATES.md), not here.
