# Upstream candidates

This document tracks helpers that should live upstream
(**raven-python** for Python, **RAVEN** for MATLAB) so yeast-GEM does not
duplicate organism-agnostic logic. The "done" section records what
already moved; the "pending" section records what remains and the API
shape we want it to land with.

**Decision #1 was revised after phase 3** (see [PORTING_PLAN.md](PORTING_PLAN.md)):
generic helpers move upstream now and yeast-GEM declares a dependency
on those toolboxes, rather than implementing locally first.

## Done — moved upstream in phase 3.5

The following were extracted from yeast-GEM into upstream branches
`feat/yeast-gem-shared` on both repos:

| yeast-GEM module (was) | Upstream home | Public API |
|---|---|---|
| `yeastgem.compare.compare_models` / `ComparisonReport` | `raven_python.comparison` | `diff_models`, `DiffReport` |
| `yeastgem.missing_fields.add_sbo_terms` | `raven_python.annotation.sbo` | `add_sbo_terms` (with `only_last_reaction_for_pseudo` legacy flag) |
| `yeastgem.missing_fields.{load,save}_delta_g` mechanism | `raven_python.annotation.delta_g` | `load_delta_g_csv`, `save_delta_g_csv` (column / note-key params) |
| `yeastgem.conditions.apply` internals (prelude / cofactor / biomass-delta / bounds) | `raven_python.conditions` | `apply_condition`, `load_condition`, `set_reaction_bounds` |
| `code/readYAML.m` | RAVEN `io/readYAML.m` | unchanged signature |
| `code/applyCondition.m` (generic core) | RAVEN `core/applyCondition.m` | takes YAML path or struct |

yeast-GEM now keeps:
- `yeastgem.compare` — re-export of the upstream `diff_models` under
  the historical names `compare_models` / `ComparisonReport`.
- `yeastgem.missing_fields` — thin wrappers passing the yeast CSV
  paths (`data/databases/model_metDeltaG.csv` etc.) and the legacy
  `only_last_reaction_for_pseudo=True` bug-compat flag.
- `yeastgem.conditions.apply` — resolves a name to `data/conditions/<name>.yml`,
  raises `NotImplementedError` on `amino_acid_ratio` (Tier 2), then
  delegates to upstream.
- `code/applyYeastCondition.m` — same shape, in MATLAB.

## Pending — not yet moved

Listing a function here does **not** create an upstream PR or a
dependency. It records the abstraction we'd want, the proposed
signature, the rationale, and what would trigger actually doing the
move. When a candidate's local implementation is touched, this
document is updated alongside so the API shape stays aligned with the
proposed upstream signature.

## How to read each entry

- **Current location** — the file(s) in yeast-GEM that hold the local
  implementation today.
- **Proposed upstream signature** — what the API would look like in
  ravengem / RAVEN, parameterised so it's not yeast-specific.
- **Why local for now** — the reason we're not upstreaming yet (usually:
  upstream not ready, API not yet validated by a second use case, or risk
  of premature abstraction).
- **Trigger to upstream** — the concrete condition under which we'd
  actually do the move.

---

## High-confidence candidates

These have clean abstractions and are clearly useful beyond yeast.

### Biomass subsystem

- **Current location (MATLAB):** `code/otherChanges/sumBioMass.m`,
  `scaleBioMass.m`, `rescalePseudoReaction.m`.
- **Current location (Python):** `code/python/yeastgem/biomass.py`.
- **Proposed upstream signature (ravengem):**
  ```python
  # ravengem.biomass
  @dataclass
  class BiomassConfig:
      biomass_rxn: str           # e.g. "r_4041"
      proton_met: str            # e.g. "s_0794"
      components: dict[str, str] # component → pseudoreaction id
                                 # e.g. {"protein": "r_4047", "carbohydrate": ...}

  def sum_biomass(model, cfg: BiomassConfig) -> dict[str, float]: ...
  def rescale_pseudoreaction(model, cfg, component, factor): ...
  def scale_biomass(model, cfg, component, new_value, balance_out=None): ...
  ```
- **Proposed upstream signature (RAVEN):** mirror as
  `sumBioMass(model, biomassConfig)` etc., where `biomassConfig` is a struct
  with the same fields.
- **Why local for now:** ravengem is pre-alpha; the `BiomassConfig` shape
  has not been exercised on a second model.
- **Trigger to upstream:** at least one other GEM project (e.g. Human-GEM,
  yeast-pcGEM) successfully uses the same `BiomassConfig` shape, or ravengem
  reaches a stable release.

### `changeGAM` → `set_gam`

- **Current location:** `code/otherChanges/changeGAM.m`,
  `code/python/yeastgem/biomass.py`.
- **Proposed upstream signature (ravengem):**
  ```python
  def set_gam(
      model,
      value: float,
      *,
      biomass_rxn: str,
      cofactor_mets: list[str],   # yeast: ["ATP","ADP","H2O","H+","phosphate"]
      ngam_rxn: str | None = None,
      ngam_value: float | None = None,
  ) -> None: ...
  ```
- **Proposed upstream signature (RAVEN):** same parameters as a MATLAB
  function.
- **Why local for now:** trivial to implement locally; abstraction is clean
  but not urgent.
- **Trigger to upstream:** bundled with the biomass subsystem move.

### `fitGAM` → `fit_gam`

- **Current location:** `code/otherChanges/fitGAM.m`, `code/python/yeastgem/biomass.py`.
- **Proposed upstream signature (ravengem):**
  ```python
  def fit_gam(
      model,
      chemostat_data: pandas.DataFrame,
      *,
      biomass_rxn: str,
      substrate_exchange_id: str,
      cofactor_mets: list[str],
      gam_bounds: tuple[float, float] = (30, 70),
      refinement: tuple[float, float, float] = (5, 1, 0.1),
  ) -> tuple[Model, FitResult]: ...
  ```
  Required `chemostat_data` schema (documented as part of the public API):
  - columns: `dilution_rate`, `<substrate>_uptake`, optionally
    `<product>_secretion` for one or more products
  - units: mmol gDW⁻¹ h⁻¹ for fluxes, h⁻¹ for dilution rate
- **Proposed upstream signature (RAVEN):** mirror with the same TSV schema.
- **Why local for now:** the data-schema contract needs at least one second
  validation against a non-yeast chemostat dataset.
- **Trigger to upstream:** schema validated against a second organism's
  chemostat data; ravengem ready.

### `simulate_chemostat` (helper extracted from `growth` / `fitGAM` / `anaerobic_flux_predictions`)

- **Current location (Python):** private helper in `yeastgem/biomass.py`
  initially; promote to `yeastgem.analysis.chemostat` once Tier 3 lands.
- **Proposed upstream signature (ravengem):**
  ```python
  def chemostat_sweep(
      model,
      dilution_rates: Iterable[float],
      *,
      biomass_rxn: str,
      substrate_exchange_id: str,
      tracked_exchanges: list[str] | None = None,
  ) -> pandas.DataFrame: ...
  ```
  Returns a DataFrame indexed by dilution rate with columns for substrate
  uptake, biomass flux, and any tracked product fluxes.
- **Why local for now:** internal helper for `fit_gam` / `growth`; we want to
  ship those first and only then extract.
- **Trigger to upstream:** `fit_gam` graduates upstream (it depends on this).

### `CheckEnergyProduction` → `check_energy_cycles`

- **Current location:** `code/modelCuration/CheckEnergyProduction.m`,
  `code/python/yeastgem/curation.py`.
- **Proposed upstream signature (ravengem):**
  ```python
  def check_energy_cycles(
      model,
      *,
      atp_id: str,
      adp_id: str,
      pi_id: str,
      h2o_id: str,
      h_id: str,
      nadh_id: str,
      nad_id: str,
      atp_threshold: float = 360,
      nadh_threshold: float = 120,
  ) -> EnergyCycleReport: ...
  ```
- **Why local for now:** the threshold defaults are yeast-tuned; needs a
  second organism's calibration before being a public API.
- **Trigger to upstream:** thresholds confirmed reasonable on a second GEM.

### `curateMetsRxnsGenes` → `batch_curate`

- **Current location:** `code/modelCuration/curateMetsRxnsGenes.m`,
  `code/python/yeastgem/curation.py`.
- **Proposed upstream signature (ravengem):**
  ```python
  def batch_curate(
      model,
      *,
      mets_tsv: str | Path | None = None,
      rxns_tsv: str | Path | None = None,
      genes_tsv: str | Path | None = None,
      rxns_coeffs_tsv: str | Path | None = None,
  ) -> Model: ...
  ```
  Each TSV has a documented schema (columns and meaning); function
  add/updates entities by matching on `name[compartment]` for metabolites,
  stoichiometry for reactions, and gene name for genes.
- **Why local for now:** TSV schemas are project conventions; we want to
  exercise them on a few yeast curation rounds before locking the columns.
- **Trigger to upstream:** TSV schemas stable for two consecutive yeast-GEM
  minor releases.

### `findDuplicatedRxns` → `find_duplicate_reactions`

- **Current location:** `code/modelTests/findDuplicatedRxns.m`,
  `code/python/yeastgem/tests/utils.py` (or similar).
- **Proposed upstream signature (ravengem):**
  ```python
  def find_duplicate_reactions(model, *, ignore_direction: bool = True) -> list[tuple[Reaction, Reaction]]: ...
  ```
  ravengem already has `remove_duplicate_reactions` — this would either
  expose a `detect_only` mode or live as a separate detector.
- **Why local for now:** trivial port; not worth coordinating across repos.
- **Trigger to upstream:** ravengem PR window opens for utils.

### Cross-language model comparator

- **Current location:** `code/python/yeastgem/tests/compare_models.py` (+ a
  MATLAB twin), built for the CI level-1 gate.
- **Proposed upstream signature (ravengem):**
  ```python
  def compare_models(
      a: cobra.Model,
      b: cobra.Model,
      *,
      stoichiometry_tol: float = 1e-9,
      ignore_annotations: Iterable[str] = (),
  ) -> ComparisonReport: ...
  ```
- **Why local for now:** comparator is shaped by our specific CI needs
  (which annotation fields are key, what "normalised GPR" means here).
- **Trigger to upstream:** at least one other project asks for cross-tool
  model diffing; comparator API stable.

---

## Boundary cases

Probably worth upstreaming eventually, but the abstraction is less obvious.

### `addSBOterms` → `add_sbo_terms`

- **Current location:** `code/missingFields/addSBOterms.m`,
  `code/python/yeastgem/missing_fields.py`.
- **Sketch:** generic SBO-by-reaction-type assignment (exchange/sink/demand/
  transport/biochemical) is organism-agnostic; the yeast-specific bit is
  identifying biomass and pseudoreactions. A `pseudoreaction_classifier`
  callback would isolate the yeast logic.
- **Why local for now:** the callback abstraction is plausible but
  speculative; no second consumer to validate it.
- **Trigger to upstream:** a second GEM wants the same SBO assignment.

### ΔG annotation⇄CSV persistence (`loadDeltaG` / `saveDeltaG`)

- **Current location:** `code/missingFields/loadDeltaG.m`, `saveDeltaG.m`,
  `code/python/yeastgem/missing_fields.py`.
- **Sketch:** generic "persist a named annotation key across mets/rxns to CSV
  and reload it." Useful for any project that stores a model-adjacent
  numeric annotation outside the SBML file.
- **Why local for now:** column schema is a yeast-GEM convention; no second
  consumer.
- **Trigger to upstream:** another project asks for the same CSV format,
  *or* this becomes a recurring pattern across multiple SysBioChalmers GEMs.

---

## Explicitly not upstream candidates

For traceability — so we don't accidentally re-litigate these:

- `addConfidenceScores` — heuristics mention "SLIME rxn", "pseudoreaction",
  "Biolog update"; the function would have to be watered down to "score by
  rule" to upstream cleanly, and that's a worse function.
- All yeast benchmark + physiology code: `growth`, `essentialGenes`,
  `anaerobic_flux_predictions`, `plotAnaerobic`, `anaerobiosis`,
  `changeAminoAcidRatio`.
- Repo orchestration: `loadYeastModel` (drop entirely or keep as
  default-path shim), `commitYeastModel`, `getEarlierModelVersion`,
  `increaseVersion`.
- Data-driven condition presets (`minimal_Y6`, `anaerobicModel`,
  `glycineNitrogenSource`, `nitrogenLimitation`) — these are *data* under
  `data/conditions/`, not functions.

## Cobrapy-direct (no upstream needed)

These are not "upstream candidates" because cobrapy already covers them and
yeastgem just calls cobrapy directly:

- gene-deletion loop in `essentialGenes` → `cobra.flux_analysis.single_gene_deletion`
- mass-balance check inside `CheckBalanceforSce` (the wrapper drops; we keep
  the yeast-specific results-table formatting) → `Reaction.check_mass_balance`

## Maintenance

When a yeastgem function listed above is touched:

1. Update the local implementation.
2. Re-read its entry here; if the **proposed upstream signature** would now
   be different (better/cleaner), update it.
3. Re-evaluate the **trigger** — has it become more or less plausible?

When ravengem (or RAVEN) absorbs a candidate:

1. Move the entry out of this document into a "Done" section (or just
   delete) and reference the upstream PR.
2. Update [PORTING_PLAN.md](PORTING_PLAN.md) if the yeast-GEM code now calls
   the upstream version instead of the local one — this is a behavioral
   change and goes through the normal lock-step parity CI.
