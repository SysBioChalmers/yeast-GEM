"""Annotation and thermodynamic-field helpers for yeast-GEM.

Currently houses:
  - `add_sbo_terms`   — port of code/missingFields/addSBOterms.m
  - `load_delta_g`    — port of code/missingFields/loadDeltaG.m
  - `save_delta_g`    — port of code/missingFields/saveDeltaG.m

`add_confidence_scores` will land alongside these later (it is not part
of the commit pipeline, so it does not block phase 3).
"""
from __future__ import annotations

import math
from pathlib import Path

import cobra
import pandas as pd

from yeastgem.io import REPO_PATH

_DELTAG_DIR = REPO_PATH / "data" / "databases"
_MET_CSV = _DELTAG_DIR / "model_metDeltaG.csv"
_RXN_CSV = _DELTAG_DIR / "model_rxnDeltaG.csv"

# Keys used to store ΔG in cobra's `notes` dict (per Decision: RAVEN-only
# fields live in `annotation`/`notes`). String storage round-trips through
# SBML cleanly.
_DELTA_G_NOTE_KEY = "deltaG"


# --- SBO terms --------------------------------------------------------

_BIOMASS_METNAMES: frozenset[str] = frozenset(
    {"biomass", "DNA", "RNA", "protein", "carbohydrate",
     "lipid", "cofactor", "ion"}
)


def add_sbo_terms(model: cobra.Model) -> cobra.Model:
    """Assign SBO terms to metabolites and reactions in-place.

    Port of `code/missingFields/addSBOterms.m`. Mirrors the legacy MATLAB
    logic, including its known idiosyncrasies (see *Limitations* below).
    SBO terms are written into ``annotation['sbo']`` only if absent
    (`'fill'` semantic from RAVEN's editMiriam).

    Metabolite SBO assignment
        - SBO:0000649 (Biomass) for the biomass / DNA / RNA / protein /
          carbohydrate / lipid / cofactor / ion pseudo-metabolites and
          their ``backbone`` / ``chain`` siblings.
        - SBO:0000247 (Simple chemical) otherwise.

    Reaction SBO assignment
        - SBO:0000176 (Metabolic reaction) default.
        - Single-reactant reactions (exchange / sink / demand):
            * extracellular → SBO:0000627 (exchange)
            * sum(S) < 0 → SBO:0000632 (sink)
            * else → SBO:0000628 (demand)
        - Transport reactions → SBO:0000655.
        - The biomass pseudoreaction → SBO:0000629.
        - The NGAM reaction → SBO:0000630.
        - Other pseudoreactions / SLIME reactions → SBO:0000395.

    Limitations (preserved for lock-step parity with MATLAB):
        - The legacy MATLAB loop for the pseudoreaction SBO assignment
          iterates only over the last reaction
          (`for i=numel(model.rxns)` rather than `1:numel(...)`), so
          pseudoreaction-specific SBOs are not actually applied. This
          port mirrors that bug. Fixing it is tracked as a separate
          behaviour-change PR; both languages must move together.
    """
    _assign_metabolite_sbo(model)
    _assign_reaction_sbo(model)
    return model


def _assign_metabolite_sbo(model: cobra.Model) -> None:
    for met in model.metabolites:
        if met.name in _BIOMASS_METNAMES or met.name.endswith(
            (" backbone", " chain")
        ):
            sbo = "SBO:0000649"
        else:
            sbo = "SBO:0000247"
        _fill_sbo(met, sbo)


def _assign_reaction_sbo(model: cobra.Model) -> None:
    transport_set = _transport_reaction_ids(model)

    rxns = list(model.reactions)
    for rxn in rxns:
        sbo = "SBO:0000176"  # default: metabolic reaction
        if len(rxn.metabolites) == 1:
            (met,) = rxn.metabolites
            coef = rxn.metabolites[met]
            if met.compartment == "e" or model.compartments.get(
                met.compartment
            ) == "extracellular":
                sbo = "SBO:0000627"
            elif coef < 0:
                sbo = "SBO:0000632"
            else:
                sbo = "SBO:0000628"
        if rxn.id in transport_set:
            sbo = "SBO:0000655"
        _fill_sbo(rxn, sbo)

    # Faithfully replicate the MATLAB bug: only the last reaction is
    # considered for the pseudoreaction overrides.
    last = rxns[-1]
    if last.name == "biomass pseudoreaction":
        _fill_sbo(last, "SBO:0000629")
    elif last.name == "non-growth associated maintenance reaction":
        _fill_sbo(last, "SBO:0000630")
    elif "pseudoreaction" in last.name or "SLIME rxn" in last.name:
        _fill_sbo(last, "SBO:0000395")


def _fill_sbo(entity, sbo: str) -> None:
    """Set ``annotation['sbo']`` only if it is missing or empty."""
    current = entity.annotation.get("sbo")
    if not current:
        entity.annotation["sbo"] = sbo


def _transport_reaction_ids(model: cobra.Model) -> set[str]:
    """Heuristic mirror of RAVEN's `getTransportRxns`.

    A reaction is treated as a transport reaction if any metabolite
    *name* appears in two or more distinct compartments among its
    participants. Yeast-GEM stores metabolite names without compartment
    suffixes, so name equality is the natural grouping key.
    """
    out: set[str] = set()
    for rxn in model.reactions:
        by_name: dict[str, set[str | None]] = {}
        for met in rxn.metabolites:
            by_name.setdefault(met.name, set()).add(met.compartment)
        if any(len(comps) >= 2 for comps in by_name.values()):
            out.add(rxn.id)
    return out


# --- ΔG persistence ---------------------------------------------------

def load_delta_g(model: cobra.Model, *,
                 met_csv: Path | str | None = None,
                 rxn_csv: Path | str | None = None) -> cobra.Model:
    """Populate ΔG annotations on the model from the project CSVs.

    Port of `code/missingFields/loadDeltaG.m`. ΔG values are stored on
    each metabolite / reaction's ``notes`` dict under
    ``notes['deltaG']`` (as a string of the float), so they survive
    SBML round-trip via the standard notes element.

    Identifiers not present in the CSV keep whatever ΔG annotation they
    already had (or none). The CSV format is two columns ``Var1, Var2``
    (id, deltaG); the legacy MATLAB sentinel value ``10000000`` is
    preserved verbatim.
    """
    met_csv = Path(met_csv) if met_csv else _MET_CSV
    rxn_csv = Path(rxn_csv) if rxn_csv else _RXN_CSV
    met_df = _read_delta_g_csv(met_csv)
    rxn_df = _read_delta_g_csv(rxn_csv)
    _stamp_delta_g(model.metabolites, met_df, "metabolite", met_csv)
    _stamp_delta_g(model.reactions, rxn_df, "reaction", rxn_csv)
    return model


def save_delta_g(model: cobra.Model, *,
                 verbose: bool = False,
                 met_csv: Path | str | None = None,
                 rxn_csv: Path | str | None = None) -> None:
    """Persist ΔG annotations to the project CSVs.

    Port of `code/missingFields/saveDeltaG.m`. Writes both CSVs
    unconditionally; entries missing a ΔG note are written as ``NaN``
    so the file stays aligned with the model's metabolite/reaction
    order (mirroring MATLAB's `array2table([ids, deltaG])`).
    """
    met_csv = Path(met_csv) if met_csv else _MET_CSV
    rxn_csv = Path(rxn_csv) if rxn_csv else _RXN_CSV
    _dump_delta_g(model.metabolites, met_csv)
    _dump_delta_g(model.reactions, rxn_csv)
    if verbose:
        print(f"Wrote {met_csv}")
        print(f"Wrote {rxn_csv}")


def _read_delta_g_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if list(df.columns) != ["Var1", "Var2"]:
        raise ValueError(
            f"{path} has unexpected columns {list(df.columns)}; "
            "expected ['Var1', 'Var2']"
        )
    return df


def _stamp_delta_g(entities, df: pd.DataFrame, kind: str, source: Path) -> None:
    lookup = dict(zip(df["Var1"], df["Var2"], strict=True))
    missing = []
    for entity in entities:
        value = lookup.get(entity.id)
        if value is None or (isinstance(value, float) and math.isnan(value)):
            missing.append(entity.id)
            continue
        entity.notes[_DELTA_G_NOTE_KEY] = str(value)
    if missing:
        print(
            f"Not all {kind} identifiers are matched to {source.name}; "
            f"{len(missing)} entries (e.g. {missing[:3]}) were left untouched."
        )


def _dump_delta_g(entities, path: Path) -> None:
    rows = []
    for entity in entities:
        raw = entity.notes.get(_DELTA_G_NOTE_KEY)
        if raw is None:
            value = math.nan
        else:
            try:
                value = float(raw)
            except ValueError:
                value = math.nan
        rows.append((entity.id, value))
    pd.DataFrame(rows, columns=["Var1", "Var2"]).to_csv(path, index=False)
