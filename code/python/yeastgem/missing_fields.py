"""Yeast-specific wrappers for raven-toolbox's annotation helpers.

The mechanism for SBO assignment lives in :mod:`raven_toolbox.annotation`.
This module configures it with the yeast-GEM data layout and the
bug-compat flag that keeps the model artifact byte-equivalent during the
migration. ΔG load/save is implemented locally (not delegated to
:func:`raven_toolbox.annotation.load_delta_g_csv`/``save_delta_g_csv``,
which are CSV-only): yeast-GEM's own ΔG tables are tab-separated, to
match ``reactions.tsv``/``metabolites.tsv``/``genes.tsv``.
"""
from __future__ import annotations

import csv
import math
from collections.abc import Iterable
from pathlib import Path

import cobra
from raven_toolbox.annotation import (
    add_sbo_terms as _ra_add_sbo_terms,
)

from yeastgem.io import REPO_PATH

_DELTAG_DIR = REPO_PATH / "data" / "databases"
_MET_TSV = _DELTAG_DIR / "model_metDeltaG.tsv"
_RXN_TSV = _DELTAG_DIR / "model_rxnDeltaG.tsv"

# Key under which the ΔG value is stored in cobra ``notes``.
_DELTA_G_NOTE_KEY = "deltaG"


# SBO:0000243 ("gene"), the same constant every gene in a committed
# model/yeast-GEM.xml already carries -- previously supplied only as a
# side effect of RAVEN's MATLAB SBML writer, never by any curation code
# (raven_toolbox.annotation.add_sbo_terms and the MATLAB addSBOterms.m it
# is ported from both only ever covered mets/rxns). Kept local to
# yeast-GEM rather than added to the shared upstream helper, mirroring
# addSBOterms.m's own local (not RAVEN-generic) gene handling.
_GENE_SBO = "SBO:0000243"


def add_sbo_terms(model: cobra.Model) -> cobra.Model:
    """Assign SBO terms with yeast-GEM defaults.

    Thin wrapper over :func:`raven_toolbox.annotation.add_sbo_terms`. The
    ``only_last_reaction_for_pseudo=True`` flag reproduces the legacy
    MATLAB ``addSBOterms.m`` typo (``for i = numel(model.rxns)``) so
    yeast-GEM stays byte-equivalent through the upstream migration.
    Fixing that bug is a future behaviour-change PR; flip this flag to
    ``False`` (the upstream default) once the change is lock-stepped
    with the MATLAB side.

    Also fills ``gene.annotation['sbo']`` with ``SBO:0000243`` for any
    gene that does not already carry one -- same "fill" semantic the
    upstream helper uses for mets/rxns (only set when missing/empty).
    """
    model = _ra_add_sbo_terms(model, only_last_reaction_for_pseudo=True)
    for gene in model.genes:
        if not gene.annotation.get("sbo"):
            gene.annotation["sbo"] = _GENE_SBO
    return model


def _load_delta_g_tsv(entities: Iterable, path: Path, *, note_key: str) -> int:
    """Record ``note_key`` on each entity from a tsv of ``id -> value``.

    Mirrors :func:`raven_toolbox.annotation.load_delta_g_csv`'s semantics
    (does not interpret values -- yeast-GEM's own "no measurement"
    placeholder, ``10000000.0``, is recorded as-is), just tab-separated
    and with ``id``/``deltaG`` headers instead of pandas' ``Var1``/
    ``Var2``.
    """
    with path.open(encoding="utf-8", newline="") as fh:
        lookup = {row["id"]: row["deltaG"] for row in csv.DictReader(fh, delimiter="\t")}

    stamped = 0
    for entity in entities:
        raw = lookup.get(entity.id)
        if raw is None:
            continue
        try:
            if math.isnan(float(raw)):
                continue
        except ValueError:
            continue
        entity.notes[note_key] = raw
        stamped += 1
    return stamped


def _save_delta_g_tsv(entities: Iterable, path: Path, *, note_key: str) -> int:
    """Dump ``entity.notes[note_key]`` for each entity to a tsv.

    Entities without ``note_key`` set get ``nan`` written, preserving
    one-row-per-entity ordering (mirrors MATLAB's ``saveDeltaG.m``).
    """
    rows: list[tuple[str, str]] = []
    for entity in entities:
        raw = entity.notes.get(note_key)
        try:
            value = str(float(raw)) if raw is not None else "nan"
        except ValueError:
            value = "nan"
        rows.append((entity.id, value))

    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "deltaG"])
        writer.writerows(rows)
    return len(rows)


def load_delta_g(model: cobra.Model, *,
                 met_tsv: Path | str | None = None,
                 rxn_tsv: Path | str | None = None) -> cobra.Model:
    """Populate ΔG annotations on the model from the project tsvs.

    The tsv paths default to ``data/databases/model_{met,rxn}DeltaG.tsv``.
    Values land in ``entity.notes['deltaG']``.

    These are estimated, not curator-verified values, so they are never
    part of ``load_yeast_yaml``'s result or ``model/yeast-GEM.yml`` --
    call this explicitly if you want them.
    """
    met_tsv = Path(met_tsv) if met_tsv else _MET_TSV
    rxn_tsv = Path(rxn_tsv) if rxn_tsv else _RXN_TSV
    _load_delta_g_tsv(model.metabolites, met_tsv, note_key=_DELTA_G_NOTE_KEY)
    _load_delta_g_tsv(model.reactions, rxn_tsv, note_key=_DELTA_G_NOTE_KEY)
    return model


def save_delta_g(model: cobra.Model, *,
                 verbose: bool = False,
                 met_tsv: Path | str | None = None,
                 rxn_tsv: Path | str | None = None) -> None:
    """Persist ΔG annotations to the project tsvs.

    These are estimated, not curator-verified values, so they are never
    part of ``save_yeast_yaml``'s or ``commit_yeast_model``'s output --
    call this explicitly if you want to persist them.
    """
    met_tsv = Path(met_tsv) if met_tsv else _MET_TSV
    rxn_tsv = Path(rxn_tsv) if rxn_tsv else _RXN_TSV
    _save_delta_g_tsv(model.metabolites, met_tsv, note_key=_DELTA_G_NOTE_KEY)
    _save_delta_g_tsv(model.reactions, rxn_tsv, note_key=_DELTA_G_NOTE_KEY)
    if verbose:
        print(f"Wrote {met_tsv}")
        print(f"Wrote {rxn_tsv}")
