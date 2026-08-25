"""Yeast-specific wrappers for raven-toolbox's annotation helpers.

The mechanism for SBO assignment and ΔG side-car CSV persistence lives
in :mod:`raven_toolbox.annotation`. This module configures those helpers
with the yeast-GEM data layout (the CSV paths under
``data/databases/``) and the bug-compat flag that keeps the model
artifact byte-equivalent during the migration.
"""
from __future__ import annotations

from pathlib import Path

import cobra
from raven_toolbox.annotation import (
    add_sbo_terms as _ra_add_sbo_terms,
)
from raven_toolbox.annotation import (
    load_delta_g_csv as _ra_load_delta_g_csv,
)
from raven_toolbox.annotation import (
    save_delta_g_csv as _ra_save_delta_g_csv,
)

from yeastgem.io import REPO_PATH

_DELTAG_DIR = REPO_PATH / "data" / "databases"
_MET_CSV = _DELTAG_DIR / "model_metDeltaG.csv"
_RXN_CSV = _DELTAG_DIR / "model_rxnDeltaG.csv"

# Key under which the ΔG value is stored in cobra ``notes``.
_DELTA_G_NOTE_KEY = "deltaG"


def add_sbo_terms(model: cobra.Model) -> cobra.Model:
    """Assign SBO terms with yeast-GEM defaults.

    Thin wrapper over :func:`raven_toolbox.annotation.add_sbo_terms`. The
    ``only_last_reaction_for_pseudo=True`` flag reproduces the legacy
    MATLAB ``addSBOterms.m`` typo (``for i = numel(model.rxns)``) so
    yeast-GEM stays byte-equivalent through the upstream migration.
    Fixing that bug is a future behaviour-change PR; flip this flag to
    ``False`` (the upstream default) once the change is lock-stepped
    with the MATLAB side.
    """
    return _ra_add_sbo_terms(model, only_last_reaction_for_pseudo=True)


def load_delta_g(model: cobra.Model, *,
                 met_csv: Path | str | None = None,
                 rxn_csv: Path | str | None = None) -> cobra.Model:
    """Populate ΔG annotations on the model from the project CSVs.

    Thin wrapper over :func:`raven_toolbox.annotation.load_delta_g_csv`.
    The CSV paths default to ``data/databases/model_{met,rxn}DeltaG.csv``.
    Values land in ``entity.notes['deltaG']``.
    """
    met_csv = Path(met_csv) if met_csv else _MET_CSV
    rxn_csv = Path(rxn_csv) if rxn_csv else _RXN_CSV
    _ra_load_delta_g_csv(model.metabolites, met_csv, note_key=_DELTA_G_NOTE_KEY)
    _ra_load_delta_g_csv(model.reactions, rxn_csv, note_key=_DELTA_G_NOTE_KEY)
    return model


def save_delta_g(model: cobra.Model, *,
                 verbose: bool = False,
                 met_csv: Path | str | None = None,
                 rxn_csv: Path | str | None = None) -> None:
    """Persist ΔG annotations to the project CSVs.

    Thin wrapper over :func:`raven_toolbox.annotation.save_delta_g_csv`.
    """
    met_csv = Path(met_csv) if met_csv else _MET_CSV
    rxn_csv = Path(rxn_csv) if rxn_csv else _RXN_CSV
    _ra_save_delta_g_csv(model.metabolites, met_csv, note_key=_DELTA_G_NOTE_KEY)
    _ra_save_delta_g_csv(model.reactions, rxn_csv, note_key=_DELTA_G_NOTE_KEY)
    if verbose:
        print(f"Wrote {met_csv}")
        print(f"Wrote {rxn_csv}")
