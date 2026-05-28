"""Read and write the yeast-GEM SBML model.

Ports `code/io.py` into the `yeastgem` package. The full
`commit_yeast_model` release pipeline (canonical state, validation
gates, multi-format export, README update) lands in a later phase —
this module currently only provides the bare read/write that the
historical `code/io.py` already offered.
"""
from __future__ import annotations

import csv
import os
from copy import copy
from pathlib import Path

import cobra
from cobra.io import read_sbml_model, write_sbml_model

try:
    from dotenv import find_dotenv  # optional, kept for backwards compat
except ImportError:  # pragma: no cover - dotenv is a soft dependency
    find_dotenv = None  # type: ignore[assignment]


def _find_repo_root() -> Path:
    """Locate the yeast-GEM repo root.

    Resolution order:
    1. The ``YEAST_GEM_PATH`` environment variable, if set.
    2. Walk up from this file looking for ``model/yeast-GEM.xml``.
    3. ``find_dotenv`` (historical convention; .env at repo root).
    4. Walk up from CWD looking for ``model/yeast-GEM.xml``.
    """
    override = os.environ.get("YEAST_GEM_PATH")
    if override:
        return Path(override).resolve()

    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "model" / "yeast-GEM.xml").exists():
            return parent

    if find_dotenv is not None:
        env = find_dotenv(usecwd=True)
        if env:
            return Path(env).parent.resolve()

    cwd = Path.cwd().resolve()
    for parent in (cwd, *cwd.parents):
        if (parent / "model" / "yeast-GEM.xml").exists():
            return parent

    raise FileNotFoundError(
        "Cannot locate the yeast-GEM repository root. "
        "Set the YEAST_GEM_PATH environment variable to the repo root, "
        "or place a .env file there."
    )


REPO_PATH: Path = _find_repo_root()
MODEL_PATH: Path = REPO_PATH / "model" / "yeast-GEM.xml"


def read_yeast_model(make_bigg_compliant: bool = False) -> cobra.Model:
    """Read the yeast-GEM SBML file via cobrapy.

    Parameters
    ----------
    make_bigg_compliant
        If ``True``, rewrite metabolite/reaction ids using the BiGG
        dictionaries under ``data/databases/``. Preserved from the
        legacy ``code/io.py`` for backwards compatibility; default
        ``False``.
    """
    model = read_sbml_model(str(MODEL_PATH))
    if make_bigg_compliant and "x" not in model.compartments:
        _make_bigg_compliant(model)
    return model


def write_yeast_model(model: cobra.Model) -> None:
    """Write the model to ``model/yeast-GEM.xml`` via cobrapy.

    NOTE: this is the bare SBML write. The full release pipeline
    (``commit_yeast_model``) — canonical state, validation gates,
    multi-format export, README/version updates — lands in a later
    phase of the Python port (see [PORTING_PLAN.md](../PORTING_PLAN.md),
    phase 3).
    """
    write_sbml_model(model, str(MODEL_PATH))


# --- BiGG compliance helper (ported verbatim from legacy code/io.py) ---

def _make_bigg_compliant(model: cobra.Model) -> None:
    data_path = REPO_PATH / "data" / "databases"
    met_bigg_dict = _load_bigg_dict(data_path / "BiGGmetDictionary_newIDs.csv")
    rxn_bigg_dict = _load_bigg_dict(data_path / "BiGGrxnDictionary_newIDs.csv")

    # Metabolite changes
    comp_dic = {"er": "r", "erm": "rm", "p": "x"}
    for met in model.metabolites:
        met.notes["Original ID"] = met.id
        comp_name = model.compartments[met.compartment]
        met.name = met.name.replace(f" [{comp_name}]", "")
        if met.compartment in comp_dic:
            met.compartment = comp_dic[met.compartment]
        if "bigg.metabolite" in met.annotation:
            _add_new_id(met, met.annotation["bigg.metabolite"])
        elif met.id in met_bigg_dict:
            _add_new_id(met, met_bigg_dict[met.id])
        else:
            met.id = met.id.replace(f"[{met.compartment}]", f"_{met.compartment}")

    # Compartment renames
    comps = model.compartments
    comps["r"] = "endoplasmic reticulum"
    comps["rm"] = "endoplasmic reticulum membrane"
    comps["x"] = "peroxisome"
    model.compartments = comps

    # Reaction changes
    for rxn in model.reactions:
        if "bigg.reaction" in rxn.annotation:
            rxn.notes["Original ID"] = rxn.id
            _add_new_id(rxn, rxn.annotation["bigg.reaction"])
        elif rxn.id in rxn_bigg_dict:
            rxn.notes["Original ID"] = rxn.id
            _add_new_id(rxn, rxn_bigg_dict[rxn.id])


def _load_bigg_dict(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(path) as f:
        for row in csv.reader(f, delimiter=","):
            out[row[0]] = row[1]
    return out


def _add_new_id(element, new_id: str) -> None:
    """Assign ``new_id`` to ``element``, appending ``_copyN`` on collision."""
    original = copy(new_id)
    copy_number = 1
    while True:
        try:
            if hasattr(element, "compartment"):  # metabolites
                element.id = f"{new_id}_{element.compartment}"
            else:  # reactions
                element.id = new_id
            return
        except ValueError:
            new_id = f"{original}_copy{copy_number}"
            copy_number += 1
