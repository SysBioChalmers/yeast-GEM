"""Read and write the yeast-GEM SBML model.

Ports `code/io.py` into the `yeastgem` package and implements
`commit_yeast_model` — the release pipeline that mirrors the MATLAB
`commitYeastModel` (renamed from `saveYeastModel` in phase 3 of the
port; see [PORTING_PLAN.md](../PORTING_PLAN.md)).
"""
from __future__ import annotations

import csv
import os
import re
import warnings
from copy import copy
from datetime import datetime
from pathlib import Path

import cobra
from cobra.io import read_sbml_model, validate_sbml_model, write_sbml_model

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
    """DEPRECATED — use :func:`commit_yeast_model` instead.

    ``write_yeast_model`` implied a casual write, but the release path
    needs the full pipeline (canonical state, validation gates, ΔG
    CSVs, README update). This shim forwards to ``commit_yeast_model``
    with its default arguments and emits a DeprecationWarning. It will
    be removed at the next minor version bump after the rename ships.
    """
    warnings.warn(
        "write_yeast_model is deprecated; use commit_yeast_model instead. "
        "See code/python/PORTING_PLAN.md (phase 3) for the rename rationale.",
        DeprecationWarning,
        stacklevel=2,
    )
    commit_yeast_model(model)


# --- the release pipeline ------------------------------------------------

# Regex matching the model-stats table row in README.md, mirroring the
# legacy MATLAB regex used by saveYeastModel.m. Captures the species
# label so the rewrite preserves it.
_README_STATS_RE = re.compile(
    r"^\| (\_Saccharomyces cerevisiae\_) \| "
    r"\d{2}-\D+-\d{4} \| \d+ \| \d+ \| \d+ \|",
    re.MULTILINE,
)


def commit_yeast_model(
    model: cobra.Model,
    *,
    update_readme: bool = True,
    allow_no_growth: bool = True,
) -> cobra.Model:
    """Prepare the yeast-GEM artifacts for a curation PR.

    NOT a casual save: this is the release pipeline. Run this *before*
    ``git commit``; it does not perform the commit itself. Mirrors
    `code/commitYeastModel.m`.

    Pipeline
    --------
    1. Apply ``minimal_Y6`` (canonical media) via :mod:`yeastgem.conditions`.
    2. Apply ``add_sbo_terms`` (canonical SBO annotations) via
       :mod:`yeastgem.missing_fields`.
    3. Validate that the model writes as valid SBML (cobrapy's
       ``validate_sbml_model``).
    4. Aerobic growth check — fail (or warn) if the model cannot grow.
    5. Anaerobic growth check — apply the ``anaerobic`` condition on a
       *copy* and confirm the resulting model still grows. (Phase 4
       activated this once the biomass / amino-acid-ratio plumbing
       landed.)
    6. Write SBML to ``model/yeast-GEM.xml``.
    7. Persist ΔG annotations via :func:`save_delta_g`.
    8. Update ``README.md`` with the current date and model size
       (if ``update_readme`` is True).

    Limitations vs. the MATLAB pipeline (will close as later phases land)
        - No ``.yml`` / ``.txt`` / ``.xlsx`` / ``.mat`` companion exports
          (RAVEN's ``exportForGit``). The canonical artifact written by
          this function is ``model/yeast-GEM.xml`` only; the companions
          must currently be regenerated by running the MATLAB
          ``commitYeastModel``.
        - No ``e-005`` → ``e-05`` exponent normalisation. Python's SBML
          writer does not produce the legacy MATLAB string; the patch
          is unnecessary here.

    Parameters
    ----------
    model
        Model to commit.
    update_readme
        Whether to rewrite the model-stats row in ``README.md``
        (default True).
    allow_no_growth
        When True (default), an aerobic-growth failure warns rather
        than raises. When False, it raises ``RuntimeError``.
    """
    # Import locally to keep the io module free of circular imports —
    # these submodules depend on REPO_PATH from this module.
    from yeastgem import conditions
    from yeastgem.missing_fields import add_sbo_terms, save_delta_g

    conditions.apply(model, "minimal_Y6")
    add_sbo_terms(model)

    _check_sbml_validity(model)
    _check_growth(model, "aerobic", allow_no_growth)
    _check_growth_anaerobic(model, allow_no_growth)

    write_sbml_model(model, str(MODEL_PATH))
    save_delta_g(model)

    if update_readme:
        _update_readme(model)

    return model


def _check_sbml_validity(model: cobra.Model) -> None:
    """Round-trip through SBML and confirm cobrapy accepts the output.

    Mirrors the MATLAB ``TranslateSBML`` round-trip gate.
    """
    tmp = MODEL_PATH.parent / ".tempModel.xml"
    try:
        write_sbml_model(model, str(tmp))
        _, errors = validate_sbml_model(str(tmp))
        fatal = errors.get("SBML_ERROR") or errors.get("SBML_FATAL")
        if fatal:
            raise RuntimeError(
                "Model is not a valid SBML structure. Fix all errors "
                f"before committing:\n{fatal[:5]}"
            )
    finally:
        tmp.unlink(missing_ok=True)


def _check_growth(model: cobra.Model, condition: str, allow_no_growth: bool) -> None:
    """Solve FBA on a copy and surface a no-growth state."""
    test = model.copy()
    try:
        solution = test.optimize()
        ok = solution.status == "optimal" and solution.objective_value > 1e-6
    except Exception:  # pragma: no cover - solver failures
        ok = False

    if ok:
        return

    msg = (
        f"The model is not able to support growth under {condition} "
        "conditions. Please ensure the model can grow before opening a PR."
    )
    if allow_no_growth:
        warnings.warn(msg, stacklevel=3)
    else:
        raise RuntimeError(msg)


def _check_growth_anaerobic(model: cobra.Model, allow_no_growth: bool) -> None:
    """Apply the anaerobic condition on a copy and confirm FBA growth.

    Mirrors the MATLAB ``checkGrowth(model, 'anaerobic', ...)`` step in
    ``commitYeastModel.m``. The copy keeps the input model intact for
    the SBML write that follows.
    """
    from yeastgem import conditions

    anaerobic = model.copy()
    conditions.apply(anaerobic, "anaerobic")
    _check_growth(anaerobic, "anaerobic", allow_no_growth)


def _update_readme(model: cobra.Model) -> None:
    """Rewrite the model-stats row in ``README.md``.

    Mirrors the MATLAB regex rewrite. The species label and column order
    are preserved; only the date and the three size counters change. No
    version is stamped -- it would go stale on ``develop`` the moment the
    next curation lands, and it duplicates the "Current release" badge on
    ``main``.
    """
    readme = REPO_PATH / "README.md"
    date = datetime.now().strftime("%d-%b-%Y")
    n_rxns = len(model.reactions)
    n_mets = len(model.metabolites)
    n_genes = len(model.genes)
    replacement = f"| \\1 | {date} | {n_rxns} | {n_mets} | {n_genes} |"
    text = readme.read_text(encoding="utf-8")
    new_text = _README_STATS_RE.sub(replacement, text)
    readme.write_text(new_text, encoding="utf-8")


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
