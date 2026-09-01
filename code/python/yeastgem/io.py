"""Read and write the yeast-GEM model.

Ports `code/io.py` into the `yeastgem` package. Implements the three
curation entry points that mirror their MATLAB counterparts
(yeast-GEM#379): `load_yeast_yaml` / `code/loadYeastYaml.m`,
`save_yeast_yaml` / `code/saveYeastYaml.m`, and `commit_yeast_model` /
`code/commitYeastModel.m` (renamed from `saveYeastModel`).
"""
from __future__ import annotations

import csv
import os
import warnings
from copy import copy
from pathlib import Path

import cobra
from cobra.io import validate_sbml_model, write_sbml_model
from raven_toolbox.io import export_for_git, read_yaml_model, write_yaml_model

try:
    from dotenv import find_dotenv  # optional, kept for backwards compat
except ImportError:  # pragma: no cover - dotenv is a soft dependency
    find_dotenv = None  # type: ignore[assignment]


def _find_repo_root() -> Path:
    """Locate the yeast-GEM repo root.

    Resolution order:
    1. The ``YEAST_GEM_PATH`` environment variable, if set.
    2. Walk up from this file looking for ``model/yeast-GEM.yml``.
    3. ``find_dotenv`` (historical convention; .env at repo root).
    4. Walk up from CWD looking for ``model/yeast-GEM.yml``.

    The marker is the ``.yml``, not ``.xml``: the yml is curators' source
    of truth and always tracked, while the xml is a generated artifact
    that may not exist at all on a fresh checkout (yeast-GEM#379 stage 2).
    """
    override = os.environ.get("YEAST_GEM_PATH")
    if override:
        return Path(override).resolve()

    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "model" / "yeast-GEM.yml").exists():
            return parent

    if find_dotenv is not None:
        env = find_dotenv(usecwd=True)
        if env:
            return Path(env).parent.resolve()

    cwd = Path.cwd().resolve()
    for parent in (cwd, *cwd.parents):
        if (parent / "model" / "yeast-GEM.yml").exists():
            return parent

    raise FileNotFoundError(
        "Cannot locate the yeast-GEM repository root. "
        "Set the YEAST_GEM_PATH environment variable to the repo root, "
        "or place a .env file there."
    )


REPO_PATH: Path = _find_repo_root()
MODEL_PATH: Path = REPO_PATH / "model" / "yeast-GEM.xml"
YAML_PATH: Path = REPO_PATH / "model" / "yeast-GEM.yml"


# --- load / save for curation --------------------------------------------


def load_yeast_yaml(*, make_bigg_compliant: bool = False) -> cobra.Model:
    """Load model/yeast-GEM.yml for curation.

    Reads via :func:`raven_toolbox.io.read_yaml_model` — not
    ``cobra.io.load_yaml_model``, which silently drops RAVEN-only fields
    (``model.id``/``name``/``version``, ``deltaG``, ``confidence_score``,
    ``notes``, ``inchis``) — merges in the reaction/metabolite/gene
    cross-reference annotation from
    ``model/{reactions,metabolites,genes}.tsv`` (yeast-GEM#379) via
    :func:`yeastgem.annotate.annotate_gem`, and restores the ΔG fields via
    :func:`yeastgem.missing_fields.load_delta_g`. Mirrors
    `code/loadYeastYaml.m`.

    Loading ``model/yeast-GEM.xml`` does not need a yeast-GEM wrapper:
    once RAVEN/raven-toolbox support ΔG in SBML, ``read_sbml_model``
    already returns the complete model.

    Parameters
    ----------
    make_bigg_compliant
        If ``True``, rewrite metabolite/reaction ids using the BiGG
        dictionaries under ``data/databases/``. Default ``False``.
    """
    from yeastgem.annotate import annotate_gem
    from yeastgem.missing_fields import load_delta_g

    model = read_yaml_model(str(YAML_PATH))
    _collapse_single_value_annotations(model)
    # Pass YAML_PATH's own directory explicitly, for the same reason as in
    # save_yeast_yaml: annotate_gem's REPO_PATH-derived default is a
    # separate binding from annotate.py's own module load time.
    annotate_gem(model, YAML_PATH.parent)
    load_delta_g(model)

    if make_bigg_compliant and "x" not in model.compartments:
        _make_bigg_compliant(model)
    return model


def _collapse_single_value_annotations(model: cobra.Model) -> None:
    """Collapse single-element list annotation values to scalars, in place.

    ``raven_toolbox.io.read_yaml_model`` represents every MIRIAM-style
    annotation value as a list, even single-valued ones (e.g. ``sbo`` —
    the only annotation key yeast-GEM's yml carries inline, everything
    else living in the tsvs, yeast-GEM#379) — unlike ``cobra``'s own SBML
    reader, which collapses a lone value to a plain string. The rest of
    this codebase (``add_sbo_terms``, every ``annotation['sbo'] ==
    '...'`` comparison) is written against the scalar form, so this keeps
    a yml-loaded model's annotation shape consistent with an xml-loaded
    one rather than leaking the reader's own representation choice.
    """
    for collection in (model.reactions, model.metabolites, model.genes):
        for entity in collection:
            for key, value in list(entity.annotation.items()):
                if isinstance(value, list) and len(value) == 1:
                    entity.annotation[key] = value[0]


def save_yeast_yaml(
    model: cobra.Model,
    *,
    allow_no_growth: bool = True,
    derive_tsvs: bool = False,
) -> cobra.Model:
    """Save a curated model back to model/yeast-GEM.yml.

    Applies ``minimal_Y6`` (canonical medium), adds SBO terms, and checks
    aerobic and anaerobic growth, then writes ``model/yeast-GEM.yml``
    (with the tsv-sourced cross-reference annotation stripped back out,
    yeast-GEM#379) and the ΔG side-car CSVs. Mirrors
    `code/saveYeastYaml.m`.

    This is the function every curation script should call after editing
    a model — for producing the ``.xml``/``.txt`` files instead, see
    :func:`commit_yeast_model`.

    Parameters
    ----------
    model
        Model to save. Mutated in place (``minimal_Y6`` + SBO terms).
    allow_no_growth
        When True (default), a growth failure warns rather than raises.
    derive_tsvs
        If ``True``, ``model/{reactions,metabolites,genes}.tsv`` are
        overwritten from the model's current annotation via
        :func:`yeastgem.annotate.derive_annotation_tsvs`, so a
        cross-reference edit made programmatically (not by hand-editing a
        tsv cell) is not silently dropped. Default ``False``: a curator's
        hand-edited tsv is never overwritten unless this is explicitly
        requested.
    """
    from yeastgem import conditions
    from yeastgem.annotate import derive_annotation_tsvs
    from yeastgem.missing_fields import add_sbo_terms, save_delta_g

    conditions.apply(model, "minimal_Y6")
    add_sbo_terms(model)

    _check_growth(model, "aerobic", allow_no_growth)
    _check_growth_anaerobic(model, allow_no_growth)

    if derive_tsvs:
        # Pass YAML_PATH's own directory explicitly, so the tsvs always
        # land next to the yml this call is about to write -- rather than
        # relying on derive_annotation_tsvs' own REPO_PATH-derived default,
        # a separate binding imported at annotate.py's module load time
        # that would not follow a later-overridden YAML_PATH/REPO_PATH.
        derive_annotation_tsvs(model, YAML_PATH.parent)

    lean_model = _strip_tsv_annotation(model)
    write_yaml_model(lean_model, str(YAML_PATH))

    save_delta_g(model)

    return model


def _strip_tsv_annotation(model: cobra.Model) -> cobra.Model:
    """Copy of ``model`` with the tsv-sourced cross-reference annotation
    removed, so it is not re-embedded in model/yeast-GEM.yml on every save
    (yeast-GEM#379). Mirrors the MATLAB ``stripTsvAnnotation``.

    Replaces each entity's ``.annotation`` with a new dict rather than
    popping keys from the existing one in place: ``cobra.Model.copy()``
    does not deep-copy per-entity ``.annotation`` dicts, so an in-place
    pop here would silently strip the same keys from the caller's own
    original model too, not just this copy.
    """
    from yeastgem.annotate import GENE_COLUMNS, MET_COLUMNS, RXN_COLUMNS

    stripped = model.copy()
    for rxn in stripped.reactions:
        rxn.annotation = {k: v for k, v in rxn.annotation.items() if k not in RXN_COLUMNS}
    for met in stripped.metabolites:
        met.annotation = {k: v for k, v in met.annotation.items() if k not in MET_COLUMNS}
    for gene in stripped.genes:
        gene.annotation = {k: v for k, v in gene.annotation.items() if k not in GENE_COLUMNS}
    return stripped


def read_yeast_model(make_bigg_compliant: bool = False) -> cobra.Model:
    """DEPRECATED — use :func:`load_yeast_yaml` instead.

    ``read_yeast_model`` read ``model/yeast-GEM.xml`` directly, without
    merging in the tsv cross-reference annotation or restoring ΔG. This
    shim forwards to ``load_yeast_yaml`` and emits a DeprecationWarning.
    """
    warnings.warn(
        "read_yeast_model is deprecated; use load_yeast_yaml instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return load_yeast_yaml(make_bigg_compliant=make_bigg_compliant)


def write_yeast_model(model: cobra.Model) -> None:
    """DEPRECATED — use :func:`commit_yeast_model` instead.

    ``write_yeast_model`` implied a casual write, but the release path
    needs the full pipeline (canonical state, validation gates, ΔG
    CSVs). This shim forwards to ``commit_yeast_model`` with its default
    arguments and emits a DeprecationWarning. It will be removed at the
    next minor version bump after the rename ships.
    """
    warnings.warn(
        "write_yeast_model is deprecated; use commit_yeast_model instead. "
        "See code/python/PORTING_PLAN.md (phase 3) for the rename rationale.",
        DeprecationWarning,
        stacklevel=2,
    )
    commit_yeast_model(model)


# --- commit / export ------------------------------------------------------

_COMMIT_FORMATS = ("xml", "txt")


def commit_yeast_model(
    model: cobra.Model,
    *,
    formats: tuple[str, ...] = _COMMIT_FORMATS,
    allow_no_growth: bool = True,
) -> cobra.Model:
    """Prepare the yeast-GEM artifacts for a curation PR or a release.

    NOT a casual save: independently re-applies ``minimal_Y6``, SBO terms
    and the aerobic/anaerobic growth checks — it does not trust that
    :func:`save_yeast_yaml` already ran — validates the model as SBML,
    and writes the requested file formats. Run this *before*
    ``git commit``; it does not perform the commit itself. Use this for a
    local build or a release; for a routine curation save, use
    :func:`save_yeast_yaml` instead, which is faster and does not require
    SBML validation. Mirrors `code/commitYeastModel.m`, except only
    ``'xml'``/``'txt'`` are available here: ``'xlsx'``/``'mat'`` need
    RAVEN, so use the MATLAB ``commitYeastModel`` for those.

    Does not write ``model/yeast-GEM.yml`` or the annotation tsvs — that
    is exclusively :func:`save_yeast_yaml`'s job.

    Pipeline
    --------
    1. Apply ``minimal_Y6`` (canonical media) via :mod:`yeastgem.conditions`.
    2. Apply ``add_sbo_terms`` (canonical SBO annotations) via
       :mod:`yeastgem.missing_fields`.
    3. Validate that the model writes as valid SBML (cobrapy's
       ``validate_sbml_model``).
    4. Aerobic growth check — fail (or warn) if the model cannot grow.
    5. Anaerobic growth check — apply the ``anaerobic`` condition on a
       *copy* and confirm the resulting model still grows.
    6. Write the requested formats to ``model/`` via
       :func:`raven_toolbox.io.export_for_git`.
    7. Persist ΔG annotations via :func:`save_delta_g`.

    Root ``README.md`` is deliberately not touched here: its model
    statistics and validation numbers are only ever current for a
    specific released version, so they are stamped once at release time
    (``code/python/release/increase_version.py``) rather than on every
    curation commit, where they would immediately start going stale on
    ``develop``.

    No ``e-005`` → ``e-05`` exponent normalisation: Python's SBML writer
    does not produce the legacy MATLAB string, so the patch is
    unnecessary here.

    Parameters
    ----------
    model
        Model to commit.
    formats
        Which file formats to write, any of ``'xml'``, ``'txt'`` (default
        both).
    allow_no_growth
        When True (default), an aerobic-growth failure warns rather
        than raises. When False, it raises ``RuntimeError``.
    """
    # Import locally to keep the io module free of circular imports —
    # these submodules depend on REPO_PATH from this module.
    from yeastgem import conditions
    from yeastgem.missing_fields import add_sbo_terms, save_delta_g

    unknown = set(formats) - set(_COMMIT_FORMATS)
    if unknown:
        raise ValueError(
            f"commit_yeast_model cannot write {sorted(unknown)}: only "
            f"{_COMMIT_FORMATS} are available in Python. 'xlsx'/'mat' need "
            "RAVEN -- use the MATLAB commitYeastModel for those."
        )

    conditions.apply(model, "minimal_Y6")
    add_sbo_terms(model)

    _check_sbml_validity(model)
    _check_growth(model, "aerobic", allow_no_growth)
    _check_growth_anaerobic(model, allow_no_growth)

    if formats:
        export_for_git(
            model, MODEL_PATH.parent, prefix="yeast-GEM",
            formats=formats, sub_dirs=False,
        )
    save_delta_g(model)

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
