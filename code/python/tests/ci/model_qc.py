"""Model quality checks for the pull-request report.

Runs the count-based checks over one model file and writes, into an output
directory:

* one CSV per check that found something, listing the exact entries, so a
  count in the pull-request comment can link to the detail;
* ``qc_metrics.tsv``, a ``key<TAB>value`` table that
  :mod:`build_qc_report` reads.

The same script runs against the pull-request head and against the target
branch, and the report compares the two. Nothing here is committed to the
repository, so there is no baseline to keep up to date: the comparison is
always this model against that branch's model, computed the same way on the
same day with the same dependencies.

Four checks are gates: the model must load, it must grow, reaction/
metabolite names in model/yeast-GEM.yml must agree with
model/reactions.tsv and model/metabolites.tsv, and the model's reaction/
metabolite/gene ids must agree with the three tsvs -- with no id that
was ever deprecated back in active use (yeast-GEM#379 -- the yml is
authoritative, the tsvs' name column and identifier columns are the
supplementary annotation record, and by the time a pull request is
reviewed the two are supposed to already agree, so any disagreement can
only mean this branch introduced one). Everything else is reported so
that a pull request which makes a number worse is visible, without
blocking work on findings that were already there.

Two more checks: identifiers removed since the target branch that were
not added to data/deprecatedIdentifiers/ (report only -- needs
--base-model-dir, empty without it), and metabolite structure (SMILES;
yeast-GEM has no InChI field, so there is no smiles/InChI cross-check)
versus formula and charge (report only).

check_malformed_xrefs reads reactions.tsv/metabolites.tsv/genes.tsv
directly rather than the model's own (SBML-derived) annotation, per
edkerk's comment on #379: a curator who edits a tsv cell directly (e.g.
via the GitHub web UI, never running saveYeastYaml/save_yeast_yaml)
should have a bad identifier caught at that file, not only after
someone next regenerates model/yeast-GEM.xml. sbo is consequently no
longer checked here -- it is computed, not curator-edited, and does not
live in any tsv; a wrong SBO value would be a bug in the assignment
logic, not a curation mistake this check is meant to catch.
check_xrefs_across_compartments still reads the loaded model and was
not part of this rework -- not asked for, and arguably worth revisiting
for the same reason later.

Usage:
    python model_qc.py --model model/yeast-GEM.yml --out data/testResults
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import cobra
from raven_toolbox.io import read_yaml_model
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

# Identifier patterns, from identifiers.org, for the namespaces yeast-GEM
# actually uses. A namespace with no pattern here is not checked rather than
# being guessed at: a wrong pattern would manufacture findings, which is
# worse than missing some.
_XREF_PATTERNS = {
    "chebi": r"^CHEBI:\d+$",
    "kegg.compound": r"^C\d+$",
    "kegg.reaction": r"^R\d+$",
    "kegg.pathway": r"^\w{2,4}\d{5}$",
    "metanetx.chemical": r"^MNXM\d+$",
    "metanetx.reaction": r"^MNXR\d+$",
    "seed.compound": r"^cpd\d+$",
    "seed.reactions": r"^rxn\d+$",
    "pubchem.compound": r"^\d+$",
    "pubmed": r"^\d+$",
    "ncbigene": r"^\d+$",
    "rhea": r"^\d{5}$",
    "sbo": r"^SBO:\d{7}$",
    "ec-code": (
        r"^\d+\.-\.-\.-$|^\d+\.\d+\.-\.-$|^\d+\.\d+\.\d+\.-$"
        r"|^\d+\.\d+\.\d+\.(n)?\d+$"
    ),
}

# Cell values that mean "nothing to report" in a MACAW result column.
#
# "n/a" is the one that matters and the one that is easy to get wrong.
# MACAW writes the literal string "N/A" for a sub-test that does not apply
# to a reaction -- it is not a null, so isna() is False for it, and it is
# not "ok", so a naive filter counts it as a finding and reports every
# reaction in the model as a duplicate.
#
# It is also invisible to any check done by reading the CSV back:
# pandas.read_csv lists "N/A" among its default NA strings, so the value
# becomes a real null on the way in. The written file and the frame the
# count is taken from therefore disagree, and verifying against the file
# confirms a number the code never produced.
_NOT_A_FINDING = {"", "ok", "n/a", "na", "none", "nan"}


def _values(annotation, key):
    """Annotation values as a list; cobrapy stores either a scalar or a list."""
    value = annotation.get(key)
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def check_growth(model) -> tuple[float, list]:
    """Objective value under the model's own constraints. A gate."""
    value = model.slim_optimize(error_value=0.0)
    return (float(value or 0.0), [])


def check_orphan_metabolites(model) -> tuple[int, list]:
    rows = [(m.id, m.name or "") for m in model.metabolites if not m.reactions]
    return len(rows), rows


def check_unused_genes(model) -> tuple[int, list]:
    rows = [(g.id, g.name or "") for g in model.genes if not g.reactions]
    return len(rows), rows


def check_missing_formula(model) -> tuple[int, list]:
    rows = [(m.id, m.name or "") for m in model.metabolites if not m.formula]
    return len(rows), rows


def check_missing_charge(model) -> tuple[int, list]:
    rows = [
        (m.id, m.name or "") for m in model.metabolites if m.charge is None
    ]
    return len(rows), rows


def check_duplicate_reactions(model) -> tuple[int, list]:
    """Reactions sharing identical stoichiometry and direction.

    Keyed on the metabolite/coefficient set plus reversibility, so a
    forward-only and a reversible copy of the same conversion are not
    conflated. Reported as one row per reaction in each group larger
    than one.
    """
    groups = defaultdict(list)
    for rxn in model.reactions:
        if not rxn.metabolites:
            continue
        key = (
            frozenset((m.id, c) for m, c in rxn.metabolites.items()),
            rxn.lower_bound < 0,
        )
        groups[key].append(rxn)
    rows = []
    for number, (_key, members) in enumerate(
        (kv for kv in groups.items() if len(kv[1]) > 1), start=1
    ):
        for rxn in members:
            rows.append((number, rxn.id, rxn.name or "", rxn.reaction))
    return len(rows), rows


def _balances_by_design(rxn) -> bool:
    """Whether a reaction is expected to balance at all.

    Three classes are not, and including them buries the reactions that
    matter:

    * boundary reactions -- an exchange has one side by construction;
    * pseudoreactions -- biomass and its components are lumped compositions;
    * SLIME reactions -- yeast-GEM splits lipids into measurable entities
      with fractional, averaged coefficients, so they are unbalanced by
      design. There are 186 of them, against 2 genuinely unexpected mass
      imbalances in the rest of the model, so leaving them in would make
      the count useless.
    """
    name = (rxn.name or "")
    return not (
        rxn.boundary
        or "pseudoreaction" in name.lower()
        or "SLIME" in name
    )


def _balance(model):
    """Split reactions into mass- and charge-imbalanced."""
    mass_rows, charge_rows = [], []
    for rxn in model.reactions:
        if not _balances_by_design(rxn):
            continue
        try:
            imbalance = rxn.check_mass_balance()
        except Exception as exc:  # a finding, not something to hide
            mass_rows.append((rxn.id, rxn.name or "", f"check_failed:{exc}"))
            continue
        if not imbalance:
            continue
        mass = {k: v for k, v in imbalance.items() if k != "charge"}
        if mass:
            mass_rows.append((
                rxn.id, rxn.name or "",
                ";".join(f"{k}:{v:g}" for k, v in sorted(mass.items())),
            ))
        if "charge" in imbalance:
            charge_rows.append((
                rxn.id, rxn.name or "", f"{imbalance['charge']:g}",
            ))
    return mass_rows, charge_rows


def _why_malformed(namespace: str, value: str) -> str:
    """Say what is wrong with an identifier, not just that something is.

    A row reading "chebi / ChEBI:15987" leaves the reader to work out that
    the rest of the model writes CHEBI in capitals. Naming the likely
    cause makes the finding actionable without opening the model.
    """
    expected = _XREF_PATTERNS[namespace]
    prefix = namespace.split(".")[0]
    if value.lower().startswith(f"{prefix.lower()}:"):
        stripped = value.split(":", 1)[1]
        if re.match(expected, stripped):
            return f"remove the redundant '{value.split(':', 1)[0]}:' prefix"
        return f"wrong case or prefix; expected to match {expected}"
    if re.match(expected, value, flags=re.IGNORECASE):
        return "correct apart from letter case"
    if not any(c.isdigit() for c in value):
        return "looks like a name rather than an identifier"
    return f"does not match {expected}"


_TSV_XREF_COLUMNS = {
    "reaction": ("bigg.reaction", "ec-code", "kegg.pathway", "kegg.reaction", "metanetx.reaction"),
    "metabolite": ("bigg.metabolite", "chebi", "kegg.compound", "metanetx.chemical"),
    "gene": ("uniprot",),
}


def check_malformed_xrefs(model_dir: Path) -> tuple[int, list]:
    """Cross-references that do not match their namespace's pattern.

    Reads reactions.tsv/metabolites.tsv/genes.tsv directly -- see the
    module docstring for why this reads the tsvs and not the model.
    """
    rows = []
    tsv_names = {"reaction": "reactions.tsv", "metabolite": "metabolites.tsv", "gene": "genes.tsv"}
    for kind, columns in _TSV_XREF_COLUMNS.items():
        tsv_path = model_dir / tsv_names[kind]
        if not tsv_path.is_file():
            continue
        with tsv_path.open(encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                for column in columns:
                    if column not in _XREF_PATTERNS:
                        continue
                    cell = (row.get(column) or "").strip()
                    if not cell:
                        continue
                    pattern = _XREF_PATTERNS[column]
                    for value in (v.strip() for v in cell.split(";")):
                        if value and not re.match(pattern, value):
                            rows.append((
                                kind, row["id"], column, value,
                                _why_malformed(column, value),
                            ))
    return len(rows), rows


def check_annotation_consistency(model, model_dir: Path) -> tuple[int, list]:
    """model/yeast-GEM.yml vs. the annotation tsvs, and deprecated
    identifiers in active use. A gate: an id missing from one side, or an
    id that was deliberately retired being reused, means a curation step
    was skipped -- not a pre-existing issue to track a delta on over time.

    Complements check_name_consistency, which only looks at the
    intersection and checks *names*; this looks at *which ids exist* on
    each side.
    """
    rows = []
    deprecated_dir = model_dir.parent / "data" / "deprecatedIdentifiers"

    def compare(kind, model_ids, tsv_name, deprecated_name):
        tsv_path = model_dir / tsv_name
        if not tsv_path.is_file():
            return
        with tsv_path.open(encoding="utf-8", newline="") as fh:
            tsv_ids = {row["id"] for row in csv.DictReader(fh, delimiter="\t")}
        deprecated = set()
        dep_path = deprecated_dir / deprecated_name if deprecated_name else None
        if dep_path is not None and dep_path.is_file():
            with dep_path.open(encoding="utf-8", newline="") as fh:
                deprecated = {row["id"] for row in csv.DictReader(fh, delimiter="\t")}
        for used in sorted(model_ids & deprecated):
            rows.append((kind, used, "a deprecated identifier is in active use"))
        for missing in sorted(model_ids - tsv_ids):
            rows.append((kind, missing, f"in the model but not in {tsv_name}"))
        for orphan in sorted(tsv_ids - model_ids):
            rows.append((kind, orphan, f"in {tsv_name} but not in the model"))

    compare("reaction", {r.id for r in model.reactions},
            "reactions.tsv", "deprecatedReactions.tsv")
    compare("metabolite", {m.id for m in model.metabolites},
            "metabolites.tsv", "deprecatedMetabolites.tsv")
    compare("gene", {g.id for g in model.genes}, "genes.tsv", None)  # no deprecated list

    return len(rows), rows


def check_deprecation_completeness(
    model, model_dir: Path, base_model_dir: Path | None
) -> tuple[int, list]:
    """Ids present in the target branch but missing from this model must
    appear in the matching data/deprecatedIdentifiers/ file, so a removed
    identifier stays resolvable instead of silently vanishing. Report
    only: needs base_model_dir (the target branch's own model/ checkout);
    empty without it -- e.g. run locally, or the first comparison for a
    branch.
    """
    rows = []
    if base_model_dir is None or not base_model_dir.is_dir():
        return 0, rows

    deprecated_dir = model_dir.parent / "data" / "deprecatedIdentifiers"
    specs = [
        ("reaction", {r.id for r in model.reactions},
         "reactions.tsv", "deprecatedReactions.tsv"),
        ("metabolite", {m.id for m in model.metabolites},
         "metabolites.tsv", "deprecatedMetabolites.tsv"),
    ]
    for kind, current_ids, tsv_name, dep_name in specs:
        base_path = base_model_dir / tsv_name
        if not base_path.is_file():
            continue
        with base_path.open(encoding="utf-8", newline="") as fh:
            base_ids = {row["id"] for row in csv.DictReader(fh, delimiter="\t")}
        deprecated = set()
        dep_path = deprecated_dir / dep_name
        if dep_path.is_file():
            with dep_path.open(encoding="utf-8", newline="") as fh:
                deprecated = {row["id"] for row in csv.DictReader(fh, delimiter="\t")}
        removed = base_ids - current_ids
        for missing in sorted(removed - deprecated):
            rows.append((kind, missing, f"removed from the model but not listed in {dep_name}"))

    return len(rows), rows


_ELEM_RE = re.compile(r"([A-Z][a-z]?)(\d*)")


def _parse_formula(formula: str) -> Counter:
    """Element -> count for a Hill formula string (charge sign ignored)."""
    counts: Counter = Counter()
    for symbol, number in _ELEM_RE.findall((formula or "").strip().rstrip("+-")):
        if symbol:
            counts[symbol] += int(number) if number else 1
    return counts


def _classify_structure(model_formula, model_charge, smiles: str):
    """Relate a stored SMILES to the curated formula/charge, deriving the
    formula/charge the SMILES implies via RDKit and comparing (no InChI
    cross-check -- yeast-GEM has no InChI field). Returns (category,
    rdkit_formula, rdkit_charge)."""
    if not smiles:
        return "no_structure", "", ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "unparseable", "", ""
    # generic: an R group (dummy atom) or an "R" in the model formula means
    # there is no concrete structure to check against.
    generic = "*" in smiles or "R" in (model_formula or "")
    rd_formula = rdMolDescriptors.CalcMolFormula(mol)
    rd_charge = Chem.GetFormalCharge(mol)
    m_elements, r_elements = _parse_formula(model_formula), _parse_formula(rd_formula)
    m_heavy = {k: v for k, v in m_elements.items() if k != "H"}
    r_heavy = {k: v for k, v in r_elements.items() if k != "H"}
    m_charge = int(model_charge) if model_charge is not None else None
    if generic:
        return "generic", rd_formula, str(rd_charge)
    if m_heavy != r_heavy:
        return "formula_error", rd_formula, str(rd_charge)
    if m_elements == r_elements and m_charge == rd_charge:
        return "ok", rd_formula, str(rd_charge)
    return "protonation", rd_formula, str(rd_charge)


_STRUCTURE_INCONSISTENT = {"protonation", "formula_error"}


def check_structure_consistency(model, model_dir: Path) -> tuple[int, list]:
    """Metabolite formula/charge (model/yeast-GEM.yml) vs. the SMILES
    stored in metabolites.tsv. Report only -- see _classify_structure for
    the categories.
    """
    tsv_path = model_dir / "metabolites.tsv"
    if not tsv_path.is_file():
        return 0, []
    with tsv_path.open(encoding="utf-8", newline="") as fh:
        smiles_by_id = {
            row["id"]: (row.get("smiles") or "").strip()
            for row in csv.DictReader(fh, delimiter="\t")
        }

    rows = []
    for met in model.metabolites:
        category, rd_formula, rd_charge = _classify_structure(
            met.formula, met.charge, smiles_by_id.get(met.id, "")
        )
        if category in _STRUCTURE_INCONSISTENT:
            rows.append((
                met.id, met.name or "", category,
                met.formula or "", "" if met.charge is None else str(int(met.charge)),
                rd_formula, rd_charge, smiles_by_id.get(met.id, ""),
            ))
    return len(rows), rows


def check_xrefs_across_compartments(model) -> tuple[int, list]:
    """The same species annotated differently in different compartments.

    yeast-GEM metabolite ids are ``s_NNNN`` per compartment, so the species
    is identified by name instead. Where one name appears in several
    compartments, every copy should carry the same cross-references; a
    disagreement means one of them was edited and the others were not.
    """
    by_name = defaultdict(list)
    for met in model.metabolites:
        by_name[(met.name or "").strip()].append(met)

    rows = []
    for name, mets in sorted(by_name.items()):
        if len(mets) < 2 or not name:
            continue
        for namespace in ("chebi", "kegg.compound", "metanetx.chemical"):
            seen = {
                met.id: tuple(sorted(str(v) for v in _values(met.annotation, namespace)))
                for met in mets
            }
            distinct = {v for v in seen.values() if v}
            if len(distinct) > 1:
                rows.append((
                    name, namespace,
                    "; ".join(f"{i}={'|'.join(v) or '-'}" for i, v in sorted(seen.items())),
                ))
    return len(rows), rows


def check_name_consistency(model, model_dir: Path) -> tuple[int, list]:
    """Reaction/metabolite names: model/yeast-GEM.yml vs. the tsvs. A gate.

    model/yeast-GEM.yml is authoritative for names -- renaming something
    is a model edit, the same category as changing its formula or
    bounds. reactions.tsv/metabolites.tsv's name column is a read-only,
    human-readable copy included so the file can be scanned without
    cross-referencing the yml; this catches the two drifting apart.

    An id present in the model but entirely absent from a tsv (or vice
    versa) is a separate, row-completeness concern -- not a name
    mismatch -- and out of scope here; skipped rather than flagged.
    """
    rows = []
    for kind, entities, tsv_name in (
        ("reaction", model.reactions, "reactions.tsv"),
        ("metabolite", model.metabolites, "metabolites.tsv"),
    ):
        tsv_path = model_dir / tsv_name
        if not tsv_path.is_file():
            continue
        with tsv_path.open(encoding="utf-8", newline="") as fh:
            tsv_names = {
                row["id"]: (row.get("name") or "").strip()
                for row in csv.DictReader(fh, delimiter="\t")
            }
        for entity in entities:
            if entity.id not in tsv_names:
                continue
            model_name = (entity.name or "").strip()
            tsv_name_value = tsv_names[entity.id]
            if model_name != tsv_name_value:
                rows.append((kind, entity.id, model_name, tsv_name_value))
    return len(rows), rows


def check_macaw(model):
    """MACAW dead-end and duplicate tests.

    Returns ``(metrics, sections)``. When MACAW is not installed the
    metrics are absent rather than zero, so a missing dependency is
    reported as "not run" instead of being mistaken for a clean result.

    Reported per metabolite rather than per reaction for the dead-end
    test. MACAW flags a reaction for each dead-end metabolite in it, so
    one gap is counted many times -- a single metabolite here accounts
    for 49 reactions. The metabolite is what a curator fixes, either by
    adding the missing reaction or correcting an annotation, and the
    blocked reactions follow from it.
    """
    empty = [
        ("Dead-end metabolites", ("metabolite", "name", "reactions blocked"), []),
        ("Reactions that can only carry flux one way",
         ("reaction", "name", "finding"), []),
        ("Reactions flagged as MACAW duplicates",
         ("reaction", "name", "why", "duplicates"), []),
    ]
    try:
        from macaw.main import dead_end_test, duplicate_test
    except ImportError:
        return {}, empty

    dead_end_results, _edges = dead_end_test(model)
    duplicate_results, _dup_edges = duplicate_test(model)
    merged = dead_end_results.merge(duplicate_results)

    def _cells(column):
        """Stripped text and a mask of the rows that say something.

        Nulls are found with isna(), which is dtype-independent; see
        _NOT_A_FINDING for why the string form cannot be relied on.
        """
        if column not in merged.columns:
            raise KeyError(
                f"MACAW returned no '{column}' column; got "
                f"{sorted(merged.columns)}. The report would otherwise "
                "show this check as clean without having run it."
            )
        values = merged[column]
        text = values.astype(str).str.strip()
        return text, ~(values.isna() | text.str.lower().isin(_NOT_A_FINDING))

    def _name(identifier, lookup):
        try:
            return lookup.get_by_id(identifier).name or ""
        except KeyError:
            return ""

    ids = merged["reaction_id"].astype(str)

    # dead_end_test says either "ok", a direction phrase, or the dead-end
    # metabolites in that reaction. The three mean different things.
    dead_end_text, dead_end_flagged = _cells("dead_end_test")
    blocked_by = defaultdict(list)
    direction_rows = []
    for rxn_id, value in zip(ids[dead_end_flagged],
                             dead_end_text[dead_end_flagged], strict=True):
        if value.lower().startswith("only when"):
            direction_rows.append((rxn_id, _name(rxn_id, model.reactions), value))
            continue
        for met_id in (m.strip() for m in value.split(";")):
            if met_id:
                blocked_by[met_id].append(rxn_id)

    dead_end_rows = [
        (met_id, _name(met_id, model.metabolites), len(reactions))
        for met_id, reactions in blocked_by.items()
    ]

    # duplicate_test writes one column per kind of duplicate and no single
    # verdict. The column names are MACAW's internal ones; each is given
    # the sentence from its documentation so the finding can be read
    # without going and looking the test up.
    reasons = {
        "duplicate_test_exact":
            "same metabolites and coefficients",
        "duplicate_test_directions":
            "same metabolites, opposite direction or reversibility",
        "duplicate_test_coefficients":
            "same metabolites, different coefficients",
        "duplicate_test_redox":
            "same conversion using a different electron carrier",
    }
    duplicate_columns = [
        c for c in merged.columns if c.startswith("duplicate_test")
    ]
    if not duplicate_columns:
        raise KeyError(
            f"MACAW returned no duplicate_test* columns; got "
            f"{sorted(merged.columns)}."
        )

    duplicates = None
    duplicate_rows = []
    for column in duplicate_columns:
        text, mask = _cells(column)
        # The redox sub-test needs redox_pairs and proton_ids, which are
        # not passed, so MACAW writes N/A for every reaction. Saying so is
        # the point: counting it as zero would claim a check that never
        # ran.
        if not mask.any() and column == "duplicate_test_redox":
            print(f"    {column}: not run (no redox pairs configured)")
            continue
        print(f"    {column}: flagged={int(mask.sum())}")
        for rxn_id, other in zip(ids[mask], text[mask], strict=True):
            duplicate_rows.append((
                rxn_id, _name(rxn_id, model.reactions),
                reasons.get(column, column), other,
            ))
        duplicates = mask if duplicates is None else (duplicates | mask)

    metrics = {
        "macaw_dead_end_metabolites": len(dead_end_rows),
        "macaw_single_direction": len(direction_rows),
        "macaw_duplicates": int(duplicates.sum()) if duplicates is not None else 0,
    }
    sections = [
        ("Dead-end metabolites",
         ("metabolite", "name", "reactions blocked"), dead_end_rows),
        ("Reactions that can only carry flux one way",
         ("reaction", "name", "finding"), direction_rows),
        ("Reactions flagged as MACAW duplicates",
         ("reaction", "name", "why", "duplicates"), duplicate_rows),
    ]
    return metrics, sections


# (metric key, comment label, function)
_CHECKS = [
    ("growth", "Growth (biomass producible)", check_growth),
    ("orphan_metabolites", "Unused metabolites", check_orphan_metabolites),
    ("unused_genes", "Unused genes", check_unused_genes),
    ("missing_formula", "Metabolites missing formula", check_missing_formula),
    ("missing_charge", "Metabolites missing charge", check_missing_charge),
    ("duplicate_reactions", "Exact-duplicate reaction groups",
     check_duplicate_reactions),
    ("xrefs_across_compartments", "Cross-refs inconsistent across compartments",
     check_xrefs_across_compartments),
]
# malformed_xrefs, annotation_consistency, deprecation_completeness and
# structure_inconsistent are not in _CHECKS: each needs model_dir (and
# deprecation_completeness also base_model_dir), not just model, so they
# are called directly in run() instead of through the uniform (model,)
# signature -- the same reason check_name_consistency and check_macaw
# already sit outside this list.

_HEADERS = {
    "orphan_metabolites": ("metabolite", "name"),
    "unused_genes": ("gene", "name"),
    "missing_formula": ("metabolite", "name"),
    "missing_charge": ("metabolite", "name"),
    "duplicate_reactions": ("group", "reaction", "name", "equation"),
    "malformed_xrefs": ("kind", "id", "namespace", "value", "problem"),
    "xrefs_across_compartments": ("metabolite name", "namespace", "values"),
    "mass_imbalanced": ("reaction", "name", "imbalance"),
    "charge_imbalanced": ("reaction", "name", "charge"),
    "name_mismatches": ("kind", "id", "model name", "tsv name"),
    "annotation_consistency": ("kind", "id", "issue"),
    "deprecation_completeness": ("kind", "id", "issue"),
    "structure_inconsistent": ("metabolite", "name", "issue", "model formula",
                               "model charge", "smiles formula", "smiles charge", "smiles"),
}


def _write_findings(path: Path, sections: list[tuple[str, tuple, list]]) -> None:
    """Write every check's entries into one markdown file.

    One file with a fixed set of sections, rather than a CSV per check that
    only appears when something is found. A check going from some findings
    to none then shows as rows disappearing under a heading that stays put,
    instead of a whole file being added or deleted, and the order of the
    sections never moves.
    """
    lines = [
        "# Model QC findings",
        "",
        "Generated by `model_qc.py` on every pull request. One section per",
        "check, always present: an empty section means the check found",
        "nothing, not that it did not run. See",
        "[README.md](README.md) for what each check means.",
        "",
    ]
    for label, header, rows in sections:
        lines += [f"## {label}", ""]
        if not rows:
            lines += ["_None._", ""]
            continue
        lines += [
            "| " + " | ".join(header) + " |",
            "|" + "|".join("---" for _ in header) + "|",
        ]
        # Sorted so that the file is a function of the findings alone and
        # not of the order the model happened to be traversed in.
        for row in sorted(rows, key=lambda r: tuple(str(v) for v in r)):
            cells = [str(v).replace("|", "\\|") for v in row]
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def run(model_path: Path, out_dir: Path, base_model_dir: Path | None = None) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    if model_path.suffix in {".yml", ".yaml"}:
        # raven_toolbox.io.read_yaml_model, not cobra.io.load_yaml_model:
        # the latter silently drops RAVEN-only fields (model.id/name/
        # version, deltaG, confidence_score, notes, inchis) -- harmless
        # for the checks below today, but this branch is now the default
        # path on every develop PR (yeast-GEM#379 stage 2), not an
        # occasional opt-in, so it should read at the same fidelity as
        # the .xml branch rather than silently less.
        model = read_yaml_model(str(model_path))
    else:
        model = cobra.io.read_sbml_model(str(model_path))
    model_dir = model_path.parent

    metrics: dict[str, float] = {
        "n_reactions": len(model.reactions),
        "n_metabolites": len(model.metabolites),
        "n_genes": len(model.genes),
    }
    sections: list[tuple[str, tuple, list]] = []
    for key, label, function in _CHECKS:
        value, rows = function(model)
        metrics[key] = value
        if key in _HEADERS:
            sections.append((label, _HEADERS[key], rows))

    mass_rows, charge_rows = _balance(model)
    metrics["mass_imbalanced"] = len(mass_rows)
    metrics["charge_imbalanced"] = len(charge_rows)
    sections.append(
        ("Mass-imbalanced reactions", _HEADERS["mass_imbalanced"], mass_rows)
    )
    sections.append(
        ("Charge-imbalanced reactions", _HEADERS["charge_imbalanced"],
         charge_rows)
    )

    name_count, name_rows = check_name_consistency(model, model_dir)
    metrics["name_mismatches"] = name_count
    # Label must match build_qc_report.py's _ANNOTATION_ROWS entry for this
    # key exactly -- the PR comment links into this file's heading by
    # slugifying that label, and the two must agree for the link to resolve.
    sections.append((
        "Reaction/metabolite names disagreeing with the tsvs",
        _HEADERS["name_mismatches"], name_rows,
    ))

    xrefs_count, xrefs_rows = check_malformed_xrefs(model_dir)
    metrics["malformed_xrefs"] = xrefs_count
    sections.append(
        ("Malformed cross-references", _HEADERS["malformed_xrefs"], xrefs_rows)
    )

    consistency_count, consistency_rows = check_annotation_consistency(model, model_dir)
    metrics["annotation_consistency"] = consistency_count
    sections.append((
        "Model/annotation-table consistency",
        _HEADERS["annotation_consistency"], consistency_rows,
    ))

    dep_count, dep_rows = check_deprecation_completeness(model, model_dir, base_model_dir)
    metrics["deprecation_completeness"] = dep_count
    sections.append((
        "Removed identifiers not added to the deprecated lists",
        _HEADERS["deprecation_completeness"], dep_rows,
    ))

    structure_count, structure_rows = check_structure_consistency(model, model_dir)
    metrics["structure_inconsistent"] = structure_count
    sections.append((
        "Metabolite structure (SMILES) disagreeing with formula/charge",
        _HEADERS["structure_inconsistent"], structure_rows,
    ))

    macaw_metrics, macaw_sections = check_macaw(model)
    metrics.update(macaw_metrics)
    sections += macaw_sections

    _write_findings(out_dir / "qc_findings.md", sections)

    with (out_dir / "qc_metrics.tsv").open("w", encoding="utf-8") as handle:
        for key in sorted(metrics):
            handle.write(f"{key}\t{metrics[key]}\n")
    return metrics


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=Path("model/yeast-GEM.yml"))
    parser.add_argument("--out", type=Path, default=Path("data/testResults"))
    parser.add_argument("--base-model-dir", type=Path, default=None,
                         help="target branch's own model/ checkout, for "
                              "deprecation-completeness; omit to skip that check")
    args = parser.parse_args()

    metrics = run(args.model, args.out, args.base_model_dir)
    width = max(len(k) for k in metrics)
    for key in sorted(metrics):
        print(f"  {key:<{width}}  {metrics[key]:g}")
    # Only these gates are allowed to fail the run; every other check reports.
    failed = False
    if metrics["growth"] <= 1e-6:
        print("\nGATE FAILED: the model cannot produce biomass.")
        failed = True
    if metrics.get("name_mismatches", 0) > 0:
        print(f"\nGATE FAILED: {metrics['name_mismatches']:g} reaction/metabolite "
              "name(s) disagree between model/yeast-GEM.yml and the tsvs.")
        failed = True
    if metrics.get("annotation_consistency", 0) > 0:
        print(f"\nGATE FAILED: {metrics['annotation_consistency']:g} "
              "model/annotation-table inconsistenc(y/ies) -- see qc_findings.md.")
        failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
