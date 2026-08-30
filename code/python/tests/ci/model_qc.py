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

Only two checks are gates -- the model must load, and it must grow.
Everything else is reported so that a pull request which makes a number
worse is visible, without blocking work on findings that were already
there.

Checks that Human-GEM runs and this does not yet: model/annotation-table
consistency, removed-identifiers-not-deprecated, and structure (SMILES /
InChI) versus formula and charge. All three need the YAML-based curation
layout with separate identifier TSVs, which yeast-GEM adopts in 9.2.0
(see #379).

Usage:
    python model_qc.py --model model/yeast-GEM.xml --out data/testResults
"""
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import cobra

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


def check_malformed_xrefs(model) -> tuple[int, list]:
    """Cross-references that do not match their namespace's pattern."""
    rows = []
    for kind, entities in (
        ("metabolite", model.metabolites),
        ("reaction", model.reactions),
        ("gene", model.genes),
    ):
        for entity in entities:
            for namespace, pattern in _XREF_PATTERNS.items():
                for value in _values(entity.annotation, namespace):
                    text = str(value)
                    if not re.match(pattern, text):
                        rows.append((
                            kind, entity.id, namespace, text,
                            _why_malformed(namespace, text),
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
    ("malformed_xrefs", "Malformed cross-references", check_malformed_xrefs),
    ("xrefs_across_compartments", "Cross-refs inconsistent across compartments",
     check_xrefs_across_compartments),
]

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


def run(model_path: Path, out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    if model_path.suffix in {".yml", ".yaml"}:
        model = cobra.io.load_yaml_model(str(model_path))
    else:
        model = cobra.io.read_sbml_model(str(model_path))

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
    parser.add_argument("--model", type=Path, default=Path("model/yeast-GEM.xml"))
    parser.add_argument("--out", type=Path, default=Path("data/testResults"))
    args = parser.parse_args()

    metrics = run(args.model, args.out)
    width = max(len(k) for k in metrics)
    for key in sorted(metrics):
        print(f"  {key:<{width}}  {metrics[key]:g}")
    # Only the gate is allowed to fail the run; every other check reports.
    if metrics["growth"] <= 1e-6:
        print("\nGATE FAILED: the model cannot produce biomass.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
