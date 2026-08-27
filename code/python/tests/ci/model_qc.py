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
import csv
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

# Namespaces expected on a Saccharomyces model. Anything outside this set is
# reported: an E. coli or human gene database on a yeast gene is a copy-paste
# error, not a deliberate cross-reference.
_EXPECTED_GENE_XREFS = {
    "sbo", "uniprot", "kegg.genes", "ncbigene", "ncbiprotein", "refseq", "sgd",
}


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


def check_empty_reactions(model) -> tuple[int, list]:
    rows = [(r.id, r.name or "") for r in model.reactions if not r.metabolites]
    return len(rows), rows


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
                    if not re.match(pattern, str(value)):
                        rows.append((kind, entity.id, namespace, str(value)))
    return len(rows), rows


def check_unexpected_gene_xrefs(model) -> tuple[int, list]:
    """Gene cross-references pointing at a database for another organism."""
    rows = [
        (gene.id, namespace, ";".join(str(v) for v in _values(gene.annotation, namespace)))
        for gene in model.genes
        for namespace in gene.annotation
        if namespace not in _EXPECTED_GENE_XREFS
    ]
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


def check_macaw(model, out_dir: Path) -> dict:
    """MACAW dead-end and duplicate tests.

    Optional: when MACAW is not installed the two rows are reported as
    unavailable rather than as zero, so a missing dependency cannot be
    mistaken for a clean result.
    """
    try:
        from macaw.main import dead_end_test, duplicate_test
    except ImportError:
        return {}

    dead_end_results, _edges = dead_end_test(model)
    duplicate_results, _dup_edges = duplicate_test(model)
    merged = dead_end_results.merge(duplicate_results)
    merged.to_csv(out_dir / "macaw_results.csv", index=False)

    def _flagged(frame, column):
        if column not in frame.columns:
            return 0
        values = frame[column].astype(str).str.strip().str.lower()
        return int((~values.isin({"", "ok", "nan", "none"})).sum())

    return {
        "macaw_dead_end": _flagged(merged, "dead_end_test"),
        "macaw_duplicates": _flagged(merged, "duplicate_test"),
    }


# (metric key, comment label, function)
_CHECKS = [
    ("growth", "Growth (biomass producible)", check_growth),
    ("empty_reactions", "Reactions with no metabolites", check_empty_reactions),
    ("orphan_metabolites", "Unused metabolites", check_orphan_metabolites),
    ("unused_genes", "Unused genes", check_unused_genes),
    ("missing_formula", "Metabolites missing formula", check_missing_formula),
    ("missing_charge", "Metabolites missing charge", check_missing_charge),
    ("duplicate_reactions", "Exact-duplicate reaction groups",
     check_duplicate_reactions),
    ("malformed_xrefs", "Malformed cross-references", check_malformed_xrefs),
    ("unexpected_gene_xrefs", "Gene cross-references from other organisms",
     check_unexpected_gene_xrefs),
    ("xrefs_across_compartments", "Cross-refs inconsistent across compartments",
     check_xrefs_across_compartments),
]

_HEADERS = {
    "empty_reactions": ("reaction", "name"),
    "orphan_metabolites": ("metabolite", "name"),
    "unused_genes": ("gene", "name"),
    "missing_formula": ("metabolite", "name"),
    "missing_charge": ("metabolite", "name"),
    "duplicate_reactions": ("group", "reaction", "name", "equation"),
    "malformed_xrefs": ("kind", "id", "namespace", "value"),
    "unexpected_gene_xrefs": ("gene", "namespace", "value"),
    "xrefs_across_compartments": ("metabolite name", "namespace", "values"),
    "mass_imbalanced": ("reaction", "name", "imbalance"),
    "charge_imbalanced": ("reaction", "name", "charge"),
}


def _write_csv(path: Path, header, rows) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def run(model_path: Path, out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    if model_path.suffix in {".yml", ".yaml"}:
        model = cobra.io.load_yaml_model(str(model_path))
    else:
        model = cobra.io.read_sbml_model(str(model_path))

    metrics: dict[str, float] = {}
    for key, _label, function in _CHECKS:
        value, rows = function(model)
        metrics[key] = value
        if rows:
            _write_csv(out_dir / f"qc_{key}.csv", _HEADERS[key], rows)

    mass_rows, charge_rows = _balance(model)
    metrics["mass_imbalanced"] = len(mass_rows)
    metrics["charge_imbalanced"] = len(charge_rows)
    if mass_rows:
        _write_csv(out_dir / "qc_mass_imbalanced.csv",
                   _HEADERS["mass_imbalanced"], mass_rows)
    if charge_rows:
        _write_csv(out_dir / "qc_charge_imbalanced.csv",
                   _HEADERS["charge_imbalanced"], charge_rows)

    metrics.update(check_macaw(model, out_dir))

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
