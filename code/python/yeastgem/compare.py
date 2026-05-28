"""Cross-language model comparator — the CI level-1 semantic-equality gate.

Compares two `cobra.Model` objects for semantic equality, ignoring
formatting differences (key ordering, whitespace, float repr). Used to
verify that the MATLAB and Python toolchains produce equivalent models.

This is one of the explicit upstream candidates (see
[UPSTREAM_CANDIDATES.md](../UPSTREAM_CANDIDATES.md)); the API shape is
chosen with future ravengem upstreaming in mind.
"""
from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field

import cobra

# Annotation keys checked by default. Add/remove via `compare_models(...,
# ignore_annotations=...)` or `extra_annotations=...`.
DEFAULT_ANNOTATION_KEYS: tuple[str, ...] = (
    "sbo",
    "ec-code",
    "kegg.reaction",
    "kegg.compound",
    "metanetx.reaction",
    "metanetx.chemical",
    "chebi",
    "bigg.reaction",
    "bigg.metabolite",
)


@dataclass
class ComparisonReport:
    """Result of :func:`compare_models`.

    Truthy if the models are semantically equal. Iterate ``differences``
    to see human-readable messages.
    """

    equal: bool
    differences: list[str] = field(default_factory=list)

    def __bool__(self) -> bool:
        return self.equal

    def __str__(self) -> str:
        if self.equal:
            return "Models are semantically equal."
        head = f"Models differ ({len(self.differences)} differences):"
        body = "\n".join(f"  - {d}" for d in self.differences)
        return f"{head}\n{body}"


def compare_models(
    a: cobra.Model,
    b: cobra.Model,
    *,
    stoichiometry_tol: float = 1e-9,
    ignore_annotations: Iterable[str] = (),
    extra_annotations: Iterable[str] = (),
    max_per_category: int = 50,
) -> ComparisonReport:
    """Compare two cobra models for semantic equality.

    Checks:

    * reaction, metabolite and gene id sets (exact)
    * stoichiometry per shared reaction (within ``stoichiometry_tol``)
    * lower/upper bounds and objective coefficients (exact)
    * GPR rules (whitespace- and case-insensitive string comparison)
    * metabolite formula, charge, compartment (exact)
    * a default set of annotation keys (see ``DEFAULT_ANNOTATION_KEYS``)
      plus any ``extra_annotations``, minus any ``ignore_annotations``

    ``max_per_category`` caps the per-category diff list so a wholesale
    mismatch produces a digestible report rather than 10 000 lines.
    """
    diffs: list[str] = []
    keys = (set(DEFAULT_ANNOTATION_KEYS) | set(extra_annotations)) - set(ignore_annotations)

    _diff_id_sets(a, b, diffs)

    common_rxns = {r.id for r in a.reactions} & {r.id for r in b.reactions}
    _diff_reactions(a, b, sorted(common_rxns), stoichiometry_tol, keys, diffs, max_per_category)

    common_mets = {m.id for m in a.metabolites} & {m.id for m in b.metabolites}
    _diff_metabolites(a, b, sorted(common_mets), keys, diffs, max_per_category)

    return ComparisonReport(equal=not diffs, differences=diffs)


# --- internals ---------------------------------------------------------

def _diff_id_sets(a: cobra.Model, b: cobra.Model, diffs: list[str]) -> None:
    for label, get in (
        ("reactions", lambda m: m.reactions),
        ("metabolites", lambda m: m.metabolites),
        ("genes", lambda m: m.genes),
    ):
        a_ids = {x.id for x in get(a)}
        b_ids = {x.id for x in get(b)}
        only_a = sorted(a_ids - b_ids)
        only_b = sorted(b_ids - a_ids)
        if only_a:
            diffs.append(f"{label} only in A ({len(only_a)}): {only_a[:10]}{_more(only_a, 10)}")
        if only_b:
            diffs.append(f"{label} only in B ({len(only_b)}): {only_b[:10]}{_more(only_b, 10)}")


def _diff_reactions(
    a: cobra.Model,
    b: cobra.Model,
    common: list[str],
    tol: float,
    anno_keys: set[str],
    diffs: list[str],
    cap: int,
) -> None:
    counters = {"stoich": 0, "bounds": 0, "objective": 0, "gpr": 0, "anno": 0}
    for rxn_id in common:
        ra = a.reactions.get_by_id(rxn_id)
        rb = b.reactions.get_by_id(rxn_id)

        sa = {m.id: c for m, c in ra.metabolites.items()}
        sb = {m.id: c for m, c in rb.metabolites.items()}
        if set(sa) != set(sb):
            _push(diffs, counters, "stoich", cap, f"{rxn_id}: stoichiometry has different mets")
        else:
            for mid in sa:
                if abs(sa[mid] - sb[mid]) > tol:
                    _push(
                        diffs, counters, "stoich", cap,
                        f"{rxn_id}: coef[{mid}] A={sa[mid]:.6g} B={sb[mid]:.6g}",
                    )

        if ra.lower_bound != rb.lower_bound or ra.upper_bound != rb.upper_bound:
            _push(
                diffs, counters, "bounds", cap,
                f"{rxn_id}: bounds A=({ra.lower_bound}, {ra.upper_bound}) "
                f"B=({rb.lower_bound}, {rb.upper_bound})",
            )

        if ra.objective_coefficient != rb.objective_coefficient:
            _push(
                diffs, counters, "objective", cap,
                f"{rxn_id}: objective A={ra.objective_coefficient} B={rb.objective_coefficient}",
            )

        ga = _normalise_gpr(ra.gene_reaction_rule)
        gb = _normalise_gpr(rb.gene_reaction_rule)
        if ga != gb:
            _push(diffs, counters, "gpr", cap, f"{rxn_id}: GPR A={ga!r} B={gb!r}")

        _diff_annotations(
            f"rxn {rxn_id}", ra.annotation, rb.annotation,
            anno_keys, diffs, counters, cap,
        )


def _diff_metabolites(
    a: cobra.Model,
    b: cobra.Model,
    common: list[str],
    anno_keys: set[str],
    diffs: list[str],
    cap: int,
) -> None:
    counters = {"formula": 0, "charge": 0, "compartment": 0, "anno": 0}
    for met_id in common:
        ma = a.metabolites.get_by_id(met_id)
        mb = b.metabolites.get_by_id(met_id)

        if (ma.formula or None) != (mb.formula or None):
            _push(
                diffs, counters, "formula", cap,
                f"met {met_id}: formula A={ma.formula} B={mb.formula}",
            )

        ca = ma.charge if ma.charge is not None else 0
        cb = mb.charge if mb.charge is not None else 0
        if ca != cb:
            _push(
                diffs, counters, "charge", cap,
                f"met {met_id}: charge A={ma.charge} B={mb.charge}",
            )

        if ma.compartment != mb.compartment:
            _push(
                diffs, counters, "compartment", cap,
                f"met {met_id}: compartment A={ma.compartment} B={mb.compartment}",
            )

        _diff_annotations(
            f"met {met_id}", ma.annotation, mb.annotation,
            anno_keys, diffs, counters, cap,
        )


def _diff_annotations(
    label: str,
    aa: dict,
    ba: dict,
    keys: set[str],
    diffs: list[str],
    counters: dict,
    cap: int,
) -> None:
    for k in keys:
        va = _normalise_annotation_value(aa.get(k))
        vb = _normalise_annotation_value(ba.get(k))
        if va != vb:
            _push(diffs, counters, "anno", cap, f"{label}.annotation[{k!r}]: A={va} B={vb}")


def _normalise_gpr(rule: str) -> str:
    """Lowercase + collapse whitespace.

    A more robust comparator would parse to an AST and compare structures;
    that is on the upstream-candidate list. For now string normalisation
    catches the formatting drift we expect from different SBML writers.
    """
    if not rule:
        return ""
    return " ".join(rule.lower().split())


def _normalise_annotation_value(v):
    if v is None:
        return None
    if isinstance(v, list):
        return tuple(sorted(str(x) for x in v))
    if isinstance(v, tuple):
        return tuple(sorted(str(x) for x in v))
    return str(v)


def _push(diffs: list[str], counters: dict, key: str, cap: int, msg: str) -> None:
    counters[key] = counters.get(key, 0) + 1
    if counters[key] <= cap:
        diffs.append(msg)
    elif counters[key] == cap + 1:
        diffs.append(f"... ({key} category truncated at {cap} entries)")


def _more(seq: list, shown: int) -> str:
    return "" if len(seq) <= shown else f" ... (+{len(seq) - shown} more)"


# --- CLI ---------------------------------------------------------------

def _main() -> int:
    import argparse

    from cobra.io import read_sbml_model

    parser = argparse.ArgumentParser(
        description="Compare two SBML models for semantic equality."
    )
    parser.add_argument("a", help="path to first model (SBML)")
    parser.add_argument("b", help="path to second model (SBML)")
    parser.add_argument(
        "--stoichiometry-tol",
        type=float,
        default=1e-9,
        help="absolute tolerance on stoichiometric coefficients (default 1e-9)",
    )
    args = parser.parse_args()

    report = compare_models(
        read_sbml_model(args.a),
        read_sbml_model(args.b),
        stoichiometry_tol=args.stoichiometry_tol,
    )
    print(report)
    return 0 if report.equal else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
