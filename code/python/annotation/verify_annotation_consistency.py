"""Check whether a metabolite's or reaction's OWN cross-reference columns
agree with each other, via MetaNetX's id graph -- independent of
verify_annotations.py's structure/smiles check, and complementary to it: a
missing or unparsable smiles leaves that check with nothing to compare
against, but two existing ids simply not sharing a MetaNetX identity is
evidence on its own (e.g. chebi and kegg.compound pointing at different
compounds).

Where possible each conflict is resolved into a concrete suggested fix, in
order: (1) the entity's own smiles/participant-structure match, when
verify_annotations.py's structure check would also confirm one side -- the
strongest evidence, independent of any curated id; (2) failing that, a
strict majority of the entity's own other columns -- weaker, but still more
than an arbitrary pick. Both compare at the *skeleton* level (see
metanetx.skeleton()): MetaNetX commonly registers a generic entry and a
specific stereoisomer/protonation microspecies of the same real compound as
separate ids, and one column citing the generic form while another cites
the specific form is normal curation, not a conflict. Neither may resolve
it -- e.g. exactly two columns disagree and there is no smiles: reported
with no suggestion, alongside MetaNetX's own name/formula for what each
column's value actually names, for a human to judge -- chemical names are
not safely comparable by string matching, so this is shown as context, not
auto-decided. A resolved conflict can still be a genuine structural
difference and not a data-entry mistake -- a ring/open-chain sugar tautomer
or a monatomic ion vs. its neutral atom register as fully distinct MetaNetX
ids even though a curator reasonably treats them as the same metabolite.

Run it on the whole model, or scoped to specific ids:

    python code/python/annotation/verify_annotation_consistency.py --all
    python code/python/annotation/verify_annotation_consistency.py --mets s_0001,s_0002

Unlike verify_annotations.py, this is NOT run by model-qc.yml: most of what
it finds predates any one pull request -- annotation drift accumulates
across the whole tsv over time, not from a single edit -- so it makes a
poor per-pull-request signal and a good occasional one. Run it yourself
after curating a batch of annotation, or before a release, and work through
the report (default --out data/annotation/consistency_report.md) by hand.
No --fix: even a resolved conflict can touch multiple columns at once,
which deserves a final look before committing, unlike
verify_annotations.py's single-column missing/drift corrections.

The MetaNetX reference tables are downloaded once to a cache directory
($MNX_CACHE, default data/databases/.metanetx), shared with
verify_annotations.py.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import cobra
from metanetx import (
    MET_DB,
    RXN_DB,
    YAML,
    _best,
    _best_mnx,
    _ids,
    _norm_met,
    _water_and_proton_ids,
    ensure_metanetx,
    load_met_rows,
    load_metanetx,
    load_metkey,
    load_reaction_stoich,
    load_rxn_rows,
    match_metabolites_by_structure,
    match_reactions_by_structure,
    reverse_xref,
    skeleton,
)


# --------------------------------------------------------------------------- #
# Verification
# --------------------------------------------------------------------------- #
def verify_metabolite_conflicts(rows, mnx, rev, best_mnx_by_id, ids):
    """Flag metabolites whose OWN cross-reference columns disagree with each
    other about which MetaNetX chemical they identify -- independent of
    structure/SMILES, via MetaNetX's id graph alone. Catches rows
    verify_annotations.py's structure check cannot: a missing or unparsable
    smiles leaves it with nothing to match against, but two existing ids
    simply not sharing any MetaNetX chemical is evidence on its own.

    Groups columns at the skeleton level (see metanetx.skeleton()) into
    clusters that agree with each other, then tries to pick a winning
    cluster two ways, in order: (1) this metabolite's own smiles, if
    match_metabolites_by_structure() matched one (best_mnx_by_id) -- the
    strongest evidence, independent of any curated id; (2) a strict
    majority of this metabolite's own columns -- weaker, but still more
    than an arbitrary pick. A losing column gets a suggested replacement
    drawn from the winning cluster's own MetaNetX cross-references. Neither
    may resolve it (e.g. exactly two columns, no smiles): reported with no
    suggestion, for a human to weigh using the name/formula context
    write_report() attaches to each column.
    """
    _, mnx2key, mnx2xref, *_ = mnx
    findings = []
    for r in rows:
        if ids and r["id"] not in ids:
            continue
        cols = {}
        mnxm_direct = _ids(r.get("metanetx.chemical", ""))
        if mnxm_direct:
            cols["metanetx.chemical"] = {"value": r["metanetx.chemical"], "mnxm": mnxm_direct}
        for col, dbs in MET_DB.items():
            have = _ids(r.get(col, ""))
            if not have:
                continue
            mnxset = set()
            for hid in have:
                bare = hid.replace("CHEBI:", "") if col == "chebi" else hid
                for db in dbs:
                    mnxset |= rev.get((db, bare), set())
            if mnxset:
                cols[col] = {"value": r[col], "mnxm": mnxset}
        if len(cols) < 2:
            continue

        clusters = []
        for col, info in cols.items():
            skel = {skeleton(mnx2key.get(m, "")) for m in info["mnxm"]} - {""}
            if not skel:
                continue
            for c in clusters:
                if c["skel"] & skel:
                    c["skel"] |= skel
                    c["mnxm"] |= info["mnxm"]
                    c["cols"].append(col)
                    break
            else:
                clusters.append({"skel": skel, "mnxm": set(info["mnxm"]), "cols": [col]})
        if len(clusters) < 2:
            continue

        winner, basis = None, None
        structure_mnxm = best_mnx_by_id.get(r["id"])
        structure_skel = skeleton(mnx2key.get(structure_mnxm, "")) if structure_mnxm else ""
        if structure_skel:
            winner = next((c for c in clusters if structure_skel in c["skel"]), None)
            basis = "structure" if winner else None
        if winner is None:
            clusters.sort(key=lambda c: -len(c["cols"]))
            if len(clusters[0]["cols"]) > len(clusters[1]["cols"]):
                winner, basis = clusters[0], "majority"

        winner_ref = _best_mnx(winner["mnxm"], mnx2xref) if winner else None
        suggestions = {}
        if winner:
            for col in cols:
                if col in winner["cols"]:
                    continue
                if col == "metanetx.chemical":
                    suggestions[col] = winner_ref
                    continue
                pool = [e for m in winner["mnxm"] for db in MET_DB[col]
                        for e in mnx2xref.get(m, {}).get(db, [])]
                suggestions[col] = _best([_norm_met(e, col) for e in pool], col) if pool else ""

        findings.append({
            "kind": "met", "id": r["id"],
            "columns": {col: {"value": info["value"], "ref": _best_mnx(info["mnxm"], mnx2xref)}
                        for col, info in cols.items()},
            "winner_cols": winner["cols"] if winner else [],
            "winner_ref": winner_ref,
            "basis": basis,
            "suggestions": suggestions,
        })
    return findings


def verify_reaction_conflicts(rows, mnxr2sig, mnxr2xref, rev, best_mnxr_by_id, ids):
    """Reaction-side equivalent of verify_metabolite_conflicts(): flags
    reactions whose kegg.reaction / bigg.reaction / metanetx.reaction
    columns don't share a common MetaNetX reaction, independent of the
    participant-structure signature check, and resolves it the same way --
    this reaction's own stoichiometry match (best_mnxr_by_id) first, a
    majority of its own columns otherwise. Compares via each MNXR's own
    participant-structure signature (mnxr2sig, the inverse of sig2mnxr)
    rather than raw MNXR ids, so two ids for what MetaNetX considers the
    same reaction (e.g. a duplicate/merged MNXR entry) do not read as a
    conflict -- the same reasoning as the metabolite side's skeleton
    comparison, one level up. No MetaNetX name/formula exists for a
    reaction the way it does for a metabolite, so a resolved reaction
    conflict's suggestion is the only extra context this can offer.
    """
    findings = []
    for r in rows:
        if ids and r["id"] not in ids:
            continue
        cols = {}
        mnxr_direct = _ids(r.get("metanetx.reaction", ""))
        if mnxr_direct:
            cols["metanetx.reaction"] = {"value": r["metanetx.reaction"], "mnxr": mnxr_direct}
        for col, dbs in RXN_DB.items():
            have = _ids(r.get(col, ""))
            if not have:
                continue
            mnxrset = set()
            for hid in have:
                for db in dbs:
                    mnxrset |= rev.get((db, hid), set())
            if mnxrset:
                cols[col] = {"value": r[col], "mnxr": mnxrset}
        if len(cols) < 2:
            continue

        clusters = []
        for col, info in cols.items():
            sigs = {mnxr2sig.get(m) for m in info["mnxr"]} - {None}
            if not sigs:
                continue
            for c in clusters:
                if c["sig"] & sigs:
                    c["sig"] |= sigs
                    c["mnxr"] |= info["mnxr"]
                    c["cols"].append(col)
                    break
            else:
                clusters.append({"sig": sigs, "mnxr": set(info["mnxr"]), "cols": [col]})
        if len(clusters) < 2:
            continue

        winner, basis = None, None
        structure_mnxr = best_mnxr_by_id.get(r["id"])
        structure_sig = mnxr2sig.get(structure_mnxr) if structure_mnxr else None
        if structure_sig is not None:
            winner = next((c for c in clusters if structure_sig in c["sig"]), None)
            basis = "structure" if winner else None
        if winner is None:
            clusters.sort(key=lambda c: -len(c["cols"]))
            if len(clusters[0]["cols"]) > len(clusters[1]["cols"]):
                winner, basis = clusters[0], "majority"

        winner_ref = _best_mnx(winner["mnxr"], mnxr2xref) if winner else None
        suggestions = {}
        if winner:
            for col in cols:
                if col in winner["cols"]:
                    continue
                if col == "metanetx.reaction":
                    suggestions[col] = winner_ref
                    continue
                pool = [e for m in winner["mnxr"] for db in RXN_DB[col]
                        for e in mnxr2xref.get(m, {}).get(db, [])]
                suggestions[col] = _best(pool, col) if pool else ""

        findings.append({
            "kind": "rxn", "id": r["id"],
            "columns": {col: {"value": info["value"], "ref": _best_mnx(info["mnxr"], mnxr2xref)}
                        for col, info in cols.items()},
            "winner_cols": winner["cols"] if winner else [],
            "winner_ref": winner_ref,
            "basis": basis,
            "suggestions": suggestions,
        })
    return findings


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
def write_report(report_path: Path, conflicts, scope: str, context: dict) -> None:
    """Write every conflict finding as a diff-friendly markdown table,
    grouped by kind (metabolite/reaction) and then by whether a suggested
    fix could be worked out."""
    names, mnx_info = context["names"], context["mnx_info"]
    resolved_n = sum(1 for c in conflicts if c["basis"])
    lines = [
        "# Annotation consistency report",
        "",
        "Generated by "
        "`code/python/annotation/verify_annotation_consistency.py`. A "
        "snapshot, not auto-refreshed by CI -- re-run to update after "
        "curating a batch of annotation. See `verify_annotations.py` for "
        "the separate structure-vs-MetaNetX check (`wrong`/`missing`/"
        "`drift`) that runs on every pull request.",
        "",
        "## Method",
        "",
        "A metabolite's or reaction's own cross-reference columns (e.g. "
        "`chebi` vs `kegg.compound`) are checked against each other via "
        "MetaNetX's id graph, whether or not a smiles/structure match "
        "exists. Where possible a disagreement is resolved into a "
        "suggested fix -- from the entity's own structure match, or "
        "otherwise a majority of its other columns -- and the value(s) it "
        "is based on are marked **bold**; where neither settles it, the "
        "MetaNetX name/formula behind each column's value is shown for a "
        "human to compare instead. See the script's own docstring for the "
        "full method, including why a structural difference (e.g. a "
        "sugar's ring/open-chain tautomer) is not always a curation "
        "mistake.",
        "",
        f"Scope: {scope}.",
        "",
        "## Summary",
        "",
        "| | count |",
        "|---|---:|",
        f"| Suggested fix available | {resolved_n} |",
        f"| No consensus -- compare names | {len(conflicts) - resolved_n} |",
        f"| **Total** | {len(conflicts)} |",
        "",
    ]

    def _ref_label(kind, ref):
        """What MetaNetX believes a column's value actually names, for a
        reader who cannot tell "CHEBI:59457" from "prbamp" apart otherwise.
        Reactions have no MetaNetX name/formula to show (see
        verify_reaction_conflicts), so this is the bare id for those."""
        if not ref:
            return "_(unknown to MetaNetX)_"
        if kind != "met":
            return ref
        info = mnx_info.get(ref, {})
        name = info.get("name", "").replace("|", "\\|")
        return f"{name} ({info.get('formula', '?')})" if name else "_(unknown to MetaNetX)_"

    def _conflict_row(c, cols_order):
        name = names.get(c["id"], "").replace("|", "\\|")
        cells = [c["kind"], c["id"], name]
        for col in cols_order:
            info = c["columns"].get(col)
            if info is None:
                cells.append("_(none)_")
                continue
            value = info["value"].replace("|", "\\|")
            cell = f"{value} -- {_ref_label(c['kind'], info['ref'])}"
            cells.append(f"**{cell}**" if col in c["winner_cols"] else cell)
        return cells

    def _conflict_fix(c):
        label = _ref_label(c["kind"], c["winner_ref"])
        parts = [f"`{col}` -> `{val}`" if val else f"`{col}`: no id known to MetaNetX"
                 for col, val in c["suggestions"].items()]
        return f"Use {label}: " + "; ".join(parts)

    def _conflict_table(kind, cols_order):
        rows = [c for c in conflicts if c["kind"] == kind]
        if not rows:
            return ["_None._"]
        resolved = sorted((c for c in rows if c["basis"]), key=lambda c: c["id"])
        unresolved = sorted((c for c in rows if not c["basis"]), key=lambda c: c["id"])
        out = []
        if resolved:
            header = ["kind", "id", "name", *cols_order, "suggested fix"]
            out += [f"### Suggested fix available ({len(resolved)})", "",
                    "| " + " | ".join(header) + " |",
                    "|" + "|".join("---" for _ in header) + "|"]
            out += [f"| {' | '.join(_conflict_row(c, cols_order))} | {_conflict_fix(c)} |"
                    for c in resolved]
            out.append("")
        if unresolved:
            header = ["kind", "id", "name", *cols_order]
            out += [f"### No consensus -- compare names ({len(unresolved)})", "",
                    "| " + " | ".join(header) + " |",
                    "|" + "|".join("---" for _ in header) + "|"]
            out += [f"| {' | '.join(_conflict_row(c, cols_order))} |" for c in unresolved]
            out.append("")
        return out

    lines += ["## Metabolites", ""]
    lines += _conflict_table("met", ("bigg.metabolite", "chebi", "kegg.compound",
                                      "metanetx.chemical"))
    lines += ["", "## Reactions", ""]
    lines += _conflict_table("rxn", ("bigg.reaction", "kegg.reaction", "metanetx.reaction"))

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true", help="check the whole model")
    ap.add_argument("--mets", help="comma-separated metabolite ids to check")
    ap.add_argument("--rxns", help="comma-separated reaction ids to check")
    ap.add_argument("--out", type=Path, default=Path("data/annotation"),
                     help="directory to write consistency_report.md into "
                          "(default: %(default)s)")
    args = ap.parse_args()
    if not (args.all or args.mets or args.rxns):
        ap.error("choose --all, --mets and/or --rxns")
    report_path = args.out / "consistency_report.md"

    ensure_metanetx()
    mnx = load_metanetx()
    _, _, mnx2xref, sig2mnxr, mnxr2xref, mnx2info = mnx
    mnxr2sig = {m: sig for sig, ms in sig2mnxr.items() for m in ms}
    met_ids = set(args.mets.split(",")) if args.mets else None
    rxn_ids = set(args.rxns.split(",")) if args.rxns else None

    met_rows = load_met_rows()
    rxn_rows = load_rxn_rows()
    names = {r["id"]: r.get("name", "").strip() for r in met_rows}
    names.update({r["id"]: r.get("name", "").strip() for r in rxn_rows})

    conflicts = []
    if args.all or args.mets:
        rev = reverse_xref(mnx2xref)
        matched = match_metabolites_by_structure(met_rows, mnx)
        best_mnx = {mid: m["best"] for mid, m in matched.items()}
        conflicts += verify_metabolite_conflicts(met_rows, mnx, rev, best_mnx, met_ids)
    if args.all or args.rxns:
        model = cobra.io.load_yaml_model(str(YAML))
        metkey = load_metkey(met_rows)
        skip_ids = _water_and_proton_ids(model)
        stoich = load_reaction_stoich(model, skip_ids)
        rxn_matched = match_reactions_by_structure(stoich, mnx, metkey)
        best_mnxr = {rid: m["best"] for rid, m in rxn_matched.items()}
        conflicts += verify_reaction_conflicts(rxn_rows, mnxr2sig, mnxr2xref,
                                                reverse_xref(mnxr2xref), best_mnxr, rxn_ids)

    if conflicts:
        resolved = sum(1 for c in conflicts if c["basis"])
        print(f"conflict ({len(conflicts)}) -- an entity's own columns disagree "
              f"with each other; {resolved} with a suggested fix, "
              f"{len(conflicts) - resolved} with no consensus")
    else:
        print("conflict (0)")

    scope = "whole model" if args.all else ", ".join(
        s for s in (f"reactions {args.rxns}" if args.rxns else "",
                    f"metabolites {args.mets}" if args.mets else "") if s)
    context = {"names": names, "mnx_info": mnx2info}
    write_report(report_path, conflicts, scope, context)
    print(f"\nFull report written to {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
