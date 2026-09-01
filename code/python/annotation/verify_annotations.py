"""Verify metabolite and reaction cross-references against chemical structure.

MetaNetX (MNXref) is the hub. chem_prop maps a structure (InChIKey) to an MNXM
and chem_xref maps that MNXM to every other database's identifier; reac_prop /
reac_xref do the same for reactions, matched by their participant-structure
signature. Each cross-reference the tsvs hold is then classified:

  confirmed  the id is among the structure-verified ids
  wrong      the id resolves to a DIFFERENT compound/reaction (a real mistake)
  missing    a structure-verified id the tsv does not have yet
  drift      (MetaNetX only) the id is deprecated; a current id exists

Run it on the whole model, or on just the metabolites / reactions a pull
request changed, to catch new annotation mistakes -- e.g. after adding a
reaction:

    python code/python/annotation/verify_annotations.py --all
    python code/python/annotation/verify_annotations.py --rxns r_0001,r_0002
    python code/python/annotation/verify_annotations.py --mets s_0001 --fix

--fix applies only the safe corrections: add missing ids and update deprecated
MetaNetX ids, one value per metabolite across compartments. "wrong" ids are
reported, never overwritten automatically, because replacing a curated id
needs a human judgement (a stereochemistry-only difference, for instance, is
not a mistake).

Every run writes a full report (annotation_report.md) listing every
finding, not just the console's first-40-per-status preview, plus a
key<TAB>value summary (annotation_metrics.tsv) that build_qc_report.py
reads for the pull-request comment -- both into --out (default
data/testResults, matching model_qc.py/memote_snapshot.py). "wrong"
findings are additionally split by how many distinct entities independently
suggest the exact same correction: many entities landing on the same
(column, suggested value) is the fingerprint of a shared
structural-representation limitation -- most often a polymer/repeat-unit or
generic (R-group) entry whose SMILES matches whatever simple molecule it
was drawn as, not the thing it represents -- so those are grouped as
"recurring" for a batch look, rather than mixed in with "isolated" findings
that are unique and more directly actionable one at a time. This is a
repetition count, not a verdict: recurring findings are not assumed to be
false positives, just worth checking as a group first.

The MetaNetX reference tables are downloaded once to a cache directory
($MNX_CACHE, default data/databases/.metanetx) and reused. Run with --all
by model-qc.yml's own annotation job on every pull request (see
data/testResults/README.md) -- unlike the other model_qc.py checks it is
not measured on the target branch in the same run (the MetaNetX download
makes that costly for little benefit) and is not a merge gate: "wrong"
needs a human decision, not a pass/fail rule, and can be driven by
MetaNetX's own reference data changing over time rather than by what a
pull request did. Also runnable directly, e.g. scoped to just the
metabolites/reactions a pull request touched with --mets/--rxns.

verify_annotation_consistency.py, a separate script, checks something this
one does not: whether an entity's OWN cross-reference columns (e.g. chebi
and kegg.compound) agree with each other, independent of structure/smiles.
That is a slower, more exploratory check -- most of what it finds predates
any one pull request, since annotation drift accumulates across the whole
tsv over time rather than from a single edit -- meant to be run
occasionally by a curator, not on every pull request; see its own
docstring.

Not implemented: gene verification (UniProt ids aren't chemical structures,
so MetaNetX's structure-matching approach doesn't apply to genes.tsv) and
ec-code/kegg.pathway (classification codes, not per-molecule identifiers a
structure match can confirm).
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import cobra
from metanetx import (
    MET_DB,
    MET_TSV,
    RXN_DB,
    RXN_TSV,
    YAML,
    _best,
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
    pragmatic,
    reverse_xref,
    skeleton,
)


# --------------------------------------------------------------------------- #
# Verification
# --------------------------------------------------------------------------- #
def verify_metabolites(rows, matched, mnx, ids, rev):
    _, mnx2key, mnx2xref, *_ = mnx
    findings = []
    best_by_id = {}
    for r in rows:
        if ids and r["id"] not in ids:
            continue
        entry = matched.get(r["id"])
        if not entry:
            continue
        ik, mnxset, best_mnx = entry["ik"], entry["mnxset"], entry["best"]
        best_by_id[r["id"]] = best_mnx
        verified = defaultdict(list)
        for m in mnxset:
            for db, es in mnx2xref.get(m, {}).items():
                verified[db].extend(es)
        # MetaNetX drift / missing
        have_mnx = _ids(r.get("metanetx.chemical", ""))
        for hid in have_mnx:
            if mnx2key.get(hid) != pragmatic(ik):
                findings.append(("met", r["id"], "metanetx.chemical", "drift", hid, best_mnx))
        if not have_mnx:
            findings.append(("met", r["id"], "metanetx.chemical", "missing", "", best_mnx))
        # other namespaces
        for col, dbs in MET_DB.items():
            pool = [e for db in dbs for e in verified.get(db, [])]
            if not pool:
                continue
            vset = {_norm_met(e, col) for e in pool}
            have = _ids(r.get(col, ""))
            if not have:
                findings.append(("met", r["id"], col, "missing", "",
                                  _best([_norm_met(e, col) for e in pool], col)))
                continue
            for hid in have:
                if hid in vset:
                    continue
                # only a real mistake if the id maps to a DIFFERENT skeleton
                old_skels = {skeleton(mnx2key.get(m, ""))
                             for db in dbs for m in rev.get((db, hid.replace("CHEBI:", "")), ())}
                if old_skels and skeleton(ik) not in old_skels:
                    findings.append(("met", r["id"], col, "wrong", hid,
                                      _best([_norm_met(e, col) for e in pool], col)))
    return findings, best_by_id


def verify_reactions(rxn_rows, matched, mnx, ids):
    _, _, _, _, mnxr2xref, _ = mnx
    by_id = {r["id"]: r for r in rxn_rows}
    findings = []
    best_by_id = {}
    for rid, entry in matched.items():
        if ids and rid not in ids:
            continue
        mnxrset, best_mnxr = entry["mnxrset"], entry["best"]
        best_by_id[rid] = best_mnxr
        verified = defaultdict(list)
        for mnxr in mnxrset:
            for db, es in mnxr2xref.get(mnxr, {}).items():
                verified[db].extend(es)
        row = by_id.get(rid, {})
        have_mnx = _ids(row.get("metanetx.reaction", ""))
        if have_mnx and not (have_mnx & mnxrset):
            findings.append(("rxn", rid, "metanetx.reaction", "drift",
                              ";".join(sorted(have_mnx)), best_mnxr))
        elif not have_mnx:
            findings.append(("rxn", rid, "metanetx.reaction", "missing", "", best_mnxr))
        for col, dbs in RXN_DB.items():
            pool = [e for db in dbs for e in verified.get(db, [])]
            if not pool:
                continue
            vset = set(pool)
            have = _ids(row.get(col, ""))
            if not have:
                findings.append(("rxn", rid, col, "missing", "", _best(pool, col)))
            elif not (have & vset):
                findings.append(("rxn", rid, col, "wrong",
                                  ";".join(sorted(have)), _best(pool, col)))
    return findings, best_by_id


# --------------------------------------------------------------------------- #
# Apply (safe operations only)
# --------------------------------------------------------------------------- #
def _apply(tsv_path, findings):
    rows = list(csv.reader(tsv_path.open(encoding="utf-8", newline=""), delimiter="\t"))
    header = rows[0]
    for col in {f[2] for f in findings if f[3] in ("missing", "drift")}:
        if col not in header:
            header.append(col)
            for r in rows[1:]:
                r.append("")
    idx = {c: i for i, c in enumerate(header)}
    by = {r[0]: r for r in rows[1:]}
    n = 0
    for _kind, ent, col, status, _old, new in findings:
        if status not in ("missing", "drift") or not new:
            continue
        r = by.get(ent)
        if not r:
            continue
        while len(r) < len(header):
            r.append("")
        r[idx[col]] = new
        n += 1
    with tsv_path.open("w", encoding="utf-8", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rows)
    return n


# --------------------------------------------------------------------------- #
# Prioritisation and reporting
# --------------------------------------------------------------------------- #
_RECURRING_THRESHOLD = 3  # distinct entities sharing a suggestion, to count as a pattern


def split_wrong_by_recurrence(wrong):
    """Split "wrong" findings into (isolated, recurring, groups).

    groups: {(col, new): [findings]} for every (column, suggested value)
    pair shared by >= _RECURRING_THRESHOLD distinct entities. isolated /
    recurring are the same findings partitioned by whether their pair is
    in groups. The threshold errs toward "isolated" (2 sharing a value is
    common enough to be coincidence) so a genuinely unusual finding is
    never buried under "probably a pattern" on weak evidence.
    """
    by_pair = defaultdict(set)
    for kind, ent, col, _status, _old, new in wrong:
        by_pair[(col, new)].add((kind, ent))
    groups = {pair: ents for pair, ents in by_pair.items() if len(ents) >= _RECURRING_THRESHOLD}
    isolated, recurring = [], []
    for f in wrong:
        pair = (f[2], f[5])
        (recurring if pair in groups else isolated).append(f)
    return isolated, recurring, groups


def write_report(report_path: Path, findings, isolated, recurring, groups,
                  scope: str, context: dict) -> None:
    """Write the full findings list -- every one, not the console's
    first-40-per-status preview -- as a diff-friendly markdown table."""
    names, best_mnx, mnx_info, model_charge = (
        context["names"], context["best_mnx"], context["mnx_info"], context["model_charge"])
    by_status = Counter(f[3] for f in findings)
    lines = [
        "# Annotation verification report",
        "",
        "Generated by `code/python/annotation/verify_annotations.py`. A "
        "snapshot, not auto-refreshed by CI -- re-run to update after "
        "curating a batch of annotation. See "
        "`verify_annotation_consistency.py` for a separate, occasional "
        "check of whether an entity's own columns agree with each other.",
        "",
        "## Method",
        "",
        "Metabolites are matched to MetaNetX by structure (`smiles` -> "
        "InChIKey); reactions by their participant-structure signature -- "
        "not by name or formula. Statuses:",
        "",
        "| status | meaning |",
        "|---|---|",
        "| `wrong` | an existing id resolves to a *different* "
        "compound/reaction -- needs a human decision |",
        "| `missing` | a structure-verified id not yet recorded -- safe "
        "to add |",
        "| `drift` | a stored MetaNetX id is deprecated -- safe to "
        "update |",
        "",
        "`wrong` findings repeated identically across 3+ entities are "
        "grouped as recurring, below, rather than mixed in with unique "
        "(isolated) ones. Metabolite tables also carry MetaNetX's own "
        "name/formula/charge for the matched compound and the model's "
        "current charge, as context for the human decision `wrong` needs "
        "-- not compared automatically, since neither is reliably judged "
        "by string/number matching alone.",
        "",
        f"Scope: {scope}.",
        "",
        "## Summary",
        "",
        "| status | count |",
        "|---|---:|",
        f"| `wrong`, isolated -- review individually | {len(isolated)} |",
        f"| `wrong`, recurring -- review as a group | {len(recurring)} |",
        f"| `drift` -- safe to `--fix` | {by_status.get('drift', 0)} |",
        f"| `missing` -- safe to `--fix` | {by_status.get('missing', 0)} |",
        "",
    ]

    def _met_context(ent):
        mnxm = best_mnx.get(ent, "")
        info = mnx_info.get(mnxm, {}) if mnxm else {}
        ref_name = info.get("name", "").replace("|", "\\|")
        ref = (f"{ref_name} ({info.get('formula', '?')}, charge {info.get('charge', '?')})"
               if ref_name else "_(n/a)_")
        charge = model_charge.get(ent, "")
        return str(charge) if charge != "" else "_(n/a)_", ref

    def _table(rows, header=("kind", "id", "name", "column", "current", "suggested"),
               enrich=False):
        if enrich:
            header = (*header, "model charge", "MetaNetX reference")
        out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
        for kind, ent, col, _status, old, new in sorted(rows, key=lambda r: (r[2], r[1])):
            name = names.get(ent, "").replace("|", "\\|")
            cells = [kind, ent, name, col, old or "_(none)_", new]
            if enrich:
                cells += list(_met_context(ent)) if kind == "met" else ["_(n/a)_", "_(n/a)_"]
            out.append("| " + " | ".join(cells) + " |")
        return out

    lines += [
        "## Needs manual curation (`wrong`)",
        "",
        "### Isolated -- unique suggested correction, likely worth a look each",
        "",
    ]
    lines += _table(isolated, enrich=True) if isolated else ["_None._"]
    lines += ["", "### Recurring -- same correction suggested for "
                    f"{_RECURRING_THRESHOLD}+ entities, check the pattern first", ""]
    if recurring:
        for (col, new), ents in sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            group_rows = [f for f in recurring if (f[2], f[5]) == (col, new)]
            lines += [f"<details><summary><code>{col}</code> &rarr; <code>{new}</code> "
                      f"({len(ents)} entities)</summary>", ""]
            lines += _table(group_rows, enrich=True)
            lines += ["", "</details>", ""]
    else:
        lines += ["_None._", ""]

    lines += ["## Safe to auto-apply with `--fix`", "", "### `drift` (deprecated MetaNetX id)", ""]
    drift_rows = [f for f in findings if f[3] == "drift"]
    lines += _table(drift_rows) if drift_rows else ["_None._"]
    lines += ["", "### `missing` (structure-verified id not yet recorded)", ""]
    missing_rows = [f for f in findings if f[3] == "missing"]
    lines += _table(missing_rows) if missing_rows else ["_None._"]

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def write_metrics(metrics_path: Path, isolated, recurring, findings) -> None:
    """Write the ``key<TAB>value`` summary CI reads (see build_qc_report.py's
    _render_annotation()) -- the same shape as qc_metrics.tsv/
    validation_metrics.tsv, kept separate from the prose report so a wording
    change there can never silently break what CI parses."""
    by_status = Counter(f[3] for f in findings)
    rows = {
        "annotation_wrong_isolated": len(isolated),
        "annotation_wrong_recurring": len(recurring),
        "annotation_drift": by_status.get("drift", 0),
        "annotation_missing": by_status.get("missing", 0),
    }
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    metrics_path.write_text(
        "".join(f"{key}\t{value}\n" for key, value in rows.items()), encoding="utf-8"
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true", help="check the whole model")
    ap.add_argument("--mets", help="comma-separated metabolite ids to check")
    ap.add_argument("--rxns", help="comma-separated reaction ids to check")
    ap.add_argument("--fix", action="store_true",
                     help="apply safe corrections (add missing, update drift)")
    ap.add_argument("--out", type=Path, default=Path("data/testResults"),
                     help="directory to write annotation_report.md and "
                          "annotation_metrics.tsv into (default: %(default)s), "
                          "matching model_qc.py/memote_snapshot.py")
    args = ap.parse_args()
    if not (args.all or args.mets or args.rxns):
        ap.error("choose --all, --mets and/or --rxns")
    report_path = args.out / "annotation_report.md"
    metrics_path = args.out / "annotation_metrics.tsv"

    ensure_metanetx()
    mnx = load_metanetx()
    _, _, mnx2xref, _, _, mnx2info = mnx
    met_ids = set(args.mets.split(",")) if args.mets else None
    rxn_ids = set(args.rxns.split(",")) if args.rxns else None

    met_rows = load_met_rows()
    rxn_rows = load_rxn_rows()
    names = {r["id"]: r.get("name", "").strip() for r in met_rows}
    names.update({r["id"]: r.get("name", "").strip() for r in rxn_rows})
    model = cobra.io.load_yaml_model(str(YAML))
    model_charge = {m.id: m.charge for m in model.metabolites}

    findings = []
    best_mnx = {}
    if args.all or args.mets:
        rev = reverse_xref(mnx2xref)
        matched = match_metabolites_by_structure(met_rows, mnx)
        met_findings, met_best = verify_metabolites(met_rows, matched, mnx, met_ids, rev)
        findings += met_findings
        best_mnx.update(met_best)
    if args.all or args.rxns:
        metkey = load_metkey(met_rows)
        skip_ids = _water_and_proton_ids(model)
        stoich = load_reaction_stoich(model, skip_ids)
        rxn_matched = match_reactions_by_structure(stoich, mnx, metkey)
        rxn_findings, rxn_best = verify_reactions(rxn_rows, rxn_matched, mnx, rxn_ids)
        findings += rxn_findings
        best_mnx.update(rxn_best)

    wrong = [f for f in findings if f[3] == "wrong"]
    isolated, recurring, groups = split_wrong_by_recurrence(wrong)

    by_status = Counter(f[3] for f in findings)
    print(f"findings: {dict(by_status)} "
          f"(wrong: {len(isolated)} isolated, {len(recurring)} recurring "
          f"in {len(groups)} group(s))")

    if isolated:
        print(f"\nwrong, isolated -- review individually ({len(isolated)}):")
        for kind, ent, col, _s, old, new in sorted(isolated, key=lambda r: (r[2], r[1]))[:40]:
            print(f"  {kind:3} {ent:12} {col:20} {old or '-'} -> {new}")
        if len(isolated) > 40:
            print(f"  ... and {len(isolated) - 40} more -- see the report")
    if recurring:
        print(f"\nwrong, recurring -- {len(groups)} group(s) of "
              f"{_RECURRING_THRESHOLD}+ entities sharing one correction, "
              f"see the report for detail ({len(recurring)} findings total)")
    for status in ("drift", "missing"):
        rows = [f for f in findings if f[3] == status]
        if rows:
            print(f"\n{status} ({len(rows)}), safe to --fix:")
            for kind, ent, col, _s, old, new in rows[:10]:
                print(f"  {kind:3} {ent:12} {col:20} {old or '-'} -> {new}")
            if len(rows) > 10:
                print(f"  ... and {len(rows) - 10} more -- see the report")

    scope = "whole model" if args.all else ", ".join(
        s for s in (f"reactions {args.rxns}" if args.rxns else "",
                    f"metabolites {args.mets}" if args.mets else "") if s)
    context = {"names": names, "best_mnx": best_mnx, "mnx_info": mnx2info,
               "model_charge": model_charge}
    write_report(report_path, findings, isolated, recurring, groups, scope, context)
    print(f"\nFull report written to {report_path}")
    write_metrics(metrics_path, isolated, recurring, findings)
    print(f"Summary metrics written to {metrics_path}")

    if args.fix:
        met_f = [f for f in findings if f[0] == "met"]
        rxn_f = [f for f in findings if f[0] == "rxn"]
        nm = _apply(MET_TSV, met_f) if met_f else 0
        nr = _apply(RXN_TSV, rxn_f) if rxn_f else 0
        print(f"\napplied safe fixes: {nm} metabolite + {nr} reaction cells "
              f"(wrong ids left for review)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
