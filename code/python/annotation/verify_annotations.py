"""Verify metabolite and reaction cross-references against chemical structure.

MetaNetX (MNXref) is the hub. chem_prop maps a structure (InChIKey) to an MNXM
and chem_xref maps that MNXM to every other database's identifier; reac_prop /
reac_xref do the same for reactions, matched by their participant-structure
signature. Each cross-reference the tsvs hold is then classified:

  confirmed  the id is among the structure-verified ids
  wrong      the id resolves to a DIFFERENT compound/reaction (a real mistake)
  conflict   this entity's OWN columns disagree about its identity, checked
             via MetaNetX's id graph directly (no structure/smiles involved)
  missing    a structure-verified id the tsv does not have yet
  drift      (MetaNetX only) the id is deprecated; a current id exists

"wrong"/"missing"/"drift" need a smiles to match against; "conflict" does
not -- it only asks whether the entity's existing columns (e.g. chebi and
kegg.compound) point MetaNetX to the same chemical/reaction as each other,
so it also catches rows a missing or unparsable smiles leaves untouched by
the other three. Both checks compare at the skeleton level (ignoring
stereochemistry/protonation, see skeleton()) so a generic-vs-specific pair
of ids for the same compound is not flagged, but a "conflict" can still be
a real structural difference and not a mistake -- a ring/open-chain sugar
tautomer or a monatomic ion vs. its neutral atom register as genuinely
distinct MetaNetX ids even though they are the same metabolite in practice.
Metabolite findings additionally carry MetaNetX's own name/formula/charge
for the matched compound and the model's current charge, as extra context
for the human call "wrong"/"conflict" need -- chemical names are not
safely comparable by string matching (synonyms, stereo-descriptors) and a
charge difference is often a legitimate protonation-state choice rather
than a mistake, so neither is used to classify a finding, only to inform
the reader.

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
structural-representation limitation
-- most often a polymer/repeat-unit or generic (R-group) entry whose SMILES
matches whatever simple molecule it was drawn as, not the thing it
represents -- so those are grouped as "recurring" for a batch look, rather
than mixed in with "isolated" findings that are unique and more directly
actionable one at a time. This is a repetition count, not a verdict:
recurring findings are not assumed to be false positives, just worth
checking as a group first.

The MetaNetX reference tables are downloaded once to a cache directory
($MNX_CACHE, default data/databases/.metanetx) and reused. Run with --all
by model-qc.yml's own annotation job on every pull request (see
data/testResults/README.md) -- unlike the other model_qc.py checks it is
not measured on the target branch in the same run (the MetaNetX download
makes that costly for little benefit) and is not a merge gate: "wrong" and
"conflict" need a human decision, not a pass/fail rule, and both can be
driven by MetaNetX's own reference data changing over time rather than by
what a pull request did. Also runnable directly, e.g. scoped to just the
metabolites/reactions a pull request touched with --mets/--rxns.

Not implemented: gene verification (UniProt ids aren't chemical structures,
so MetaNetX's structure-matching approach doesn't apply to genes.tsv) and
ec-code/kegg.pathway (classification codes, not per-molecule identifiers a
structure match can confirm).
"""

import argparse
import csv
import os
import re
import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import cobra
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

MET_TSV = Path("model/metabolites.tsv")
RXN_TSV = Path("model/reactions.tsv")
YAML = Path("model/yeast-GEM.yml")
CACHE = Path(os.environ.get("MNX_CACHE", "data/databases/.metanetx"))
MNX_URL = "https://www.metanetx.org/ftp/latest"

# tsv column -> MetaNetX chem_xref / reac_xref source prefixes. Narrower than
# Human-GEM's own MET_DB/RXN_DB: yeast-GEM's tsvs (yeast-GEM#379) carry fewer
# cross-reference columns (no rhea, hmdb, lipidmaps, seed.compound).
MET_DB = {
    "kegg.compound": ("kegg.compound", "keggC", "kegg.drug", "keggD", "kegg.glycan", "keggG"),
    "chebi": ("chebi", "CHEBI"),
    "bigg.metabolite": ("bigg.metabolite",),
}
RXN_DB = {
    "kegg.reaction": ("kegg.reaction",),
    "bigg.reaction": ("bigg.reaction",),
}
SKIP_MNX = {"MNXM01", "MNXM1", "WATER", "MNXM2"}  # H+, H2O
MNXM_OK = re.compile(r"^MNXM\d+$")


# --------------------------------------------------------------------------- #
# MetaNetX reference data (download + compact extract, cached)
# --------------------------------------------------------------------------- #
def _extract(url, out_path, keep):
    """Stream a MetaNetX tsv and write only the columns `keep(fields)` returns."""
    tmp = out_path.with_suffix(".tmp")
    with urllib.request.urlopen(url) as resp, tmp.open("w", encoding="utf-8", newline="") as fh:
        for raw in resp:
            line = raw.decode("utf-8", "replace")
            if line.startswith("#"):
                continue
            row = keep(line.rstrip("\n").split("\t"))
            if row:
                fh.write("\t".join(row) + "\n")
    tmp.replace(out_path)


def ensure_metanetx():
    CACHE.mkdir(parents=True, exist_ok=True)
    jobs = {
        "chem_prop.tsv": ("mnx_chem.tsv",
                          lambda f: [f[0], f[1], f[3], f[4], f[7]] if len(f) > 7 and f[7]
                          else None),
        "chem_xref.tsv": ("mnx_xref.tsv",
                          lambda f: [f[0], f[1]] if len(f) > 1 and ":" in f[0] else None),
        "reac_prop.tsv": ("mnxr_eq.tsv",
                          lambda f: [f[0], f[1]] if len(f) > 1 and f[1] else None),
        "reac_xref.tsv": ("mnxr_xref.tsv",
                          lambda f: [f[0], f[1]] if len(f) > 1 and ":" in f[0] else None),
    }
    for src, (dst, keep) in jobs.items():
        out = CACHE / dst
        if out.exists():
            continue
        print(f"downloading {src} -> {out} (first run only) ...", file=sys.stderr, flush=True)
        _extract(f"{MNX_URL}/{src}", out, keep)


def pragmatic(ik):
    """InChIKey without the final protonation character (protonation-invariant)."""
    return ik.rsplit("-", 1)[0] if ik else ""


def skeleton(key):
    """A pragmatic() key's connectivity block alone, ignoring stereochemistry
    too -- MetaNetX frequently registers a generic entry and a specific
    stereoisomer/protonation microspecies of the same real compound as
    separate MNXM/MNXR ids (adjacent-numbered, same formula and charge,
    InChIKey differing only in the stereo layer); treating those as
    "different identities" would be noise, not a finding."""
    return key.split("-")[0] if key else ""


def load_metanetx():
    key2mnx, mnx2key, mnx2info = defaultdict(set), {}, {}
    for line in (CACHE / "mnx_chem.tsv").open(encoding="utf-8"):
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue
        mnxm, name, formula, charge, ik = parts
        mnx2info[mnxm] = {"name": name, "formula": formula, "charge": charge}
        if ik:
            key2mnx[pragmatic(ik)].add(mnxm)
            mnx2key[mnxm] = pragmatic(ik)
    mnx2xref = defaultdict(lambda: defaultdict(list))
    for line in (CACHE / "mnx_xref.tsv").open(encoding="utf-8"):
        src, _, mnxm = line.rstrip("\n").partition("\t")
        db, _, extid = src.partition(":")
        if mnxm and db and extid:
            mnx2xref[mnxm][db].append(extid)
    sig2mnxr = defaultdict(set)
    for line in (CACHE / "mnxr_eq.tsv").open(encoding="utf-8"):
        mnxr, _, eq = line.rstrip("\n").partition("\t")
        sig = _mnx_signature(eq, mnx2key)
        if sig:
            sig2mnxr[sig].add(mnxr)
    mnxr2xref = defaultdict(lambda: defaultdict(list))
    for line in (CACHE / "mnxr_xref.tsv").open(encoding="utf-8"):
        src, _, mnxr = line.rstrip("\n").partition("\t")
        db, _, extid = src.partition(":")
        if mnxr and db and extid:
            mnxr2xref[mnxr][db].append(extid)
    return key2mnx, mnx2key, mnx2xref, sig2mnxr, mnxr2xref, mnx2info


def reverse_xref(xref: dict) -> dict:
    """(namespace, external id) -> set of MetaNetX ids -- the inverse of
    mnx2xref / mnxr2xref, shared by the structure-based and the
    cross-annotation-conflict checks so it is only built once per run."""
    rev = defaultdict(set)
    for m, dbs in xref.items():
        for db, es in dbs.items():
            for e in es:
                rev[(db, e)].add(m)
    return rev


# --------------------------------------------------------------------------- #
# Structural signatures
# --------------------------------------------------------------------------- #
def inchikey(smiles):
    if not smiles:
        return ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    try:
        return Chem.MolToInchiKey(mol)
    except Exception:
        return ""


def _signature(sub, prod):
    s, p = frozenset(sub.items()), frozenset(prod.items())
    return frozenset((s, p)) if s and p and s != p else None


_EQ_TOKEN = re.compile(r"([\d.]+)\s+(\S+?)@\S+")


def _mnx_signature(eq, mnx2key):
    if " = " not in eq:
        return None
    sides = []
    for side in eq.split(" = ", 1):
        d = defaultdict(float)
        for coeff, mnxm in _EQ_TOKEN.findall(side):
            if mnxm in SKIP_MNX:
                continue
            k = mnx2key.get(mnxm)
            if not k:
                return None
            d[k] += float(coeff)
        sides.append(dict(d))
    return _signature(*sides)


# --------------------------------------------------------------------------- #
# Model loading
# --------------------------------------------------------------------------- #
def load_met_rows():
    with MET_TSV.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_rxn_rows():
    with RXN_TSV.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_reaction_stoich(model, skip_met_ids):
    """reaction id -> {met_id: coeff}, from the already-loaded model.

    Unlike Human-GEM's equivalent (a hand-rolled yml line-scanner), this
    reads the cobra Model directly -- yeast-GEM's yml no longer carries
    cross-reference annotation to work around (yeast-GEM#379), so there is
    nothing a generic yml reader can't already give here.
    """
    return {
        rxn.id: {met.id: coeff for met, coeff in rxn.metabolites.items()
                 if met.id not in skip_met_ids}
        for rxn in model.reactions
    }


def _water_and_proton_ids(model) -> set[str]:
    """Every compartment's H+ and H2O id, by name -- not hardcoded, since
    yeast-GEM ids carry no shared prefix across compartments the way
    Human-GEM's MAMnnnnn[c/m/...] ids do."""
    names = {"h+", "h2o"}
    return {m.id for m in model.metabolites if (m.name or "").strip().lower() in names}


def _ids(val):
    return {x.strip() for x in (val or "").split(";") if x.strip()}


def _norm_met(idv, col):
    if col == "chebi":
        return idv if idv.startswith("CHEBI:") else "CHEBI:" + idv
    return idv


def _best(candidates, col):
    c = Counter(candidates)
    if not c:
        return ""
    if col == "kegg.compound":
        return sorted(c, key=lambda i: (-c[i], {"C": 0, "D": 1, "G": 2}.get(i[:1], 3), i))[0]
    return sorted(c, key=lambda i: (-c[i], i))[0]


# --------------------------------------------------------------------------- #
# Verification
# --------------------------------------------------------------------------- #
def verify_metabolites(rows, mnx, ids, rev):
    key2mnx, mnx2key, mnx2xref, *_ = mnx
    findings = []
    best_by_id = {}
    for r in rows:
        if ids and r["id"] not in ids:
            continue
        ik = inchikey(r.get("smiles", "").strip())
        if not ik:
            continue
        mnxset = key2mnx.get(pragmatic(ik), set())
        if not mnxset:
            continue
        verified = defaultdict(list)
        for m in mnxset:
            for db, es in mnx2xref.get(m, {}).items():
                verified[db].extend(es)
        # MetaNetX drift / missing
        have_mnx = _ids(r.get("metanetx.chemical", ""))
        best_mnx = sorted(mnxset, key=lambda m: (-len(mnx2xref.get(m, {})), m))[0]
        best_by_id[r["id"]] = best_mnx
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


def _best_mnx(mnxset, mnx2xref):
    return sorted(mnxset, key=lambda m: (-len(mnx2xref.get(m, {})), m))[0]


def verify_metabolite_conflicts(rows, mnx, rev, best_mnx_by_id, ids):
    """Flag metabolites whose OWN cross-reference columns disagree with each
    other about which MetaNetX chemical they identify -- independent of
    structure/SMILES, via MetaNetX's id graph alone. Catches rows a
    structure check cannot: a missing or unparsable smiles leaves
    verify_metabolites() with nothing to match against, but two existing
    ids simply not sharing any MetaNetX chemical is evidence on its own.

    Groups columns at the skeleton level (see skeleton()) into clusters that
    agree with each other, then tries to pick a winning cluster two ways,
    in order: (1) this metabolite's own smiles, if verify_metabolites()
    matched one (best_mnx_by_id) -- the strongest evidence, independent of
    any curated id; (2) a strict majority of this metabolite's own columns
    -- weaker, but still more than an arbitrary pick. A losing column gets
    a suggested replacement drawn from the winning cluster's own MetaNetX
    cross-references. Neither may resolve it (e.g. exactly two columns,
    no smiles): reported with no suggestion, for a human to weigh using
    the name/formula context write_report() attaches to each column.
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


def verify_reactions(rxn_rows, stoich, mnx, metkey, ids):
    _, _, _, sig2mnxr, mnxr2xref, _ = mnx
    by_id = {r["id"]: r for r in rxn_rows}
    findings = []
    best_by_id = {}
    for rid, mets in stoich.items():
        if ids and rid not in ids:
            continue
        sub, prod, ok = defaultdict(float), defaultdict(float), True
        for met, coeff in mets.items():
            k = metkey.get(met)
            if not k:
                ok = False
                break
            (sub if coeff < 0 else prod)[k] += abs(coeff)
        sig = _signature(dict(sub), dict(prod)) if ok else None
        mnxrset = sig2mnxr.get(sig) if sig else None
        if not mnxrset:
            continue
        verified = defaultdict(list)
        for mnxr in mnxrset:
            for db, es in mnxr2xref.get(mnxr, {}).items():
                verified[db].extend(es)
        row = by_id.get(rid, {})
        best_mnxr = sorted(mnxrset, key=lambda m: (-len(mnxr2xref.get(m, {})), m))[0]
        best_by_id[rid] = best_mnxr
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


def write_report(report_path: Path, findings, isolated, recurring, groups, conflicts,
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
        "curating a batch of annotation.",
        "",
        "## Method",
        "",
        "Metabolites are matched to MetaNetX by structure (`smiles` -> "
        "InChIKey); reactions by their participant-structure signature -- "
        "not by name or formula. Separately, an entity's own columns (e.g. "
        "`chebi` vs `kegg.compound`) are checked against each other via "
        "MetaNetX's id graph, even without a smiles. Statuses:",
        "",
        "| status | meaning |",
        "|---|---|",
        "| `wrong` | an existing id resolves to a *different* "
        "compound/reaction |",
        "| `conflict` | this entity's own columns disagree about its "
        "identity |",
        "| `missing` | a structure-verified id not yet recorded -- safe "
        "to add |",
        "| `drift` | a stored MetaNetX id is deprecated -- safe to "
        "update |",
        "",
        "`wrong`/`conflict` need a human decision; `missing`/`drift` are "
        "safe with `--fix`. `wrong` findings repeated identically across "
        "3+ entities are grouped as recurring, below. Where possible a "
        "`conflict` also gets a suggested fix, from this entity's own "
        "structure match or a majority of its other columns; where "
        "neither settles it, the names/formulas shown are the best "
        "evidence available. See the script's own docstring for the full "
        "method, including why a structural difference (e.g. a sugar's "
        "ring/open-chain tautomer) is not always a curation mistake.",
        "",
        f"Scope: {scope}.",
        "",
        "## Summary",
        "",
        "| status | count |",
        "|---|---:|",
        f"| `wrong`, isolated -- review individually | {len(isolated)} |",
        f"| `wrong`, recurring -- review as a group | {len(recurring)} |",
        f"| `conflict` -- review individually | {len(conflicts)} |",
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
            out += [f"#### Suggested fix available ({len(resolved)})", "",
                    "| " + " | ".join(header) + " |",
                    "|" + "|".join("---" for _ in header) + "|"]
            out += [f"| {' | '.join(_conflict_row(c, cols_order))} | {_conflict_fix(c)} |"
                    for c in resolved]
            out.append("")
        if unresolved:
            header = ["kind", "id", "name", *cols_order]
            out += [f"#### No consensus -- compare names ({len(unresolved)})", "",
                    "| " + " | ".join(header) + " |",
                    "|" + "|".join("---" for _ in header) + "|"]
            out += [f"| {' | '.join(_conflict_row(c, cols_order))} |" for c in unresolved]
            out.append("")
        return out

    lines += [
        "## Cross-annotation conflicts (`conflict`)",
        "",
        "This entity's own columns imply different MetaNetX identities "
        "from each other -- see Method above. **Bold** marks the value(s) "
        "the suggested fix is based on.",
        "",
        "### Metabolites",
        "",
    ]
    lines += _conflict_table("met", ("bigg.metabolite", "chebi", "kegg.compound",
                                      "metanetx.chemical"))
    lines += ["", "### Reactions", ""]
    lines += _conflict_table("rxn", ("bigg.reaction", "kegg.reaction", "metanetx.reaction"))
    lines += [""]

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


def write_metrics(metrics_path: Path, isolated, recurring, conflicts, findings) -> None:
    """Write the ``key<TAB>value`` summary CI reads (see build_qc_report.py's
    _render_annotation()) -- the same shape as qc_metrics.tsv/
    validation_metrics.tsv, kept separate from the prose report so a wording
    change there can never silently break what CI parses."""
    by_status = Counter(f[3] for f in findings)
    rows = {
        "annotation_wrong_isolated": len(isolated),
        "annotation_wrong_recurring": len(recurring),
        "annotation_conflict": len(conflicts),
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
    _, _, mnx2xref, sig2mnxr, mnxr2xref, mnx2info = mnx
    mnxr2sig = {m: sig for sig, ms in sig2mnxr.items() for m in ms}
    met_ids = set(args.mets.split(",")) if args.mets else None
    rxn_ids = set(args.rxns.split(",")) if args.rxns else None

    met_rows = load_met_rows()
    rxn_rows = load_rxn_rows()
    names = {r["id"]: r.get("name", "").strip() for r in met_rows}
    names.update({r["id"]: r.get("name", "").strip() for r in rxn_rows})
    model = cobra.io.load_yaml_model(str(YAML))
    model_charge = {m.id: m.charge for m in model.metabolites}

    findings = []
    conflicts = []
    best_mnx = {}
    if args.all or args.mets:
        rev = reverse_xref(mnx2xref)
        met_findings, met_best = verify_metabolites(met_rows, mnx, met_ids, rev)
        findings += met_findings
        conflicts += verify_metabolite_conflicts(met_rows, mnx, rev, met_best, met_ids)
        best_mnx.update(met_best)
    if args.all or args.rxns:
        smiles_by_id = {r["id"]: r.get("smiles", "").strip() for r in met_rows}
        metkey = {mid: pragmatic(inchikey(smiles)) for mid, smiles in smiles_by_id.items()}
        metkey = {k: v for k, v in metkey.items() if v}
        skip_ids = _water_and_proton_ids(model)
        stoich = load_reaction_stoich(model, skip_ids)
        rxn_findings, rxn_best = verify_reactions(rxn_rows, stoich, mnx, metkey, rxn_ids)
        findings += rxn_findings
        conflicts += verify_reaction_conflicts(rxn_rows, mnxr2sig, mnxr2xref,
                                                reverse_xref(mnxr2xref), rxn_best, rxn_ids)
        best_mnx.update(rxn_best)

    wrong = [f for f in findings if f[3] == "wrong"]
    isolated, recurring, groups = split_wrong_by_recurrence(wrong)

    by_status = Counter(f[3] for f in findings)
    print(f"findings: {dict(by_status)}, conflict: {len(conflicts)} "
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
    if conflicts:
        resolved = sum(1 for c in conflicts if c["basis"])
        print(f"\nconflict ({len(conflicts)}) -- an entity's own columns disagree "
              f"with each other; {resolved} with a suggested fix, "
              f"{len(conflicts) - resolved} with no consensus -- see the report")
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
    write_report(report_path, findings, isolated, recurring, groups, conflicts, scope, context)
    print(f"\nFull report written to {report_path}")
    write_metrics(metrics_path, isolated, recurring, conflicts, findings)
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
