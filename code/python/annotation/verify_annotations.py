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

The MetaNetX reference tables are downloaded once to a cache directory
($MNX_CACHE, default data/databases/.metanetx) and reused. Not part of
model_qc.py / the pull-request gate -- this needs a first-run download of
MetaNetX's full reference tables (tens of MB) and is meant to be run by a
curator on demand, the same way Human-GEM's equivalent tool is.

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
        "chem_prop.tsv": ("mnx_struct.tsv",
                          lambda f: [f[0], f[7]] if len(f) > 7 and f[7] else None),
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


def load_metanetx():
    key2mnx, mnx2key = defaultdict(set), {}
    for line in (CACHE / "mnx_struct.tsv").open(encoding="utf-8"):
        mnxm, _, ik = line.rstrip("\n").partition("\t")
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
    return key2mnx, mnx2key, mnx2xref, sig2mnxr, mnxr2xref


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
def verify_metabolites(rows, mnx, ids):
    key2mnx, mnx2key, mnx2xref, *_ = mnx
    rev = defaultdict(set)
    for m, dbs in mnx2xref.items():
        for db, es in dbs.items():
            for e in es:
                rev[(db, e)].add(m)
    findings = []
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
                old_skels = {mnx2key.get(m, "").split("-")[0]
                             for db in dbs for m in rev.get((db, hid.replace("CHEBI:", "")), ())}
                if old_skels and pragmatic(ik).split("-")[0] not in old_skels:
                    findings.append(("met", r["id"], col, "wrong", hid,
                                      _best([_norm_met(e, col) for e in pool], col)))
    return findings


def verify_reactions(rxn_rows, stoich, mnx, metkey, ids):
    _, _, _, sig2mnxr, mnxr2xref = mnx
    by_id = {r["id"]: r for r in rxn_rows}
    findings = []
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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true", help="check the whole model")
    ap.add_argument("--mets", help="comma-separated metabolite ids to check")
    ap.add_argument("--rxns", help="comma-separated reaction ids to check")
    ap.add_argument("--fix", action="store_true",
                     help="apply safe corrections (add missing, update drift)")
    args = ap.parse_args()
    if not (args.all or args.mets or args.rxns):
        ap.error("choose --all, --mets and/or --rxns")

    ensure_metanetx()
    mnx = load_metanetx()
    met_ids = set(args.mets.split(",")) if args.mets else None
    rxn_ids = set(args.rxns.split(",")) if args.rxns else None

    findings = []
    if args.all or args.mets:
        findings += verify_metabolites(load_met_rows(), mnx, met_ids)
    if args.all or args.rxns:
        model = cobra.io.load_yaml_model(str(YAML))
        met_rows = load_met_rows()
        smiles_by_id = {r["id"]: r.get("smiles", "").strip() for r in met_rows}
        metkey = {mid: pragmatic(inchikey(smiles)) for mid, smiles in smiles_by_id.items()}
        metkey = {k: v for k, v in metkey.items() if v}
        skip_ids = _water_and_proton_ids(model)
        stoich = load_reaction_stoich(model, skip_ids)
        findings += verify_reactions(load_rxn_rows(), stoich, mnx, metkey, rxn_ids)

    by_status = Counter(f[3] for f in findings)
    print(f"findings: {dict(by_status)}")
    for status in ("wrong", "drift", "missing"):
        rows = [f for f in findings if f[3] == status]
        if rows:
            print(f"\n{status} ({len(rows)}):")
            for kind, ent, col, _s, old, new in rows[:40]:
                print(f"  {kind:3} {ent:12} {col:20} {old or '-'} -> {new}")
            if len(rows) > 40:
                print(f"  ... and {len(rows) - 40} more")

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
