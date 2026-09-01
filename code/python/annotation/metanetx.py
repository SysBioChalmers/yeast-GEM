"""Shared MetaNetX plumbing for verify_annotations.py and
verify_annotation_consistency.py: downloading/caching MetaNetX's reference
tables, the structure-matching primitives both scripts compare against, and
the small tsv/id-handling helpers neither needs to reimplement.
"""
from __future__ import annotations

import csv
import os
import re
import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

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
    mnx2xref / mnxr2xref, shared by both scripts so it is only built once
    per run."""
    rev = defaultdict(set)
    for m, dbs in xref.items():
        for db, es in dbs.items():
            for e in es:
                rev[(db, e)].add(m)
    return rev


# --------------------------------------------------------------------------- #
# Structural signatures and matching
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


def _best_mnx(mnxset, xref):
    return sorted(mnxset, key=lambda m: (-len(xref.get(m, {})), m))[0]


def match_metabolites_by_structure(rows, mnx):
    """metabolite id -> {"ik", "mnxset", "best"} for every row whose smiles
    resolves to at least one MetaNetX chemical -- the structure-match step
    shared by verify_annotations.py (compares each tsv column against it)
    and verify_annotation_consistency.py ("best" is the strongest tie-break
    when an entity's own columns disagree with each other)."""
    key2mnx, _, mnx2xref, *_ = mnx
    result = {}
    for r in rows:
        ik = inchikey(r.get("smiles", "").strip())
        if not ik:
            continue
        mnxset = key2mnx.get(pragmatic(ik), set())
        if not mnxset:
            continue
        result[r["id"]] = {"ik": ik, "mnxset": mnxset, "best": _best_mnx(mnxset, mnx2xref)}
    return result


def match_reactions_by_structure(stoich, mnx, metkey):
    """reaction id -> {"mnxrset", "best"} for every reaction whose
    participant-structure signature resolves to at least one MetaNetX
    reaction -- the reaction-side equivalent of
    match_metabolites_by_structure, one level up."""
    _, _, _, sig2mnxr, mnxr2xref, _ = mnx
    result = {}
    for rid, mets in stoich.items():
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
        result[rid] = {"mnxrset": mnxrset, "best": _best_mnx(mnxrset, mnxr2xref)}
    return result


# --------------------------------------------------------------------------- #
# tsv / model loading
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


def load_metkey(met_rows):
    """metabolite id -> pragmatic InChIKey, for every row with a valid
    smiles -- the per-metabolite structure key match_reactions_by_structure
    needs to key each participant."""
    smiles_by_id = {r["id"]: r.get("smiles", "").strip() for r in met_rows}
    metkey = {mid: pragmatic(inchikey(smiles)) for mid, smiles in smiles_by_id.items()}
    return {k: v for k, v in metkey.items() if v}


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
