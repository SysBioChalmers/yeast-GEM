"""Supplementary tables for the revision, from the study's result files.

Writes results/supplementary/*.tsv. Numbers are for the published model
versions (Yeast9.0.2 / 9.1.0, original / updated anaerobic constraints);
rerun after any model change.
"""
from __future__ import annotations

import csv
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
R = Path(os.environ.get("ANAEROBIC_RESULTS", Path(__file__).resolve().parents[1] / "results"))
OUT = R / "supplementary"
FLUX = REPO / "data" / "physiology" / "flux_data_anaerobic.tsv"
VARIANT = {"vdev_new": "Yeast9.1.2, updated constraints", "vdev_old": "Yeast9.1.2, original constraints plus thiamine uptake",
           "v902_new": "Yeast9.0.2, updated constraints", "v902_old": "Yeast9.0.2, original constraints"}


def read(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(name: str, rows: list[dict], fields: list[str] | None = None) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fields = fields or list(rows[0])
    with open(OUT / name, "w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"{name}: {len(rows)} rows")


def rxn_names() -> dict[str, str]:
    names = {}
    with open(FLUX, newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) > 7 and row[6].startswith("r_"):
                names[row[6]] = row[7]
    return names


def fmt(x: str, d: int = 3) -> str:
    try:
        return f"{float(x):.{d}f}"
    except ValueError:
        return x


# SX: constraints of the original and updated anaerobic functions
SX = [
    ("O2 uptake", "r_1992", "closed", "closed", "anaerobic", "", "function"),
    ("ergosterol, lanosterol, zymosterol, 14-demethyllanosterol uptake", "r_1757; r_1915; r_2106; r_2134", "open", "open",
     "O2 needed for sterol synthesis", "", "function"),
    ("oleate, palmitoleate uptake", "r_2189; r_1994", "open", "open",
     "O2 needed by Ole1; note: yeast can grow anaerobically without (Dekker et al. 2019), the model cannot", "10.1093/femsyr/foz060", "function"),
    ("ergosta-5,7,22,24(28)-tetraen-3beta-ol uptake", "r_2137", "open", "closed",
     "not a medium component; its reduction to ergosterol by Erg4 consumes NADPH", "", "function"),
    ("nicotinate, (R)-pantothenate uptake", "r_1967; r_1548", "closed", "open",
     "O2 needed for their synthesis (Perli et al. 2020; 2021)", "10.1002/yea.3461; 10.1128/mbio.00967-21", "function"),
    ("heme a in cofactor pseudoreaction", "r_4598", "removed", "removed", "heme synthesis needs O2", "", "function"),
    ("NAD, NADH, NADP, NADPH, CoA in cofactor pseudoreaction", "r_4598", "removed", "kept",
     "made anaerobically when nicotinate and pantothenate are supplied", "", "function"),
    ("GAM / NGAM / protein", "r_4041; r_4046", "30.49 / 0 / 0.461 (Nissen et al. 1997)", "55.3 / 0.7 / 0.46 (yeast-GEM)",
     "same maintenance for all conditions", "", "function"),
    ("amino-acid ratio of protein", "r_4047", "aerobic", "anaerobic", "measured (Bjorkeroth et al. 2020)", "", "function"),
    ("FADH2 production in biomass", "r_4041", "none", "0.08 mmol/gDW", "Ero1 disulphide formation, FAD reoxidised by fumarate reductase",
     "10.1038/s41467-018-07285-9", "function"),
    ("Mdh1, mitochondrial malate dehydrogenase", "r_0713", "malate -> OAA only", "reversible", "no evidence found for the restriction", "", "-"),
    ("Mdh2, cytosolic malate dehydrogenase", "r_0714", "malate -> OAA only", "blocked",
     "glucose-induced inactivation (Hung et al. 2004); 260-fold lower anaerobically (Table 1)", "10.1074/jbc.M404544200", "function"),
    ("Idp2, cytosolic NADP-isocitrate dehydrogenase", "r_0659", "reversible", "blocked",
     "glucose-repressed, Cat8-activated (Bojunga & Entian 1999); not detected (Table 1)", "", "function"),
    ("Gcy1", "r_0487", "blocked (written glycerol + NADP -> DHA + NADPH)", "DHA + NADPH -> glycerol, irreversible",
     "thermodynamics (reply to R1 comment 6)", "10.1038/srep12846; 10.1074/jbc.M203030200", "model"),
    ("Glt1, glutamate synthase", "r_0472", "blocked", "blocked",
     "ancillary with ammonium (Valenzuela et al. 1998); Gdh1 activity 50-100x higher (Nissen et al. 2000); enzyme present (Table 1)",
     "10.1128/JB.180.14.3533-3540.1998; 10.1006/mben.1999.0140", "model; condition file nitrogen_limitation"),
    ("glycine cleavage system", "r_0501; r_0507; r_0509", "open", "blocked",
     "nitrogen catabolite repression, induced by glycine; 13C activity reported (Maaheimo et al. 2001; Jouhten et al. 2008)",
     "10.1111/j.1567-1364.2002.tb00069.x", "model; condition files glycine_nitrogen, nitrogen_limitation"),
    ("carnitine O-acetyltransferase (cytosolic)", "r_0252", "open", "blocked",
     "no carnitine in the medium; no synthesis (van Rossum et al. 2016)", "10.1128/mBio.00520-16", "model; condition file carnitine"),
]


def main() -> None:
    names = rxn_names()
    # glycerol when the single constraint is undone: function_ablation, sectionB, curation_ablation
    def sjoberg(folder):
        out = {}
        for r in read(R / folder / "fig1_sjoberg.tsv"):
            out.setdefault(r["variant"], {})[r["quantity"]] = float(r["simulated"])
        return out
    abl, opened = sjoberg("function_ablation"), sjoberg("sectionB")
    cur = {r["set"]: r for r in read(R / "curation_ablation.tsv")}
    carn = next(r for k, r in cur.items() if k.startswith("all but 15 "))

    def show(glycerol, growth, ref=abl["vdev_new"]["growth"]):
        text = f"{glycerol:.2f}"
        return text + (f" (growth {growth:.2f})" if abs(growth - ref) > 0.02 else "")
    undone = {"r_2137": "vdev_new+undoTetraenol", "r_4598": "vdev_new+oldCofactor", "r_4041; r_4046": "vdev_new+oldGAM",
              "r_4047": "vdev_new+undoAA", "r_4041": "vdev_new+fadh2x0", "r_0713": "vdev_new+mdh1old", "r_0714": "vdev_new+undoMdh2",
              "r_0659": "vdev_new+undoIdp2"}
    effect = {r: show(abl[v]["glycerol"], abl[v]["growth"]) for r, v in undone.items()}
    for r, v in {"r_0487": "vdev_new+gcy1rev", "r_0472": "vdev_new+glt1open", "r_0501; r_0507; r_0509": "vdev_new+gcvopen"}.items():
        effect[r] = show(opened[v]["glycerol"], opened[v]["growth"])
    effect["r_0252"] = show(float(carn["glycerol"]), float(carn["growth"]))
    write("SX_anaerobic_constraints.tsv",
          [{"constraint": c, "reactions": r, "original": o, "updated": u, "basis": b, "doi": d, "where": w,
            "glycerol_if_undone": effect.get(r, "")} for c, r, o, u, b, d, w in SX])

    # SY: opening GCV, Glt1, Gcy1 against the 13C data
    pts = read(R / "sectionB" / "fig3cd_points.tsv")
    cols = {"vdev_new": "published", "vdev_new+gcvopen": "GCV_open", "vdev_new+glt1open": "Glt1_open",
            "vdev_new+gcy1rev": "Gcy1_reversible"}
    by = {}
    for p in pts:
        if p["variant"] in cols:
            by.setdefault((p["dataset"], p["rxn"]), {"dataset": p["dataset"], "rxn": p["rxn"],
                                                    "name": names.get(p["rxn"], ""), "measured": fmt(p["measured"], 2)})
            by[(p["dataset"], p["rxn"])][cols[p["variant"]]] = fmt(p["simulated"], 2)
    write("SY_opening_constraints_13C.tsv", list(by.values()),
          ["dataset", "rxn", "name", "measured", "published", "GCV_open", "Glt1_open", "Gcy1_reversible"])

    # SZ1: 2x2
    sj = read(R / "repro" / "fig1_sjoberg.tsv")
    gr = {r["variant"]: r for r in read(R / "repro" / "fig3ab_summary.tsv")}
    fx = {r["variant"]: r for r in read(R / "repro" / "fig3cd_summary.tsv")}
    rows = []
    for v in ("v902_old", "v902_new", "vdev_old", "vdev_new"):
        s = {r["quantity"]: r["simulated"] for r in sj if r["variant"] == v}
        rows.append({"model": VARIANT[v], "growth": fmt(s["growth"]), "glycerol": fmt(s["glycerol"]),
                     "ethanol": fmt(s["ethanol"], 2), "CO2": fmt(s["CO2"], 2),
                     "growth_R2": fmt(gr[v]["R2_identity"]), "growth_MRE": fmt(gr[v]["MRE"]),
                     "MRE_C_limited_anaerobic": fmt(gr[v]["MRE_C_limited_anaerobic"]),
                     "MRE_N_limited_anaerobic": fmt(gr[v]["MRE_N_limited_anaerobic"]),
                     "flux_R2": fmt(fx[v]["R2_identity"]), "flux_MRE": fmt(fx[v]["MRE"])})
    write("SZ1_model_x_constraints.tsv", rows)

    # SZ2: anaerobic constraints undone one at a time
    label = {"vdev_new": "none (updated constraints)", "vdev_new+undoTetraenol": "tetraenol uptake allowed",
             "vdev_new+undoMdh2": "Mdh2 unblocked", "vdev_new+undoIdp2": "Idp2 unblocked",
             "vdev_new+fadh2x0": "no FADH2 in biomass", "vdev_new+undoAA": "aerobic amino-acid ratio",
             "vdev_new+oldGAM": "GAM 30.49, NGAM 0, protein 0.461 (original)",
             "vdev_new+oldCofactor": "NAD(P)(H), CoA removed from cofactors (original)",
             "vdev_new+mdh1old": "Mdh1 malate -> OAA only (original)"}
    sj = read(R / "function_ablation" / "fig1_sjoberg.tsv")
    gr = {r["variant"]: r for r in read(R / "function_ablation" / "fig3ab_summary.tsv")}
    fx = {r["variant"]: r for r in read(R / "function_ablation" / "fig3cd_summary.tsv")}
    rows = []
    for v, lab in label.items():
        s = {r["quantity"]: r["simulated"] for r in sj if r["variant"] == v}
        rows.append({"change on Yeast9.1.2 with updated constraints": lab, "growth": fmt(s["growth"]),
                     "glycerol": fmt(s["glycerol"]), "ethanol": fmt(s["ethanol"], 2),
                     "growth_R2": fmt(gr[v]["R2_identity"]), "MRE_C_limited_anaerobic": fmt(gr[v]["MRE_C_limited_anaerobic"]),
                     "flux_MRE": fmt(fx[v]["MRE"])})
    write("SZ2_constraints_undone.tsv", rows)

    # SZ3: curation of Yeast9.1.2 with one part left out or added; SZ4: two parts left out
    METRICS = ["growth", "glycerol", "ethanol", "CO2", "growth_R2", "growth_MRE", "MRE_Clim_anaerobic",
               "MRE_Nlim_anaerobic", "flux_R2", "flux_MRE"]
    cur = read(R / "curation_ablation.tsv")
    ref = next(r for r in cur if r["set"] == "all sections")
    delta = {"glycerol": 3, "flux_MRE": 3, "growth_R2": 3}
    rows = []
    for r in cur:
        row = {"curation of Yeast9.1.2, one part left out or added (replay of v9_1_0.m on Yeast9.0.2 plus the later changes; updated constraints)": r["set"]}
        row.update({m: fmt(r[m]) for m in METRICS})
        row.update({f"change_in_{m}": fmt(float(r[m]) - float(ref[m])) for m in delta})
        rows.append(row)
    write("SZ3_curation_parts.tsv", rows)

    single = {}
    for r in cur:
        if r["set"].startswith("all but ") and " and " not in r["set"]:
            single[r["set"].split(" ", 3)[2] if r["set"].split(" ", 3)[2].isdigit() else r["set"]] = r
    pairs = read(R / "curation_pairs.tsv")
    rows = []
    for r in pairs:
        row = {"parts left out": r["set"]}
        row.update({m: fmt(r[m]) for m in METRICS})
        d = {m: float(r[m]) - float(ref[m]) for m in delta}
        row.update({f"change_in_{m}": fmt(d[m]) for m in delta})
        label = r["set"]
        if label.startswith("all but ") and not label.startswith("all but categories"):
            i, j = label.split()[2], label.split()[4]
            si, sj = single.get(i), single.get(j)
            if si and sj:
                row.update({f"interaction_{m}": fmt(d[m] - (float(si[m]) - float(ref[m])) - (float(sj[m]) - float(ref[m]))) for m in delta})
        rows.append(row)
    write("SZ4_curation_pairs.tsv", rows)

    # SZ5: two anaerobic constraints undone
    sj = read(R / "function_pairs" / "fig1_sjoberg.tsv")
    gr = {r["variant"]: r for r in read(R / "function_pairs" / "fig3ab_summary.tsv")}
    fx = {r["variant"]: r for r in read(R / "function_pairs" / "fig3cd_summary.tsv")}
    rows = []
    for v in gr:
        s = {r["quantity"]: r["simulated"] for r in sj if r["variant"] == v}
        rows.append({"constraints undone on Yeast9.1.2": v.replace("vdev_new+", "").replace("+", " and ") or "none",
                     "growth": fmt(s["growth"]), "glycerol": fmt(s["glycerol"]), "ethanol": fmt(s["ethanol"], 2),
                     "CO2": fmt(s["CO2"], 2), "growth_R2": fmt(gr[v]["R2_identity"]), "growth_MRE": fmt(gr[v]["MRE"]),
                     "MRE_C_limited_anaerobic": fmt(gr[v]["MRE_C_limited_anaerobic"]),
                     "MRE_N_limited_anaerobic": fmt(gr[v]["MRE_N_limited_anaerobic"]),
                     "flux_R2": fmt(fx[v]["R2_identity"]), "flux_MRE": fmt(fx[v]["MRE"])})
    write("SZ5_constraints_pairs.tsv", rows)

    # S10: every compared 13C flux, pFBA and sampled interval
    rows = []
    for p in read(R / "analysis_opt0.99999" / "fluxomics_points.tsv"):
        rows.append({"model": VARIANT.get(p["variant"], p["variant"]), "dataset": p["dataset"], "rxn": p["rxn"],
                     "name": names.get(p["rxn"], ""), "measured": fmt(p["measured"], 2), "pFBA": fmt(p["pfba"], 2),
                     "sampled_2.5%": fmt(p["q2_5"], 2), "sampled_median": fmt(p["median"], 2),
                     "sampled_97.5%": fmt(p["q97_5"], 2), "glucose_point": p["is_glucose"]})
    write("S10_fluxomics_points.tsv", rows)


if __name__ == "__main__":
    main()
