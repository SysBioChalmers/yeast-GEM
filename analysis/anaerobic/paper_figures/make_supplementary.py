"""Supplementary tables for the revision, from the study's result files.

Writes results/supplementary/*.tsv. Numbers are for the published model
versions (Yeast9.0.2 / 9.1.0, original / updated anaerobic constraints);
rerun after any model change.
"""
from __future__ import annotations

import csv
from pathlib import Path

R = Path(__file__).resolve().parents[1] / "results"
OUT = R / "supplementary"
FLUX = R.parents[1] / "data" / "physiology" / "flux_data_anaerobic.tsv"
VARIANT = {"v910_new": "Yeast9.1.0, updated constraints", "v910_old": "Yeast9.1.0, original constraints",
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
    # glycerol when the single constraint is undone (function_ablation, curation_ablation_fine)
    effect = {"r_2137": "4.35", "r_4598": "4.35", "r_4041; r_4046": "5.38 (growth 0.47)", "r_4047": "4.36",
              "r_0713": "4.48 if restricted", "r_0714": "4.35", "r_0659": "3.99", "r_0487": "0.00",
              "r_0472": "2.30", "r_0501; r_0507; r_0509": "3.41", "r_0252": "3.81"}
    write("SX_anaerobic_constraints.tsv",
          [{"constraint": c, "reactions": r, "original": o, "updated": u, "basis": b, "doi": d, "where": w,
            "glycerol_if_undone": effect.get(r, "")} for c, r, o, u, b, d, w in SX])

    # SY: opening GCV, Glt1, Gcy1 against the 13C data
    pts = read(R / "sectionB" / "fig3cd_points.tsv")
    cols = {"v910_new": "published", "v910_new+gcvopen": "GCV_open", "v910_new+glt1open": "Glt1_open",
            "v910_new+gcy1rev": "Gcy1_reversible"}
    by = {}
    for p in pts:
        if p["variant"] in cols:
            by.setdefault((p["dataset"], p["rxn"]), {"dataset": p["dataset"], "rxn": p["rxn"],
                                                    "name": names.get(p["rxn"], ""), "measured": fmt(p["measured"], 2)})
            by[(p["dataset"], p["rxn"])][cols[p["variant"]]] = fmt(p["simulated"], 2)
    write("SY_opening_constraints_13C.tsv", list(by.values()),
          ["dataset", "rxn", "name", "measured", "published", "GCV_open", "Glt1_open", "Gcy1_reversible"])

    # SZ1: 2x2
    sj = read(R / "repro_3.0.0-rc1" / "fig1_sjoberg.tsv")
    gr = {r["variant"]: r for r in read(R / "repro_3.0.0-rc1" / "fig3ab_summary.tsv")}
    fx = {r["variant"]: r for r in read(R / "repro_3.0.0-rc1" / "fig3cd_summary.tsv")}
    rows = []
    for v in ("v902_old", "v902_new", "v910_old", "v910_new"):
        s = {r["quantity"]: r["simulated"] for r in sj if r["variant"] == v}
        rows.append({"model": VARIANT[v], "growth": fmt(s["growth"]), "glycerol": fmt(s["glycerol"]),
                     "ethanol": fmt(s["ethanol"], 2), "CO2": fmt(s["CO2"], 2),
                     "growth_R2": fmt(gr[v]["R2_identity"]), "growth_MRE": fmt(gr[v]["MRE"]),
                     "MRE_C_limited_anaerobic": fmt(gr[v]["MRE_C_limited_anaerobic"]),
                     "MRE_N_limited_anaerobic": fmt(gr[v]["MRE_N_limited_anaerobic"]),
                     "flux_R2": fmt(fx[v]["R2_identity"]), "flux_MRE": fmt(fx[v]["MRE"])})
    write("SZ1_model_x_constraints.tsv", rows)

    # SZ2: anaerobic constraints undone one at a time
    label = {"v910_new": "none (updated constraints)", "v910_new+undoTetraenol": "tetraenol uptake allowed",
             "v910_new+undoMdh2": "Mdh2 unblocked", "v910_new+undoIdp2": "Idp2 unblocked",
             "v910_new+fadh2x0": "no FADH2 in biomass", "v910_new+undoAA": "aerobic amino-acid ratio",
             "v910_new+oldGAM": "GAM 30.49, NGAM 0, protein 0.461 (original)",
             "v910_new+oldCofactor": "NAD(P)(H), CoA removed from cofactors (original)",
             "v910_new+mdh1old": "Mdh1 malate -> OAA only (original)"}
    sj = read(R / "function_ablation" / "fig1_sjoberg.tsv")
    gr = {r["variant"]: r for r in read(R / "function_ablation" / "fig3ab_summary.tsv")}
    fx = {r["variant"]: r for r in read(R / "function_ablation" / "fig3cd_summary.tsv")}
    rows = []
    for v, lab in label.items():
        s = {r["quantity"]: r["simulated"] for r in sj if r["variant"] == v}
        rows.append({"change on Yeast9.1.0 with updated constraints": lab, "growth": fmt(s["growth"]),
                     "glycerol": fmt(s["glycerol"]), "ethanol": fmt(s["ethanol"], 2),
                     "growth_R2": fmt(gr[v]["R2_identity"]), "MRE_C_limited_anaerobic": fmt(gr[v]["MRE_C_limited_anaerobic"]),
                     "flux_MRE": fmt(fx[v]["MRE"])})
    write("SZ2_constraints_undone.tsv", rows)

    # SZ3: curation replayed, one section left out
    rows = []
    for r in read(R / "curation_ablation_fine.tsv"):
        rows.append({"curation (v9_1_0.m) on Yeast9.0.2, updated constraints": r["set"],
                     "growth": fmt(r["growth"]), "glycerol": fmt(r["glycerol"]), "ethanol": fmt(r["ethanol"], 2),
                     "growth_R2": fmt(r["growth_R2"]), "flux_MRE": fmt(r["flux_MRE"])})
    write("SZ3_curation_sections.tsv", rows)

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
