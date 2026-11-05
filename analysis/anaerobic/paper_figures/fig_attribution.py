"""Effect of each curation and constraint change on the predictions of Yeast9.1.2.

figAttribution: one row for each change that is left out of the finished model,
grouped by the categories of reviewer 1. Panels: change in glycerol production, in the mean relative
error of the 13C fluxes, and in the R2 of the chemostat growth rates.
figCategories: for each category the model without it and the model with only it
(Yeast9.0.2 plus the category). Writes results/figures/figAttribution.{png,pdf}
and figCategories.{png,pdf}.
"""
import csv
import os
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = Path(os.environ.get("ANAEROBIC_RESULTS", Path(__file__).resolve().parents[1] / "results"))
OUT = R / "figures"
OUT.mkdir(parents=True, exist_ok=True)


def read(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


cur = {r["set"]: r for r in read(R / "curation_ablation.tsv")}
ref = cur["all sections"]
CAT = {  # section number -> category (reviewer 1, major 3)
    **{k: "charge and element balance" for k in (3, 12)},
    **{k: "biomass equation" for k in (16, 17, 18, 19)},
    **{k: "mitochondrial and nitrogen metabolism" for k in (1, 5, 13, 14, 15, 20)},
    **{k: "reaction directions and cofactors" for k in (2, 4, 6, 7, 8, 9, 10, 11)},
}
COLOUR = {"charge and element balance": "#2a6fdb", "biomass equation": "#eb6834",
          "mitochondrial and nitrogen metabolism": "#1bae7a", "reaction directions and cofactors": "#eda100",
          "changes of 9.1.2 after 9.1.0": "#8a8a8a", "anaerobic constraints": "#e87bb0"}
METRICS = [("glycerol", "glycerol (mmol gDW$^{-1}$ h$^{-1}$)"), ("flux_MRE", "13C flux error"),
           ("growth_R2", "growth R$^2$")]

rows = []   # (category, label, values, added-alone values or None)
sections = []
for name, r in cur.items():
    m = re.match(r"all but (\d+) \((.*)\)$", name)
    if m:
        sections.append((CAT[int(m.group(1))], m.group(2), r))
for cat in list(COLOUR)[:4]:
    for c, lab, r in sections:
        if c == cat:
            rows.append((c, lab, r, None))
    rows.append((cat, "ALL " + cat.upper(), cur["all but category: " + cat], cur["only category: " + cat]))
GROUP = {"thiamine": "thiamine (HMP-P synthase, #423)", "fumarate": "fumarate reductases (#423)",
         "tyr1": "TYR1 cofactor (#426)", "acetate_rnr": "acetate transport, ribonucleotide reductase (#389)"}
for g in GROUP:
    rows.append(("changes of 9.1.2 after 9.1.0", GROUP[g], cur["all but 9.1.2 change: " + g], None))

fa = {}
for sets in (read(R / "function_ablation" / "fig1_sjoberg.tsv"),):
    for r in sets:
        fa.setdefault(r["variant"], {})[r["quantity"]] = float(r["simulated"])
gr = {r["variant"]: r for r in read(R / "function_ablation" / "fig3ab_summary.tsv")}
fx = {r["variant"]: r for r in read(R / "function_ablation" / "fig3cd_summary.tsv")}
LAB = {"undoTetraenol": "tetraenol uptake allowed", "undoMdh2": "Mdh2 unblocked", "undoIdp2": "Idp2 unblocked",
       "undoAA": "aerobic amino-acid ratio", "oldGAM": "original GAM, NGAM, protein", "oldCofactor": "original cofactors",
       "mdh1old": "Mdh1 as originally", "fadh2x0": "no FADH2 in biomass"}
base = "vdev_new"
const_ref = {"glycerol": fa[base]["glycerol"], "flux_MRE": float(fx[base]["MRE"]), "growth_R2": float(gr[base]["R2_identity"])}
const_rows = []
for key, label in LAB.items():
    v = f"{base}+{key}"
    const_rows.append(("anaerobic constraints", label,
                       {"glycerol": fa[v]["glycerol"], "flux_MRE": float(fx[v]["MRE"]), "growth_R2": float(gr[v]["R2_identity"])}))

fig, axes = plt.subplots(1, 3, figsize=(11, 9), sharey=True)
fig.subplots_adjust(left=0.27, right=0.98, top=0.98, bottom=0.14, wspace=0.08)
labels = [r[1] for r in rows] + [r[1] for r in const_rows]
y = {i: len(labels) - i for i in range(len(labels))}
for ax, (metric, title) in zip(axes, METRICS):
    ax.axvline(0, color="#999999", lw=0.8)
    for i, (cat, lab, v, alone) in enumerate(rows):
        d = float(v[metric]) - float(ref[metric])
        ax.plot(d, y[i], "o", color=COLOUR[cat], ms=6)
    for j, (cat, lab, v) in enumerate(const_rows):
        ax.plot(float(v[metric]) - const_ref[metric], y[len(rows) + j], "o", color=COLOUR[cat], ms=6)
    if metric == "growth_R2":
        ax.set_xlim(-0.3, 0.03)
        ax.annotate("glycine export: -0.84", xy=(-0.3, y[next(i for i, r in enumerate(rows) if "glycine" in r[1])]),
                    xytext=(-0.29, y[next(i for i, r in enumerate(rows) if "glycine" in r[1])] - 1.3), fontsize=7,
                    arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.set_xlabel("change when left out: " + title, fontsize=8)
    ax.tick_params(labelsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
axes[0].set_yticks([y[i] for i in range(len(labels))])
axes[0].set_yticklabels(labels, fontsize=7)
for tick, r in zip(axes[0].get_yticklabels(), rows + [(None, "", None)] * len(const_rows)):
    if r[1].startswith("ALL "):
        tick.set_fontweight("bold")
handles = [plt.Line2D([], [], marker="o", ls="", color=c, label=k) for k, c in COLOUR.items()]
fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=7.5, frameon=False, bbox_to_anchor=(0.5, 0.0))
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"figAttribution.{ext}", dpi=300)
print("written", OUT / "figAttribution.png")


# Categories: model without the category and model with only the category
import numpy as np
refs = [("Yeast9.1.2 (all changes)", cur["all sections"]), ("Yeast9.0.2 (no curation)", cur["9.0.2 (no section)"])]
cats = list(COLOUR)[:4]
table = [(c, cur["all but category: " + c], cur["only category: " + c]) for c in cats]
fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
fig.subplots_adjust(left=0.30, right=0.99, top=0.84, bottom=0.12, wspace=0.08)
names = [r[0] for r in refs] + [f"{c}" for c in cats]
for ax, (metric, title) in zip(axes, METRICS):
    vals = [float(r[1][metric]) for r in refs]
    ax.axvline(vals[0], color="#999999", lw=0.8, ls=":")
    ax.axvline(vals[1], color="#999999", lw=0.8, ls=":")
    for i, (c, without, only) in enumerate(table):
        yy = len(names) - 3 - i
        ax.plot([float(without[metric]), float(only[metric])], [yy, yy], "-", color="#cccccc", lw=2, zorder=1)
        ax.plot(float(without[metric]), yy, "s", color=COLOUR[c], ms=7, zorder=3)
        ax.plot(float(only[metric]), yy, "o", color=COLOUR[c], ms=7, zorder=3)
    for i, (n, r) in enumerate(refs):
        ax.plot(float(r[metric]), len(names) - 1 - i, "D", color="#333333", ms=6)
    ax.set_title(title, fontsize=8)
    ax.tick_params(labelsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
axes[0].set_yticks([len(names) - 1 - i for i in range(len(names))])
axes[0].set_yticklabels(names, fontsize=7)
fig.legend(handles=[plt.Line2D([], [], marker="s", ls="", color="#555555", label="Yeast9.1.2 without the category"),
                    plt.Line2D([], [], marker="o", ls="", color="#555555", label="Yeast9.0.2 with only the category"),
                    plt.Line2D([], [], marker="D", ls="", color="#333333", label="reference models")],
           loc="upper center", ncol=3, fontsize=7.5, frameon=False, bbox_to_anchor=(0.6, 1.0))
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"figCategories.{ext}", dpi=300)
print("written", OUT / "figCategories.png")
