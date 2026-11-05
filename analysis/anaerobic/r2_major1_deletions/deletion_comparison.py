"""Deletion strains: measured phenotype against the prediction of Yeast9.1.2.

Reads data/knockout_data_extended.tsv (anaerobic physiology of deletion strains
with their same-paper references) and results/ko_panel.tsv (pFBA predictions of
anaerobic batch growth at a glucose uptake of 20 mmol/gDW/h, updated
constraints). A strain is matched when its deleted genes are the genes removed
in the prediction.

* viability: the growth statement of the paper (growth_qualitative) is
  classified as normal, impaired or none, and the prediction from the ratio of
  the predicted growth to that of the wild type (below 0.05 none, below 0.9
  impaired);
* yields: ratio of the mutant to the same-paper wild type for the yields of
  biomass, ethanol and glycerol (the units cancel in the ratio).

Writes results/deletion_comparison.tsv and results/figures/figDeletions.{png,pdf}.
Usage: python deletion_comparison.py [variant]   (default vdev_new)
"""
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
R = Path(os.environ.get("ANAEROBIC_RESULTS", STUDY / "results"))
VARIANT = sys.argv[1] if len(sys.argv) > 1 else "vdev_new"


def read(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


data = read(STUDY / "data" / "knockout_data_extended.tsv")
panel = {r["strain"]: r for r in read(R / "ko_panel.tsv") if r["variant"] == VARIANT}
by_genes = {frozenset(r["genes_in_model"].split()): (name, r) for name, r in panel.items() if r["genes_in_model"] or name == "wild type"}
WT = panel["wild type"]


def predicted(genes):
    return by_genes.get(frozenset(genes.split()))


def classify(text):
    t = text.lower()
    if any(k in t for k in ("rescued", "restored", "relieved", "alleviated", "resumes", "not rescued")):
        return "rescue"
    if any(k in t for k in ("no growth", "unable", "abolished", "completely inhibited")):
        return "none"
    if any(k in t for k in ("not observed", "no effect", "not affected", "similar to wild", "wild-type growth",
                            "wild-type-like", "grows like", "grew well", "viable")):
        return "normal"
    if any(k in t for k in ("reduced", "poorly", "slow", "decreased", "defect", "lag")):
        return "impaired"
    return "unclear"


def predicted_class(match):
    ratio = float(match["mu"]) / float(WT["mu"])
    return "none" if ratio < 0.05 else "impaired" if ratio < 0.9 else "normal"


rows, viability = [], []
ref = defaultdict(dict)   # (source, quantity, unit) -> wild type value
for r in data:
    if r["genotype"].startswith("wild type"):
        try:
            ref[(r["source"], r["quantity"], r["unit"])][r["cultivation"]] = float(r["value"])
        except ValueError:
            pass
for r in data:
    m = predicted(r["deleted_genes_systematic"]) if r["deleted_genes_systematic"] else None
    if m is None or r["genotype"].startswith("wild type"):
        continue
    name, pr = m
    if r["quantity"] == "growth_qualitative":
        measured = classify(r["value"])
        if measured in ("rescue", "unclear"):
            continue
        pc = predicted_class(pr)
        viability.append((r["source"], name, measured, pc))
        rows.append({"source": r["source"], "strain": name, "genotype": r["genotype"], "quantity": "growth",
                     "measured": r["value"], "measured_class": measured, "predicted": f"{float(pr['mu']):.3f}",
                     "predicted_class": pc, "agree": int(measured == pc)})
    elif r["quantity"] in ("Y_biomass", "Y_ethanol", "Y_glycerol"):
        try:
            value = float(r["value"])
            wt = ref[(r["source"], r["quantity"], r["unit"])][r["cultivation"]]
        except (ValueError, KeyError):
            continue
        measured_ratio = value / wt
        predicted_ratio = float(pr[r["quantity"]]) / float(WT[r["quantity"]])
        rows.append({"source": r["source"], "strain": name, "genotype": r["genotype"], "quantity": r["quantity"],
                     "measured": f"{measured_ratio:.3f}", "measured_class": "", "predicted": f"{predicted_ratio:.3f}",
                     "predicted_class": "", "agree": ""})

out = R / "deletion_comparison.tsv"
with open(out, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
n = len(viability)
unique = set(viability)
print(f"{len(rows)} rows; growth class agrees in {sum(v[2] == v[3] for v in unique)} of {len(unique)} different strain, paper and class combinations")

# figure: yields and viability
fig, (axa, axb) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1.3, 1]})
colour = {"Y_biomass": "#2a6fdb", "Y_ethanol": "#eb6834", "Y_glycerol": "#1bae7a"}
for q, c in colour.items():
    pts = [(float(r["measured"]), float(r["predicted"]), r["strain"]) for r in rows if r["quantity"] == q]
    axa.scatter([p[0] for p in pts], [p[1] for p in pts], s=26, color=c, label=q.replace("Y_", "yield of "), zorder=3)
lim = [0, max(1.6, max(float(r["measured"]) for r in rows if r["quantity"] in colour))]
axa.plot(lim, lim, color="#999999", lw=0.8, zorder=1)
axa.set_xlim(lim); axa.set_ylim(0, 1.6)
axa.set_xlabel("measured, mutant / wild type", fontsize=9)
axa.set_ylabel("predicted, mutant / wild type", fontsize=9)
axa.legend(fontsize=8, frameon=False)
axa.set_title("a   Yields", fontsize=10, loc="left")
# viability: one row per strain and source
seen = {}
for s, name, m, p in viability:
    seen[(name, s)] = (m, p)
keys = sorted(seen)
cmap = {"normal": "#1bae7a", "impaired": "#eda100", "none": "#d6403f"}
for i, k in enumerate(keys):
    m, p = seen[k]
    y = len(keys) - i
    axb.scatter(0, y, marker="s", s=90, color=cmap[m])
    axb.scatter(1, y, marker="s", s=90, color=cmap[p])
axb.set_yticks([len(keys) - i for i in range(len(keys))])
axb.set_yticklabels([f"{k[0]} ({k[1]})" for k in keys], fontsize=7)
axb.set_xticks([0, 1]); axb.set_xticklabels(["measured", "predicted"], fontsize=9)
axb.set_xlim(-0.5, 1.5)
axb.set_title("b   Growth", fontsize=10, loc="left")
axb.legend(handles=[plt.Line2D([], [], marker="s", ls="", color=c, label=k) for k, c in cmap.items()], fontsize=8, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.2), ncol=3)
for ax in (axa, axb):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(R / "figures" / f"figDeletions.{ext}", dpi=300)
print("written", R / "figures" / "figDeletions.png")
