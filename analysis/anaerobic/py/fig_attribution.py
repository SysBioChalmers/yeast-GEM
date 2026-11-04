"""Supplementary figure: effect of leaving out one section of the 9.0.2 -> 9.1.0
curation (results/curation_ablation_fine.tsv) on Sjöberg glycerol production
and on the mean relative error against the 13C flux data. Two panels sharing
the section axis (no dual axis). Writes results/supplementary/FigS_attribution.{png,pdf}."""
from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = Path(__file__).resolve().parents[1] / "results"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
SERIES = "#2a78d6"          # categorical slot 1 (reference palette, light)
BAND = "#ecebe7"


def main() -> None:
    rows = list(csv.DictReader(open(R / "curation_ablation_fine.tsv"), delimiter="\t"))
    full = next(r for r in rows if r["set"] == "all sections")
    data = []
    for r in rows:
        m = re.match(r"all but \d+ \((.*)\)$", r["set"])
        if not m or "glycine export" in m.group(1) or "renamed" in m.group(1):   # without glycine export the GCV-blocked model barely grows
            continue
        data.append((m.group(1), float(r["glycerol"]), 100 * float(r["flux_MRE"])))
    data.sort(key=lambda d: d[1])
    labels = [d[0] for d in data]
    y = range(len(data))

    plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.edgecolor": MUTED,
                         "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK})
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 4.6), sharey=True, gridspec_kw={"wspace": 0.08})
    for ax in (a, b):
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)

    a.axvspan(4.1, 4.9, color=BAND, zorder=0)
    a.text(4.05, len(data) - 0.45, "measured 4.5 ± 0.4", ha="right", va="bottom", color=INK2, fontsize=7)
    a.barh(y, [d[1] for d in data], height=0.55, color=SERIES, zorder=2)
    a.axvline(float(full["glycerol"]), color=INK, linewidth=1, linestyle=(0, (3, 2)), zorder=3)
    a.set_xlabel("glycerol production (mmol gDW$^{-1}$ h$^{-1}$)")
    a.set_yticks(list(y), labels)
    a.set_xlim(0, 5.2)
    for i, d in enumerate(data):
        if abs(d[1] - float(full["glycerol"])) > 0.2:
            if 1 < d[1] < float(full["glycerol"]):
                a.text(d[1] - 0.06, i, f"{d[1]:.2f}", va="center", ha="right", color="white", fontsize=7, zorder=4)
            else:
                a.text(d[1] + 0.06, i, f"{d[1]:.2f}", va="center", color=INK2, fontsize=7)

    b.barh(y, [d[2] for d in data], height=0.55, color=SERIES, zorder=2)
    b.axvline(100 * float(full["flux_MRE"]), color=INK, linewidth=1, linestyle=(0, (3, 2)), zorder=3)
    b.set_xlabel("mean relative error vs. 13C fluxes (%)")
    for i, d in enumerate(data):
        if abs(d[2] - 100 * float(full["flux_MRE"])) > 3:
            b.text(d[2] + 0.8, i, f"{d[2]:.0f}%", va="center", color=INK2, fontsize=7)
    b.set_xlim(0, 85)

    fig.suptitle("Yeast9.0.2 → 9.1.0 curation with one section left out (updated anaerobic constraints)",
                 x=0.02, ha="left", fontsize=9, color=INK)
    fig.text(0.02, 0.005, "Dashed line: all sections (Yeast9.1.0). Not shown: glycine export, without which the "
             "model with glycine cleavage blocked barely grows.", fontsize=6.5, color=MUTED)
    out = R / "supplementary"
    out.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"FigS_attribution.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    print("written", out / "FigS_attribution.png")


if __name__ == "__main__":
    main()
