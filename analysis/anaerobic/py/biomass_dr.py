"""Degree of reduction of yeast-GEM biomass, computed directly from the
biomass reaction: the pseudo-reactions (protein, carbohydrate, RNA, DNA,
lipid backbone/chain, cofactor, ion) are expanded recursively until every
input has a chemical formula; the net elemental composition and charge of
what goes into 1 gDW is summed. Degree of reduction per C-mol with C 4, H 1,
O -2, N -3, S 6, P 5, minus the charge (NH3, H2SO4, H3PO4, H2O, CO2 = 0).
GAM (ATP + H2O -> ADP + Pi + H+) is excluded, as it is element-neutral.
Usage: python biomass_dr.py <model.xml> [label]"""
from __future__ import annotations

import re
import sys
from collections import defaultdict

import cobra

DRW = {"C": 4, "H": 1, "O": -2, "N": -3, "S": 6, "P": 5}


def parse(f: str) -> dict[str, float]:
    out: dict[str, float] = defaultdict(float)
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", f or ""):
        out[el] += float(n) if n else 1.0
    return out


def expand(model, rxn, scale, acc, path, log):
    """Add the formula-bearing inputs (negative) and outputs of rxn*scale to acc."""
    for met, coef in rxn.metabolites.items():
        if met.id == "s_0450":          # biomass itself
            continue
        c = coef * scale
        if met.formula:
            acc[met.id] += c
            continue
        if c > 0:                        # a formula-less by-product: report
            log.append(f"formula-less product {met.id} {met.name} {c:.4g}")
            continue
        producers = [r for r in met.reactions if r.metabolites.get(met, 0) > 0 and r.id != rxn.id]
        if len(producers) != 1:
            log.append(f"{met.id} {met.name}: {len(producers)} producers, not expanded")
            continue
        p = producers[0]
        expand(model, p, -c / p.metabolites[met], acc, path + [p.id], log)


def main(xml: str, label: str) -> None:
    model = cobra.io.read_sbml_model(xml)
    bio = model.reactions.get_by_id("r_4041")
    acc: dict[str, float] = defaultdict(float)
    log: list[str] = []
    expand(model, bio, 1.0, acc, ["r_4041"], log)
    # remove GAM: ATP + H2O -> ADP + Pi + H+ (element- and charge-neutral pairs)
    els: dict[str, float] = defaultdict(float)
    charge = 0.0
    contrib = []
    for mid, c in acc.items():
        m = model.metabolites.get_by_id(mid)
        f = parse(m.formula)
        for el, n in f.items():
            els[el] += -c * n           # consumed (c<0) goes into biomass
        ch = m.charge if m.charge is not None else 0
        charge += -c * ch
        dr = sum(DRW.get(el, 0) * n for el, n in f.items()) - ch
        contrib.append((-c * dr, -c * f.get("C", 0), mid, m.name, c))
    C = els["C"]
    DR = sum(DRW.get(el, 0) * n for el, n in els.items()) - charge
    print(f"== {label} ({xml})")
    for l in log:
        print("  note:", l)
    comp = " ".join(f"{el}{els[el] / C:.4f}" for el in ["C", "H", "O", "N", "S", "P"] if els.get(el))
    print(f"  per gDW: C {C:.3f} H {els['H']:.3f} O {els['O']:.3f} N {els['N']:.3f} S {els['S']:.4f} P {els['P']:.4f} mmol; charge {charge:.3f} mmol")
    print(f"  per C-mol: {comp}; charge {charge / C:.4f}")
    print(f"  electrons (DR units) {DR:.2f} mmol/gDW over {C:.2f} C-mmol/gDW -> DR {DR / C:.4f} per C-mol")
    others = {k: v for k, v in els.items() if k not in DRW}
    if others:
        print("  other elements (DR 0):", {k: round(v, 4) for k, v in others.items()})
    contrib.sort(key=lambda t: -abs(t[0]))
    print("  largest contributions (electrons, C-mmol, met, coefficient):")
    for e, c, mid, name, coef in contrib[:12]:
        print(f"    {e:8.2f} {c:8.3f}  {mid} {name} {coef:.4g}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "")
