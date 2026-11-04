"""Elemental composition and degree of reduction of the biomass made by a
yeast-GEM model, from the exchange balance: at steady state everything that
ends up in biomass enters or leaves through exchange reactions, so
biomass elements = -sum(v_ex * elements_ex) / growth. Growth fixed at 0.1 h-1
on the model's default (aerobic minimal glucose) medium, pFBA. Also checks
that every reaction carrying flux is element- and charge-balanced.
DR weights: C 4, H 1, O -2, N -3, S 6, P 5, minus charge.
Usage: python biomass_dr_exchange.py <model.xml> <label> [remove_nadph]"""
from __future__ import annotations

import re
import sys
from collections import defaultdict

import cobra
from cobra.flux_analysis import pfba

DRW = {"C": 4, "H": 1, "O": -2, "N": -3, "S": 6, "P": 5}


def parse(f):
    out = defaultdict(float)
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", f or ""):
        out[el] += float(n) if n else 1.0
    return out


def main(xml, label, remove_nadph=False):
    m = cobra.io.read_sbml_model(xml)
    bio = m.reactions.r_4041
    if remove_nadph:
        g = m.metabolites.get_by_id
        bio.add_metabolites({g("s_1212"): 3, g("s_1207"): -3, g("s_0794"): 3})  # undo -3 NADPH, +3 NADP, -3 H+
    m.reactions.r_4041.bounds = (0.1, 0.1)
    m.reactions.r_1714.lower_bound = -10   # glucose uptake up to 10 (default medium allows 1)
    m.objective = "r_1714"           # glucose exchange (negative = uptake); maximise = least uptake
    m.objective_direction = "max"
    sol = pfba(m)
    els = defaultdict(float)
    charge = 0.0
    lines = []
    for r in m.exchanges:
        v = sol.fluxes[r.id]
        if abs(v) < 1e-9 or r.id == "r_2111":
            continue
        (met, c), = r.metabolites.items()
        f = parse(met.formula)
        for el, n in f.items():
            els[el] += -v * c * n / 0.1
        charge += -v * c * (met.charge or 0) / 0.1
        lines.append((r.id, met.name, met.formula, met.charge, v))
    unbalanced = []
    for r in m.reactions:
        if r in m.exchanges or abs(sol.fluxes[r.id]) < 1e-9 or "pseudoreaction" in r.name or r.id == "r_4041" or "SLIME" in r.name:
            continue
        bal = r.check_mass_balance()
        if bal:
            unbalanced.append((r.id, r.name, sol.fluxes[r.id], dict(bal)))
    C = els["C"]
    DR = sum(DRW.get(el, 0) * n for el, n in els.items()) - charge
    print(f"== {label}")
    for l in lines:
        print(f"   exchange {l[0]} {l[1]} ({l[2]}, charge {l[3]}): {l[4]:.4f}")
    print(f"   per gDW: " + " ".join(f"{el} {els[el]:.3f}" for el in ['C','H','O','N','S','P']) + f" mmol; charge {charge:.3f}")
    x = {el: els[el] / C for el in ['H','O','N','S','P']}
    print(f"   per C-mol: CH{x['H']:.4f} O{x['O']:.4f} N{x['N']:.4f} S{x['S']:.4f} P{x['P']:.4f}, charge {charge / C:.4f}")
    print(f"   DR = 4 + {x['H']:.4f} - 2*{x['O']:.4f} - 3*{x['N']:.4f} + 6*{x['S']:.4f} + 5*{x['P']:.4f} - ({charge / C:.4f}) = {DR / C:.4f}")
    print(f"   unbalanced reactions carrying flux: {len(unbalanced)}")
    for u in unbalanced[:10]:
        print("     ", u)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], len(sys.argv) > 3)
