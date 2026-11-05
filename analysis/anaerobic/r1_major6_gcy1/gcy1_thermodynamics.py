"""Gcy1 direction: reaction energy of glycerol + NADP+ = dihydroxyacetone + NADPH.

dG'm (all reactants at 1 mM) from eQuilibrator at the cytosolic pH, ionic
strength and pMg of the model checks, and the largest dihydroxyacetone to
glycerol ratio at which the oxidative direction can run for a measured
NADPH/NADP+ ratio. Writes results/gcy1_thermodynamics.tsv.

Usage: python gcy1_thermodynamics.py [results_dir]
Needs equilibrator-api.
"""
import math
import os
import sys
from pathlib import Path

from equilibrator_api import ComponentContribution, Q_

RT = 8.314e-3 * 298.15  # kJ/mol
NADPH_NADP = {"Zhang 2015, low": 16, "Zhang 2015, high": 22}

cc = ComponentContribution()
cc.p_h = Q_(7.2)
cc.ionic_strength = Q_("0.25 M")
cc.p_mg = Q_(3.0)
reaction = cc.parse_reaction_formula("kegg:C00116 + kegg:C00006 = kegg:C00184 + kegg:C00005")
estimate = cc.physiological_dg_prime(reaction)
dg = estimate.value.m_as("kJ/mol")
sd = estimate.error.m_as("kJ/mol")

out = Path(sys.argv[1] if len(sys.argv) > 1 else os.environ.get("ANAEROBIC_RESULTS", Path(__file__).resolve().parents[1] / "results"))
out.mkdir(parents=True, exist_ok=True)
rows = ["case\tdG_m_kJ_mol\tsd_kJ_mol\tNADPH_NADP\tmax_DHA_per_glycerol\tinverse"]
for name, ratio in NADPH_NADP.items():
    for label, g in (("dG-sd", dg - sd), ("dG", dg), ("dG+sd", dg + sd)):
        bound = math.exp(-g / RT) / ratio
        rows.append(f"{name} {label}\t{g:.1f}\t{sd:.1f}\t{ratio}\t{bound:.2e}\t{1 / bound:.0f}")
(out / "gcy1_thermodynamics.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
print("\n".join(rows))
