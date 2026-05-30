"""Relative fermentation-product bar plot for the anaerobic model.

Port of ``code/modelTests/plotAnaerobic.m``. The model must already be
in anaerobic state. Glucose uptake is fixed to 23 mmol/gDW/h; FBA is
run, and predicted glycerol / ethanol / CO2 / biomass fluxes are
compared against experimental measurements (4.5 ± 0.4 mmol gly,
31 ± 2 mmol eth, 38 ± 10 mmol CO2, 0.36 ± 0.02 1/h biomass).
"""
from __future__ import annotations

import cobra
import numpy as np

# Reaction ids the plot drives. Comments mirror the labels in plotAnaerobic.m.
_GLC_EX_ID = "r_1714"
_ETHANOL_EX_ID = "r_1761"
_CO2_EX_ID = "r_1672"
_GLYCEROL_EX_ID = "r_1808"
_BIOMASS_RXN_ID = "r_4041"

# Experimental measurements: gly, eth, CO2 (mmol/gDW/h) and biomass (1/h).
_DATA = np.array([4.5, 31.0, 38.0, 0.36])
_ERROR = np.array([0.4, 2.0, 10.0, 0.02])
_LABELS = ("Glycerol", "Ethanol", "CO2", "Biomass")


def plot_anaerobic(
    model_anaerobic: cobra.Model,
    *,
    glucose_uptake: float = 23.0,
    write_output: bool = False,
):
    """Plot relative fermentation-product fluxes and return the
    predicted vector ``(gly, eth, CO2, biomass)``.

    Parameters
    ----------
    model_anaerobic
        Model with the anaerobic condition already applied.
    glucose_uptake
        Glucose uptake rate (mmol/gDW/h). Default 23, matching the
        legacy MATLAB script.
    write_output
        When True, save the figure + a markdown summary under
        ``data/testResults/``.
    """
    import matplotlib

    if write_output:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with model_anaerobic:
        model_anaerobic.reactions.get_by_id(_GLC_EX_ID).bounds = (
            -glucose_uptake, -glucose_uptake,
        )
        sol = model_anaerobic.optimize()
        if sol.status != "optimal":
            raise RuntimeError(
                f"anaerobic FBA returned status {sol.status!r}; cannot plot."
            )
        sim = np.array([
            sol.fluxes[_GLYCEROL_EX_ID],
            sol.fluxes[_ETHANOL_EX_ID],
            sol.fluxes[_CO2_EX_ID],
            sol.fluxes[_BIOMASS_RXN_ID],
        ])

    x = np.arange(len(_LABELS))
    fig, ax = plt.subplots()
    ax.bar(x, _DATA / _DATA, alpha=0.5, label="data")
    ax.bar(x, sim / _DATA, alpha=0.5, label="simulation")
    ax.errorbar(x, _DATA / _DATA, yerr=_ERROR / _DATA,
                fmt="none", color="black")
    ax.set_xticks(x)
    ax.set_xticklabels(_LABELS)
    ax.set_ylabel("Relative value (data / data)")
    ax.legend()

    if write_output:
        from yeastgem.io import REPO_PATH
        out_dir = REPO_PATH / "data" / "testResults"
        out_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_dir / "anaerobic_products.png", bbox_inches="tight")
        plt.close(fig)
    return sim
