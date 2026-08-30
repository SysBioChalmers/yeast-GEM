"""Exchange-rate comparison for the anaerobic model.

Port of ``code/modelTests/plotAnaerobic.m``. The model must already be in
the anaerobic state. The measurements come from Sjoberg et al. (2024),
doi:10.1016/j.ymben.2024.01.007, and are read from
``data/physiology/exchange_data_anaerobic.tsv`` together with the glucose
uptake rate they were measured at, so both languages compare against the
same numbers.

No single R2 is reported. The product rates are in mmol/gDW/h and growth
is in 1/h, so pooling them into one coefficient of determination would not
mean anything; the relative deviation per measurement is comparable across
units.
"""
from __future__ import annotations

from dataclasses import dataclass

import cobra
import pandas as pd

from yeastgem.io import REPO_PATH

_GLC_EX_ID = "r_1714"
# Ammonium uptake has to run against the proton gradient the plasma
# membrane ATPase maintains, so the two are coupled. The measured ratio is
# close to 1.
_AMMONIUM_EX_ID = "r_1115"
_ATPASE_ID = "r_0227"

_MEASUREMENTS = REPO_PATH / "data" / "physiology" / "exchange_data_anaerobic.tsv"


@dataclass(frozen=True)
class AnaerobicExchangeResult:
    """Predicted against measured anaerobic exchange rates.

    ``results`` holds one row per measurement, with the measured value,
    its error, the prediction, the relative deviation and whether the
    prediction falls within the experimental error.
    """

    mean_relative_error: float
    max_relative_error: float
    fraction_within_error: float
    n_measurements: int
    ammonium_exchange: float
    atpase: float
    ammonium_per_atpase: float
    results: pd.DataFrame


def plot_anaerobic(
    model_anaerobic: cobra.Model,
    *,
    plot: bool = True,
    write_output: bool = False,
) -> AnaerobicExchangeResult:
    """Compare predicted anaerobic exchange rates against measurements.

    Parameters
    ----------
    model_anaerobic
        Model with the anaerobic condition already applied.
    plot
        Draw the bar plot. Set False to use this as a metric alone.
    write_output
        When True, save the figure under ``data/testResults/``. Implies
        ``plot``.
    """
    measurements = pd.read_csv(_MEASUREMENTS, sep="\t")

    # All rows share one glucose uptake rate; check rather than silently
    # using the first of several different values.
    uptake = measurements["glucoseUptake"].unique()
    if len(uptake) != 1:
        raise ValueError(
            f"{_MEASUREMENTS} gives more than one glucose uptake rate: {uptake}"
        )
    glucose_uptake = float(uptake[0])

    with model_anaerobic:
        model_anaerobic.reactions.get_by_id(_GLC_EX_ID).bounds = (
            -glucose_uptake, -glucose_uptake,
        )
        sol = model_anaerobic.optimize()
        if sol.status != "optimal":
            raise RuntimeError(
                f"anaerobic FBA returned status {sol.status!r} at a glucose "
                f"uptake rate of {glucose_uptake}."
            )
        predicted = [abs(sol.fluxes[rxn]) for rxn in measurements["rxnID"]]
        ammonium = float(sol.fluxes[_AMMONIUM_EX_ID])
        atpase = float(sol.fluxes[_ATPASE_ID])

    results = measurements.copy()
    results["predicted"] = predicted
    results["relativeError"] = (
        (results["predicted"] - results["measured"]) / results["measured"]
    )
    results["withinError"] = (
        (results["predicted"] - results["measured"]).abs() <= results["stdev"]
    )

    if plot or write_output:
        _plot(results, write_output=write_output)

    return AnaerobicExchangeResult(
        mean_relative_error=float(results["relativeError"].abs().mean()),
        max_relative_error=float(results["relativeError"].abs().max()),
        fraction_within_error=float(results["withinError"].mean()),
        n_measurements=len(results),
        ammonium_exchange=ammonium,
        atpase=atpase,
        ammonium_per_atpase=ammonium / atpase,
        results=results,
    )


def _plot(results: pd.DataFrame, *, write_output: bool) -> None:
    """Bar plot of each prediction relative to its measured value.

    Takes the computed table, so the plot shows exactly the numbers that
    were reported rather than solving again.
    """
    import matplotlib

    if write_output:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    measured = results["measured"].to_numpy()
    x = np.arange(len(results))
    fig, ax = plt.subplots()
    ax.bar(x, measured / measured, alpha=0.5, label="data")
    ax.bar(x, results["predicted"].to_numpy() / measured, alpha=0.5,
           label="simulation")
    ax.errorbar(x, measured / measured,
                yerr=results["stdev"].to_numpy() / measured,
                fmt="none", color="black")
    ax.set_xticks(x)
    ax.set_xticklabels(results["rxnName"].str.replace(" exchange", "",
                                                      regex=False))
    ax.set_ylabel("Relative value")
    ax.legend()

    if write_output:
        out_dir = REPO_PATH / "data" / "testResults"
        out_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_dir / "anaerobic_products.png", bbox_inches="tight")
        plt.close(fig)
