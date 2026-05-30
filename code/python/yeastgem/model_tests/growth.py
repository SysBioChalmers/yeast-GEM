"""Chemostat growth-rate validation against Tobias 2013.

Port of ``code/modelTests/growth.m``. Runs FBA at four conditions
(N-/C-limited, aerobic/anaerobic) with the substrate uptake rates
fixed to experimental values, then compares the predicted growth
against the experimental dilution rates.

Returns the R² (coefficient of determination) across all 32 data
points; an R² ≥ 0.9 is typical for a healthy yeast-GEM release.
"""
from __future__ import annotations

import cobra
import numpy as np
import pandas as pd

from yeastgem import biomass, conditions
from yeastgem.io import REPO_PATH, read_yeast_model

# Yeast-GEM reaction ids the chemostat sweep needs to drive.
_GLC_EX_ID = "r_1714"
_O2_EX_ID = "r_1992"
_NH3_EX_ID = "r_1654"
_GROWTH_RXN_ID = "r_2111"

# N-limited nitrogen-source-derepression: glutamine synthase + glycine
# cleavage system are derepressed under N limitation.
_N_DEREPRESS_UB_RXNS = ("r_0472", "r_0501", "r_0507", "r_0509")

# Slice of the Tobias data per condition. Row indices are 0-based,
# upper bound exclusive — same partition as the legacy MATLAB.
_CONDITIONS = (
    ("N-limited aerobic",   slice(0, 9),   "aerobic",   "N"),
    ("C-limited aerobic",   slice(9, 20),  "aerobic",   "C"),
    ("C-limited anaerobic", slice(20, 26), "anaerobic", "C"),
    ("N-limited anaerobic", slice(26, 32), "anaerobic", "N"),
)


def growth(
    model: cobra.Model | None = None,
    *,
    write_output: bool = False,
    plot: bool = False,
) -> float:
    """Return R² of predicted vs experimental growth rate.

    Parameters
    ----------
    model
        Model to test. Defaults to a fresh :func:`read_yeast_model` load.
    write_output
        When True, save the scatter plot to
        ``data/testResults/growth.png`` and a markdown report to
        ``growth.md`` alongside it. Implies ``plot=True``.
    plot
        Draw the figure into the active matplotlib axes (does not save).
    """
    if model is None:
        model = read_yeast_model()

    exp = _load_tobias_chemostat_data()
    predicted = np.zeros(len(exp))
    for _label, rows, oxygen_mode, lim_mode in _CONDITIONS:
        sub = exp.iloc[rows]
        preds = _simulate_chemostat(model, sub, oxygen_mode, lim_mode)
        predicted[rows] = preds["growth"]

    exp_growth = exp["growth"].to_numpy(dtype=float)
    r2 = float(np.corrcoef(exp_growth, predicted)[0, 1] ** 2)

    if plot or write_output:
        _plot(exp, predicted, r2, write_output=write_output)
    return r2


# --- chemostat sweep ---------------------------------------------------

def _simulate_chemostat(
    base_model: cobra.Model,
    exp: pd.DataFrame,
    oxygen_mode: str,
    lim_mode: str,
) -> pd.DataFrame:
    """Solve FBA at each row of ``exp`` and return the matching uptake
    rates + growth.

    Mirrors the inner ``simulateChemostat`` of ``growth.m``: switch the
    base model into anaerobic / N-limited mode if requested, then for
    each data row fix the substrate uptake rates and maximise growth.
    """
    model = base_model.copy()

    if oxygen_mode == "anaerobic":
        conditions.apply(model, "anaerobic")

    if lim_mode == "N":
        # Protein content under NH3-lim 0.1/h chemostat (Lahtvee et al.
        # 2017, doi:10.1016/j.femsyr.2005.04.003).
        biomass.scale_biomass(model, "protein", 0.28)
        # RNA decreased by the same ~40%, balanced into carbohydrate.
        biomass.scale_biomass(model, "RNA", 0.0329, balance_out="carbohydrate")
        # Derepress glutamate synthase + glycine cleavage system.
        for rxn_id in _N_DEREPRESS_UB_RXNS:
            model.reactions.get_by_id(rxn_id).upper_bound = 1000

    glc = model.reactions.get_by_id(_GLC_EX_ID)
    o2 = model.reactions.get_by_id(_O2_EX_ID)
    nh3 = model.reactions.get_by_id(_NH3_EX_ID)
    growth_rxn = model.reactions.get_by_id(_GROWTH_RXN_ID)
    model.objective = growth_rxn

    out_rows: list[dict[str, float]] = []
    for _, row in exp.iterrows():
        with model:
            _fix_uptake(glc, row["GLCxtI"])
            _fix_uptake(o2, row["O2xtI"])
            _fix_uptake(nh3, row["NH3xtI"])
            try:
                sol = model.optimize()
                if sol.status != "optimal":
                    raise RuntimeError(sol.status)
                out_rows.append({
                    "GLCxtI": abs(sol.fluxes[glc.id]),
                    "O2xtI": abs(sol.fluxes[o2.id]),
                    "NH3xtI": abs(sol.fluxes[nh3.id]),
                    "growth": abs(sol.fluxes[growth_rxn.id]),
                })
            except Exception:  # pragma: no cover - infeasible path
                out_rows.append({"GLCxtI": 0, "O2xtI": 0, "NH3xtI": 0, "growth": 0})
    return pd.DataFrame(out_rows, index=exp.index)


def _fix_uptake(rxn: cobra.Reaction, exp_value: float) -> None:
    """Mirror ``setParam(model,'eq'|'lb',...,-exp_value)``.

    The MATLAB convention treats an experimental value of exactly 1000
    as "open uptake" (set lb to -1000, ub stays at default); any other
    value is a hard equality constraint (lb = ub = -value).
    """
    if abs(exp_value) == 1000:
        rxn.lower_bound = -float(exp_value)
    else:
        rxn.bounds = (-float(exp_value), -float(exp_value))


# --- data + plotting ---------------------------------------------------

def _load_tobias_chemostat_data() -> pd.DataFrame:
    """Read Tobias 2013 chemostat TSV. 32 rows, columns GLCxtI / O2xtI /
    NH3xtI / experimental growth (the latter renamed to ``growth``)."""
    path = REPO_PATH / "data" / "physiology" / "chemostatData_Tobias2013.tsv"
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"experimental growth": "growth"})
    if len(df) != 32:
        raise RuntimeError(f"Expected 32 rows in {path.name}, got {len(df)}")
    return df


def _plot(exp: pd.DataFrame, predicted: np.ndarray, r2: float, *,
          write_output: bool) -> None:
    """Scatter plot of experimental vs predicted growth, colour-coded
    per condition. Saves to ``data/testResults/growth.png`` when
    ``write_output`` is True."""
    import matplotlib

    if write_output:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = np.array([
        [215, 25, 28], [253, 174, 97], [171, 217, 233], [44, 123, 182]
    ]) / 256
    markers = ("o", "s", "d", ">")
    fig, ax = plt.subplots()
    for (label, rows, _oxy, _lim), color, marker in zip(
        _CONDITIONS, colors, markers, strict=True
    ):
        ax.plot(
            exp["growth"].iloc[rows], predicted[rows],
            linestyle="", marker=marker, markersize=8,
            markeredgecolor="k", markerfacecolor=color, label=label,
        )
    lim = max(exp["growth"].max(), predicted.max()) + 0.05
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.plot([0, lim], [0, lim], "--", color=np.array([64, 64, 64]) / 256)
    ax.set_xlabel("Experimental growth rate [1/h]")
    ax.set_ylabel("In silico growth rate [1/h]")
    ax.legend(loc="upper left")
    ax.text(0.25 * lim, 0.1 * lim, f"R² = {r2:.4g}")

    if write_output:
        out_dir = REPO_PATH / "data" / "testResults"
        out_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_dir / "growth.png", bbox_inches="tight")
        (out_dir / "growth.md").write_text(
            f"## R2 of growth rate prediction\n{r2:.4g}\n\n"
            "![Growth curve](growth.png)\n"
        )
        plt.close(fig)
