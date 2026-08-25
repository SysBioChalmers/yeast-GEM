"""Anaerobic intracellular flux validation.

Port of ``code/modelTests/anaerobic_flux_predictions.m``. For each
dataset in ``data/physiology/flux_data_anaerobic.tsv``: fix glucose
uptake to the experimental rate, maximise growth, and compare
predicted fluxes (scaled to 100 * v_i / v_glc) against the
experimental values. Returns the coefficient of determination about
the line of identity across all data points.
"""
from __future__ import annotations

from collections.abc import Iterable

import cobra
import numpy as np
import pandas as pd

from yeastgem.io import REPO_PATH


def anaerobic_flux_predictions(
    model: cobra.Model,
    *,
    threshold: float = 30.0,
    plot: bool = False,
    write_output: bool = False,
) -> tuple[float, float]:
    """Return ``(R², mean_relative_error)`` for the anaerobic flux fit.

    The model must already be in anaerobic state (caller's
    responsibility — typically obtained via
    ``conditions.apply(model, 'anaerobic')``).

    Parameters
    ----------
    threshold
        Upper limit of the plot axes, in 100 · v_i / v_glc units.
        Default 30, matching the MATLAB script. It no longer filters
        the data entering the R² and MRE.
    plot
        Draw a scatter plot into the active matplotlib axes.
    write_output
        Save the plot + a markdown summary under
        ``data/testResults/``. Implies ``plot=True``.
    """
    flux_df = _load_flux_data()
    glc_rxn_id = "r_1714"

    merged_data: list[float] = []
    merged_sim: list[float] = []
    merged_names: list[str] = []
    per_dataset: list[tuple[str, np.ndarray, np.ndarray]] = []

    for dataset_name, sub in flux_df.groupby("dataset", sort=False):
        glc_row = sub.loc[sub["rxn_id"] == glc_rxn_id]
        if glc_row.empty:
            continue
        glc_value = float(glc_row["target_flux"].iloc[0])

        with model:
            model.reactions.get_by_id(glc_rxn_id).bounds = (-glc_value, -glc_value)
            sol = model.optimize()
            if sol.status != "optimal":
                continue
            glc_predicted = sol.fluxes[glc_rxn_id]
            if glc_predicted == 0:
                continue
            sim_for_set: list[float] = []
            data_for_set: list[float] = []
            names_for_set: list[str] = []
            for _, row in sub.iterrows():
                rxn_id = row["rxn_id"]
                if rxn_id not in model.reactions:
                    continue
                sim = abs(-100.0 * sol.fluxes[rxn_id] / glc_predicted)
                sim_for_set.append(sim)
                data_for_set.append(abs(float(row["experimental_flux"])))
                names_for_set.append(rxn_id)
            merged_data.extend(data_for_set)
            merged_sim.extend(sim_for_set)
            merged_names.extend(names_for_set)
            per_dataset.append((
                str(dataset_name),
                np.array(data_for_set),
                np.array(sim_for_set),
            ))

    data = np.array(merged_data)
    sim = np.array(merged_sim)
    if data.size < 2:
        raise RuntimeError(
            f"Need >=2 data points to compute R2; got {data.size}."
        )
    # Coefficient of determination about the line of identity, with the
    # experimental values as reference; all data points, no thresholding.
    # Matches code/modelTests/anaerobic_flux_predictions.m, which reports
    # the CoD rather than Pearson's r squared. ``threshold`` now only
    # bounds the plot axes.
    r2 = float(
        1 - np.sum((data - sim) ** 2) / np.sum((data - data.mean()) ** 2)
    )
    mre_terms = np.abs(sim - data) / np.where(data == 0, np.nan, data)
    mre = float(np.nanmean(mre_terms))

    if plot or write_output:
        _plot(per_dataset, r2, mre, threshold, write_output=write_output)
    return r2, mre


def _load_flux_data() -> pd.DataFrame:
    """Parse the anaerobic flux TSV (headerless, 9 columns)."""
    path = REPO_PATH / "data" / "physiology" / "flux_data_anaerobic.tsv"
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=[
            "name_in", "name_out",
            "exp_low", "target_flux", "experimental_flux",
            "dataset", "rxn_id", "rxn_name", "equation",
        ],
    )
    return df


def _plot(
    per_dataset: Iterable[tuple[str, np.ndarray, np.ndarray]],
    r2: float,
    mre: float,
    threshold: float,
    *,
    write_output: bool,
) -> None:
    import matplotlib

    if write_output:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import colormaps

    cmap = colormaps.get_cmap("tab10")
    fig, ax = plt.subplots()
    for i, (name, data, sim) in enumerate(per_dataset):
        ax.plot(data, sim, "^", color=cmap(i % cmap.N), label=name)
    ax.plot([0, threshold], [0, threshold], "--",
            color=np.array([64, 64, 64]) / 256)
    ax.set_xlim(0, threshold)
    ax.set_ylim(0, threshold)
    ax.text(0.4 * threshold, 0.15 * threshold,
            f"mean relative error: {mre:.4g}\nR² = {r2:.4g}",
            verticalalignment="top")
    ax.set_xlabel("Experimental 100 · v_i / v_glc")
    ax.set_ylabel("In silico 100 · v_i / v_glc")
    ax.legend(loc="upper left", fontsize="small")
    if write_output:
        out_dir = REPO_PATH / "data" / "testResults"
        out_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_dir / "anaerobic_flux_predictions.png",
                    bbox_inches="tight")
        (out_dir / "anaerobic_flux_predictions.md").write_text(
            f"## Anaerobic flux R²\n{r2:.4g}\n\n"
            f"Mean relative error: {mre:.4g}\n\n"
            "![Anaerobic fluxes](anaerobic_flux_predictions.png)\n"
        )
        plt.close(fig)
