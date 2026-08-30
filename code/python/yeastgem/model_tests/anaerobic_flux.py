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
from dataclasses import dataclass

import cobra
import numpy as np
import pandas as pd

from yeastgem.io import REPO_PATH

# Below this, in 100 * v_i / v_glucose units, a flux is treated as zero:
# the measurements are not reported to that precision and a ratio taken
# against such a value is noise.
_FLUX_TOL = 1e-9


@dataclass(frozen=True)
class AnaerobicFluxResult:
    """Predicted against measured intracellular fluxes.

    ``median_fold_error`` of 1.2 means the typical prediction is 1.2x out,
    in either direction. ``n_unpredicted`` counts pairs with no comparable
    ratio -- a measured flux predicted as zero, or predicted in the
    opposite direction -- which are errors that a fold error cannot
    express.
    """

    median_fold_error: float
    mean_fold_error: float
    fraction_within_two_fold: float
    n_unpredicted: int
    n_measurements: int
    # One row per measurement: (reaction, measured, predicted, fold error
    # or the reason there is none). Carried out so the pull-request report
    # can show what the median is a summary of.
    measurements: list[tuple]
    # Retained for comparison with earlier releases; see the note in the
    # function body for why it is not the headline number.
    r2: float
    mean_relative_error: float


def anaerobic_flux_predictions(
    model: cobra.Model,
    *,
    threshold: float = 30.0,
    plot: bool = False,
    write_output: bool = False,
) -> AnaerobicFluxResult:
    """Compare predicted against measured intracellular fluxes.

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
            f"Need >=2 data points to compare; got {data.size}."
        )

    # Coefficient of determination about the line of identity, kept so that
    # a release stays comparable with earlier ones. It is a poor summary of
    # this particular fit and is not the headline number: the measured
    # fluxes span three orders of magnitude, so the sum of squares is
    # dominated by a handful of large glycolytic fluxes whose values follow
    # from the glucose constraint and stoichiometry almost regardless of
    # the model. On the committed model it reads 0.997 while the smaller
    # half of the data is out by about 50%.
    r2 = float(
        1 - np.sum((data - sim) ** 2) / np.sum((data - data.mean()) ** 2)
    )
    mre_terms = np.abs(sim - data) / np.where(data == 0, np.nan, data)
    mre = float(np.nanmean(mre_terms))

    # Fold error, which is what the headline numbers are built from. On a
    # quantity spanning 0.1 to 163 the meaningful question is how far out a
    # prediction is *proportionally*: being 2x out matters as much at a
    # flux of 1 as at a flux of 100, and a ratio says that where a
    # difference does not.
    #
    # A ratio needs both values non-zero and pointing the same way. Pairs
    # that fail either test are counted rather than dropped: predicting no
    # flux where flux was measured, or predicting the opposite direction,
    # is a worse error than any finite ratio, and letting it vanish from
    # the denominator would make the score improve as the model got worse.
    comparable = (
        (np.abs(data) > _FLUX_TOL)
        & (np.abs(sim) > _FLUX_TOL)
        & (np.sign(sim) == np.sign(data))
    )
    n_unpredicted = int(data.size - comparable.sum())
    log_ratio = np.abs(
        np.log10(np.abs(sim[comparable]) / np.abs(data[comparable]))
    )

    rows = []
    for name, measured, predicted, ok in zip(
        merged_names, data, sim, comparable, strict=True
    ):
        if ok:
            fold = 10 ** abs(np.log10(abs(predicted) / abs(measured)))
            rows.append((name, f"{measured:.4g}", f"{predicted:.4g}",
                         f"{fold:.3f}"))
        elif abs(predicted) <= _FLUX_TOL:
            rows.append((name, f"{measured:.4g}", f"{predicted:.4g}",
                         "no flux predicted"))
        elif abs(measured) <= _FLUX_TOL:
            rows.append((name, f"{measured:.4g}", f"{predicted:.4g}",
                         "no flux measured"))
        else:
            rows.append((name, f"{measured:.4g}", f"{predicted:.4g}",
                         "opposite direction"))

    result = AnaerobicFluxResult(
        median_fold_error=float(10 ** np.median(log_ratio)),
        mean_fold_error=float(10 ** np.mean(log_ratio)),
        fraction_within_two_fold=float(np.mean(log_ratio < np.log10(2))),
        n_unpredicted=n_unpredicted,
        n_measurements=int(data.size),
        r2=r2,
        mean_relative_error=mre,
        measurements=rows,
    )

    if plot or write_output:
        _plot(per_dataset, r2, mre, threshold, write_output=write_output)
    return result


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
