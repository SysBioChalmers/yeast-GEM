"""ATP maintenance (GAM and NGAM) of yeast-GEM: data, simulation and fit.

GAM is the ATP hydrolysed per gram of biomass formed (mmol/gDW), held in
the ATP, water, ADP, phosphate and proton coefficients of the biomass
pseudoreaction. NGAM is the ATP hydrolysed per hour and gram of biomass
(mmol/gDW/h), held as the fixed flux of the non-growth associated
maintenance reaction.

The two are determined separately, because the data that fix them differ:

* NGAM from the maintenance coefficient measured in retentostats and in
  chemostats at very low growth rates
  (``data/physiology/maintenanceRetentostat.tsv``). The glucose uptake
  that the model needs at zero growth is proportional to NGAM, so NGAM
  is the value at which the model reproduces the measured coefficient
  (:func:`ngam_from_retentostat`). The ATP per glucose comes from the
  model, not from the assumption the authors used to convert to ATP.
* GAM from the glucose uptake of C-limited chemostats at the measured
  dilution rate (``data/physiology/maintenanceChemostats.tsv``), with
  NGAM held at that value (:func:`fit_gam`).

Aerobic and anaerobic growth are fitted separately. The chemostat data
come from many studies; the fit weights every study equally (each unit
counts 1/sqrt of the number of units of its study, and units read from
figures or from compilations count half) and uses a Huber loss on the
log ratio of model to measured uptake, so that one study that disagrees
with the others does not move the result. The spread between studies is
estimated by resampling the studies (:func:`bootstrap_gam`).

The uptake at a fixed growth rate does not depend on enzyme constraints
as long as the cell is not limited by enzyme capacity, which holds for
all fitted chemostats (D below 0.3 1/h), so the fit applies to the
enzyme-constrained model as well.

Command line::

    python -m yeastgem.maintenance fit [--bootstrap N] [--write]
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass

import cobra
import numpy as np
import pandas as pd

from yeastgem import biomass, conditions
from yeastgem.io import REPO_PATH, load_yeast_yaml

_DATA = REPO_PATH / "data" / "physiology"
CHEMOSTAT_FILE = _DATA / "maintenanceChemostats.tsv"
RETENTOSTAT_FILE = _DATA / "maintenanceRetentostat.tsv"
FIT_FILE = _DATA / "maintenanceFit.tsv"

_GLC_EX_ID = "r_1714"
_O2_EX_ID = "r_1992"
_GROWTH_RXN_ID = "r_2111"
_BIOMASS_RXN_ID = "r_4041"
_NGAM_RXN_ID = "r_4046"
_ATP_ID = "s_0434"
# Metabolites of the GAM term of the biomass pseudoreaction:
# ATP + H2O -> ADP + phosphate + H+, with the sign they have there.
_GAM_TERMS = {"s_0434": -1, "s_0803": -1, "s_0394": 1, "s_1322": 1, "s_0794": 1}

# Studies whose units count half in the fit: values read from figures,
# and values taken from a compilation of other studies.
_HALF_WEIGHT = ("figure", "compiled")
_HUBER_DELTA = 0.15
# Penalty (log ratio) for a dilution rate that the model cannot sustain.
_INFEASIBLE_LOG_RATIO = 0.7


# --- GAM and NGAM on the model ---------------------------------------------


def get_gam(model: cobra.Model) -> float:
    """GAM of the model: the ATP coefficient of the biomass pseudoreaction."""
    rxn = model.reactions.get_by_id(_BIOMASS_RXN_ID)
    return -float(rxn.metabolites[model.metabolites.get_by_id(_ATP_ID)])


def set_maintenance(model: cobra.Model, gam: float, ngam: float) -> None:
    """Set GAM and NGAM in place.

    The five GAM terms of the biomass pseudoreaction move by the same
    amount, so the proton coefficient keeps its offset from the ATP
    coefficient (the offset balances the charge of the other biomass
    terms). Inside ``with model:`` both changes revert on exit.
    """
    biomass_rxn = model.reactions.get_by_id(_BIOMASS_RXN_ID)
    change = gam - get_gam(model)
    biomass_rxn.add_metabolites(
        {model.metabolites.get_by_id(met): sign * change
         for met, sign in _GAM_TERMS.items()},
        combine=True,
    )
    model.reactions.get_by_id(_NGAM_RXN_ID).bounds = (ngam, ngam)


# --- data ----------------------------------------------------------------


def load_chemostat_data(path=CHEMOSTAT_FILE) -> pd.DataFrame:
    """C-limited chemostats: one row per dilution rate and study.

    Uptakes and productions are in mmol/gDW/h; ``protein`` is g/gDW where
    measured.
    """
    return pd.read_csv(path, sep="\t")


def load_retentostat_data(path=RETENTOSTAT_FILE) -> pd.DataFrame:
    """Maintenance coefficients measured at near-zero growth rates."""
    return pd.read_csv(path, sep="\t")


def study_weights(data: pd.DataFrame) -> np.ndarray:
    """Weight of every row: 1/sqrt(rows of its study), halved for figure
    and compiled values."""
    counts = data.groupby(["reference", "oxygen"])["D"].transform("size")
    weight = 1.0 / np.sqrt(counts.to_numpy(dtype=float))
    return weight * np.where(data["source_type"].isin(_HALF_WEIGHT), 0.5, 1.0)


# --- simulation ------------------------------------------------------------


def regime_models(model: cobra.Model) -> dict[str, cobra.Model]:
    """The model as it is used for aerobic and for anaerobic growth."""
    aerobic = model.copy()
    anaerobic = model.copy()
    conditions.apply(anaerobic, "anaerobic")
    return {"aerobic": aerobic, "anaerobic": anaerobic}


def min_glucose_uptake(
    model: cobra.Model, growth: float, *, protein: float | None = None,
) -> float:
    """Smallest glucose uptake (mmol/gDW/h) at a fixed growth rate.

    ``protein`` (g/gDW) sets the protein content of the biomass, balanced
    with the carbohydrate content. NaN when the growth rate cannot be
    reached.
    """
    with model:
        if protein is not None and not np.isnan(protein):
            biomass.scale_biomass(
                model, "protein", protein, balance_out="carbohydrate",
            )
        glucose = model.reactions.get_by_id(_GLC_EX_ID)
        glucose.lower_bound = -1000.0
        growth_rxn = model.reactions.get_by_id(_GROWTH_RXN_ID)
        growth_rxn.bounds = (growth, growth)
        model.objective = {glucose: 1.0}  # fluxes of uptake are negative
        value = model.slim_optimize(error_value=np.nan)
    return -value


def uptake_grid(
    model: cobra.Model, data: pd.DataFrame, gams, ngam: float,
) -> np.ndarray:
    """Model uptake for every GAM in ``gams`` (rows) and data row (columns)."""
    out = np.full((len(gams), len(data)), np.nan)
    protein = data["protein"].to_numpy(dtype=float)
    growth = data["D"].to_numpy(dtype=float)
    for i, gam in enumerate(gams):
        with model:
            set_maintenance(model, float(gam), ngam)
            for j in range(len(data)):
                out[i, j] = min_glucose_uptake(
                    model, growth[j], protein=protein[j],
                )
    return out


# --- fit -------------------------------------------------------------------


def _huber(r: np.ndarray, delta: float = _HUBER_DELTA) -> np.ndarray:
    a = np.abs(r)
    return np.where(a <= delta, 0.5 * r * r, delta * (a - 0.5 * delta))


def _log_ratio(grid: np.ndarray, measured: np.ndarray) -> np.ndarray:
    with np.errstate(all="ignore"):
        r = np.log(grid / measured[None, :])
    return np.where(np.isfinite(r), r, _INFEASIBLE_LOG_RATIO)


def _best_gam(gams, loss: np.ndarray) -> float:
    return float(gams[int(np.argmin(loss))])


def ngam_from_retentostat(
    model: cobra.Model, ms_glucose: float, *, gam: float | None = None,
) -> float:
    """NGAM at which the model needs ``ms_glucose`` mmol glucose/gDW/h at
    zero growth.

    The glucose uptake at zero growth is proportional to NGAM, so one
    solve is enough.
    """
    with model:
        if gam is not None:
            set_maintenance(model, gam, 1.0)
        else:
            model.reactions.get_by_id(_NGAM_RXN_ID).bounds = (1.0, 1.0)
        per_unit = min_glucose_uptake(model, 0.0)
    return float(ms_glucose / per_unit)


@dataclass
class MaintenanceFit:
    """Fitted values for one oxygen regime."""

    oxygen: str
    ngam: float
    gam: float
    gam_p16: float
    gam_p84: float
    units: int
    studies: int
    median_error: float


def fit_gam(
    model: cobra.Model, data: pd.DataFrame, ngam: float, gams=None,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Best GAM for ``data`` (one oxygen regime) at a given NGAM.

    Returns the GAM, the GAM grid and the matrix of log ratios (GAM by
    row of ``data``) so the caller can resample without new solves.
    """
    if gams is None:
        gams = np.arange(25.0, 70.01, 1.25)
    measured = data["glucose_uptake"].to_numpy(dtype=float)
    ratios = _log_ratio(uptake_grid(model, data, gams, ngam), measured)
    weight = study_weights(data)
    loss = (_huber(ratios) * weight[None, :]).sum(axis=1) / weight.sum()
    return _best_gam(gams, loss), np.asarray(gams), ratios


def bootstrap_gam(
    gams, ratios: np.ndarray, data: pd.DataFrame, n: int = 400, seed: int = 0,
) -> tuple[float, float]:
    """16th and 84th percentile of the best GAM when the studies are
    resampled with replacement."""
    rng = np.random.default_rng(seed)
    weight = study_weights(data)
    study = (data["reference"] + data["oxygen"]).to_numpy()
    names = sorted(set(study))
    best = []
    for _ in range(n):
        drawn = rng.choice(len(names), len(names))
        count = np.zeros(len(data))
        for k in drawn:
            count[study == names[k]] += 1
        w = weight * count
        if w.sum() == 0:
            continue
        loss = (_huber(ratios) * w[None, :]).sum(axis=1) / w.sum()
        best.append(_best_gam(gams, loss))
    return float(np.percentile(best, 16)), float(np.percentile(best, 84))


def fit(
    model: cobra.Model | None = None, *, bootstrap: int = 400, gams=None,
) -> list[MaintenanceFit]:
    """NGAM from the retentostats and GAM from the chemostats, for aerobic
    and anaerobic growth."""
    model = load_yeast_yaml() if model is None else model
    chemostats = load_chemostat_data()
    retentostats = load_retentostat_data()
    models = regime_models(model)
    results = []
    for oxygen, regime_model in models.items():
        ms = float(retentostats.loc[retentostats["oxygen"] == oxygen,
                                    "ms_glucose"].iloc[0])
        ngam = ngam_from_retentostat(regime_model, ms)
        data = chemostats[chemostats["oxygen"] == oxygen].reset_index(drop=True)
        gam, grid, ratios = fit_gam(regime_model, data, ngam, gams)
        lo, hi = bootstrap_gam(grid, ratios, data, bootstrap)
        best = ratios[int(np.argmin(np.abs(grid - gam)))]
        results.append(MaintenanceFit(
            oxygen=oxygen, ngam=ngam, gam=gam, gam_p16=lo, gam_p84=hi,
            units=len(data), studies=data["reference"].nunique(),
            median_error=float(np.median(np.abs(np.exp(best) - 1.0))),
        ))
    return results


def write_fit(results: list[MaintenanceFit], path=FIT_FILE) -> None:
    pd.DataFrame([r.__dict__ for r in results]).round(3).to_csv(
        path, sep="\t", index=False,
    )


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("command", choices=["fit"])
    parser.add_argument("--bootstrap", type=int, default=400)
    parser.add_argument("--write", action="store_true",
                        help=f"write the result to {FIT_FILE.name}")
    args = parser.parse_args(argv)
    results = fit(bootstrap=args.bootstrap)
    for r in results:
        print(f"{r.oxygen:9s} NGAM {r.ngam:.2f}  GAM {r.gam:.2f} "
              f"({r.gam_p16:.1f}-{r.gam_p84:.1f})  {r.units} units, "
              f"{r.studies} studies, median error {r.median_error * 100:.1f}%")
    if args.write:
        write_fit(results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
