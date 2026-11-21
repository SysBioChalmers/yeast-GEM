"""ATP maintenance (GAM and NGAM) of yeast-GEM: data files and the fit.

GAM is the ATP hydrolysed per gram of biomass formed (mmol/gDW), held in the ATP, water, ADP,
phosphate and proton coefficients of the biomass pseudoreaction. NGAM is the ATP hydrolysed per
hour and gram of biomass (mmol/gDW/h), held as the fixed flux of the non-growth associated
maintenance reaction. The method that fits them, with its loss, study weights and resampling of
studies, is :func:`raven_toolbox.biomass.fit_maintenance`; this module supplies what is specific to
yeast-GEM:

* the data files in ``data/physiology``:
  ``maintenanceRetentostat.tsv`` has the maintenance coefficient of glucose (``ms_glucose``) in
  retentostats and in chemostats at very low growth rates, and ``maintenanceChemostats.tsv`` the
  glucose uptake of C-limited chemostats at the measured dilution rate ``D``, with the study of
  every row (``reference``), the origin of the value (``source_type``) and the protein content of
  the biomass where it was measured;
* the ids of the biomass pseudoreaction, the NGAM reaction and the cofactors (``CONFIG``, from
  ``data/yeastgem/ids.yml``);
* the protein content of each chemostat, which is set in the biomass before the uptake is computed
  (:func:`set_chemostat_protein`);
* the two regimes, aerobic and anaerobic growth, which are fitted separately with the model as it
  is used in each (:func:`regime_models`).

The uptake at a fixed growth rate does not depend on enzyme constraints as long as the cell is not
limited by enzyme capacity, which holds for all fitted chemostats (D below 0.3 1/h), so the fit
applies to the enzyme-constrained model as well.

Command line::

    python -m yeastgem.maintenance fit [--bootstrap N] [--write]
"""
from __future__ import annotations

import argparse
from dataclasses import asdict

import cobra
import numpy as np
import pandas as pd
from raven_toolbox.biomass import (
    MaintenanceConfig,
    MaintenanceFit,
    fit_maintenance,
)
from raven_toolbox.biomass import (
    get_gam as _get_gam,
)
from raven_toolbox.biomass import (
    set_gam as _set_gam,
)
from raven_toolbox.biomass import (
    study_weights as _study_weights,
)

from yeastgem import biomass, conditions
from yeastgem.config import load_ids
from yeastgem.io import REPO_PATH, load_yeast_yaml

_DATA = REPO_PATH / "data" / "physiology"
CHEMOSTAT_FILE = _DATA / "maintenanceChemostats.tsv"
RETENTOSTAT_FILE = _DATA / "maintenanceRetentostat.tsv"
FIT_FILE = _DATA / "maintenanceFit.tsv"

_ids = load_ids()
CONFIG = MaintenanceConfig(
    biomass_rxn=_ids.biomass_rxn,
    ngam_rxn="r_4046",
    growth_rxn="r_2111",
    substrate_exchange="r_1714",  # D-glucose exchange
    cofactor_met_names=tuple(_ids.gam_cofactors),
)

# Studies whose units count half in the fit: values read from figures,
# and values taken from a compilation of other studies.
_HALF_WEIGHT = ("figure", "compiled")


# --- GAM and NGAM on the model ---------------------------------------------


def get_gam(model: cobra.Model) -> float:
    """GAM of the model (mmol/gDW): the size of the ATP coefficient of the biomass
    pseudoreaction. See :func:`raven_toolbox.biomass.get_gam`."""
    return _get_gam(model, biomass_rxn=CONFIG.biomass_rxn)


def set_maintenance(model: cobra.Model, gam: float, ngam: float) -> None:
    """Set GAM (mmol/gDW) and NGAM (mmol/gDW/h) of the model in place.

    This is :func:`raven_toolbox.biomass.set_gam` with ``keep_offsets_of`` set to ATP: the other
    GAM terms of the biomass pseudoreaction move by the same amount as ATP, so the proton
    coefficient keeps the offset that balances the charge of the other biomass terms. Inside
    ``with model:`` both changes revert on exit.
    """
    _set_gam(
        model, gam, biomass_rxn=CONFIG.biomass_rxn,
        cofactor_met_names=CONFIG.cofactor_met_names,
        ngam_rxn=CONFIG.ngam_rxn, ngam_value=ngam,
        keep_offsets_of=CONFIG.atp_met_name,
    )


# --- data ----------------------------------------------------------------


def load_chemostat_data(path=CHEMOSTAT_FILE) -> pd.DataFrame:
    """Read the C-limited chemostats: one row per dilution rate ``D`` (1/h) and study.

    Uptakes and productions are in mmol/gDW/h; ``protein`` is g/gDW where measured.
    """
    return pd.read_csv(path, sep="\t")


def load_retentostat_data(path=RETENTOSTAT_FILE) -> pd.DataFrame:
    """Read the maintenance coefficients of glucose measured at near-zero growth rates
    (``ms_glucose``, mmol/gDW/h), one row per oxygen regime."""
    return pd.read_csv(path, sep="\t")


def study_weights(data: pd.DataFrame) -> np.ndarray:
    """Weight of every row of the chemostat table: 1/sqrt(rows of its study), halved for
    values read from a figure or taken from a compilation. See
    :func:`raven_toolbox.biomass.study_weights`."""
    return _study_weights(
        data["reference"] + data["oxygen"],
        discounted=data["source_type"].isin(_HALF_WEIGHT),
    )


def set_chemostat_protein(model: cobra.Model, row: pd.Series) -> None:
    """Set the protein content of the biomass to the one measured in the
    chemostat of ``row`` (g/gDW), balanced with the carbohydrate content.
    Nothing changes where the content was not measured."""
    protein = row["protein"]
    if protein is not None and not np.isnan(protein):
        biomass.scale_biomass(model, "protein", protein, balance_out="carbohydrate")


# --- fit -------------------------------------------------------------------


def regime_models(model: cobra.Model) -> dict[str, cobra.Model]:
    """The model as it is used for aerobic and for anaerobic growth."""
    aerobic = model.copy()
    anaerobic = model.copy()
    conditions.apply(anaerobic, "anaerobic")
    return {"aerobic": aerobic, "anaerobic": anaerobic}


def fit(
    model: cobra.Model | None = None, *, bootstrap: int = 400, gams=None,
) -> dict[str, MaintenanceFit]:
    """Fit NGAM from the retentostats and GAM from the chemostats, for aerobic and anaerobic
    growth, with :func:`raven_toolbox.biomass.fit_maintenance`.

    ``model`` defaults to the yeast-GEM of the repository. ``bootstrap`` is the number of
    resamplings of the studies that give the spread of GAM. ``gams`` is the grid of GAM values
    (mmol/gDW) on which the uptake is computed. The result is keyed by the oxygen regime.
    """
    model = load_yeast_yaml() if model is None else model
    chemostats = load_chemostat_data()
    retentostats = load_retentostat_data()
    results = {}
    for oxygen, regime_model in regime_models(model).items():
        ms = float(retentostats.loc[retentostats["oxygen"] == oxygen,
                                    "ms_glucose"].iloc[0])
        data = chemostats[chemostats["oxygen"] == oxygen].reset_index(drop=True)
        results[oxygen] = fit_maintenance(
            regime_model, data, ms, CONFIG,
            study_col="reference", growth_col="D", uptake_col="glucose_uptake",
            discounted=data["source_type"].isin(_HALF_WEIGHT),
            gams=gams, bootstrap=bootstrap,
            condition_setup=set_chemostat_protein,
        )
    return results


def write_fit(results: dict[str, MaintenanceFit], path=FIT_FILE) -> None:
    """Write the fit of each oxygen regime, one row each, to a tab-separated file."""
    rows = [{"oxygen": oxygen, **asdict(r)} for oxygen, r in results.items()]
    pd.DataFrame(rows).round(3).to_csv(path, sep="\t", index=False)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("command", choices=["fit"])
    parser.add_argument("--bootstrap", type=int, default=400)
    parser.add_argument("--write", action="store_true",
                        help=f"write the result to {FIT_FILE.name}")
    args = parser.parse_args(argv)
    results = fit(bootstrap=args.bootstrap)
    for oxygen, r in results.items():
        print(f"{oxygen:9s} NGAM {r.ngam:.2f}  GAM {r.gam:.2f} "
              f"({r.gam_p16:.1f}-{r.gam_p84:.1f})  {r.units} units, "
              f"{r.studies} studies, median error {r.median_error * 100:.1f}%")
    if args.write:
        write_fit(results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
