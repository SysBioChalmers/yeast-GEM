"""Tests for ``yeastgem.maintenance`` (GAM and NGAM: data, model values, fit)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from raven_toolbox.biomass import (
    huber_loss,
    log_ratio,
    min_uptake,
    ngam_from_maintenance_coefficient,
    uptake_grid,
)

from yeastgem import conditions, maintenance

# GAM and NGAM that the model and the anaerobic environment carry: the
# fitted values of data/physiology/maintenanceFit.tsv, rounded.
AEROBIC_GAM, AEROBIC_NGAM = 39.0, 0.75
ANAEROBIC_GAM, ANAEROBIC_NGAM = 50.0, 1.0


def _charge(model, rxn_id):
    rxn = model.reactions.get_by_id(rxn_id)
    return sum((met.charge or 0) * coef for met, coef in rxn.metabolites.items())


# --- data files (no model needed) ------------------------------------------


def test_chemostat_data_shape():
    data = maintenance.load_chemostat_data()
    assert set(data["oxygen"]) == {"aerobic", "anaerobic"}
    assert set(data["source_type"]) <= {"table", "derived", "figure", "compiled"}
    assert data["glucose_uptake"].notna().all() and (data["glucose_uptake"] > 0).all()
    assert (data["D"] > 0).all() and (data["D"] < 0.4).all()
    assert data["reference"].notna().all()
    # Many studies, none dominating.
    per_study = data.groupby("reference").size()
    assert len(per_study) >= 20 and per_study.max() <= 0.25 * len(data)


def test_retentostat_data_shape():
    data = maintenance.load_retentostat_data()
    assert set(data["oxygen"]) == {"aerobic", "anaerobic"}
    assert (data["ms_glucose"] > 0).all()
    anaerobic = data[data["oxygen"] == "anaerobic"].iloc[0]
    assert anaerobic["m_atp"] == pytest.approx(anaerobic["ms_glucose"] * anaerobic["atp_per_glucose"])


def test_study_weights():
    data = pd.DataFrame({
        "reference": ["a"] * 4 + ["b"] * 4 + ["c"],
        "oxygen": ["aerobic"] * 9,
        "D": [0.1] * 9,
        "source_type": ["table"] * 4 + ["figure"] * 4 + ["table"],
    })
    w = maintenance.study_weights(data)
    assert w[:4] == pytest.approx(0.5)      # 1/sqrt(4)
    assert w[4:8] == pytest.approx(0.25)    # halved for figure values
    assert w[8] == pytest.approx(1.0)


# --- GAM and NGAM on the model ---------------------------------------------


def test_model_carries_the_fitted_aerobic_values(model):
    assert maintenance.get_gam(model) == pytest.approx(AEROBIC_GAM)
    assert model.reactions.get_by_id("r_4046").bounds == (AEROBIC_NGAM, AEROBIC_NGAM)


def test_fit_file_matches_the_model():
    fit = pd.read_csv(maintenance.FIT_FILE, sep="\t").set_index("oxygen")
    for oxygen, gam, ngam in (("aerobic", AEROBIC_GAM, AEROBIC_NGAM),
                              ("anaerobic", ANAEROBIC_GAM, ANAEROBIC_NGAM)):
        row = fit.loc[oxygen]
        assert row["gam_p16"] <= gam <= row["gam_p84"]
        assert abs(row["gam"] - gam) <= 1.5
        assert abs(row["ngam"] - ngam) <= 0.01


def test_anaerobic_environment_carries_the_anaerobic_values(model):
    anaerobic = model.copy()
    conditions.apply(anaerobic, "anaerobic")
    assert maintenance.get_gam(anaerobic) == pytest.approx(ANAEROBIC_GAM)
    assert anaerobic.reactions.get_by_id("r_4046").bounds == (ANAEROBIC_NGAM, ANAEROBIC_NGAM)


def test_set_maintenance_reverts_and_keeps_the_charge_balance(model):
    charge = _charge(model, "r_4041")
    with model:
        maintenance.set_maintenance(model, 30.0, 0.2)
        assert maintenance.get_gam(model) == pytest.approx(30.0)
        assert model.reactions.get_by_id("r_4046").bounds == (0.2, 0.2)
        assert _charge(model, "r_4041") == pytest.approx(charge)
    assert maintenance.get_gam(model) == pytest.approx(AEROBIC_GAM)
    assert model.reactions.get_by_id("r_4046").bounds == (AEROBIC_NGAM, AEROBIC_NGAM)


# --- simulation and fit ------------------------------------------------------


def test_uptake_at_zero_growth_is_proportional_to_ngam(model):
    with model:
        values = []
        for ngam in (0.5, 1.0):
            maintenance.set_maintenance(model, AEROBIC_GAM, ngam)
            values.append(min_uptake(model, 0.0, maintenance.CONFIG))
    assert values[1] == pytest.approx(2 * values[0], rel=1e-6)


def test_ngam_reproduces_the_retentostat_coefficients(model):
    retentostat = maintenance.load_retentostat_data().set_index("oxygen")
    models = maintenance.regime_models(model)
    for oxygen, ngam in (("aerobic", AEROBIC_NGAM), ("anaerobic", ANAEROBIC_NGAM)):
        ms = float(retentostat.loc[oxygen, "ms_glucose"])
        assert ngam_from_maintenance_coefficient(
            models[oxygen], ms, maintenance.CONFIG) == pytest.approx(
            ngam, abs=0.01)


def test_chemostat_uptake_is_lowest_near_the_fitted_gam(model):
    """With NGAM fixed, the weighted loss over a few studies is smaller at the
    fitted GAM than at a GAM 15 mmol/gDW lower or higher."""
    data = maintenance.load_chemostat_data()
    data = data[(data["oxygen"] == "aerobic")
                & data["study"].isin(["Van Hoek 1998", "Overkamp 2000"])]
    data = data.reset_index(drop=True)
    gams = np.array([AEROBIC_GAM - 15, AEROBIC_GAM, AEROBIC_GAM + 15])
    grid = uptake_grid(
        model, data, gams, AEROBIC_NGAM, maintenance.CONFIG, growth_col="D",
        condition_setup=maintenance.set_chemostat_protein,
    )
    loss = huber_loss(log_ratio(
        grid, data["glucose_uptake"].to_numpy(dtype=float))).sum(axis=1)
    assert loss[1] < loss[0] and loss[1] < loss[2]
