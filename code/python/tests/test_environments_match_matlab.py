"""Python ``conditions.apply`` makes the changes that MATLAB ``applyEnvironment`` makes.

``reference/environments/<name>.tsv`` holds, per environment in
``data/conditions``, the bounds and stoichiometric coefficients that MATLAB
changed (written by ``reference/dumpEnvironments.m``). A new or changed
environment file needs its reference regenerated with MATLAB.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

from yeastgem import conditions
from yeastgem.paths import find_repo_root

_REFERENCE = Path(__file__).parent / "reference" / "environments"
_NAMES = sorted(p.stem for p in (find_repo_root() / "data" / "conditions").glob("*.yml"))


def _changes(before, after):
    bounds = {r.id: r.bounds for r in after.reactions if r.bounds != before.reactions.get_by_id(
        r.id).bounds}
    coefficients = {}
    for rxn in after.reactions:
        old = before.reactions.get_by_id(rxn.id).metabolites
        new = rxn.metabolites
        for met in {m.id for m in old} | {m.id for m in new}:
            was = next((c for m, c in old.items() if m.id == met), 0.0)
            now = next((c for m, c in new.items() if m.id == met), 0.0)
            if abs(now - was) > 1e-12:
                coefficients[(met, rxn.id)] = now
    return bounds, coefficients


def _reference(name):
    path = _REFERENCE / f"{name}.tsv"
    assert path.is_file(), f"No MATLAB reference for {name}: run reference/dumpEnvironments.m"
    bounds, coefficients = {}, {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["kind"] == "bounds":
                bounds[row["id"]] = (float(row["value"]), float(row["value2"]))
            else:
                coefficients[(row["id"], row["id2"])] = float(row["value"])
    return bounds, coefficients


def _open_exchanges(model):
    """The model with every exchange reaction opened, as in dumpEnvironments.m: the
    model file already holds the minimal medium."""
    opened = model.copy()
    for rxn in opened.reactions:
        if not any(c > 0 for c in rxn.metabolites.values()) \
                or not any(c < 0 for c in rxn.metabolites.values()):
            rxn.bounds = (-1000.0, 1000.0)
    return opened


@pytest.mark.parametrize("name", _NAMES)
def test_environment_changes_equal_the_matlab_changes(model, name):
    base = _open_exchanges(model)
    applied = conditions.apply(base.copy(), name)
    bounds, coefficients = _changes(base, applied)
    expected_bounds, expected_coefficients = _reference(name)
    assert bounds == expected_bounds
    assert coefficients.keys() == expected_coefficients.keys()
    for key, value in coefficients.items():
        assert value == pytest.approx(expected_coefficients[key], abs=1e-9)


def test_every_reference_belongs_to_an_environment():
    assert {p.stem for p in _REFERENCE.glob("*.tsv")} <= set(_NAMES)
