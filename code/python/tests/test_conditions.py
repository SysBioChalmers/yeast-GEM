"""Tests for ``yeastgem.conditions`` (data-driven condition presets)."""
from __future__ import annotations

from pathlib import Path

import cobra
import pytest

from yeastgem import compare_models, conditions

# --- data-file shape checks (no model load needed) --------------------

def test_minimal_Y6_loads():
    cfg = conditions.load_condition("minimal_Y6")
    assert cfg["name"] == "minimal_Y6"
    assert cfg["prelude"]["reset_exchanges"] == "out"
    assert cfg["expected_uptake_count"] == 15


def test_anaerobic_loads():
    cfg = conditions.load_condition("anaerobic")
    assert cfg["name"] == "anaerobic"
    assert cfg["amino_acid_ratio"] == "anaerobic"
    assert cfg["cofactor_pseudoreaction"]["rxn_id"] == "r_4598"
    assert cfg["biomass_stoichiometry_delta"]["rxn_id"] == "r_4041"


def test_glycine_nitrogen_loads():
    cfg = conditions.load_condition("glycine_nitrogen")
    assert cfg["name"] == "glycine_nitrogen"
    assert {b["rxn"] for b in cfg["bounds"]} == {"r_0501", "r_0507", "r_0509"}


def test_nitrogen_limitation_loads():
    cfg = conditions.load_condition("nitrogen_limitation")
    assert cfg["name"] == "nitrogen_limitation"


def test_unknown_condition_raises():
    with pytest.raises(FileNotFoundError):
        conditions.load_condition("does_not_exist")


# --- application checks (need the model) ------------------------------

def test_apply_glycine_nitrogen_sets_bounds(model):
    mutated = model.copy()
    conditions.apply(mutated, "glycine_nitrogen")
    for rxn_id in ("r_0501", "r_0507", "r_0509"):
        rxn = mutated.reactions.get_by_id(rxn_id)
        assert rxn.lower_bound == 0
        assert rxn.upper_bound == 1000


def test_apply_nitrogen_limitation_sets_bounds(model):
    mutated = model.copy()
    conditions.apply(mutated, "nitrogen_limitation")
    assert mutated.reactions.get_by_id("r_0472").upper_bound == 1000
    for rxn_id in ("r_0501", "r_0507", "r_0509"):
        assert mutated.reactions.get_by_id(rxn_id).lower_bound == 0
        assert mutated.reactions.get_by_id(rxn_id).upper_bound == 1000


def test_apply_carnitine_opens_shuttle(model):
    mutated = model.copy()
    assert mutated.reactions.get_by_id("r_0252").upper_bound == 0
    conditions.apply(mutated, "carnitine")
    assert mutated.reactions.get_by_id("r_0252").upper_bound == 1000
    assert mutated.reactions.get_by_id("r_1545").lower_bound == -1000


@pytest.mark.parametrize("preset", ["glycine_nitrogen", "nitrogen_limitation"])
@pytest.mark.parametrize("environment", [None, "anaerobic"])
def test_nitrogen_presets_keep_the_model_growing(model, preset, environment):
    """The presets open reactions, so growth is feasible and does not fall."""
    mutated = model.copy()
    if environment:
        conditions.apply(mutated, environment)
    before = mutated.slim_optimize()
    conditions.apply(mutated, preset)
    after = mutated.slim_optimize()
    assert before > 1e-6
    assert after >= before - 1e-9


def test_apply_minimal_Y6_caps_glucose_and_zeros_bicarbonate(model):
    mutated = model.copy()
    conditions.apply(mutated, "minimal_Y6")
    glucose = mutated.reactions.get_by_id("r_1714")
    assert glucose.lower_bound == -1
    bicarbonate = mutated.reactions.get_by_id("r_1663")
    assert bicarbonate.lower_bound == 0
    assert bicarbonate.upper_bound == 0
    # Allowed uptakes (sample a few)
    for rxn_id in ("r_1654", "r_1992", "r_2005", "r_2060"):
        assert mutated.reactions.get_by_id(rxn_id).lower_bound == -1000


def test_apply_minimal_Y6_resets_all_exchanges(model):
    """Prelude should set all "out" exchanges to (lb=0, ub=1000) before
    the targeted overrides — confirmed by the bicarbonate/oxygen path."""
    mutated = model.copy()
    # Pick an exchange that the condition does NOT touch and verify the
    # prelude zeroed its lb (uptake blocked) and capped its ub at 1000.
    untouched_exchange = next(
        r for r in mutated.exchanges
        if r.id not in {b["rxn"] for b in conditions.load_condition("minimal_Y6")["bounds"]}
    )
    conditions.apply(mutated, "minimal_Y6")
    assert untouched_exchange.lower_bound == 0
    assert untouched_exchange.upper_bound == 1000


def test_apply_anaerobic_runs_end_to_end(model):
    """The anaerobic environment as a whole. The resulting model must have O2 uptake
    blocked, ergosterol uptake allowed, MDH2 blocked, and the cofactor
    pseudoreaction's heme coefficient set to zero."""
    mutated = model.copy()
    conditions.apply(mutated, "anaerobic")
    assert mutated.reactions.get_by_id("r_1992").lower_bound == 0   # O2 blocked
    assert mutated.reactions.get_by_id("r_1757").lower_bound == -1000  # ergosterol
    assert mutated.reactions.get_by_id("r_0714").bounds == (0, 0)   # MDH2
    cofac = mutated.reactions.get_by_id("r_4598")
    heme = mutated.metabolites.get_by_id("s_3714")
    assert cofac.metabolites.get(heme, 0) == 0


def test_apply_is_idempotent_for_glycine(model):
    """Applying glycine_nitrogen twice must produce the same model."""
    once = model.copy()
    conditions.apply(once, "glycine_nitrogen")
    twice = model.copy()
    conditions.apply(twice, "glycine_nitrogen")
    conditions.apply(twice, "glycine_nitrogen")
    report = compare_models(once, twice)
    assert report.equal, report


# --- partial-anaerobic checks on the real model ---------------------
#
# Single steps of the anaerobic environment, to catch ID drift (heme a,
# FADH2 / FAD / H+, biomass reaction).


def test_anaerobic_cofactor_step_removes_heme_on_real_model(model):
    """Apply only the cofactor step. The cofactor pseudoreaction (r_4598) should lose heme a
    (s_3714)."""
    mutated = model.copy()
    cofac = mutated.reactions.get_by_id("r_4598")
    heme = mutated.metabolites.get_by_id("s_3714")
    assert heme in cofac.metabolites

    full_cfg = conditions.load_condition("anaerobic")
    sub_cfg = {"cofactor_pseudoreaction": full_cfg["cofactor_pseudoreaction"]}
    conditions.apply_condition(mutated, sub_cfg)
    assert cofac.metabolites.get(heme, 0) == 0


def test_anaerobic_biomass_step_adds_fadh2_on_real_model(model):
    """Same idea for the biomass stoichiometry delta block."""
    mutated = model.copy()
    bio = mutated.reactions.get_by_id("r_4041")
    fadh2 = mutated.metabolites.get_by_id("s_0689")
    fad = mutated.metabolites.get_by_id("s_0687")
    proton = mutated.metabolites.get_by_id("s_0794")

    before = {
        fadh2.id: bio.metabolites.get(fadh2, 0),
        fad.id: bio.metabolites.get(fad, 0),
        proton.id: bio.metabolites.get(proton, 0),
    }

    full_cfg = conditions.load_condition("anaerobic")
    sub_cfg = {"biomass_stoichiometry_delta": full_cfg["biomass_stoichiometry_delta"]}
    conditions.apply_condition(mutated, sub_cfg)

    after = bio.metabolites
    assert after[fadh2] == pytest.approx(before[fadh2.id] + 0.08)
    assert after[fad] == pytest.approx(before[fad.id] - 0.08)
    assert after[proton] == pytest.approx(before[proton.id] - 0.16)


def _top_level_imports(path: Path) -> set[str]:
    import ast

    tree = ast.parse(path.read_text())
    return {
        (n.module if isinstance(n, ast.ImportFrom) else a.name)
        for n in ast.walk(tree) if isinstance(n, (ast.Import, ast.ImportFrom))
        for a in (n.names if isinstance(n, ast.Import) else [n])
    }


def test_conditions_module_needs_no_toolbox():
    """Environments depend on no toolbox: yeastgem.conditions imports only the
    standard library, cobra, yaml and yeastgem.paths (standard library only)."""
    import sys

    allowed = set(sys.stdlib_module_names) | {"__future__"}
    imported = _top_level_imports(Path(conditions.__file__))
    assert {m.split(".")[0] for m in imported} <= allowed | {"cobra", "yaml", "yeastgem"}
    assert {m for m in imported if m.startswith("yeastgem")} == {"yeastgem.paths"}
    paths = Path(conditions.__file__).with_name("paths.py")
    assert {m.split(".")[0] for m in _top_level_imports(paths)} <= allowed | {"dotenv"}


def test_load_condition_follows_yeast_gem_path(tmp_path, monkeypatch):
    (tmp_path / "data" / "conditions").mkdir(parents=True)
    (tmp_path / "data" / "conditions" / "test_env.yml").write_text("name: test_env\n")
    monkeypatch.setenv("YEAST_GEM_PATH", str(tmp_path))
    assert conditions.load_condition("test_env")["name"] == "test_env"


def test_load_condition_from_path(tmp_path):
    path = tmp_path / "my_env.yml"
    path.write_text("name: my_env\nbounds:\n  - { rxn: r_1992, lb: 0 }\n")
    assert conditions.load_condition(str(path))["bounds"][0]["rxn"] == "r_1992"


@pytest.mark.parametrize("cfg", [
    {"bounds": [{"rxn": "r_1992", "lb": -5}, {"rxn": "r_0714", "lb": 5, "ub": 1}]},
    {"bounds": [{"rxn": "r_1992", "lb": -5}, {"rxn": "r_0714", "lb": None}]},
    {"bounds": [{"rxn": "r_1992", "lb": -5}, {"rxn": "r_0714", "ub": "high"}]},
    {"prelude": {"reset_exchanges": "outt"}, "bounds": [{"rxn": "r_1992", "lb": -5}]},
    {"amino_acid_ratio": "anoxic", "bounds": [{"rxn": "r_1992", "lb": -5}]},
])
def test_invalid_environment_leaves_model_unchanged(model, cfg):
    mutated = model.copy()
    before = {r.id: r.bounds for r in mutated.reactions}
    protein = dict(mutated.reactions.get_by_id("r_4047").metabolites)
    with pytest.raises(ValueError):
        conditions.apply_condition(mutated, cfg)
    assert {r.id: r.bounds for r in mutated.reactions} == before
    assert dict(mutated.reactions.get_by_id("r_4047").metabolites) == protein


def test_empty_lists_are_allowed(model):
    mutated = model.copy()
    conditions.apply_condition(mutated, {
        "cofactor_pseudoreaction": {"rxn_id": "r_4598", "remove_mets": None},
        "biomass_stoichiometry_delta": {"rxn_id": "r_4041", "add": None},
        "bounds": None,
    })


# --- environments on a small model ------------------------------------

def _toy_model() -> cobra.Model:
    """A -> B with an exchange of A (A =>), an import of B (=> B) and H+."""
    a = cobra.Metabolite("a", compartment="c", charge=0)
    b = cobra.Metabolite("b", compartment="c", charge=-1)
    x = cobra.Metabolite("x", compartment="c")  # without charge
    h = cobra.Metabolite("h", compartment="c", charge=1)
    model = cobra.Model("toy")
    ex_a = cobra.Reaction("EX_a", lower_bound=-5, upper_bound=7)
    ex_a.add_metabolites({a: -1})
    in_b = cobra.Reaction("IN_b", lower_bound=-3, upper_bound=4)
    in_b.add_metabolites({b: 1})
    r = cobra.Reaction("R", lower_bound=-10, upper_bound=10)
    r.add_metabolites({a: -1, b: 1})
    pseudo = cobra.Reaction("P", lower_bound=0, upper_bound=1000)
    pseudo.add_metabolites({a: -1, x: -1, b: 1, h: 0.5})
    model.add_reactions([ex_a, in_b, r, pseudo])
    return model


def _state(model):
    return ({r.id: r.bounds for r in model.reactions},
            {r.id: dict(r.metabolites) for r in model.reactions})


@pytest.mark.parametrize("reset, changed", [
    ("out", {"EX_a"}), ("in", {"IN_b"}), ("all", {"EX_a", "IN_b"}),
])
def test_reset_exchanges_selects_by_the_side_without_products_or_substrates(reset, changed):
    model = _toy_model()
    before = {r.id: r.bounds for r in model.reactions}
    conditions.apply_condition(model, {"prelude": {"reset_exchanges": reset}})
    after = {r.id: r.bounds for r in model.reactions}
    assert {k for k in after if after[k] != before[k]} == changed
    assert all(after[k] == (0, 1000) for k in changed)


@pytest.mark.parametrize("cfg", [
    {"bounds": [{"rxn": "R", "lb": -1}, {"rxn": "nope", "ub": 1}]},
    {"cofactor_pseudoreaction": {"rxn_id": "nope"}},
    {"cofactor_pseudoreaction": {"rxn_id": "P", "remove_mets": [{"met": "nope"}]}},
    {"cofactor_pseudoreaction": {"rxn_id": "P", "charge_balance_met": "nope"}},
    {"biomass_stoichiometry_delta": {"rxn_id": "nope", "add": []}},
    {"biomass_stoichiometry_delta": {"rxn_id": "P", "add": [{"met": "nope", "coef": 1}]}},
])
def test_an_identifier_that_is_not_in_the_model_is_an_error_and_changes_nothing(cfg):
    model = _toy_model()
    before = _state(model)
    cfg = {"prelude": {"reset_exchanges": "all"}, **cfg}
    with pytest.raises(ValueError, match="not in the model"):
        conditions.apply_condition(model, cfg)
    assert _state(model) == before


def test_bound_entries_of_one_reaction_are_combined_in_order():
    model = _toy_model()
    conditions.apply_condition(model, {"bounds": [
        {"rxn": "R", "lb": -1}, {"rxn": "R", "ub": 3}, {"rxn": "R", "lb": -2},
    ]})
    assert model.reactions.R.bounds == (-2, 3)


def test_uptake_count_warns_only_on_a_difference():
    cfg = {"bounds": [{"rxn": "R", "lb": -1000}], "expected_uptake_count": 1}
    conditions.apply_condition(_toy_model(), cfg)
    with pytest.warns(UserWarning, match="Expected 2 uptake reactions, applied 1"):
        conditions.apply_condition(_toy_model(), {**cfg, "expected_uptake_count": 2})


def test_charge_balance_counts_a_metabolite_without_charge_as_zero_and_warns():
    model = _toy_model()
    cfg = {"cofactor_pseudoreaction": {"rxn_id": "P", "charge_balance_met": "h"}}
    with pytest.warns(UserWarning, match="no charge for x"):
        conditions.apply_condition(model, cfg)
    # b has charge -1 with coefficient 1, a and x count 0: h balances it
    assert model.reactions.P.metabolites[model.metabolites.h] == pytest.approx(1.0)


def test_numbers_and_words_are_read_as_in_the_matlab_reader(tmp_path):
    path = tmp_path / "env.yml"
    path.write_text("a: 1e3\nb: 010\nc: yes\nd: true\ne: 'no'\nf: -.5\n")
    cfg = conditions.load_condition(str(path))
    assert cfg == {"a": 1000.0, "b": 10.0, "c": "yes", "d": True, "e": "no", "f": -0.5}


def test_importing_conditions_does_not_import_raven_toolbox():
    import subprocess
    import sys

    code = (
        "import sys; sys.modules['raven_toolbox'] = None\n"
        "import yeastgem.conditions as c\n"
        "import cobra\n"
        "m = cobra.Model('t'); r = cobra.Reaction('r', lower_bound=-1, upper_bound=1)\n"
        "m.add_reactions([r]); c.apply_condition(m, {'bounds': [{'rxn': 'r', 'lb': 0}]})\n"
        "assert r.bounds == (0, 1)\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                            check=False)
    assert result.returncode == 0, result.stderr
