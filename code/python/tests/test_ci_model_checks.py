"""Tests for the structural checks of the pull-request report."""
from __future__ import annotations

import csv

import cobra

from tests.ci import model_checks as mc


def _metabolite(identifier: str, name: str, compartment: str = "c") -> cobra.Metabolite:
    return cobra.Metabolite(identifier, name=name, compartment=compartment)


def _reaction(identifier: str, stoichiometry: dict, lower: float = -1000, upper: float = 1000,
              name: str = "") -> cobra.Reaction:
    reaction = cobra.Reaction(identifier, name=name, lower_bound=lower, upper_bound=upper)
    reaction.add_metabolites(stoichiometry)
    return reaction


def _write_tsv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def test_duplicate_keys_names_the_entry_and_the_lines(tmp_path):
    model_file = tmp_path / "model.yml"
    model_file.write_text(
        "!!omap\n"
        "- metabolites:\n"
        "  - !!omap\n"
        "    - id: s_0001\n"
        "    - name: first\n"
        "    - name: second\n"
        "  - !!omap\n"
        "    - id: s_0002\n"
        "    - name: only\n",
        encoding="utf-8",
    )
    count, rows = mc.check_duplicate_keys(model_file)
    assert count == 1
    assert rows == [("s_0001", "field", "name", 5, 6)]


def test_duplicate_keys_clean_file(tmp_path):
    model_file = tmp_path / "model.yml"
    model_file.write_text(
        "!!omap\n- metabolites:\n  - !!omap\n    - id: s_0001\n    - name: first\n",
        encoding="utf-8",
    )
    assert mc.check_duplicate_keys(model_file) == (0, [])


def test_empty_reactions():
    model = cobra.Model("toy")
    a, b = _metabolite("s_1", "A"), _metabolite("s_2", "B")
    model.add_reactions([cobra.Reaction("r_1", name="nothing"), _reaction("r_2", {a: -1, b: 1})])
    assert mc.check_empty_reactions(model) == (1, [("r_1", "nothing")])


def test_reaction_sanity(tmp_path):
    _write_tsv(tmp_path / "genes.tsv", ["id", "uniprot"], [["G1", "P1"]])
    model = cobra.Model("toy")
    a, b = _metabolite("s_1", "A"), _metabolite("s_2", "B")
    good = _reaction("r_good", {a: -1, b: 1})
    good.gene_reaction_rule = "G1"
    inverted = _reaction("r_inverted", {a: -1, b: 1})
    wide = _reaction("r_wide", {a: -1, b: 1}, lower=-1e6, upper=1000)
    unknown = _reaction("r_unknown", {a: -1, b: 1})
    unknown.gene_reaction_rule = "G1 and G2"
    exchange = _reaction("EX_a", {a: -1})
    exchange.gene_reaction_rule = "G1"
    model.add_reactions([good, inverted, wide, unknown, exchange])
    inverted._lower_bound, inverted._upper_bound = 5, 1  # the solver refuses these bounds

    count, rows = mc.check_reaction_sanity(model, tmp_path)

    issues = {row[0]: row[2] for row in rows}
    assert count == 4
    assert "r_good" not in issues
    assert issues["r_inverted"] == "lb>ub"
    assert issues["r_wide"] == "|bound|>1000"
    assert issues["r_unknown"] == "undocumented_genes:G2"
    assert issues["EX_a"] == "boundary_with_gpr"


def test_split_compartments_flags_chemistry_but_not_transport():
    model = cobra.Model("toy")
    atp_c, atp_e = _metabolite("s_1", "ATP", "c"), _metabolite("s_2", "ATP", "e")
    glucose_c, glucose_m = _metabolite("s_3", "glucose", "c"), _metabolite("s_4", "glucose", "m")
    x_c, y_m, z_c = (_metabolite("s_5", "X", "c"), _metabolite("s_6", "Y", "m"),
                     _metabolite("s_7", "Z", "c"))
    transport = _reaction("r_transport", {atp_c: -1, atp_e: 1})
    split = _reaction("r_split", {x_c: -1, y_m: -1, z_c: 1})
    within = _reaction("r_within", {glucose_c: -1, z_c: 1})
    slime = _reaction("r_slime", {x_c: -1, y_m: 1}, name="lipid SLIME rxn")
    lumen, membrane = _metabolite("s_8", "ceramide", "g"), _metabolite("s_9", "lipid", "gm")
    organelle = _reaction("r_organelle", {lumen: -1, membrane: 1})
    model.add_reactions([transport, split, within, slime, organelle,
                         _reaction("r_mixed", {glucose_c: -1, glucose_m: 1})])

    count, rows = mc.check_split_compartments(model)

    assert count == 1
    assert rows[0][0] == "r_split"
    assert rows[0][2] == "c,m"
