"""Constrain yeast-GEM to anaerobic conditions (pure cobrapy).

This is a self-contained port of ``code/otherChanges/anaerobicModel.m`` (and
the ``changeAminoAcidRatio`` / ``sumBioMass`` / ``rescalePseudoReaction``
helpers it relies on) to cobrapy. It has no dependency on the RAVEN toolbox
or raven-python; the only third-party requirement is ``cobra``.

By default yeast-GEM represents aerobic metabolism. :func:`anaerobic_model`
edits exchange reactions and a few intracellular reactions in place so the
model reaches the exchange rates measured under anaerobic batch growth on
minimal glucose media.

Usage::

    import cobra
    from anaerobic_model import anaerobic_model

    model = cobra.io.read_sbml_model("model/yeast-GEM.xml")
    anaerobic_model(model)
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

import cobra

# Repo paths: this file lives in ``code/``; the repo root is its parent.
_REPO_ROOT = Path(__file__).resolve().parent.parent
_AA_TSV = _REPO_ROOT / "data" / "physiology" / "aminoacid_Bjorkeroth2020.tsv"

# Metabolite/reaction IDs reused across helpers.
_PROTON_ID = "s_0794"        # H+ [cytoplasm]
_PROTEIN_RXN = "r_4047"      # protein pseudoreaction
_PROTEIN_MET_NAME = "protein"  # the protein pseudoreaction product (s_3717)

# Atomic weights, identical to the table in ``sumBioMass.m``
# (``parseChemicalFormula``). ``R`` is the generic residue placeholder and
# carries no mass, matching the MATLAB convention.
_ELEMENT_WEIGHTS = {
    "C": 12.01, "H": 1.008, "N": 14.007, "O": 15.999, "P": 30.974,
    "S": 32.06, "R": 0.0, "Fe": 55.845, "K": 39.098, "Na": 22.99,
    "Cl": 35.45, "Mn": 54.938, "Zn": 65.38, "Ca": 40.078, "Mg": 24.305,
    "Cu": 63.546,
}

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]*)(\d*)")


def anaerobic_model(model: cobra.Model) -> cobra.Model:
    """Constrain ``model`` to anaerobic conditions in place.

    Parameters
    ----------
    model : cobra.Model
        yeast-GEM model, aerobic by default.

    Returns
    -------
    cobra.Model
        The same model object, modified to match anaerobic conditions.
    """
    # 1. Cofactor pseudoreaction (r_4598): heme a synthesis requires O2, so
    #    remove heme a from the biomass cofactor pool and rebalance protons.
    cofactor_rxn = model.reactions.get_by_id("r_4598")
    _set_coefficient(cofactor_rxn, model.metabolites.get_by_id("s_3714"), 0.0)
    _rebalance_proton(model, cofactor_rxn)

    # 2. Switch the protein pseudoreaction to the anaerobic amino-acid ratio.
    _change_amino_acid_ratio(model, aerobic=False)

    # 3. Exchange reactions: block O2 uptake, and open the sterol / fatty-acid
    #    / vitamin uptakes that are essential supplements for anaerobic growth.
    lower_bounds = {
        "r_1992": 0.0,        # oxygen (block uptake)
        "r_1757": -1000.0,    # ergosterol
        "r_1915": -1000.0,    # lanosterol
        "r_2106": -1000.0,    # zymosterol
        "r_2134": -1000.0,    # 14-demethyllanosterol
        "r_1994": -1000.0,    # palmitoleate
        "r_2189": -1000.0,    # oleate
        "r_2137": 0.0,        # ergosta-5,7,22,24(28)-tetraen-3beta-ol (block)
        "r_1967": -1000.0,    # nicotinate
        "r_1548": -1000.0,    # (R)-pantothenate
    }
    for rxn_id, lower_bound in lower_bounds.items():
        model.reactions.get_by_id(rxn_id).lower_bound = lower_bound

    # 4. Block MDH2 (r_0714) and IDP2 (r_0659): repressed/undetected under
    #    anaerobic glucose growth.
    for rxn_id in ("r_0714", "r_0659"):
        model.reactions.get_by_id(rxn_id).bounds = (0.0, 0.0)

    # 5. Fumarate reductase recycles FADH2 from Ero1-driven disulphide-bond
    #    formation; add the net FADH2/FAD/H+ turnover to the biomass reaction.
    biomass_rxn = model.reactions.get_by_id("r_4041")
    fadh2_prod = 0.08
    biomass_rxn.add_metabolites(
        {
            model.metabolites.get_by_id("s_0689"): fadh2_prod,        # FADH2
            model.metabolites.get_by_id("s_0687"): -fadh2_prod,       # FAD
            model.metabolites.get_by_id("s_0794"): -2.0 * fadh2_prod,  # H+
        },
        combine=True,
    )
    return model


# --- helpers ----------------------------------------------------------------

def _set_coefficient(
    reaction: cobra.Reaction, metabolite: cobra.Metabolite, value: float
) -> None:
    """Set an *absolute* stoichiometric coefficient (the ``S(i,j) = x``
    semantics of the MATLAB code). Setting ``value`` to 0 removes the
    metabolite from the reaction.
    """
    current = reaction.metabolites.get(metabolite, 0.0)
    delta = float(value) - current
    if delta != 0.0:
        reaction.add_metabolites({metabolite: delta}, combine=True)


def _rebalance_proton(model: cobra.Model, reaction: cobra.Reaction) -> None:
    """Set the H+ coefficient so ``reaction`` is charge balanced.

    Mirrors the H+ correction repeated in ``anaerobicModel.m`` and
    ``rescalePseudoReaction.m``: zero the proton, then set it to the negative
    of the remaining charge imbalance. Missing charges (cobra ``None``) are
    treated as zero, matching MATLAB's ``'omitnan'``.
    """
    proton = model.metabolites.get_by_id(_PROTON_ID)
    _set_coefficient(reaction, proton, 0.0)
    imbalance = sum(
        (met.charge or 0.0) * coeff for met, coeff in reaction.metabolites.items()
    )
    _set_coefficient(reaction, proton, -imbalance)


def _change_amino_acid_ratio(model: cobra.Model, *, aerobic: bool = True) -> None:
    """Rewrite the protein pseudoreaction's amino-acid ratios.

    Port of ``changeAminoAcidRatio.m``. The protein mass is snapshotted, the
    tRNA stoichiometries in ``r_4047`` are replaced with the aerobic or
    anaerobic ratios from ``aminoacid_Bjorkeroth2020.tsv``, and the protein
    pseudoreaction is rescaled so the protein fraction returns to its
    pre-switch value.
    """
    target_protein = _protein_mass(model)

    protein_rxn = model.reactions.get_by_id(_PROTEIN_RXN)
    # TSV columns (0-indexed): 0 aa, 1 substrate met, 2 product met, 3 MW,
    # 4 aerobic ratio, 5 anaerobic ratio.
    ratio_column = 4 if aerobic else 5
    for row in _read_amino_acid_ratios():
        ratio = float(row[ratio_column])
        substrate = model.metabolites.get_by_id(row[1])
        product = model.metabolites.get_by_id(row[2])
        _set_coefficient(protein_rxn, substrate, -ratio)
        _set_coefficient(protein_rxn, product, ratio)

    factor = target_protein / _protein_mass(model)
    _rescale_protein(model, factor)


def _protein_mass(model: cobra.Model) -> float:
    """Protein fraction [g/gDW], port of ``sumBioMass.m`` ``getFraction('P')``.

    Sums the molecular weight of the protein pseudoreaction substrates, with
    two protons (2.016) removed from each charged-tRNA formula.
    """
    protein_rxn = model.reactions.get_by_id(_PROTEIN_RXN)
    total = 0.0
    for met, coeff in protein_rxn.metabolites.items():
        if coeff < 0:  # substrate
            total += -coeff * (_formula_weight(met) - 2.016)
    return total / 1000.0


def _rescale_protein(model: cobra.Model, factor: float) -> None:
    """Multiply every coefficient in the protein pseudoreaction by ``factor``
    except the ``protein`` product, then rebalance H+.

    Port of ``rescalePseudoReaction.m`` for the protein component.
    """
    protein_rxn = model.reactions.get_by_id(_PROTEIN_RXN)
    deltas = {
        met: (factor - 1.0) * coeff
        for met, coeff in protein_rxn.metabolites.items()
        if met.name != _PROTEIN_MET_NAME
    }
    if deltas:
        protein_rxn.add_metabolites(deltas, combine=True)
    _rebalance_proton(model, protein_rxn)


def _formula_weight(metabolite: cobra.Metabolite) -> float:
    """Molecular weight from a metabolite formula, using the same atomic
    weights as ``sumBioMass.m`` (so the generic residue ``R`` weighs 0)."""
    formula = metabolite.formula
    if not formula:
        raise ValueError(
            f"Biomass metabolite {metabolite.id} has an empty formula field."
        )
    weight = 0.0
    for element, count in _FORMULA_TOKEN.findall(formula):
        if not element:
            continue
        try:
            weight += (int(count) if count else 1) * _ELEMENT_WEIGHTS[element]
        except KeyError as exc:
            raise ValueError(
                f"Unknown element '{element}' in formula '{formula}'."
            ) from exc
    return weight


def _read_amino_acid_ratios() -> list[list[str]]:
    """Read ``aminoacid_Bjorkeroth2020.tsv`` (one header line, tab separated).

    Returns the data rows verbatim; columns are
    ``[aa, substrate_met, product_met, MW, aerobic, anaerobic]``.
    """
    rows: list[list[str]] = []
    with open(_AA_TSV, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)  # skip header
        for record in reader:
            if record and record[0].strip():
                rows.append(record)
    return rows
