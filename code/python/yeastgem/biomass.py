"""Yeast-specific biomass helpers — wrappers over raven_toolbox.biomass.

The generic mechanism (sum / scale / rescale / set_gam) lives upstream.
This module configures it with the yeast layout (``data/yeastgem/ids.yml``)
and adds the one yeast-only operation: :func:`change_amino_acid_ratio`,
which rewrites the protein pseudoreaction from
``data/physiology/aminoAcid_Bjorkeroth2020.tsv``.

The legacy MATLAB names (``sumBioMass`` etc.) are intentionally not
mirrored — the upstream API names (``sum_biomass`` etc.) are cleaner
and the local wrappers just hand them the yeast :class:`BiomassConfig`.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import cobra
import pandas as pd
from raven_toolbox.biomass import (
    BiomassComponent,
    BiomassConfig,
)
from raven_toolbox.biomass import (
    rescale_pseudoreaction as _ra_rescale_pseudoreaction,
)
from raven_toolbox.biomass import (
    scale_biomass as _ra_scale_biomass,
)
from raven_toolbox.biomass import (
    set_gam as _ra_set_gam,
)
from raven_toolbox.biomass import (
    sum_biomass as _ra_sum_biomass,
)

from yeastgem.config import load_ids
from yeastgem.io import REPO_PATH

_AA_TSV = REPO_PATH / "data" / "physiology" / "aminoAcid_Bjorkeroth2020.tsv"


@lru_cache(maxsize=1)
def yeast_biomass_config() -> BiomassConfig:
    """Build a :class:`BiomassConfig` from ``data/yeastgem/ids.yml``."""
    ids = load_ids()
    components = tuple(
        BiomassComponent(
            name=c.name,
            pseudoreaction_name=ids.pseudoreaction_names[c.name],
            mass_strategy=c.mass_strategy,  # type: ignore[arg-type]
        )
        for c in ids.biomass_components
    )
    return BiomassConfig(
        biomass_rxn=ids.biomass_rxn,
        proton_met=ids.proton_met,
        components=components,
    )


def sum_biomass(model: cobra.Model) -> dict[str, float]:
    """Mass fraction (g/gDW) per yeast-GEM biomass component, plus total."""
    return _ra_sum_biomass(model, yeast_biomass_config())


def scale_biomass(
    model: cobra.Model,
    component: str,
    new_value: float,
    *,
    balance_out: str | None = None,
) -> None:
    """Scale a biomass component to a target g/gDW.

    With ``balance_out`` set, the second component is adjusted so the
    biomass total stays at 1 g/gDW.
    """
    _ra_scale_biomass(
        model, yeast_biomass_config(), component, new_value,
        balance_out=balance_out,
    )


def rescale_pseudoreaction(
    model: cobra.Model,
    component: str,
    factor: float,
) -> None:
    """Multiply the substrate coefficients of a component pseudoreaction
    by ``factor`` and rebalance H+.

    Yeast-specific aggregation: ``component='lipid'`` rescales both
    ``lipid_backbone`` and ``lipid_chain`` together (mirroring the
    legacy MATLAB ``rescalePseudoReaction``).
    """
    cfg = yeast_biomass_config()
    if component == "lipid":
        _ra_rescale_pseudoreaction(model, cfg, "lipid_backbone", factor)
        # lipid_chain is not in the default config (it doesn't contribute
        # to the mass total). Look it up from ids.yml + apply by name.
        _rescale_named_pseudoreaction(model, "lipid chain pseudoreaction", factor)
    else:
        _ra_rescale_pseudoreaction(model, cfg, component, factor)


def set_gam(
    model: cobra.Model,
    value: float,
    *,
    ngam: float | None = None,
) -> None:
    """Set GAM (and optionally NGAM) on the yeast-GEM biomass pseudoreaction.

    NGAM is the reaction whose name is "non-growth associated maintenance
    reaction" (yeast-GEM convention); pass ``ngam`` to fix its bounds at
    ``(ngam, ngam)``. The cofactor metabolite set comes from
    ``ids.yml::gam_cofactors``.
    """
    ids = load_ids()
    ngam_rxn = None
    if ngam is not None:
        # NGAM rxn is looked up by name to mirror the legacy MATLAB.
        for rxn in model.reactions:
            if rxn.name == "non-growth associated maintenance reaction":
                ngam_rxn = rxn.id
                break
        if ngam_rxn is None:
            raise ValueError(
                "Could not find a reaction named "
                "'non-growth associated maintenance reaction' in the model."
            )
    _ra_set_gam(
        model, value,
        biomass_rxn=ids.biomass_rxn,
        cofactor_met_names=tuple(ids.gam_cofactors),
        ngam_rxn=ngam_rxn,
        ngam_value=ngam,
    )


def change_amino_acid_ratio(
    model: cobra.Model,
    *,
    aerobic: bool = True,
    aa_tsv: Path | str | None = None,
) -> cobra.Model:
    """Switch the protein pseudoreaction's amino-acid ratios.

    Ports yeast-GEM's ``changeAminoAcidRatio.m``. Reads
    ``data/physiology/aminoAcid_Bjorkeroth2020.tsv`` (20 rows; columns:
    aa name, tRNA substrate id, charged-tRNA product id, MW, aerobic
    fraction, anaerobic fraction). Replaces the tRNA stoichiometries
    in the protein pseudoreaction and rescales protein back to its
    pre-switch mass via :func:`scale_biomass`.
    """
    path = Path(aa_tsv) if aa_tsv else _AA_TSV
    aa_df = _read_aa_ratio_tsv(path)
    column = "aerobic" if aerobic else "anaerobic"

    # Snapshot current protein mass so we can rescale after replacing
    # the stoichiometry.
    cfg = yeast_biomass_config()
    fractions_before = _ra_sum_biomass(model, cfg)
    protein_target = fractions_before["protein"]

    ids = load_ids()
    rxn = model.reactions.get_by_id(ids.protein_rxn)
    for _i, row in enumerate(aa_df.itertuples(index=False)):
        sub = model.metabolites.get_by_id(row.tRNA_substrate)
        prod = model.metabolites.get_by_id(row.tRNA_product)
        ratio = float(row[aa_df.columns.get_loc(column)])  # type: ignore[index]
        _set_coefficient(rxn, sub, -ratio)
        _set_coefficient(rxn, prod,  ratio)

    # Rescale protein content back to its pre-switch mass to keep the
    # biomass equation summing to 1 g/gDW.
    scale_biomass(model, "protein", protein_target)
    return model


# --- helpers ----------------------------------------------------------

def _read_aa_ratio_tsv(path: Path) -> pd.DataFrame:
    """Parse the AA-ratio TSV (columns shared with MATLAB's textscan).

    The TSV header line has the layout:
        <tab><tab><tab>MW<tab>aerobic<tab>anaerobic
    so pandas can't autodetect column names. We name them explicitly.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=0,
        names=["aa", "tRNA_substrate", "tRNA_product", "MW", "aerobic", "anaerobic"],
    )
    return df


def _rescale_named_pseudoreaction(
    model: cobra.Model,
    pseudoreaction_name: str,
    factor: float,
) -> None:
    """Rescale a pseudoreaction located by ``model.reactions[*].name``.

    Used for the yeast lipid_chain aggregation case where the
    component isn't in the BiomassConfig (lipid_chain doesn't
    contribute mass) but still needs to be rescaled in lock-step with
    lipid_backbone. Mirrors :func:`raven_toolbox.biomass.rescale_pseudoreaction`
    in shape, but identifies the rxn by name and treats every
    metabolite as a "substrate" (the matching product check is
    elided because the lipid-chain product is unique to the rxn).
    """
    cfg = yeast_biomass_config()
    proton_met = model.metabolites.get_by_id(cfg.proton_met)

    rxn = next((r for r in model.reactions if r.name == pseudoreaction_name), None)
    if rxn is None:
        return  # no-op if the pseudoreaction is absent

    # Treat "the metabolite whose name appears after 'lipid '" as the
    # product to mirror rescale_pseudoreaction's logic; for any other
    # use case the caller would go through the BiomassConfig path.
    product_name = pseudoreaction_name.removesuffix(" pseudoreaction")
    deltas = {}
    for met, coef in rxn.metabolites.items():
        if met.name == product_name:
            continue
        deltas[met] = (factor - 1.0) * coef
    if deltas:
        rxn.add_metabolites(deltas, combine=True)

    _set_coefficient(rxn, proton_met, 0.0)
    total_charge = sum((m.charge or 0) * c for m, c in rxn.metabolites.items())
    _set_coefficient(rxn, proton_met, -total_charge)


def _set_coefficient(rxn: cobra.Reaction, met: cobra.Metabolite, value: float) -> None:
    current = rxn.metabolites.get(met, 0.0)
    delta = float(value) - current
    if delta != 0:
        rxn.add_metabolites({met: delta}, combine=True)
