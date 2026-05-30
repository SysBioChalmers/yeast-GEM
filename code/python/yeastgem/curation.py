"""Yeast-GEM batch curation entry point.

Thin wrapper over :func:`raven_python.curation.batch_curate` (and its
``from_tsv`` companion) that pins the yeast-GEM id prefixes
(``'s_'`` for new metabolites, ``'r_'`` for new reactions). All
schema details (column names, MIRIAM-auto-detection, match keys) live
upstream — see ``raven_python.curation`` for the full reference.

The MATLAB counterpart is ``code/modelCuration/curateMetsRxnsGenes.m``
(also a shim over RAVEN's ``curateModelFromTables``).
"""
from __future__ import annotations

from pathlib import Path

import cobra
import pandas as pd
from raven_python.curation import (
    CurationResult,
)
from raven_python.curation import (
    batch_curate as _ra_batch_curate,
)
from raven_python.curation import (
    batch_curate_from_tsv as _ra_batch_curate_from_tsv,
)

# Yeast-GEM id prefixes — frozen for both Python and MATLAB callers.
_MET_ID_PREFIX = "s_"
_RXN_ID_PREFIX = "r_"


def curate_mets_rxns_genes(
    model: cobra.Model,
    *,
    mets_df: pd.DataFrame | None = None,
    genes_df: pd.DataFrame | None = None,
    rxns_df: pd.DataFrame | None = None,
    rxns_coeffs_df: pd.DataFrame | None = None,
) -> CurationResult:
    """Add or update metabolites / reactions / genes from DataFrames.

    Yeast-GEM-specific id prefixes are applied automatically (``s_`` /
    ``r_``); everything else is delegated to
    :func:`raven_python.curation.batch_curate`. See its docstring for
    the schema, match-key rules and the MIRIAM-auto-detection
    convention.
    """
    return _ra_batch_curate(
        model,
        mets_df=mets_df,
        genes_df=genes_df,
        rxns_df=rxns_df,
        rxns_coeffs_df=rxns_coeffs_df,
        met_id_prefix=_MET_ID_PREFIX,
        rxn_id_prefix=_RXN_ID_PREFIX,
    )


def curate_mets_rxns_genes_from_tsv(
    model: cobra.Model,
    *,
    mets_tsv: str | Path | None = None,
    genes_tsv: str | Path | None = None,
    rxns_tsv: str | Path | None = None,
    rxns_coeffs_tsv: str | Path | None = None,
) -> CurationResult:
    """File-path convenience wrapper — same shape as the MATLAB
    ``curateMetsRxnsGenes(model, metsInfo, genesInfo, rxnsCoeffs,
    rxnsInfo)``."""
    return _ra_batch_curate_from_tsv(
        model,
        mets_tsv=mets_tsv,
        genes_tsv=genes_tsv,
        rxns_tsv=rxns_tsv,
        rxns_coeffs_tsv=rxns_coeffs_tsv,
        met_id_prefix=_MET_ID_PREFIX,
        rxn_id_prefix=_RXN_ID_PREFIX,
    )
