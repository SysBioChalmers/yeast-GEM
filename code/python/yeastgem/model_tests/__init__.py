"""yeast-GEM model-validation tests — Python counterparts of `code/modelTests/`.

These are *integration* tests against experimental data (the
``increaseVersion`` PR gate, plus standalone benchmarks):

* :mod:`growth` — chemostat growth across 4 conditions vs Tobias 2013.
* :mod:`essential_genes` — single-gene-knockout vs Stanford KO collection.
* :mod:`anaerobic_flux` — anaerobic intracellular flux vs Jouhten 2008
  + Frick & Wittmann 2005.
* :mod:`plot_anaerobic` — relative fermentation product bar plot.
* :mod:`find_duplicated_rxns` — print duplicate stoichiometry pairs.

They reuse the yeast-specific bits from :mod:`yeastgem.conditions` and
:mod:`yeastgem.biomass`, and delegate generic FBA / deletion / duplicate
detection to cobrapy / raven-python.
"""
from yeastgem.model_tests.anaerobic_flux import anaerobic_flux_predictions
from yeastgem.model_tests.essential_genes import EssentialGeneResult, essential_genes
from yeastgem.model_tests.find_duplicated_rxns import find_duplicated_rxns
from yeastgem.model_tests.growth import growth
from yeastgem.model_tests.plot_anaerobic import plot_anaerobic

__all__ = [
    "EssentialGeneResult",
    "anaerobic_flux_predictions",
    "essential_genes",
    "find_duplicated_rxns",
    "growth",
    "plot_anaerobic",
]
