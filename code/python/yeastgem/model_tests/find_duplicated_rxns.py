"""Print duplicate-stoichiometry reaction pairs.

Port of ``code/modelTests/findDuplicatedRxns.m``. Identifies reactions
with identical (or reversed) stoichiometry and prints their names,
GPRs and bounds — the same shape of output the legacy MATLAB function
produced. Detection itself is delegated to
:func:`raven_toolbox.manipulation.find_duplicate_reactions`.
"""
from __future__ import annotations

import cobra
from raven_toolbox.manipulation import find_duplicate_reactions


def find_duplicated_rxns(model: cobra.Model) -> list[list[cobra.Reaction]]:
    """Print duplicate-stoichiometry reaction groups and return them.

    Each group's reactions are listed with ``name``, ``gene_reaction_rule``,
    ``lower_bound`` and ``upper_bound`` — the legacy MATLAB output
    format. yeast-GEM's convention treats A→B and B→A as duplicates
    (``ignore_direction=True``).
    """
    groups = find_duplicate_reactions(model, ignore_direction=True)
    for group in groups:
        for rxn in group:
            print(
                f"Name: {rxn.name} - GPR: {rxn.gene_reaction_rule} - "
                f"LB={rxn.lower_bound} - UB={rxn.upper_bound}"
            )
        print()  # blank line between groups
    return groups
