"""Load canonical yeast-GEM identifiers from data/yeastgem/ids.yml."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from yeastgem.io import REPO_PATH

_IDS_PATH = REPO_PATH / "data" / "yeastgem" / "ids.yml"


@dataclass(frozen=True)
class BiomassComponentConfig:
    """One entry under ``biomass_components`` in ids.yml."""

    name: str
    mass_strategy: str  # see raven_python.biomass.config.MassStrategy


@dataclass(frozen=True)
class YeastIDs:
    """Canonical yeast-GEM identifiers consumed by generic algorithms.

    Mirrors the MATLAB `applyIDs()` struct exactly: the YAML file is the
    single source of truth for both languages.
    """

    biomass_rxn: str
    protein_rxn: str
    cofactor_rxn: str
    proton_met: str
    pseudoreaction_names: dict[str, str]
    gam_cofactors: list[str]
    biomass_components: tuple[BiomassComponentConfig, ...]


def load_ids(path: Path | str | None = None) -> YeastIDs:
    """Load and return the canonical yeast IDs from ids.yml."""
    path = Path(path) if path else _IDS_PATH
    with open(path) as f:
        data = yaml.safe_load(f)
    components = tuple(
        BiomassComponentConfig(name=c["name"], mass_strategy=c["mass_strategy"])
        for c in data.get("biomass_components", [])
    )
    return YeastIDs(
        biomass_rxn=data["biomass_rxn"],
        protein_rxn=data["protein_rxn"],
        cofactor_rxn=data["cofactor_rxn"],
        proton_met=data["proton_met"],
        pseudoreaction_names=dict(data["pseudoreaction_names"]),
        gam_cofactors=list(data["gam_cofactors"]),
        biomass_components=components,
    )
