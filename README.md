# yeast-GEM: The consensus genome-scale metabolic model of _Saccharomyces cerevisiae_

[![DOI](https://zenodo.org/badge/52777598.svg)](https://zenodo.org/badge/latestdoi/52777598) [![GitHub version](https://badge.fury.io/gh/sysbiochalmers%2Fyeast-gem.svg)](https://badge.fury.io/gh/sysbiochalmers%2Fyeast-gem) [![Join the chat at https://gitter.im/SysBioChalmers/yeast-GEM](https://badges.gitter.im/SysBioChalmers/yeast-GEM.svg)](https://gitter.im/SysBioChalmers/yeast-GEM?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)[![Memote history](https://github.com/SysBioChalmers/yeast-GEM/workflows/Memote%20history/badge.svg)](https://sysbiochalmers.github.io/yeast-GEM/history_report.html)

# Description

This repository contains the current consensus genome-scale metabolic model of _Saccharomyces cerevisiae_. It is the continuation of the legacy project [yeastnet](https://sourceforge.net/projects/yeast/). For the latest release please [click here](https://github.com/SysBioChalmers/yeast-GEM/releases).

# Citation

* If you use yeast-GEM please cite the yeast8 paper:
  > Lu, H. et al. _A consensus S. cerevisiae metabolic model Yeast8 and its ecosystem for comprehensively probing cellular metabolism._ Nature Communications 10, 3586 (2019). [doi:10.1038/s41467-019-11581-3](https://doi.org/10.1038/s41467-019-11581-3)
* Additionally, all yeast-GEM releases are archived in [Zenodo](https://zenodo.org/badge/latestdoi/52777598), for you to cite the specific version of yeast-GEM that you used in your study, to ensure reproducibility. You should always cite the original publication + the specific version, for instance:
  > _The yeast consensus genome-scale model [Lu et al. 2019], version 8.3.4 [Sánchez et al. 2019], was used._

  Find the citation details for your specific version [here](https://zenodo.org/search?page=1&size=20&q=conceptrecid:%221494182%22&sort=-publication_date&all_versions=True).
  
# Keywords

**Utilisation:** experimental data reconstruction; multi-omics integrative analysis; _in silico_ strain design; model template  
**Field:** metabolic-network reconstruction  
**Type of model:** reconstruction; curated  
**Model source:** [YeastMetabolicNetwork](http://doi.org/10.1038/nbt1492)  
**Omic source:** genomics; metabolomics  
**Taxonomic name:** _Saccharomyces cerevisiae_  
**Taxonomy ID:** [taxonomy:559292](https://identifiers.org/taxonomy:559292)  
**Genome ID:** [insdc.gca:GCA_000146045.2](https://identifiers.org/insdc.gca:GCA_000146045.2)  
**Metabolic system:** general metabolism  
**Strain:** S288C  
**Condition:** aerobic, glucose-limited, defined media  

## Model overview

| Taxonomy | Latest update | Version | Reactions | Metabolites | Genes |
|:-------|:--------------|:------|:------|:----------|:-----|
| _Saccharomyces cerevisiae_ | 08-Jul-2023 | develop | 4131 | 2806 | 1163 |

# Installation & usage

## Obtain model

You can obtained the model by any of the following methods:
1. If you have a Git client installed on your computer, you can clone the [`main`](https://github.com/SysBioChalmers/yeast-GEM) branch of the yeast-GEM repository.
2. You can directly download [the latest release](https://github.com/SysBioChalmers/yeast-GEM/releases) as a ZIP file.
3. If you want to contribute to the development of yeast-GEM (see [below](#below)), it is best to [fork](https://github.com/SysBioChalmers/yeast-GEM/fork) the yeast-GEM repository to your own Github account.

## Required software

### Basic user

If you want to use the model for your own model simulations, you can use **any software** that accepts SBML L3V1 FBCv3 formatted model files. This includes any of the following:
* MATLAB-based
  * [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN) version 2.7.1 or later (recommended)  
  * [COBRA Toolbox](https://github.com/opencobra/cobratoolbox)

* Python-based
  * [cobrapy](https://github.com/opencobra/cobrapy)  

Please see the installation instructions for each software package.

### Developer

* MATLAB-based  
  If you want to contribute to the development of yeast-GEM, or otherwise want to run any of the [provided](https://github.com/SysBioChalmers/yeast-GEM/tree/main/code) MATLAB functions, then the following software is required:
  * [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN) version 2.7.1 or later

* Python-based  
  Contribution via python (cobrapy) is not yet functional. In essence, if you can retain the same format of the model files, you can still contribute to the development of yeast-GEM. However, you cannot use the MATLAB functions.

  If you want to use any of the [provided](https://github.com/SysBioChalmers/yeast-GEM/tree/main/code) Python functions, you may create an environment with all requirements:
  ```bash
  pip install -r code/requirements/requirements.txt  # installs all dependencies
  touch .env                                    # creates a .env file for locating the root
  ```

If you want to locally run `memote run` or `memote report history`, you should also install [git lfs](https://git-lfs.github.com/), as `results.db` (the database that stores all memote results) is tracked with git lfs.

## Model usage

Make sure to load/save the model with the corresponding wrapper functions
* In Matlab:
  ```matlab
  cd ./code
  model = loadYeastModel(); % loading
  saveYeastModel(model);    % saving
  ```
  * If RAVEN is not installed, you can also use COBRA-native functions (`readCbModel`, `writeCbModel`), but these model-files cannot be committed back to the GitHub repository.
* In Python:
  ```python
  import code.io as io
  model = io.read_yeast_model() # loading
  io.write_yeast_model(model)   # saving
  ```

### Online visualization/simulation

* You can visualize selected pathways of yeast-GEM and perform online constraint-based simulations using [Caffeine](https://caffeine.dd-decaf.eu/interactive-map), by creating a simulation with the latest yeast-GEM version available, and choosing any _S. cerevisiae_ map (currently only `iMM904` maps are available). Learn more [about Caffeine](https://caffeine.dd-decaf.eu).
* Additionally, you can interactively navigate model components and visualize 3D representations of all compartments and subsystems of yeast-GEM at [Metabolic Atlas](https://metabolicatlas.org/explore?selected=Yeast-GEM). Learn more [about Metabolic Atlas](https://www.metabolicatlas.org/about).

# Contributing

Contributions are always welcome! Please read the [contributions guideline](https://github.com/SysBioChalmers/yeast-GEM/blob/main/.github/CONTRIBUTING.md) to get started.

## Contributors

Code contributors are reported automatically by GitHub under [Contributors](https://github.com/SysBioChalmers/yeast-GEM/graphs/contributors), while other contributions come in as [Issues](https://github.com/SysBioChalmers/yeast-GEM/issues).
