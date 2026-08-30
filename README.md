# yeast-GEM: The consensus genome-scale metabolic model of _Saccharomyces cerevisiae_

[![Current release](https://img.shields.io/github/release/SysBioChalmers/yeast-GEM/all.svg)](https://github.com/SysBioChalmers/yeast-GEM/releases/)
[![GitHub Discussions](https://img.shields.io/github/discussions/sysbiochalmers/yeast-gem)](https://github.com/SysBioChalmers/yeast-GEM/discussions)
[![Memote report](https://github.com/SysBioChalmers/yeast-GEM/actions/workflows/memote-release.yml/badge.svg)](https://sysbiochalmers.github.io/yeast-GEM/release_report.html)
[![DOI](https://zenodo.org/badge/52777598.svg)](https://zenodo.org/badge/latestdoi/52777598)

# Description

This repository contains the current consensus genome-scale metabolic model of _Saccharomyces cerevisiae_. It is the continuation of the legacy project [yeastnet](https://sourceforge.net/projects/yeast/).

# Citation

* If you use yeast-GEM please cite the yeast9 paper:
  > Zhang, C. et al. _Yeast9: a consensus genome-scale metabolic model for S. cerevisiae curated by the community._ Molecular Systems Biology (2024) doi:[10.1038/s44320-024-00060-7](https://doi.org/10.1038/s44320-024-00060-7)
* For pre-yeast9 versions:
  > Lu, H. et al. _A consensus S. cerevisiae metabolic model Yeast8 and its ecosystem for comprehensively probing cellular metabolism._ Nature Communications 10, 3586 (2019). doi:[10.1038/s41467-019-11581-3](https://doi.org/10.1038/s41467-019-11581-3)
* Additionally, every yeast-GEM release is archived in [Zenodo](https://zenodo.org/badge/latestdoi/52777598). To ensure reproducibility, always cite both the original publication and the specific version you used, for instance:
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

# Model overview

Model statistics (reaction, metabolite and gene counts) and validation
results against experimental data (gene-essentiality and growth-rate
prediction accuracy) are stamped here at release time, since `develop`
itself is a moving target that would make those numbers stale the moment
the next curation lands. See the
[latest release](https://github.com/SysBioChalmers/yeast-GEM/releases/latest)
for those numbers, or
[data/testResults/README.md](data/testResults/README.md) for this
branch's current validation results.

# Installation & usage

## Obtain model

You can obtain the model by any of the following methods:
1. If you have a Git client installed on your computer, you can clone the [`main`](https://github.com/SysBioChalmers/yeast-GEM) branch of the yeast-GEM repository.
2. You can directly download [the latest release](https://github.com/SysBioChalmers/yeast-GEM/releases) as a ZIP file.
3. If you want to contribute to the development of yeast-GEM (see [Contributing](#contributing)), it is best to [fork](https://github.com/SysBioChalmers/yeast-GEM/fork) the yeast-GEM repository to your own Github account.

## Required software

### Basic user

If you want to use the model for your own model simulations, you can use **any software** that accepts SBML L3V1 FBCv3 formatted model files. This includes any of the following:
* MATLAB-based
  * [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN) version 2.8.3 or later (recommended) — except for loading `model/yeast-GEM.yml` itself, which `loadYeastModel` does by default; see the developer note below  
  * [COBRA Toolbox](https://github.com/opencobra/cobratoolbox)

* Python-based
  * [cobrapy](https://github.com/opencobra/cobrapy)  

Please see the installation instructions for each software package.

### Developer

* MATLAB-based  
  If you want to contribute to the development of yeast-GEM, or otherwise want to run any of the [provided](https://github.com/SysBioChalmers/yeast-GEM/tree/main/code) MATLAB functions, then the following software is required:
  * [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN) version 2.8.3 or later, **except** for reading or writing `model/yeast-GEM.yml` (what `loadYeastModel`/`saveYeastModel` do by default): that needs RAVEN's `develop3` branch, until its changes reach a numbered release

* Python-based  
  The `yeastgem` package (`code/python/`) is the actively maintained Python
  interface to the model: a cobrapy-based API built on
  [raven-toolbox](https://github.com/SysBioChalmers/raven-python) for the generic,
  organism-agnostic parts. It covers loading and committing the model, applying
  growth conditions, biomass composition, and the validation test suite. See
  [code/python/README.md](code/python/README.md) for the API map.
  ```bash
  pip install -e code/python/
  touch .env # create a .env file for locating the root
  ```

  The standalone functions in
  [code/](https://github.com/SysBioChalmers/yeast-GEM/tree/main/code) (`code/io.py`
  and its neighbours) predate `yeastgem` and are being phased out; prefer
  `yeastgem` above for anything new. If you still need them:
  ```bash
  pip install -r code/requirements/requirements.txt  # install all dependencies
  touch .env # create a .env file for locating the root
  ```

MEMOTE scores are reported automatically — a fast score on every pull request,
and the full suite at release, published [here](https://sysbiochalmers.github.io/yeast-GEM/release_report.html).
To run MEMOTE locally instead, see [`memote_snapshot.py`](code/python/tests/ci/memote_snapshot.py) or run `memote` directly.

## Model usage

Make sure to load/save the model with the corresponding wrapper functions:
* In Matlab:
  ```matlab
  cd ./code
  model = loadYeastModel(); % loading
  saveYeastModel(model);    % saving
  ```
  * If RAVEN is not installed, you can also use COBRA-native functions (`readCbModel`, `writeCbModel`), but these model-files cannot be committed back to the GitHub repository.
* In Python, with the [`yeastgem`](code/python/README.md) package:  
Before opening Python, the following command only needs to be run once, in the yeast-GEM root folder:  
  ```bash
  touch .env # create a .env file for locating the root
  ```
  Afterwards, the model can be loaded in Python with:
  ```python
  from yeastgem import read_yeast_model, commit_yeast_model
  model = read_yeast_model()      # loading
  commit_yeast_model(model)       # saving
  ```

# Contributing

Contributions are always welcome! Please read the [contributions guideline](https://github.com/SysBioChalmers/yeast-GEM/blob/main/.github/CONTRIBUTING.md) to get started.

## Contributors

Code contributors are reported automatically by GitHub under [Contributors](https://github.com/SysBioChalmers/yeast-GEM/graphs/contributors), while other contributions come in as [Issues](https://github.com/SysBioChalmers/yeast-GEM/issues).
