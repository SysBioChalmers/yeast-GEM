# yeast-GEM: The consensus genome-scale metabolic model of _Saccharomyces cerevisiae_

[![GitHub version](https://badge.fury.io/gh/sysbiochalmers%2Fyeast-gem.svg)](https://badge.fury.io/gh/sysbiochalmers%2Fyeast-gem) [![Join the chat at https://gitter.im/SysBioChalmers/yeast-GEM](https://badges.gitter.im/SysBioChalmers/yeast-GEM.svg)](https://gitter.im/SysBioChalmers/yeast-GEM?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

* Brief Model Description:

This repository contains the current consensus genome-scale metabolic model of _Saccharomyces cerevisiae_. It is the continuation of the legacy project [yeastnet](http://yeast.sourceforge.net/). For the latest release please [click here](https://github.com/SysBioChalmers/yeast-GEM/releases).

* Model KeyWords:

**GEM Category:** Species; **Utilisation:** predictive simulation, multi-omics integrative analysis, _in silico_ strain design, model template; **Field:** metabolic-network reconstruction; **Type of Model:** curated, reconstruction; **Model Source:** [Yeast 7.6](https://sourceforge.net/projects/yeast/); **Taxonomy:** _Saccharomyces cerevisiae_; **Metabolic System:** General Metabolism; **Condition:** aerobic, glucose-limited, defined media, maximization of growth.

* Last update: 2018-05-17

* Main Model Descriptors:

|Taxonomy | Template Model | Reactions | Metabolites| Genes |
|:-------:|:--------------:|:---------:|:----------:|:-----:|
|_Saccharomyces cerevisiae_|[Yeast 7.6](https://sourceforge.net/projects/yeast/)|3496|2224|922|

This repository is administered by Benjamín J. Sánchez ([@BenjaSanchez](https://github.com/benjasanchez)), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.

## Installation

### Required Software - User:

* Matlab user:
  * A functional Matlab installation (MATLAB 7.3 or higher).
  * The [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox).
* Python user:
  * Python 2.7, 3.4, 3.5 or 3.6
  * [cobrapy](https://github.com/opencobra/cobrapy)

### Required Software - Contributor:

* Both of the previous Matlab requirements.
* The [RAVEN toolbox for MATLAB](https://github.com/SysBioChalmers/RAVEN).
* A [git wrapper](https://github.com/manur/MATLAB-git) added to the search path.

### Dependencies - Recommended Software:
* For Matlab, libSBML MATLAB API (version [5.16.0](https://sourceforge.net/projects/sbml/files/libsbml/5.16.0/stable/MATLAB%20interface/) is recommended).
* [Gurobi Optimizer](http://www.gurobi.com/registration/download-reg) for any simulations.

### Installation Instructions
* For users: Clone it from [`master`](https://github.com/SysBioChalmers/yeast-GEM) in the Github repo, or just download [the latest release](https://github.com/SysBioChalmers/yeast-GEM/releases).
* For contributors: Fork it to your Github account, and create a new branch from [`devel`](https://github.com/SysBioChalmers/yeast-GEM/tree/devel).

## Usage

Make sure to load/save the model with the corresponding wrapper functions!
* In Matlab:
  * Loading: `complementaryScripts/loadYeastModel.m`
  * Saving: `complementaryScripts/saveYeastModel.m`
* In Python:
  * Loading: `complementaryScripts/loadYeastModel.py`
  * Saving: currently unavailable


## Complementary Scripts

* `missingFields`: Folder with functions for adding missing fields to the model.
   * `addGeneNames.m`: Adds the field `geneNames` by extracting the data from KEGG.
   * `getConfidenceScores.m`: Assigns confidence scores based in a basic automatic criteria.
   * `getMissingFields.m`: Retrieves missing information (`rxnECNumbers` and `subSystems`) from KEGG & Swissprot. It uses `changeRules.m` for properly reading the gene-reaction rules, and `findInDB.m`, `getAllPath.m`, `findSubSystem.m` & `deleteRepeated.m` for reading the databases. The latter 4 functions are adapted versions of functions from the [GECKO toolbox](https://github.com/SysBioChalmers/GECKO).
* `modelCuration`: Folder with curation functions.
   * `addiSce926changes.m`: Updates the model to include curation from [the iSce926 model](http://www.maranasgroup.com/submission_models/iSce926.htm).
   * `calculateContent.m`: Calculates the protein and carb fraction in the biomass pseudo-rxn.
   * `changeBiomass.m`: Rescales the biomass composition for varying protein content in anaerobic case. Also changes GAM and NGAM.
   * `checkMetBalance.m`: Shows rxns that consume/produce a given metabolite in the model.
   * `makeFormulasCompliant.m`: Fixes the compliance problems of all metabolite formulas.
   * `modelCorrections.m`: Corrects various issues in yeast7 (biomass composition, proton balance, Ox.Pho., GAM and NGAM).
   * `takeOutFromFormula.m`: Takes away from formula each of the elements specified.
   * `updateMetaboliteAnnotation.m`: Reads `metabolite_manual_curation.tsv` and updates the model with it.
* `otherChanges`: Folder with other types of changes.
   * `anaerobicModel.m`: Transforms the model to anaerobic conditions.
   * `clusterBiomass.m`: Separates the biomass in 5 main components: protein, carbohydrate, lipid, RNA and DNA.
   * `convertYmn2FBC2.m`: Converts yeast7 from COBRA-compatible SBML2 to FBC v2, thereby adding the missing annotation data, which could not be retained with the older COBRA versions.
   * `getNewIndex.m`: Finds the highest index available in either metabolites or rxns, and then adds one to it, for creating any new species.
* `dependencies.txt`: Tracks SBML versions and levels used for saving the model.
* `boundaryMets.txt`: Contains a list of all boundary metabolites in model, listing the id and name.
* `increaseVersion.m`: Updates the version of the model in `version.txt` and as metaid in the `.xml` file. Saves the model as `.mat` and as `.xlsx`
* `loadYeastModel.m`: Loads the yeast model from the `.xml` file for Matlab.
* `loadYeastModel.py`: Loads the yeast model from the `.xml` file for Python.
* `saveYeastModel.m`: Saves yeast model as a `.xml`, `.yml` and `.txt` file, and updates `boundaryMets.txt` and `dependencies.txt`.

## Complementary Data

* `iSce926curatedGeneRules.tsv`: Manually curated gene rules added to the model based on [the iSce926 model](http://www.maranasgroup.com/submission_models/iSce926.htm).
* `iSce926newGenes.tsv`: New genes added to the model based on [the iSce926 model](http://www.maranasgroup.com/submission_models/iSce926.htm).
* `metabolite_manual_curation.tsv`: All manually curated data added to metabolites.
* `SGDgeneNames.tsv`: Short gene names for each gene ID.
* `kegg.tsv`: KEGG data for _S. cerevisiae_.
* `swissprot.tsv`: SWISSPROT data for _S. cerevisiae_.

## Contributors

* [Eduard J. Kerkhoven](https://www.chalmers.se/en/staff/Pages/Eduard-Kerkhoven.aspx) ([@edkerk](https://github.com/edkerk)), Chalmers University of Technology, Gothenburg Sweden
* [Feiran Li](https://www.chalmers.se/en/staff/Pages/feiranl.aspx) ([@feiranl](https://github.com/feiranl)), Chalmers University of Technology, Gothenburg Sweden
* [Hongzhong Lu](https://www.chalmers.se/en/Staff/Pages/luho.aspx) ([@hongzhonglu](https://github.com/hongzhonglu)), Chalmers University of Technology, Gothenburg Sweden
* [Simonas Marcišauskas](https://www.chalmers.se/en/Staff/Pages/simmarc.aspx) ([@simas232](https://github.com/simas232)), Chalmers University of Technology, Gothenburg Sweden
* [Benjamín J. Sánchez](https://www.chalmers.se/en/staff/Pages/bensan.aspx) ([@BenjaSanchez](https://github.com/benjasanchez)), Chalmers University of Technology, Gothenburg Sweden
