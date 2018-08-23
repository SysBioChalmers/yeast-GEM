# History

### yeast 8.2.0:
* Features:
  * Fixes #38: Added 183 new reactions, 277 new metabolites and 163 new genes based on the latest genome annotation in SGD, uniprot, KEGG, Biocyc & Reactome (PR #142).
* Fixes:
  * `grRules` deleted from pseudoreactions, removing with this 49 genes (PR #145).
* Chores:
  * Updated COBRA, which changed the number of decimals in some stoichiometric coefficients in `.txt` (PR #143)

### yeast 8.1.3:
* Features:
  * Added SBO terms for all metabolites and reactions, based on an automatic script now part of `saveYeastModel.m` (PR #132).
  * `increaseVersion.m` now avoids conflicts between `devel` and `master` by erroring before releasing and guiding the admin to change first `devel` (PR #133).
  * Website now available in `gh-pages` branch: http://sysbiochalmers.github.io/yeast-GEM/
* Fixes:
  * Standardize naming of pseudo-metabolites "lipid backbone" & "lipid chain" (PR #130).
* Chores:
  * Updated COBRA, which swapped around the order of the `bqbiol:is` and `bqbiol:isDescribedBy` qualifiers in the `.xml` file (PR #131).

### yeast 8.1.2:
* New features:
  * `saveYeastModel.m` now checks if the model is a valid SBML structure; if it isn't it will error (PR #126).
  * Date + model size in `README.md` updates automatically when saving the model (PR #123).
  * Added `modelName` and `modelID`; the latter which will now store the version number (PR #127).
* Fixes:
  * Fixes #60: New GPR relations for existing reactions were added according to new annotation from 5 different databases (PR #124).
  * Various fixes in `README.md` (PR #123).

### yeast 8.1.1:
* Fixes:
  * Fixes #96: regardless if the model is saved with a windows or a MAC machine, the `.xml` file is now stored with the same scientific format.
  * Fixes #108: No CHEBI or KEGG ids are now shared by different metabolites. Also, updated the metabolites that were skipped in the previous manual curation (PR #74).
  * Remade function for defining confidence scores, which fixed 38 scores in `rxnConfidenceScores` (most of them from pseudoreactions).
  * `loadYeastModel` and `saveYeastModel` were improved to allow their use also when outside of the actual folder.

### yeast 8.1.0:
* New features:
  * SLIME reactions added to the model using [SLIMEr](https://github.com/SysBioChalmers/SLIMEr), to properly account for constraints on lipid metabolism (fixes #21):
    * SLIME rxns replace old ISA rxns for lumping lipids. They create 2 types of lipid pseudometabolites: backbones and acyl chains.
    * There are now 3 lipid pseudoreactions: 1 constrains backbones, 1 constrains acyl chains, 1 merges both.
* Fixes:
  * All metabolite formulas made compliant with SBML (fixes #19). Model is now a valid SBML object.
  * Biomass composition was rescaled to experimental data from [Lahtvee et al. 2017](https://www.sciencedirect.com/science/article/pii/S2405471217300881), including protein and RNA content, trehalose and glycogen concentrations, lipid profile and FAME data. Biomass was fitted to add up to 1 g/gDW by rescaling total carbohydrate content (unmeasured).
* Refactoring:
  * Organized all files in `ComplementaryData`

### yeast 8.0.2:
* New features:
  * Model can now be used with cobrapy by running `loadYeastModel.py`
  * `loadYeastModel.m` now adds the `rxnGeneMat` field to the model
* Refactoring:
  * Moved `pmids` of model from `rxnNotes` to `rxnReferences` (COBRA-compliant)
  * `yeastGEM.yml` and `dependencies.txt` are now updated by RAVEN (a few dependencies added)
  * Moved `boundaryMets.txt` and `dependencies.txt` to the `ModelFiles` folder
* Documentation:
  * Added badges and adapted README ro reflect new features

### yeast 8.0.1:
* `.yml` format included for easier visualization of model changes
* Empty notes removed from model
* Issue and PR templates included
* `README.md` updated to comply with new repo's name

### yeast 8.0.0:
First version of the yeast8 model, to separate it from previous versions:

* Manual curation project:
  * All metabolite information manually curated (names, charges, kegg IDs, chebi IDs)
  * Reaction gene rules updated with curation from [the iSce926 model](http://www.maranasgroup.com/submission_models/iSce926.htm). 13 genes added in this process
* Format changes:
  * Folder `ComplementaryData` introduced
  * All data is stored in `.tsv` format now (can be navigated in Github)
  * Releases now come in `.xlsx` as well
* Other new features:
  * Added `loadYeastModel.m`
  * A much smarter `increaseVersion.m`
  * Lots of refactoring

### yeast 7.8.3:
* curated tRNA's formulas
* started tracking COBRA and RAVEN versions
* dropped SBML toolbox as requirement
* reorganized `complementaryScripts`
* switched to a CC-BY-4.0 license

### yeast 7.8.2:
* fixed subSystems bug: now they are saved as individual groups
* solved inter-OS issues
* remade license to follow GitHub format
* added `history.md` and made it a requirement to update when increasing version

### yeast 7.8.1:
* started following dependencies
* started keeping track of the version in the repo (`version.txt`)
* included `.gitignore`
* dropped `.mat` storage for `devel` + feature branches (but kept it in `master`)

### yeast 7.8.0:
* Added information:
  * `metFormulas` added for all lipids
  * `rxnKEGGID` added from old version
  * `rxnNotes` enriched with Pubmed ids (`pmid`) from old version
  * `rxnConfidenceScores` added based on  automatic script (available in [`ComplementaryScripts`](https://github.com/SysBioChalmers/yeast-GEM/blob/master/ComplementaryScripts))
* Format changes:
  * Biomass clustered by 5 main groups: protein, carbohydrate, lipid, RNA and DNA

### yeast 7.7.0:
* Format changes:
  * FBCv2 compliant
  * Compatible with latest COBRA and RAVEN parsers
  * Created main structure of repository
* Added information:
  * `geneNames` added to genes based on [KEGG](http://www.genome.jp/kegg/) data
  * `subSystems` and `rxnECnumbers` added to reactions based on [KEGG](http://www.genome.jp/kegg/) & [Swissprot](http://www.uniprot.org/uniprot/?query=*&fil=organism%3A%22Saccharomyces+cerevisiae+%28strain+ATCC+204508+%2F+S288c%29+%28Baker%27s+yeast%29+%5B559292%5D%22+AND+reviewed%3Ayes) data
  * Boundary metabolites tracked (available in [`ComplementaryScripts`](https://github.com/SysBioChalmers/yeast-GEM/blob/master/ComplementaryScripts))
* Simulation improvements:
  * Glucan composition fixed in biomass pseudo-rxn
  * Proton balance in membrane restored
  * Ox.Pho. stoichiometry fixed
  * NGAM rxn introduced
  * GAM in biomass pseudo-rxn fixed and refitted to chemostat data
  
  ### yeast 7.6.0:
First release of the yeast model in GitHub, identical to the last model available at https://sourceforge.net/projects/yeast/
