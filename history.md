# History

### yeast 9.1.0:
- Fixes:
  - Curations related to reactions involved in anaerobic metabolism, as well as other errorenous reactions that we encountered during the curation process. (PR #386)
  - Updated anaerobicModel.m and new anaerobic_model.py functions. (PR #386)
  - Mitochondrial carnitine shuttle (`r_0252`) is disabled (LB=UB=0) by default, to represent that defined medium does not (typically) contain carnitine. The UB should be set to 1000 when simulating carnitine in the growtn medium. (PR #386)
  - Glutamate synthase (`r_0472`) is disabled (LB=UB=0) by default, as the expression of this gene is severely repressed in excess nitrogen. The UB should be set to 1000 when simulating nitrogen-limiting conditions. (PR #386)
  - Glycine cleavage ('r_0501', 'r_0507', 'r_0509') is disabled (LB=UB=0) by default, as the expression is repressed when nitrogen is in excess The UB should be set to 1000 when glycine is used as nitrogen source, or nitrogen limiting conditions are simulated. (PR #386)
  - Various bugs that gave misleading growth predictions under nitrogen-limiting anaerobic conditions (closes #376)

### yeast 9.0.2:
- Fixes:
  - Curation based on auxotrophy predictions. (closes #371)
  - Correct UniProt ID for YCR024C. (closes #372)

### yeast 9.0.1:
- Features:
  - Provide reference metabolite concentration from YMDB, which can be loaded with `code/missingFields/addYMDBconcentrations.m`. (PR #364)
- Documentation:
  - Update citation and `README.md`. (PR #368)
- Chore:
  - Use latest RAVEN (2.9.3) for import/export (resolved bug https://github.com/SysBioChalmers/RAVEN/issues/529 results in higher accuracy in gene essentiality prediction).

### yeast 9.0.0:
- Fixes:
  - Correct unbalanced reactions `r_4196` (closes #358), `r_4240`, `r_4722`. (PR #360)
  - Correct metFormula of various dolichol-related metabolites. (PR #360)
  - `increaseVersion` uses YAML-model as input, to ensure it contains the deltaG annotations. (PR #360)

### yeast 8.7.1:
- Fix:
  - Failing GitHub Action workflows. (PR #350)
  - Correct 'enxyme' in r_1024 reaction name. (PR #355)
  - Make `.standard-GEM.md` hidden. (PR #357)
- Features:
  - UniProt annotations for all genes. (PR #349)

### yeast 8.7.0:
- Features:
  - Add reactions based on KEGG and MetaCyc annotations. (PR #304)
  - Add `metDeltaG` and `rxnDeltaG` fields. (PR #330)
  - Add `metSmiles` field. (PR #344)
- Fixes:
  - Fructose transporter `r_1134` as facilitated diffusion, not symport. (closes #332)
  - Remove genes from spontaneous half-reactions (`r_0029` and `r_0030`) of 3-isopropylmalate dehydrogenase. (closes #331)

### yeast 8.6.3:
- Features:
  - New reactions and metabolites related to volatile esters & polyphosphates. (PR #337)

### yeast 8.6.2:
- Fixes:
  - Correct ATP synthase mitochondrial complex gene associations. (PR #323)

### yeast 8.6.1:
- Fixes:
  - Manual curation of gene associations of transport reactions, based on various databases. (PR #306, closes #160)
  - Correct annotation of gene associations of enzyme complex, based on Complex Portal, Uniprot and SGD. (PR #305)
  - Curate 19 new GPR and consolidate curations between model releases. (PR #313)
- Features:
  - Assignment of single `subSystem` per reactions. (PR #307, closes #11) 
- Refactor:
  - Reduce software dependencies of `modelTests`. (PR #305, closes #309)

### yeast 8.6.0:
- Fixes:
  - Closes #265: Make `r_0446` (formate-tetrahydrofolate ligase) irreversible, to prevent non-zero flux through TCA cycle. (PR #290)
- Features:
  - Add pathways responsible for the formation of hydrogen sulfide as well as other volatile sulfur compounds during fermentation (PR #300)
  - Closes #302: Simplify model curation with `curateRxnsGenesMets` function (PR #300)
  - Remove COBRA Toolbox and MATLAB-git dependencies for the MATLAB-based curation pipeline (PR #303)
  - Closes #308: Distribute `yeast-GEM.mat` in RAVEN's format, to include `grRules` and `metComps` fields (PR #301)
- Refactor:
  - Change format of `yeast-GEM.txt` file to include metabolite names and compartments, instead of metabolite identifiers, to simplify `diff`-ing (metabolite identifiers are already trackable in the `yml`-file) (PR #312)

### yeast 8.5.0:
- Features:
  - Set up memote as GitHub Action for pull requests (PR #162)
  - Moved old subSystems to reaction annotations (in `rxnMiriams` or `rxnKEGGpathways`) (PR #253)
- Fixes:
  - Combine glycolysis + gluconeogenesis as single KEGG pathway annotations (PR #251).
  - Closes #252: Correct grRule of r_4590 (PR #255).
  - Closes #254: Corrects name of s_1218 (PR #255).
  - `saveYeastModel.m` now correctly handles Unicode characters (PR #255).
  - Closes #238: Correct indentation of `yeastGEM.yml` (PR #236 and #255).
- Chore:
  - Update dependencies in `*requirements.txt` (PR #256).
  - Minor changes in model file formatting due to updates in COBRA+RAVEN toolboxes (PR #253).
- Refactor:
  - Closes #232: Follow `standard-GEM` specifications (PR #257).
  - Closes #258: Rename git branches `master` and `devel` to `main` and `develop` (PR #261).

### yeast 8.4.2:
* Features:
  * `saveYeastModel.m` now checks if the model can grow and, based on the `allowNoGrowth` flag, returns either warnings or errors if not (PR #244).
  * Added several fatty acid ester producing reactions to the model, for improved simulation of alcoholic fermentation conditions (PRs #190 and #248).
* Fixes:
  * Closes #242: Fixed a bug that prevented the model from growing (PR #243).
  * Corrected directionality of 23 reactions and removed a generic reaction (PR #228).

### yeast 8.4.1:
* Features:
  * Switched to `pip-tools` for managing python dependencies, distinguishing between user requirements `/requirements/requirements.txt` and developer requirements `/requirements/dev-requirements.txt` (PR #235).
* Fixes:
  * Closes #201: Changed generic protein name to avoid confusion (PR #237).
  * Closes #205: Finished correcting reactions' stoichiometry based on KEGG data (PR #237).
  * Closes #215: Corrected wrong gene rule in reaction (PR #237).
  * Closes #225: Moved MNX rxn ids from notes to the proper annotation field (PR #226).
* Documentation:
  * Closes #223: Clarified releasing steps, including authorship criteria for Zenodo releases, in contributing guidelines (PR #233).
  * Closes #227: Removed authorships/dates from all scripts, as it is redundant information (PR #230).
  * Added admin guidelines for managing python dependencies (PR #235).
  * Included links for model visualization in README file (PR #240).

### yeast 8.4.0:
* Features:
  * New functions `mapKEGGID.m ` and `mapMNXMID.m` for adding ids in model. Used them to add missing KEGG and MetaNetX ids for both metabolites and reactions (PR #220).
  * Solves #197: Added missing MetaNetX ids using KEGG ids and ChEBI ids (PR #220).
  * Added BiGG ids for all matched metabolites/reactions using MetaNetX + manual curation, together with lists containing new BiGG ids for the unmatched ones (PR #188).
  * New functions `read_yeast_model` and `write_yeast_model` for easier usage in python (PR #224).
  * Solves #172: Model can now be loaded with BiGG ids as main ids, for better compliance with cobrapy (PR #224).
* Fixes:
  * Solves #102: Every component of the model is now preserved when the model is opened with cobrapy, including gene names (PR #216).
  * Manual curation of MetaNetX, KEGG and ChEBI ids for metabolites/reactions (PRs #188 and #220).
  * Solves #187: Removed some duplicate reactions in the model (PR #188).
  * Mass/charge balanced most unbalanced reactions in model using `checkSmatrixMNX.m`, bringing the number down to 17 reactions (PR #222).
* Others:
  * Configured repo to ensure that files always use `LF` as EOL character (PR #221).
  * Gene SBO terms are now recorded, after update in COBRA toolbox (PR #188).

### yeast 8.3.5:
* Fixes:
  * Closes #129: Removed non-S288C genes (PR #211).
  * Closes #198: Fixes function for converting model to anaerobic (PR #199).
* Tests:
  * Added growth tests for carbon & nitrogen limitation (#199).
  * Added test for computing gene essentiality (PR #200).
* Documentation/Others:
  * Clarified with README's the purpose of each script/data folder (#209).
  * Closes #206: Updated citation guidelines (PR #210).
  * Updated contribution guidelines + issue/PR templates (PR #210).
  * Created folders with deprecated scripts (PR #209).

### yeast 8.3.4:
* Features:
  * Fixes #171: Added 101 GPR rules to transport rxns according to TCDB database (PR #178).
  * Added 18 met formulas from manual search (PR #155).
  * Performed gap-filling for connecting 29 dead-end mets, by adding 28 transport rxns (PR #185). Added documentation to the gap-filling functions in PR #195.
* Fixes:
  * Corrected typo in gene ID (PR #186).

### yeast 8.3.3:
* Features:
  * Fixes #107: Two new pseudoreactions (`cofactor pseudoreaction` and `ion pseudoreaction`) added to the model as extra requirements to the biomass pseudoreaction (PRs #174 & #183).
* Fixes:
  * `addSBOterms.m` adapted to identify new pseudoreactions (PR #180).
  * Removed non-compliant symbol from a reaction name to avoid parsing errors (PR #179).
* Documentation:
  * Model keywords modified to comply with the sysbio rulebook (PR #173).
  * Added citation guidelines (PR #181).

### yeast 8.3.2:
* Features:
  * Fixes #154: MetaNetX IDs added from the yeast7.6 [MetaNetX](https://www.metanetx.org) model & from existing ChEBI and KEGG IDs in the model (PR #167).
  * Introduced contributing guidelines + code of conduct (PR #175).
* Fixes:
  * Fixes #161: Added as `rxnNotes` and `metNotes` the corresponding PR number (#112, #142, #149 or #156) in which each rxn and met was introduced (PR #170).
  * Fixes #169: Compartment error for `r_4238` (PR #170).
  * Corrected confidence score of rxns from PR #142 (PR #170).

### yeast 8.3.1:
* Features:
  * Added 21 reactions & 14 metabolites based on metabolomics data (PR #156).
  * Added metadata to the excel version of the model (PR #163).
  * Added `ComplementaryData/physiology` with biological data of yeast (PR #159).
* Fixes/Others:
  * Fixed bug that underestimated the biomass content (PR #159).
  * Fitted GAM to chemostat data (PR #159).

### yeast 8.3.0:
* Features:
  * Added 225 new reactions and 148 new metabolites, based on growth data from a Biolog substrate usage experiment on carbon, nitrogen, sulfur and phosphorus substrates (PR #149).
* Fixes/Others:
  * Removed verbose details from `README.md` (PR #150).
  * Updated RAVEN, which added extra annotation to the `.yml` file (PR #151).
  * Minor changes to `saveYeastModel.m` (PR #152).
  * Model is now stored simulating minimal media conditions (PR #157).

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