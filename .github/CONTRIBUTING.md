## Contributor guidelines

First of all, thank you for contributing to yeast-GEM! Anybody is welcome to contribute, but please abide by the following guidelines.

You can contribute in 2 main ways: by creating issues, and by sending pull requests (PRs) with additions, deletions, corrections, etc. to the model.

### Reporting issues in the model

Report an issue at https://github.com/SysBioChalmers/yeast-GEM/issues if you note any of the following:

* Incorrect annotation for any model component.
* Missing feature or field you would like the model to have.
* Bug/weird simulation results.
* Lacking documentation.
* Any type of feedback.

If you are unsure about the issue, consider asking first in our [Gitter chat room](https://gitter.im/SysBioChalmers/yeast-GEM) or [GitHub Discussion forum](https://github.com/SysBioChalmers/yeast-GEM/discussions).

When creating the issue, please make sure:

* You tested your code (if any) with all requirements for running the model.
* You did your analysis in the `main` branch of the repository.
* You provide any necessary files/links needed for understanding the issue.
* You checked that a similar issue does not exist already

Feel free to also comment on any of the [open issues](https://github.com/SysBioChalmers/yeast-GEM/issues). When doing so, please comply with our [code of conduct](https://github.com/SysBioChalmers/yeast-GEM/blob/main/.github/CODE_OF_CONDUCT.md).

Finally, if you like yeast-GEM please remember to 'star' our Github page (click on the star at top right corner), that way we also have an idea of who is using yeast-GEM!

### Contributing to the model

Do you want to contribute to the model with some additions or improvements? Consider starting by raising an issue and assign it to yourself to describe what you want to achieve. This way, we reduce the risk of duplicated efforts and you may also get suggestions on how to best proceed, e.g. there may be half-finished work in some branch that you could start with. Also, feel free to browse our [open issues](https://github.com/SysBioChalmers/yeast-GEM/issues) and our [ongoing projects](https://github.com/SysBioChalmers/yeast-GEM/projects): Anything tagged with "help wanted" is open to whoever wants to implement it!

Here is how to set up yeast-GEM for local development to contribute smaller features or changes that you can implement yourself:

1. First of all, make sure that you have all [requirements](https://github.com/SysBioChalmers/yeast-GEM#developer) for contributing to yeast-GEM.

2. Fork the yeast-GEM repository on GitHub (go to https://github.com/SysBioChalmers/yeast-GEM & click on the upper right corner).

3. Clone your fork locally:
    ```
    git clone https://github.com/<your Github name>/yeast-GEM.git
    ```

4. Check out the branch that you want to contribute to. Most likely that will be ``develop``:
    ```
    git checkout develop
    ```

5. Create a branch for local development based on the previously checked out branch ([see below](#branching-model) for details on the branching model and how to name your branch):
    ```
    git checkout -b name-of-your-branch
    ```

6. Now you can make your changes locally!
    * Always make your changes in MATLAB and never directly editing the model files. For loading the model use `loadYeastYaml.m`, and for saving it back use `saveYeastYaml.m` — it applies the canonical medium, adds SBO terms, checks aerobic/anaerobic growth, and writes `model/yeast-GEM.yml`. Only when you need a local `.xml`/`.txt`/`.xlsx`/`.mat` build (e.g. to test against another tool that needs SBML) run `commitYeastModel.m` separately; most curation never needs it, since routine PRs only change `model/yeast-GEM.yml` and the annotation tsvs. Python curators have the same three functions: `load_yeast_yaml`/`save_yeast_yaml`/`commit_yeast_model` in the `yeastgem` package.
    * Reaction, metabolite and gene cross-reference identifiers (KEGG, BiGG, ChEBI, MetaNetX, EC numbers, UniProt) live in `model/reactions.tsv`, `model/metabolites.tsv` and `model/genes.tsv`, not in `model/yeast-GEM.yml` itself (see [model/README.md](../model/README.md)). Curate those directly by editing the matching column — no RAVEN or cobrapy needed. `loadYeastYaml.m`/`load_yeast_yaml` merges them back into the model automatically, so this is transparent for any other kind of curation. A newly added reaction, metabolite or gene needs a row added to the matching tsv too, and CI fails the build if the model and the tsvs disagree about which ids exist. `reactions.tsv`/`metabolites.tsv` also carry a `name` column: the yml stays authoritative, so if you rename a reaction or metabolite, update the tsv row's `name` to match in the same change — CI fails the build otherwise. Removing a reaction or metabolite also needs an entry in `data/deprecatedIdentifiers/deprecatedReactions.tsv` or `deprecatedMetabolites.tsv` (CI reports, but does not block, if this is missed) — and reusing a deprecated id for something new fails the build.
    * `python code/python/annotation/verify_annotations.py --all` (or `--mets`/`--rxns` with a comma-separated id list) cross-checks `reactions.tsv`/`metabolites.tsv` identifiers against MetaNetX's structure-matched reference tables and reports any that are wrong, missing, or drifted to a deprecated MetaNetX id. Runs on every pull request (a dedicated `annotation` job in `model-qc.yml`, since it needs a ~200MB MetaNetX download — cached across runs) and reports into the same pull-request comment as the rest of `model-qc.yml`, compared against the target branch's own last-committed summary rather than a fresh run on it, the same way the MEMOTE score is — see `data/testResults/README.md`. Not a gate: `wrong` needs a human decision, not a pass/fail rule — a real structural difference (e.g. a sugar's ring vs. open-chain tautomer, or a monatomic ion vs. its neutral atom) can register as a mismatch without being a curation mistake, so check the metabolite's MetaNetX reference name/formula/charge (shown alongside the finding) before deciding. You can also run it yourself locally, e.g. scoped to just the metabolites/reactions you are curating with `--mets`/`--rxns`. `--fix` applies the safe corrections (missing/drift) directly to the tsvs; it rewrites the whole file with `\n` line endings, so the diff will show every row as touched even though only the fixed cells' values actually changed — check the *content* diff before committing, not the line count. Every run also writes a full markdown report and a `key<TAB>value` summary into `--out` (default `data/testResults`, matching `model_qc.py`) with every finding, not just the console preview. `wrong` findings are further split into isolated — a unique suggested correction, worth a look on its own — and recurring — the same correction shared by 3+ entities, grouped into a collapsible section so you can sanity-check the pattern once instead of row by row.
    * `python code/python/annotation/verify_annotation_consistency.py --all` (or `--mets`/`--rxns`) checks something the tool above does not: whether an entity's own columns (e.g. `chebi` and `kegg.compound`) agree with each other via MetaNetX's id graph, independent of structure — resolved, where possible, using the entity's own structure match or a majority of its other columns, giving a suggested fix and the MetaNetX name/formula behind each column's value; where nothing settles it, those same names/formulas are shown for you to compare by eye. Not run in CI — most of what it finds predates any one pull request, since annotation drift accumulates across the whole tsv over time rather than from a single edit — so run it yourself after curating a batch of annotation, or before a release, and work through the report (default `--out data/annotation`) by hand. No `--fix`: even a resolved conflict can touch multiple columns at once, which deserves a final look before committing.
    * Curations to the model are consolidated in scripts that can convert a model from one yeast-GEM release (=version) to the next. Therefore, use the `/code/modelCuration/TEMPLATEcuration.m` file, or contribute directly to an existing `v$VERSION.m` script, where `$VERSION` matches the most recent yeast-GEM release. If the most recent release is 8.6.0, there will be a `/code/modelCuration/v8_6_0.m` script that consolidates curations to be applied to release 8.6.0.
    * See the `TEMPLATEcuration.m` script with more description on how to use this file. Look at previous curations files (e.g. `v8_6_0.m`) to get an idea of how this works.
    * As much as possible follow the model conventions:
      - For metabolite and/or reaction names, please avoid any unconventional characters (e.g. Greek letters).
      - For new genes, please use as id the [systematic names from SGD](http://seq.yeastgenome.org/help/community/nomenclature-conventions).
    * Store generic scripts in the appropriate folder in `/code` and data (as `.tsv` files) in the appropriate folder in `/data`. If a data file is just to apply one round of curation (for instance a table with selected gene associations to be corrected), it should be stored in the appropriate subfolder in `/code/modelCuration`. If you think no folder is adequate for your script/data, feel free to create your own folder. Note that binary data such as `.mat` structures or `.xls` tables cannot be stored in the repo (as they cannot be version-controlled, and they increment too much the size of the repo).
    * When you are done making changes, review locally your changes with `git diff` or any git client, to make sure you are modifying the model as you intended.

7. Commit your changes and push your branch to GitHub.
    ```
    git add .
    git commit -m "Title of your commit"
    git push origin name-of-your-branch
    ```
    [See below](#semantic-commits) for recommendations on how to name your commits. In case of larger updates, you can of course make several commits on a single contribution. However, if you need to do too many commits, consider if your contribution could be instead split into separate branches (making it easier for reviewing later).

8. Submit a pull request through the GitHub website (https://help.github.com/articles/creating-a-pull-request-from-a-fork/) to the `develop` branch of the original SysBioChalmers repo (not your fork). We recommend ticking the box "Allow edits from maintainers" if you wish for us to be able to contribute directly to your branch (speeding-up the reviewing process).

Note that steps 3, 4, 5 and 7 can be done, if you prefer, with any git client, such as [Github Desktop](https://desktop.github.com/).

Finally, and for larger features that you want to work on collaboratively, you may consider to first request to join our development team to get write access to the repository so that you can create a branch directly in the main repository (or simply ask the administrator to create a branch for you). Once you have a new branch, you can push your changes directly to the main repository and when finished, submit a pull request from that branch to ``develop``. [See below](#development-team-guidelines) for more details.

Thank you very much for contributing to yeast-GEM!

#### Branching model

* `develop`: Is the branch all pull-requests should be based on.

* `main`: Is only touched by the administrator and is the branch with the tested & reviewed model that is released or ready for the next release.

* `gh-pages`: Is only touched by the administrator and for maintaining the [landing page](http://sysbiochalmers.github.io/yeast-GEM/) of yeast-GEM.

* `{chore, doc, feat, fix, refactor, style}/descriptive-name`: Any other branch created in the model. If you work on a fix, start the branch name with `fix/`, if you work on a feature, start the branch name with `feat/`. Examples: `fix/format_reactions` or `feat/new_algorithms`. [See below](#semantic-commits) for more details on the possible actions you can use.

* `release/X.Y.Z`: Created and merged automatically by the [Release workflow](#releasing-a-new-version). Never create or push to one of these by hand.

Never merge `main` into `develop` directly, and never push to `develop` from `main`. The binary model files (`.mat`, `.xlsx`) are tracked on `main` only; merging `main` into `develop` by hand would carry them across, which `branch-hygiene.yml` then blocks on every subsequent pull request until they are removed again. The release workflow's sync-back pull request (see [below](#releasing-a-new-version)) is the only sanctioned way for `main`'s content to reach `develop`.

#### Semantic commits

Please use concise descriptive commit messages. Ideally, use semantic commit messages to make it easier to show what you are aiming to do:

`action-object: brief description`
* `action` = {`feat`, `fix`, `refactor`, `style`, `doc`, `chore`, `test`}
* `object` (optional field) = {`rxn`, `rxn.annot`, `rxn.prop`, `met`, `met.annot`, `met.prop`, `gene`, `gene.annot`, `comp`, `comp.annot`, `data`}

`action` refers to what exactly are you doing in the commit, following a [standard definition](http://karma-runner.github.io/2.0/dev/git-commit-msg.html) in software development:
* `chore`: updating toolbox, data files, etc.
* `doc`: updating documentation or explanatory comments in functions.
* `feat`: new feature added, e.g. new reaction / metabolite / function / etc.
* `fix`: something that was incorrect in the model and now has been corrected.
* `refactor`: see [code refactoring](https://en.wikipedia.org/wiki/Code_refactoring).
* `style`: minor format changes of model, functions or data (spaces, semi-colons, etc., no code change).

`object` may or may not be used, and refers to the receiver of the action, i.e. what part of the model/data are you modifying:
* `object`: Used when adding/removing/modifying the whole object (`rxn`, `met`, `gene`, `comp` or `data`).
* `object.annot`: Annotation of the object (e.g. KEGG IDs, EC numbers, etc). Applies for `rxn`, `met`, `gene` or `comp`.
* `object.prop`: Any other property of the object (e.g. formulas, charge, stoich coeffs, LB, UB, etc). Applies for `rxn` or `met`.

Some examples:

|commit|commit message|
|:---:|:---:|
|Add new rxns|`feat-rxn: methanol pathway`|
|Remove a metabolite|`fix-met: duplicated citrate`|
|Add metabolite formula|`feat-met.prop: carbohydrate formulas`|
|Fix rxn stoichiometry|`fix-rxn.prop: complex V stoich coeffs`|
|Update gene IDs|`fix-gene.annot: update IDs from SWISSPROT`|
|Format name of compartment|`style-comp.annot: remove uppercases`|
|Split a rxn in 2|`refactor-rxn: split isomerase in 2 steps`|
|Add some data|`feat-data: metabolomics data`|
|Update documentation of function|`doc: addDBnewRxn.m`|
|Update toolbox|`chore: update RAVEN version`|

More examples [here](https://github.com/SysBioChalmers/yeast-GEM/commits/main). A more detailed explanation or comments is encouraged to be left in the commit description.

## Development team guidelines

This section is meant for the development team of yeast-GEM. As a member of the development team, you should comply with all [previous contributor guidelines](#contributor-guidelines) as well. Besides that, please also consider the following guidelines.

### Creating pull requests

Changes to the model _cannot_ be directly committed to the `main` or `develop` branches (in fact they are protected). Commits are made to side-branches, after which pull requests are made for merging with `develop`. For this, follow the [instructions](#contributing-to-the-model) for contributors, but consider that as members of the development team have write access to the repository, you can create a branch directly in the main repository without needing to fork, for your convenience. This means that you can:

* Skip step 2 of the contribution process.
* In step 3 of the contribution process, clone directly the original repo:
  ```
  git clone https://github.com/SysBioChalmers/yeast-GEM.git
  ```

Follow all other steps in the same way. Also, when creating your pull request (or after you have created it):
* Choose 1 or more members of the team (ideally with knowledge on the pull request) as reviewers. Note that the person making the pull request and the reviewer _cannot_ be the same person.
* Assign appropriate [labels](https://github.com/SysBioChalmers/yeast-GEM/issues/labels).
* Assign the pull request to one of the available [projects](https://github.com/SysBioChalmers/yeast-GEM/projects), if applicable.

### Reviewing pull requests

Every pull request must be approved by at least one reviewer before it can be merged. When reviewing someone else's pull request, keep in mind the following aspects:
* **Compatibility:** First of all, make sure that the model is still compatible with the loading/saving wrappers (`loadYeastYaml.m` & `saveYeastYaml.m`) and that no errors appear. Check also that [`dependencies.txt`](https://github.com/SysBioChalmers/yeast-GEM/blob/main/model/dependencies.txt) does not change in any unexpected ways (e.g. an "unknown" toolbox version). Finally, ensure that the SBML fields `model metaid`, `model id` and `model name` never change, as if they change it would create a conflict in the next release.
* **Documentation:** Every change should be justified with a reference/link/argument. This can be provided as data in `/data`, or directly as a comment in the pull request.
* **Reproducibility:** If there are any added scripts, make sure that if you run them, the model gets updated from how it was in `develop` to how it is in the pull request. For this, you may _locally_ switch to the corresponding branch, replace the `.xml` file with the one from `develop` before the changes, and run the associated scripts. Remember to stash any undesired changes afterwards.
* **Style:** Ensure that the changes to the model are compliant with the model's rxn/met/gene naming conventions (when unsure, take a look at a similar field in the model). Also, make sure that scripts have a compliant style, and datasets are straight-forward to understand.
* When commenting in the review, please comply with our [code of conduct](https://github.com/SysBioChalmers/yeast-GEM/blob/main/.github/CODE_OF_CONDUCT.md).
* Avoid vague comments and try to be as explicit as possible (e.g.: _"Please include X here"_ instead of _"X could be included here"_).
* As much as possible, try to keep the review process in the pull request discussion, and not in separate private emails.

## Administrator guidelines

This section is meant for the administrator of this repo. The main duties of the administrator are:
* To make sure conventions and standards in the model are kept.
* To keep the repository clean and organized, i.e. avoid redundancy in functions and/or data, and keep coherency in naming of files.
* To manage package dependencies and regularly update them.
* To help in the reviewing process of external pull requests by assigning reviewers, [labels](https://github.com/SysBioChalmers/yeast-GEM/issues/labels) and [projects](https://github.com/SysBioChalmers/yeast-GEM/projects), if applicable.
* To keep [issues](https://github.com/SysBioChalmers/yeast-GEM/issues) with the proper labels, and to close them once they are fixed in the `main` branch.
* In cases of disagreement between contributors, to decide how to resolve the issue.
* To merge open pull requests into `develop` regularly (see [below](#merging-contributions)).
* To generate new releases of the model regularly (see [below](#releasing-a-new-version)).

### Merging contributions
The following points should be considered when merging branches to `develop`:
* Make sure the branch gets accepted by at least one developer with writing access.
* Wait at least a day before merging, to allow other developers to inspect the pull request.
* As soon as the branch is merged, check if `develop` is still possible to merge to `main` (this can be checked [here](https://github.com/SysBioChalmers/yeast-GEM/compare/develop)). If conflicts appear (which should rarely happen and only if the `.xml` file was changed in an unexpected way by a toolbox update), fix the conflict _locally_ as soon as possible in `develop` and then push it (note, **DO NOT** pull any other changes from `main` to `develop`, just the single file that is creating the conflict).

### Managing python dependencies

We use [pip-tools](https://github.com/jazzband/pip-tools) for managing dependencies:

* If a new dependency is needed for users/developers:
  1. Add the dependency to `/requirements/requirements.in` or `/requirements/dev-requirements.in`, respectively (note that by default any requirement in `requirements.in` will also be in `dev-requirements.in`).
  2. From `/requirements`, run `pip-compile requirements.in` and/or `pip-compile dev-requirements.in` as needed.

* If dependencies need to be upgraded, from `/requirements` run:
  ```
  pip-compile --upgrade requirements.in
  pip-compile --upgrade dev-requirements.in
  ```
  Dependencies should be upgraded regularly, but always first tested in separate branches.

### Releasing a new version

* A merge of `develop` with `main` invokes a new release.
* A new release should be made as soon as there is substantial new work in `develop` (as rule of thumb, after around 3 pull request merges).

yeast-GEM follows [semantic versioning](https://semver.org/), adapted to GEMs:
* A `major` release is seldom used and only meant for a new publication. Backwards compatibility should be, ideally, always preserved.
* A `minor` release involves a substantial change in the model (several new reactions/metabolites/genes), such as:
  * Addition of genes/reactions/metabolites from a whole genome annotation.
  * Addition of several annotation fields.
  * Inclusion of a major new formalism in the model.
  * Addition of a plurality of pathways.
* A `patch` release is the most common one and is done when only few things have changed in the model, or there are only changes that have to do with format, such as:
  * Adding a single new annotation field.
  * Fixing some chemical formulas/charges.
  * Updating toolboxes.
  * Re-organization of data
  * Refactoring of code.

The release itself is automated by [`.github/workflows/release.yml`](../.github/workflows/release.yml). The administrator's job is to trigger it, review what it produces, and approve the resulting pull requests -- not to run any of the steps by hand.

1. Update `history.md` by hand, putting the description of everything to be released at the top under a new `### yeast X.Y.Z:` heading. This is the only manual editing step; everything else reads this heading to know what version is being released.
2. Trigger the [Release workflow](https://github.com/SysBioChalmers/yeast-GEM/actions/workflows/release.yml) with "Run workflow", giving it the version number (`X.Y.Z`, matching the heading from step 1). It will:
   * validate that `X.Y.Z` is a legal bump from `version.txt` and that `history.md` has a matching heading, and cut `release/X.Y.Z` from `develop`;
   * run `code/python/release/increase_version.py`, the Python port of `increaseVersion.m`: stamp the version and date into `model/yeast-GEM.yml`, run the model tests, and refresh `README.md`, `data/testResults/README.md`, `growth.md`/`growth.png` and `essentialGenes.tsv`;
   * export `model/yeast-GEM.xml`, `.txt`, `.xlsx` and `.mat` from the now-stamped `.yml`, via RAVEN -- the one step that still needs MATLAB: there is no Python writer for `.txt`/`.xlsx`/`.mat`, and `.xml` is produced alongside them here so all four come from the same, single RAVEN export rather than `.xml` from one writer and the rest from another;
   * take a memote snapshot and open a pull request from `release/X.Y.Z` into `main`.
3. Review the pull request like any other: check the model-quality and validation-metrics comment for regressions, and the linked memote snapshot. Wait at least a day, get at least one approving review, then merge.
4. On merge, the workflow tags `vX.Y.Z`, publishes the GitHub release, and opens a pull request syncing `main` back into `develop` -- with the binary model files removed again, so they never reach `develop` (see `branch-hygiene.yml`, which would otherwise flag them). Merge that pull request too.
5. Review the [Zenodo](https://zenodo.org) release: every new release from GitHub (step 4) automatically triggers a new release in Zenodo. However, so far it is not possible to fully customize this release, and some manual curation is needed. This includes:
    * Ensuring the title of the release has the format `SysBioChalmers/yeast-GEM: yeast X.Y.Z`.
    * Correcting author names to include all commit authors and PR reviewers that contributed to the release.
    * Ensuring the version of the release has the format `vX.Y.Z`.
    * Setting the language to English.
    * Adding any grant IDs (if applicable).

   Make sure to both save & publish your edits. You will find the new release at the top of [all Zenodo releases](https://zenodo.org/search?page=1&size=20&q=conceptrecid:%221494182%22&sort=-publication_date&all_versions=True). Note that it might take some minutes for the Zenodo release to appear after you create the release in GitHub.

## Acknowledgments

These contribution guidelines were written based on the contribution guidelines of [opencobra/cobrapy](https://github.com/opencobra/cobrapy/blob/develop/.github/CONTRIBUTING.rst) and [SysBioChalmers/RAVEN](https://github.com/SysBioChalmers/RAVEN/wiki/DevGuidelines).
