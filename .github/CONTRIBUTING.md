## Contributor guidelines

First of all, thank you for contributing to yeast-GEM!! Anybody is welcome to contribute, but please abide by the following guidelines.

You can contribute in 2 main ways: by creating issues, and by sending pull requests (PRs) with additions, deletions, corrections, etc. to the model.

### Reporting issues in the model

Report an issue at https://github.com/SysBioChalmers/yeast-GEM/issues if you note any of the following:

* Incorrect annotation for any model component.
* Missing feature or field you would like the model to have.
* Bug/weird simulation results.
* Lacking documentation.
* Any type of feedback.

If you are unsure about the issue, consider asking first in our [Gitter chat room](https://gitter.im/SysBioChalmers/yeast-GEM).

When creating the issue, please make sure:

* You tested your code (if any) with all requirements for running the model.
* You did your analysis in the `master` branch of the repository.
* You provide any necessary files/links needed for understanding the issue.
* You checked that a similar issue does not exist already

Feel free to also comment on any of the [open issues](https://github.com/SysBioChalmers/yeast-GEM/issues). When doing so, please comply with our [code of conduct](https://github.com/SysBioChalmers/yeast-GEM/blob/master/.github/CODE_OF_CONDUCT.md).

Finally, if you like yeast-GEM please remember to 'star' our Github page (click on the star at top right corner), that way we also have an idea of who is using yeast-GEM!

### Contributing to the model

Want to contribute to the model with some additions or improvements? Consider starting by raising an issue and assign it to yourself to describe what you want to achieve. This way, we reduce the risk of duplicated efforts and you may also get suggestions on how to best proceed, e.g. there may be half-finished work in some branch that you could start with. Also, feel free to browse our [open issues](https://github.com/SysBioChalmers/yeast-GEM/issues) and our [ongoing projects](https://github.com/SysBioChalmers/yeast-GEM/projects): Anything tagged with "help wanted" is open to whoever wants to implement it!

Here's how to set up yeast-GEM for local development to contribute smaller features or changes that you can implement yourself:

1. First of all, make sure that you have all [requirements](https://github.com/SysBioChalmers/yeast-GEM#required-software---contributor) for contributing to yeast-GEM. Note that COBRA and RAVEN should be cloned in your computer and not just downloaded.

2. Fork the yeast-GEM repository on GitHub (go to https://github.com/SysBioChalmers/yeast-GEM & click on the upper right corner).

3. Clone your fork locally:
    ```
    git clone https://github.com/<your Github name>/yeast-GEM.git
    ```

4. Check out the branch that you want to contribute to. Most likely that will be ``devel``:
    ```
    git checkout devel
    ```

5. Create a branch for local development based on the previously checked out branch ([see below](#branching-model) for details on the branching model and how to name your branch):
    ```
    git checkout -b name-of-your-branch
    ```

6. Now you can make your changes locally!
    * Always make your changes in Matlab and never directly editing the model files. For loading the model use `loadYeastModel.m`, and for saving the model use `saveYeastModel.m`. The latter will ensure all 3 versions of the model (`.xml`, `.yml` & `.txt`) get updated in the same way.
    * If your changes are minor (e.g. a single chemical formula you wish to correct), you can do it directly from the command line.
    * If your changes are not so small and require several steps, create a script that loads the model, reads data (if applicable), changes the model accordingly, and saves the model back.
    * Each script should start with a commented section describing the script, explaining the parameters, and indicating your name and the date it was written. Existing functions can clarify what style should be used.
    * As much as possible follow the model conventions:
      - For new metabolites and/or reactions, please use `/code/otherChanges/getNewIndex.m` for obtaining new ids.
      - For metabolite and/or reaction names, please avoid any unconventional characters (e.g. Greek letters).
      - For new genes, please use as id the [systematic names from SGD](http://seq.yeastgenome.org/help/community/nomenclature-conventions).
    * Store scripts in the appropriate folder in `/code` and data (as `.tsv` files) in the appropriate folder in `/data`. If you think no folder is adequate for your script/data, feel free to create your own folder. Note that binary data such as `.mat` structures or `.xls` tables cannot be stored in the repo (as they cannot be version-controlled, and they increment too much the size of the repo).
    * When you are done making changes, review locally your changes with `git diff` or any git client, to make sure you are modifying the model as you intended.

7. Commit your changes and push your branch to GitHub.
    ```
    git add .
    git commit -m "Title of your commit"
    git push origin name-of-your-branch
    ```
    [See below](#semantic-commits) for recommendations on how to name your commits. In case of larger updates, you can of course make several commits on a single contribution. However, if you need to do too many commits, consider if your contribution could be instead split into separate branches (making it easier for reviewing later).

8. Submit a pull request through the GitHub website (https://help.github.com/articles/creating-a-pull-request-from-a-fork/) to the `devel` branch of the original SysBioChalmers repo (not your fork). We recommend ticking the box "Allow edits from maintainers" if you wish for us to be able to contribute directly to your branch (speeding-up the reviewing process).

Note that steps 3, 4, 5 and 7 can be done, if you prefer, with any git client, such as [Github Desktop](https://desktop.github.com/).

Finally, and for larger features that you want to work on collaboratively, you may consider to first request to join our development team to get write access to the repository so that you can create a branch directly in the main repository (or simply ask the administrator to create a branch for you). Once you have a new branch, you can push your changes directly to the main repository and when finished, submit a pull request from that branch to ``devel``. [See below](#development-team-guidelines) for more details.

Thank you very much for contributing to yeast-GEM!

#### Branching model

* `devel`: Is the branch all pull-requests should be based on.

* `master`: Is only touched by the administrator and is the branch with the tested & reviewed model that is released or ready for the next release.

* `gh-pages`: Is only touched by the administrator and for maintaining the [landing page](http://sysbiochalmers.github.io/yeast-GEM/) of yeast-GEM.

* `{chore, doc, feat, fix, refactor, style}/descriptive-name`: Any other branch created in the model. If you work on a fix, start the branch name with `fix/`, if you work on a feature, start the branch name with `feat/`. Examples: `fix/format_reactions` or `feat/new_algorithms`. [See below](#semantic-commits) for more details on the possible actions you can use.

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

More examples [here](https://github.com/SysBioChalmers/yeast-GEM/commits/master). A more detailed explanation or comments is encouraged to be left in the commit description.

## Development team guidelines

This section is meant for the development team of yeast-GEM. As a member of the development team, you should comply with all [previous contributor guidelines](#contributor-guidelines) as well. Besides that, please also consider the following guidelines.

### Creating pull requests

Changes to the model _cannot_ be directly committed to the `master` or `devel` branches (in fact they are protected). Commits are made to side-branches, after which pull requests are made for merging with `devel`. For this, follow the [instructions](#contributing-to-the-model) for contributors, but consider that as members of the development team have write access to the repository, you can create a branch directly in the main repository without needing to fork, for your convenience. This means that you can:

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
* **Compatibility:** First of all, make sure that the model is still compatible with the loading/saving wrappers (`loadYeastModel.m` & `saveYeastModel.m`) and that no errors appear. Check also that [`dependencies.txt`](https://github.com/SysBioChalmers/yeast-GEM/blob/master/model/dependencies.txt) does not change in any unexpected ways (e.g. an "unknown" toolbox version). Finally, ensure that the SBML fields `model metaid`, `model id` and `model name` never change, as if they change it would create a conflict in the next release.
* **Documentation:** Every change should be justified with a reference/link/argument. This can be provided as data in `/data`, or directly as a comment in the pull request.
* **Reproducibility:** If there are any added scripts, make sure that if you run them, the model gets updated from how it was in `devel` to how it is in the pull request. For this, you may _locally_ switch to the corresponding branch, replace the `.xml` file with the one from `devel` before the changes, and run the associated scripts. Remember to stash any undesired changes afterwards.
* **Style:** Ensure that the changes to the model are compliant with the model's rxn/met/gene naming conventions (when unsure, take a look at a similar field in the model). Also, make sure that scripts have a compliant style, and datasets are straight-forward to understand.
* When commenting in the review, please comply with our [code of conduct](https://github.com/SysBioChalmers/yeast-GEM/blob/master/.github/CODE_OF_CONDUCT.md).
* Avoid vague comments and try to be as explicit as possible (e.g.: _"Please include X here"_ instead of _"X could be included here"_).
* As much as possible, try to keep the review process in the pull request discussion, and not in separate private emails.

## Administrator guidelines

This section is meant for the administrator of this repo. The main duties of the administrator are:
* To make sure conventions and standards in the model are kept.
* To keep the repository clean and organized, i.e. avoid redundancy in functions and/or data, and keep coherency in naming of files.
* To help in the reviewing process of external pull requests by assigning reviewers, [labels](https://github.com/SysBioChalmers/yeast-GEM/issues/labels) and [projects](https://github.com/SysBioChalmers/yeast-GEM/projects), if applicable.
* To keep [issues](https://github.com/SysBioChalmers/yeast-GEM/issues) with the proper labels, and to close them once they are fixed in the `master` branch.
* In cases of disagreement between contributors, to decide how to resolve the issue.
* To merge open pull requests into `devel` regularly (see [below](#merging-contributions)).
* To generate new releases of the model regularly (see [below](#releasing-a-new-version)).

### Merging contributions
The following points should be considered when merging branches to `devel`:
* Make sure the branch gets accepted by at least one developer with writing access.
* Wait at least a day before merging, to allow other developers to inspect the pull request.
* As soon as the branch is merged, check if `devel` is still possible to merge to `master` (this can be checked [here](https://github.com/SysBioChalmers/yeast-GEM/compare/devel)). If conflicts appear (which should rarely happen and only if the `.xml` file was changed in an unexpected way by a toolbox update), fix the conflict _locally_ as soon as possible in `devel` and then push it (note, **DO NOT** pull any other changes from `master` to `devel`, just the single file that is creating the conflict).

### Releasing a new version

* A merge of `devel` with `master` invokes a new release.
* A new release should be made as soon as there is substantial new work in `devel` (as rule of thumb, after around 3 pull request merges).

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

When releasing, please follow these steps:
  1. Make sure all dependencies in `devel` correspond to the setup from the local computer from which the release will be made. If not, make a single commit in `devel` updating this with a `loadYeastModel`/`saveYeastModel` cycle.
  2. Create a pull request from `devel` to `master`, indicating all new features/fixes/etc. and referencing every previous pull request included (examples [here](https://github.com/SysBioChalmers/yeast-GEM/releases)). Tip: if any [issue](https://github.com/SysBioChalmers/yeast-GEM/issues) gets solved in the release, write in the pull request description "Fixes #X", where "X" is the issue number. That way the issue will be automatically closed after merge.
  3. Merge at least a day after (having at least one accepted review).
  4. Switch locally to `master` and update `history.md`, by putting at the top the same description of the corresponding PR from step 2.
  5. Bump version with `increaseVersion.m`. **NOTE:** The function will error if unexpected changes are occurring. If this happens, probably step 1 was done incorrectly. To fix it, commit in `devel` any necessary changes and make a new pull request.
  6. Commit changes from steps 4 and 5 with the message `chore: version X.Y.Z`.
  7. Make the new release at GitHub [here](https://github.com/SysBioChalmers/yeast-GEM/releases/new), using the proper tag "vX.Y.Z" and with the same description as the corresponding PR from step 2.
  8. Pull from `master` to `gh-pages` to update the landing page.
  9. Review the [Zenodo](https://zenodo.org) release: Every new release from Github (step 7) automatically triggers a new release in Zenodo. However, so far it is not possible to fully customize this release, and some manual curation is needed. This includes:
      * Ensuring the title of the release has the format `SysBioChalmers/yeast-GEM: yeast X.Y.Z`.
      * Correcting author names.
      * Ensuring the version of the release has the format `vX.Y.Z`.
      * Setting the language to English.
      * Adding any Grant IDs (if applicable).

     Make sure to both save & publish your edits. You will find the new release at the top of [all Zenodo releases](https://zenodo.org/search?page=1&size=20&q=conceptrecid:%221494182%22&sort=-publication_date&all_versions=True). Note that it might take some minutes for the Zenodo release to appear after you create the release in Github.

## Acknowledgments

These contribution guidelines were written based on the contribution guidelines of [opencobra/cobrapy](https://github.com/opencobra/cobrapy/blob/devel/.github/CONTRIBUTING.rst) and [SysBioChalmers/RAVEN](https://github.com/SysBioChalmers/RAVEN/wiki/DevGuidelines).
