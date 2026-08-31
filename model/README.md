# yeast-GEM model and annotation files

This directory contains the yeast-GEM model and annotation files.

## Model

The model is available as `.xml`, `.txt`, `.yml`, and (on `main` only) `.xlsx` and `.mat`.

`model/yeast-GEM.yml` is the file to edit directly when curating the model
(yeast-GEM#379) — reactions, metabolites, genes, stoichiometry, bounds,
gene-reaction rules, subsystems, SBO terms, ΔG values and confidence scores
all live there. It is intentionally **not** the place to look for
cross-reference identifiers to external databases (KEGG, BiGG, ChEBI,
MetaNetX, EC numbers, UniProt): those live in the three TSV files below
instead, kept separate so the yml stays lean, diffable, and reviewable
directly on GitHub.

`loadYeastModel.m` merges the TSV annotation back in automatically after
reading the yml (via `annotateGEM.m`), so day-to-day curation in MATLAB sees
exactly the same fully-annotated model as before this split — the merge is
transparent. `saveYeastModel.m` writes `model/yeast-GEM.yml` from a stripped
copy (so it stays lean) and writes `.xml`/`.txt`/`.xlsx`/`.mat` from the full,
unstripped model (so those keep carrying every identifier, exactly as
before). Python's `code/python/yeastgem` reads/writes `.xml` directly and is
unaffected by any of this.

## Reaction, Metabolite, and Gene Annotation

External identifiers for yeast-GEM reactions, metabolites, and genes are
provided as `tsv` files, one row per entity, keyed by `id`. A multi-valued
column (e.g. a reaction with more than one EC number) is `;`-joined within
the cell.

`reactions.tsv` and `metabolites.tsv` also carry a `name` column, right
after `id`, so the file is readable without cross-referencing the yml.
**`model/yeast-GEM.yml` is authoritative for names** — renaming a reaction
or metabolite is a model edit, the same category as changing its formula or
bounds, so it's still done in the yml, never in the tsv directly. The tsv's
`name` column is a read-only copy, kept in sync by convention; CI fails the
build if the two disagree (`model_qc.py`'s `name_mismatches` check, a gate —
see [data/testResults/README.md](../data/testResults/README.md)), so drift
gets caught rather than silently accumulating. If you rename something in
the yml, update the matching tsv row's `name` column in the same change.

* `reactions.tsv`:

| column            | content                                    |
|-------------------|---------------------------------------------|
| `id`              | identical to the reaction's `id` in the yml |
| `name`            | identical to the reaction's `name` in the yml (read-only copy — see above) |
| `bigg.reaction`   | BiGG reaction ID                            |
| `ec-code`         | EC number(s)                                |
| `kegg.pathway`    | KEGG pathway ID(s)                          |
| `kegg.reaction`   | KEGG reaction ID                            |
| `metanetx.reaction` | MetaNetX reaction ID                      |

* `metabolites.tsv`:

| column              | content                                       |
|---------------------|------------------------------------------------|
| `id`                | identical to the metabolite's `id` in the yml  |
| `name`              | identical to the metabolite's `name` in the yml (read-only copy — see above) |
| `bigg.metabolite`   | BiGG metabolite ID                             |
| `chebi`             | ChEBI ID                                       |
| `kegg.compound`     | KEGG compound ID                               |
| `metanetx.chemical` | MetaNetX chemical ID                           |
| `smiles`            | SMILES structure                               |

* `genes.tsv`:

| column    | content     |
|-----------|-------------|
| `id`      | identical to the gene's `id` in the yml |
| `uniprot` | UniProt ID  |

To curate a cross-reference identifier, edit the matching column directly —
no RAVEN, cobrapy, or other package required, the same way the yml itself is
meant to be hand-edited (yeast-GEM#379). A newly added reaction, metabolite,
or gene needs a (possibly mostly-blank) row added to the matching TSV too;
CI's `model_qc.py` gates the build on the model and the tsvs agreeing exactly
(id-for-id, and name-for-name — see [data/testResults/README.md](../data/testResults/README.md)).

Removing a reaction or metabolite needs its row deleted here *and* an entry
added to `data/deprecatedIdentifiers/deprecatedReactions.tsv` or
`deprecatedMetabolites.tsv`, so the retired id stays resolvable instead of
silently vanishing — CI reports (does not block) when this is missed.
