"""Independent verification for the annotation-tsv migration (yeast-GEM#379).

The migration itself runs through real MATLAB
(``code/modelCuration/migrateAnnotationTsvs.m``): it extracts reaction,
metabolite and gene cross-reference annotation (KEGG, BiGG, ChEBI,
MetaNetX, EC codes, UniProt) out of ``model/yeast-GEM.yml`` into
``model/reactions.tsv``, ``model/metabolites.tsv`` and
``model/genes.tsv``, rewrites the yml without those fields, and verifies
itself by merging the tsvs back in (via ``code/annotateGEM.m``) and
comparing against the original load. reactions.tsv/metabolites.tsv also
carry a read-only ``name`` column (the yml stays authoritative for
names); CI's ``model_qc.py`` gates on the two staying in agreement going
forward, but that only catches *future* drift -- this script checks the
migration itself got the initial copy right.

This script re-derives the same result independently, in Python, via
``raven_toolbox`` (a separate reader implementation from RAVEN's MATLAB
``readYAMLmodel.m``) and ``yeastgem.annotate.annotate_gem`` (a separate
merge implementation from ``annotateGEM.m``), and reports any
disagreement. Two independently-written implementations of the same
extraction/merge agreeing is much stronger evidence than either one
checking itself.

Kept as a historical record of how the migration was verified, like
``code/modelCuration/vX_Y_Z.m`` -- not meant to be run again after the
migration lands, though it will still pass if run against an unmodified
model/ (nothing here writes anything).

Run locally, after the MATLAB migration has produced the tsvs and the
lean yml:

    python code/python/release/migrate_annotation_tsvs.py \\
        --original-yml <path to the pre-migration yeast-GEM.yml>

``--original-yml`` is required and not defaulted: the whole point is
comparing against the model *before* the migration touched it, which
this repository's working tree no longer has once the migration has
been applied and committed. Get one with, e.g.:

    git show <commit-before-migration>:model/yeast-GEM.yml > /tmp/original.yml
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import raven_toolbox.io as rio

from yeastgem.annotate import GENE_COLUMNS, MET_COLUMNS, RXN_COLUMNS, annotate_gem
from yeastgem.io import REPO_PATH

_MODEL_DIR = REPO_PATH / "model"
_LEAN_YML = _MODEL_DIR / "yeast-GEM.yml"


def _annotation_values(annotation: dict, key: str) -> set:
    if key not in annotation:
        return set()
    value = annotation[key]
    return set(value) if isinstance(value, list) else {value}


def _check_entities(entities_a, entities_b, columns: tuple[str, ...], label: str) -> list[str]:
    mismatches = []
    for entity_a in entities_a:
        entity_b = entities_b.get_by_id(entity_a.id)
        for column in (*columns, "sbo"):
            values_a = _annotation_values(entity_a.annotation, column)
            values_b = _annotation_values(entity_b.annotation, column)
            if values_a != values_b:
                mismatches.append(
                    f"{label} {entity_a.id} {column}: original={sorted(values_a)} "
                    f"merged={sorted(values_b)}"
                )
    return mismatches


def _check_names(entities, tsv_path: Path, label: str) -> list[str]:
    """model/yeast-GEM.yml is authoritative for names; the tsv's name column
    is a read-only copy that must match it exactly (yeast-GEM#379)."""
    with tsv_path.open(encoding="utf-8", newline="") as fh:
        tsv_names = {row["id"]: row.get("name", "") for row in csv.DictReader(fh, delimiter="\t")}
    mismatches = []
    for entity in entities:
        if entity.id not in tsv_names:
            continue
        model_name = entity.name or ""
        tsv_name = tsv_names[entity.id]
        if model_name != tsv_name:
            mismatches.append(
                f"{label} {entity.id} name: model={model_name!r} tsv={tsv_name!r}"
            )
    return mismatches


def verify(original_yml: Path) -> list[str]:
    """Compare the original model against the lean yml + tsvs, merged back
    together via annotate_gem(). Returns the list of mismatches found (empty
    means the migration is verified lossless)."""
    print(f"Loading original model from {original_yml}")
    original = rio.read_yaml_model(original_yml)
    print(f"  {len(original.reactions)} reactions, {len(original.metabolites)} metabolites, "
          f"{len(original.genes)} genes")

    print(f"Loading lean model from {_LEAN_YML} and merging {_MODEL_DIR}/*.tsv via annotate_gem()")
    merged = annotate_gem(rio.read_yaml_model(_LEAN_YML), _MODEL_DIR)

    mismatches: list[str] = []
    mismatches += _check_entities(original.reactions, merged.reactions, RXN_COLUMNS, "reaction")
    mismatches += _check_entities(
        original.metabolites, merged.metabolites, MET_COLUMNS, "metabolite"
    )
    mismatches += _check_entities(original.genes, merged.genes, GENE_COLUMNS, "gene")

    mismatches += _check_names(original.reactions, _MODEL_DIR / "reactions.tsv", "reaction")
    mismatches += _check_names(original.metabolites, _MODEL_DIR / "metabolites.tsv", "metabolite")

    # Row-count/id-set sanity: every tsv id set matches the model exactly.
    for tsv_name, entities in (
        ("reactions.tsv", original.reactions),
        ("metabolites.tsv", original.metabolites),
        ("genes.tsv", original.genes),
    ):
        with (_MODEL_DIR / tsv_name).open(encoding="utf-8", newline="") as fh:
            tsv_ids = {row["id"] for row in csv.DictReader(fh, delimiter="\t")}
        model_ids = {e.id for e in entities}
        if tsv_ids != model_ids:
            mismatches.append(
                f"{tsv_name} id set mismatch: {tsv_ids - model_ids} extra, "
                f"{model_ids - tsv_ids} missing"
            )

    return mismatches


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--original-yml", required=True, type=Path,
                         help="path to model/yeast-GEM.yml as it was before the migration")
    args = parser.parse_args()

    mismatches = verify(args.original_yml)
    if mismatches:
        print(f"\n{len(mismatches)} MISMATCHES:")
        for m in mismatches[:50]:
            print(" ", m)
        return 1

    print("\nVerified: model/reactions.tsv, metabolites.tsv and genes.tsv, merged back into "
          "the lean model/yeast-GEM.yml, reproduce the original model's annotation exactly "
          "-- independently of the MATLAB migration's own self-check.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
