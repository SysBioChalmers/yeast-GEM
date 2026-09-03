"""Merge reaction, metabolite and gene cross-reference annotation from
model/reactions.tsv, model/metabolites.tsv and model/genes.tsv into a
cobra Model.

Python port of `code/annotateGEM.m` (yeast-GEM#379): the yml stores only
the inline fields (sbo, deltaG, confidence_score, notes); the cross-
reference identifiers (KEGG, BiGG, ChEBI, MetaNetX, EC codes, UniProt)
live in the three tsv tables and are merged in here. Used to
independently cross-check the MATLAB migration, and available for
Stage 3's annotation-consistency checks.
"""
from __future__ import annotations

import csv
from pathlib import Path

import cobra

from yeastgem.io import REPO_PATH

RXN_COLUMNS = ("bigg.reaction", "ec-code", "kegg.pathway", "kegg.reaction", "metanetx.reaction")
MET_COLUMNS = ("bigg.metabolite", "chebi", "kegg.compound", "metanetx.chemical", "smiles")
GENE_COLUMNS = ("uniprot",)


def _read_tsv(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return {row["id"]: row for row in csv.DictReader(fh, delimiter="\t")}


def _merge_row(annotation: dict, row: dict, columns: tuple[str, ...]) -> None:
    """Merge each ';'-joined tsv cell in ``row`` into ``annotation`` by column
    name, appending to (not overwriting) whatever is already there under
    that key -- e.g. an existing 'sbo' entry, or another column's earlier
    write, are left alone."""
    for column in columns:
        cell = (row.get(column) or "").strip()
        if not cell:
            continue
        parts = [p.strip() for p in cell.split(";") if p.strip()]
        if not parts:
            continue
        existing = annotation.get(column)
        if existing is None:
            merged = parts
        elif isinstance(existing, list):
            merged = list(dict.fromkeys(existing + parts))
        else:
            merged = list(dict.fromkeys([existing, *parts]))
        annotation[column] = merged if len(merged) > 1 else merged[0]


def annotate_gem(
    model: cobra.Model,
    model_dir: str | Path | None = None,
    *,
    types: tuple[str, ...] = ("rxn", "met", "gene"),
) -> cobra.Model:
    """Merge the tsv cross-references into ``model``'s annotation, in place.

    Parameters
    ----------
    model
        Model whose reactions/metabolites/genes carry yeast-GEM ids.
    model_dir
        Directory holding reactions.tsv / metabolites.tsv / genes.tsv.
        Defaults to the repository's ``model/`` folder.
    types
        Which annotation classes to merge (``"rxn"``, ``"met"``, ``"gene"``).

    Returns
    -------
    cobra.Model
        The same ``model`` object, now annotated. Pass a copy if the
        caller needs to keep an un-annotated version -- this mirrors
        annotateGEM.m, which is likewise called on a copy in
        saveYeastYaml.m so that model/yeast-GEM.yml stays lean.
    """
    model_dir = Path(model_dir) if model_dir is not None else REPO_PATH / "model"

    if "rxn" in types:
        rows = _read_tsv(model_dir / "reactions.tsv")
        for rxn in model.reactions:
            row = rows.get(rxn.id)
            if row is not None:
                _merge_row(rxn.annotation, row, RXN_COLUMNS)

    if "met" in types:
        rows = _read_tsv(model_dir / "metabolites.tsv")
        for met in model.metabolites:
            row = rows.get(met.id)
            if row is not None:
                _merge_row(met.annotation, row, MET_COLUMNS)

    if "gene" in types:
        rows = _read_tsv(model_dir / "genes.tsv")
        for gene in model.genes:
            row = rows.get(gene.id)
            if row is not None:
                _merge_row(gene.annotation, row, GENE_COLUMNS)

    return model


def _cell(value: object) -> str:
    """Render an annotation value back into a tsv cell: ``None`` -> '',
    a list -> ';'-joined, anything else -> its string form."""
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return ";".join(str(v) for v in value)
    return str(value)


def _write_tsv(
    path: Path, entities, columns: tuple[str, ...], *, include_name: bool
) -> None:
    fieldnames = ["id", *(["name"] if include_name else []), *columns]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for entity in entities:
            row = {"id": entity.id}
            if include_name:
                row["name"] = entity.name
            for column in columns:
                row[column] = _cell(entity.annotation.get(column))
            writer.writerow(row)


def derive_annotation_tsvs(
    model: cobra.Model,
    model_dir: str | Path | None = None,
    *,
    types: tuple[str, ...] = ("rxn", "met", "gene"),
) -> None:
    """Write reactions.tsv / metabolites.tsv / genes.tsv from ``model``'s
    current annotation -- the reverse of :func:`annotate_gem`.

    Each requested tsv is rewritten from scratch, one row per id, in the
    model's current order. Closes a gap where a cross-reference edit made
    programmatically (e.g. by assigning ``rxn.annotation[...]`` directly),
    rather than by hand-editing a tsv cell, would otherwise be silently
    discarded the next time the model is saved -- the annotation strip
    that keeps model/yeast-GEM.yml lean (yeast-GEM#379) drops these
    fields and never writes them anywhere else.

    Parameters
    ----------
    model
        Model whose reactions/metabolites/genes carry the cross-reference
        annotation to persist.
    model_dir
        Directory to write reactions.tsv / metabolites.tsv / genes.tsv
        into. Defaults to the repository's ``model/`` folder.
    types
        Which annotation classes to write (``"rxn"``, ``"met"``,
        ``"gene"``).
    """
    model_dir = Path(model_dir) if model_dir is not None else REPO_PATH / "model"

    if "rxn" in types:
        _write_tsv(model_dir / "reactions.tsv", model.reactions, RXN_COLUMNS,
                   include_name=True)
    if "met" in types:
        _write_tsv(model_dir / "metabolites.tsv", model.metabolites, MET_COLUMNS,
                   include_name=True)
    if "gene" in types:
        _write_tsv(model_dir / "genes.tsv", model.genes, GENE_COLUMNS,
                   include_name=False)
