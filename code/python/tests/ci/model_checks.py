"""Structural checks for the pull-request report.

Each function takes the loaded model (or the model file) and returns
``(count, rows)``, the shape ``model_qc.py`` writes into ``qc_findings.md`` and
``qc_metrics.tsv``.

Checks:

* duplicate keys inside a ``!!omap`` entry of ``model/yeast-GEM.yml``. A
  repeated key (two ``name`` fields, or the same metabolite twice in a
  stoichiometry) is accepted by RAVEN but makes ``cobra.io.load_yaml_model``
  raise a bare ``AssertionError``; the scan names the entry, the key and the
  line numbers;
* reactions without metabolites;
* reaction bounds and gene rules that cannot be right (lower bound above upper
  bound, a bound beyond +-1000, a gene that is not in ``model/genes.tsv``, a
  gene rule on a boundary reaction);
* reactions whose chemistry is split over compartments.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import yaml

# A membrane compartment and its lumen belong to the same organelle.
ORGANELLE = {"erm": "er", "gm": "g", "mm": "m", "vm": "v"}


def is_pseudoreaction(rxn) -> bool:
    """Biomass and its components are lumped compositions, written as pseudoreactions."""
    return "pseudoreaction" in (rxn.name or "").lower()


def in_pseudoreaction(met) -> bool:
    """Whether the metabolite is a substrate or product of a pseudoreaction: the lumped
    pools (biomass, protein, lipid, ...) and what they are made of."""
    return any(is_pseudoreaction(rxn) for rxn in met.reactions)


def balances_by_design(rxn) -> bool:
    """Whether a reaction is expected to balance and to keep its chemistry in one place.

    Three classes are not:

    * boundary reactions, an exchange has one side by construction;
    * pseudoreactions, biomass and its components are lumped compositions;
    * SLIME reactions, yeast-GEM splits lipids into measurable entities with
      fractional, averaged coefficients, so they are unbalanced by design.
      There are 186 of them, against 2 genuinely unexpected mass imbalances in
      the rest of the model, so leaving them in would make the count useless.
    """
    name = rxn.name or ""
    return not (rxn.boundary or is_pseudoreaction(rxn) or "SLIME" in name)


def check_duplicate_keys(model_file: Path) -> tuple[int, list]:
    """Rows ``(entry, scope, key, first_line, duplicate_line)`` for repeated !!omap keys."""
    stack: list[dict] = []
    expect: list[str] = []
    firsts: list[list] = []
    duplicates: list[tuple] = []

    def owner() -> str:
        for context in reversed(stack):
            if context["type"] == "omap" and context["id"]:
                return context["id"]
        return "(top level)"

    with open(model_file, encoding="utf-8") as handle:
        for event in yaml.parse(handle):
            kind = type(event).__name__
            line = event.start_mark.line + 1
            if kind == "SequenceStartEvent":
                is_omap = event.tag == "tag:yaml.org,2002:omap"
                stack.append({"type": "omap" if is_omap else "seq", "keys": {}, "id": None})
            elif kind == "SequenceEndEvent":
                stack.pop()
            elif kind == "MappingStartEvent":
                stack.append({"type": "map"})
                expect.append("key")
                firsts.append([None, None, line])
            elif kind == "MappingEndEvent":
                stack.pop()
                first_key, first_value, key_line = firsts.pop()
                expect.pop()
                if stack and stack[-1]["type"] == "omap" and first_key is not None:
                    parent = stack[-1]
                    if first_key in parent["keys"]:
                        scope = "field" if parent["id"] else "member"
                        duplicates.append(
                            (owner(), scope, first_key, parent["keys"][first_key], key_line))
                    else:
                        parent["keys"][first_key] = key_line
                    if first_key == "id" and parent["id"] is None:
                        parent["id"] = first_value
            elif kind == "ScalarEvent" and stack and stack[-1]["type"] == "map":
                position = len(expect) - 1
                if expect[position] == "key":
                    if firsts[position][0] is None:
                        firsts[position][0] = event.value
                    expect[position] = "value"
                else:
                    if firsts[position][1] is None:
                        firsts[position][1] = event.value
                    expect[position] = "key"
    return len(duplicates), sorted(duplicates)


def check_empty_reactions(model) -> tuple[int, list]:
    rows = sorted((r.id, r.name or "") for r in model.reactions if not r.metabolites)
    return len(rows), rows


def check_reaction_sanity(model, model_dir: Path) -> tuple[int, list]:
    """Rows ``(reaction, name, issues)`` for bounds and gene rules that cannot be right."""
    genes_tsv = model_dir / "genes.tsv"
    with genes_tsv.open(encoding="utf-8", newline="") as handle:
        documented = {row["id"] for row in csv.DictReader(handle, delimiter="\t")}
    rows = []
    for rxn in model.reactions:
        issues = []
        if rxn.lower_bound > rxn.upper_bound:
            issues.append("lb>ub")
        if rxn.lower_bound < -1000 or rxn.upper_bound > 1000:
            issues.append("|bound|>1000")
        undocumented = sorted(g.id for g in rxn.genes if g.id not in documented)
        if undocumented:
            issues.append("undocumented_genes:" + ";".join(undocumented))
        if rxn.boundary and rxn.gene_reaction_rule.strip():
            issues.append("boundary_with_gpr")
        if issues:
            rows.append((rxn.id, rxn.name or "", ";".join(issues)))
    rows.sort()
    return len(rows), rows


def dead_end_kind(met) -> str:
    """Why a metabolite is a dead end, from the bounds of its reactions.

    ``never consumed`` or ``never produced`` when no reaction can run in the direction
    that would consume or produce it; ``no open reaction`` when none can carry it either
    way; ``after other dead ends`` when each direction is open in some reaction, so the
    metabolite is a dead end only once the reactions blocked by other dead ends are
    removed."""
    produced = consumed = False
    for rxn in met.reactions:
        coefficient = rxn.metabolites[met]
        forward, backward = rxn.upper_bound > 0, rxn.lower_bound < 0
        produced |= (coefficient > 0 and forward) or (coefficient < 0 and backward)
        consumed |= (coefficient < 0 and forward) or (coefficient > 0 and backward)
    if produced and consumed:
        return "after other dead ends"
    if produced:
        return "never consumed"
    if consumed:
        return "never produced"
    return "no open reaction"


def equation_with_compartments(rxn) -> str:
    """The reaction equation with the metabolite names and their compartments, as
    ``2 ATP [c] + H2O [c] --> ADP [m]``."""
    def side(items):
        return " + ".join(
            (f"{abs(c):g} " if abs(c) != 1 else "") + f"{m.name} [{m.compartment}]"
            for m, c in sorted(items, key=lambda item: (item[0].name, item[0].compartment)))

    arrow = "<=>" if rxn.lower_bound < 0 < rxn.upper_bound else (
        "<--" if rxn.upper_bound <= 0 else "-->")
    reactants = [(m, c) for m, c in rxn.metabolites.items() if c < 0]
    products = [(m, c) for m, c in rxn.metabolites.items() if c > 0]
    return f"{side(reactants)} {arrow} {side(products)}"


def check_split_compartments(model) -> tuple[int, list]:
    """Reactions whose own chemistry, apart from what they transport, spans compartments.

    A metabolite that appears in two compartments of the same reaction is what
    the reaction transports and is ignored. The metabolites left over are the
    reaction's own chemistry, which normally takes place in one compartment; one
    in another compartment (a cytosolic substrate of a mitochondrial reaction)
    is usually a curation error. Some enzymes do work across a membrane, so this
    is a report, and only a rise in the count is flagged. Metabolites of
    different compartments are the same chemical when they share a name, and an
    organelle membrane counts as the organelle it surrounds.
    """
    rows = []
    for rxn in model.reactions:
        if not balances_by_design(rxn):
            continue
        compartments_of = defaultdict(set)
        for met in rxn.metabolites:
            compartments_of[met.name].add(met.compartment)
        chemistry = [m for m in rxn.metabolites if len(compartments_of[m.name]) == 1]
        compartments = sorted({ORGANELLE.get(m.compartment, m.compartment) for m in chemistry})
        if len(compartments) > 1:
            rows.append((rxn.id, rxn.name or "", ",".join(compartments),
                         equation_with_compartments(rxn)))
    rows.sort()
    return len(rows), rows
