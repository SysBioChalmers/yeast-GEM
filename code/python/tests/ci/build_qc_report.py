"""Render the model-quality report for the pull-request comment.

Compares the checks computed on this pull request against the same checks
computed on the target branch, and renders one markdown table per group.
The comparison is what makes the report readable: a model of this size
carries findings that predate any one pull request, so an absolute count
says little, while a count that *rose* is a regression this branch
introduced.

Icon rule, per row:

* growth -- :white_check_mark: if the model grows, :x: if it cannot. This
  is a gate and blocks the merge.
* name_mismatches and annotation_consistency (both a "gate_count") --
  :white_check_mark: if zero, :x: if non-zero. Also gates: by the time a
  pull request is reviewed, the yml and the tsvs are supposed to already
  agree, so any non-zero count can only mean this branch introduced one
  -- unlike an ordinary count, a pre-existing non-zero value here is not
  a legitimate "no worse than before" state to warn about instead of
  blocking.
* model-file checks (qc_status.tsv) -- :white_check_mark: for ``pass``,
  :x: and fatal for ``fail``, :grey_question: if the check did not run.
* every other count -- :x: if it rose against the target branch (a
  regression), :warning: if it is non-zero but no worse (a pre-existing
  finding, not blocking), :white_check_mark: if it is zero.
* :hourglass_flowing_sand: -- the group is still being computed on this run.

One exception: the model-size line (reaction/metabolite/gene counts) at
the top of the comment carries no icon at all. Unlike every other count
here, a curation pull request is expected to change it, so treating a
rise like a regression would flag normal work.

MEMOTE and the MetaNetX annotation-verification sections are rendered
separately (_render_memote, _render_annotation) rather than through the
icon rule above: both are compared against the target branch's own
last-committed result instead of a fresh run on it (each has its own
reason -- MEMOTE's cost, MetaNetX's download), and neither counts a rise
as a regression, since annotation `wrong` findings need a human decision
and MetaNetX's own reference data can move them on its own.

Only the gates (growth, duplicate keys, reactions without metabolites,
name_mismatches, annotation_consistency, the model-file checks) fail the
build. Everything else is reported.

Two renderings come from the same rows: the condensed comment (one row per
section and three lists) and the full report, committed as
model_qc_summary.md.

Usage:
    python build_qc_report.py --current data/testResults --base /tmp/base \\
        --base-ref develop --url-base https://github.com/.../data/testResults \\
        --out report.md --summary-out data/testResults/model_qc_summary.md
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

# (metric key, label)
#
# No icon, unlike every other row in this file: a curation PR is expected
# to add or remove reactions, metabolites and genes, so a changed count is
# not itself a finding -- only a place for the reader to notice a size
# change they were not expecting.
_SIZE_ROWS = [
    ("n_reactions", "reactions"),
    ("n_metabolites", "metabolites"),
    ("n_genes", "genes"),
]

# (metric key, comment label, kind)
_MODEL_ROWS = [
    ("growth", "Growth (biomass producible)", "growth"),
    ("duplicate_keys", "Duplicate keys in the model file", "gate_count"),
    ("empty_reactions", "Reactions without metabolites", "gate_count"),
    ("duplicate_reactions", "Exact-duplicate reaction groups", "count"),
    ("orphan_metabolites", "Unused metabolites", "count"),
    ("unused_genes", "Unused genes", "count"),
    ("missing_formula", "Metabolites missing formula", "count"),
    ("missing_charge", "Metabolites missing charge", "count"),
    ("reaction_sanity", "Reaction bound and gene rule issues", "count"),
    ("split_compartments", "Reactions split across compartments", "count"),
]

_BALANCE_ROWS = [
    ("mass_imbalanced", "Mass-imbalanced reactions", "count"),
    ("charge_imbalanced", "Charge-imbalanced reactions", "count"),
    ("macaw_dead_end_metabolites", "Dead-end metabolites", "count"),
    ("macaw_single_direction",
     "Reactions that can only carry flux one way", "count"),
    ("macaw_duplicates", "Reactions flagged as MACAW duplicates", "count"),
]

_ANNOTATION_ROWS = [
    ("malformed_xrefs", "Malformed cross-references", "count"),
    ("xrefs_across_compartments", "Cross-refs inconsistent across compartments",
     "count"),
    ("name_mismatches", "Reaction/metabolite names disagreeing with the tsvs",
     "gate_count"),
    ("annotation_consistency", "Model/annotation-table consistency", "gate_count"),
    ("deprecation_completeness", "Removed identifiers not added to the deprecated lists",
     "count"),
    ("structure_inconsistent", "Metabolite structure (SMILES) vs. formula/charge",
     "count"),
    ("delta_g_consistency", "deltaG tsvs out of sync with the model", "count"),
]

# (status key, comment label). Read from qc_status.tsv, where each check
# records "pass" or "fail"; all three are gates.
_STATUS_ROWS = [
    ("roundtrip_yaml", "YAML round-trip"),
    ("roundtrip_sbml", "SBML round-trip"),
    ("yamllint", "YAML lint"),
]

# (metric key -> the reason clause used in the "Merge blocked" verdict)
_FATAL_MESSAGES = {
    "growth": "the model cannot grow",
    "duplicate_keys": "the model file has duplicate keys",
    "empty_reactions": "a reaction has no metabolites",
    "roundtrip_yaml": "the YAML round-trip failed",
    "roundtrip_sbml": "the SBML round-trip failed",
    "yamllint": "the YAML lint failed",
    "name_mismatches": "reaction/metabolite names disagree between the model "
                        "and the tsvs",
    "annotation_consistency": "the model and its annotation tables disagree, "
                               "or a deprecated identifier is in active use",
}


# (metric key, comment label, direction, tolerance)
#
# "direction" says which way is an improvement, which is what makes a
# delta on a continuous metric readable: a rising R2 is a gain, a rising
# error or false-call count is not. A move within tolerance still shows
# its real value -- only the icon (not the number) treats it as clean, so
# solver noise gets a checkmark without hiding what actually moved.
# Validation metrics whose ideal value is zero. These follow the same icon
# rule as the QC counts -- non-zero is a standing finding worth a warning
# even when this branch did not make it worse -- rather than showing a
# green tick merely because nothing moved. The others have no meaningful
# zero: a growth R2 of 0 would be a disaster, not a clean result.
_ZERO_IS_IDEAL = {
    "geneEssentialityFalseNonEssential",
    "geneEssentialityFalseEssential",
}

# Metrics whose detail lives under another metric's heading. The four
# gene counts are the confusion matrix behind the accuracy, so they share
# its section rather than repeating the same table four times.
_DETAIL_ANCHOR = {
    "geneEssentialityTrueNonEssential": "Gene essentiality accuracy",
    "geneEssentialityTrueEssential": "Gene essentiality accuracy",
    "geneEssentialityFalseNonEssential": "Gene essentiality accuracy",
    "geneEssentialityFalseEssential": "Gene essentiality accuracy",
}

_VALIDATION_ROWS = [
    ("growthPredictionR2", "Growth prediction R2", "higher", 5e-3),
    ("anaerobicFluxMedianFoldError", "Anaerobic flux median fold error",
     "lower", 5e-2),
    ("anaerobicExchangeMeanRelativeError",
     "Anaerobic exchange mean relative error", "lower", 1e-2),
    ("anaerobicAmmoniumPerATPase", "Ammonium per ATPase", "none", 5e-2),
    ("geneEssentialityAccuracy", "Gene essentiality accuracy", "higher", 5e-3),
    ("geneEssentialityTrueNonEssential", "True non-essential genes",
     "higher", 2),
    ("geneEssentialityTrueEssential", "True essential genes", "higher", 2),
    ("geneEssentialityFalseNonEssential", "False non-essential genes",
     "lower", 2),
    ("geneEssentialityFalseEssential", "False essential genes", "lower", 2),
]


def read_metrics(directory: Path, name: str = "qc_metrics.tsv") -> dict[str, float]:
    """Read a ``key<TAB>value`` metrics file, if it exists."""
    path = directory / name
    if not path.is_file():
        return {}
    metrics = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        key, _, value = line.partition("\t")
        try:
            metrics[key.strip()] = float(value)
        except ValueError:
            continue
    return metrics


def _slug(label: str) -> str:
    """GitHub heading anchor for a label.

    Mirrors GitHub's algorithm -- lowercase, drop punctuation, spaces to
    hyphens -- so a row label links to the identically-named heading in
    the testResults README.
    """
    text = label.lower().replace("`", "")
    text = "".join(c for c in text if c.isalnum() or c in " -_")
    return text.strip().replace(" ", "-")


def _icon(value: float, base: float | None, kind: str):
    """Return ``(delta_text, icon, is_regression, is_fatal)``."""
    if kind == "growth":
        grows = value > 1e-6
        icon = ":white_check_mark:" if grows else ":x:"
        if base is None:
            return "new", icon, False, not grows
        change = value - base
        text = f"{change:+.3g}" if abs(change) > 1e-9 else "0"
        return text, icon, False, not grows

    if kind == "gate_count":
        # Like "count" for the delta/icon shown, but any non-zero value is
        # fatal outright -- not just a rise vs. base -- because by the time
        # this exists the two sides being compared are supposed to already
        # agree, so a non-zero count can only mean this branch broke it.
        is_fatal = value > 0
        icon = ":x:" if is_fatal else ":white_check_mark:"
        if base is None:
            return "new", icon, False, is_fatal
        change = int(value) - int(base)
        text = f"+{change}" if change > 0 else ("0" if change == 0 else str(change))
        return text, icon, change > 0, is_fatal

    if base is None:
        return "new", (":warning:" if value else ":white_check_mark:"), False, False
    change = int(value) - int(base)
    if change > 0:
        return f"+{change}", ":x:", True, False
    text = "0" if change == 0 else str(change)
    if value > 0:
        return text, ":warning:", False, False
    return text, ":white_check_mark:", False, False


def _render_size(current, base, base_ref, running) -> str:
    """One-line model-size summary: no icon, no gate.

    A changed reaction/metabolite/gene count is expected on essentially
    every curation pull request, so treating it like the other counts --
    green unless it rose -- would flag normal work as a regression. This
    is context for the reader, not a check.
    """
    if running:
        return "_running_"
    if not current:
        return "_not run._"
    parts = []
    for key, label in _SIZE_ROWS:
        if key not in current:
            parts.append(f"{label}: _not run_")
            continue
        value = int(current[key])
        if key in base:
            change = int(value - int(base[key]))
            delta = f"{change:+d}" if change else "0"
        else:
            delta = "new"
        parts.append(f"{value} {label} ({delta})")
    return ", ".join(parts) + f" vs `{base_ref}`"


def _row_records(rows, current, base, running):
    """One record per row: the value, delta, icon and how it counts in the verdict.

    ``state`` is "running" (the group is still being computed), "skipped" (the
    check produced no value, so its result is unknown) or "measured".
    """
    records = []
    for key, label, kind in rows:
        record = {"key": key, "label": label, "kind": kind, "state": "measured",
                  "value": None, "text": "", "delta": "", "icon": "",
                  "regression": False, "fatal": False, "improved": False}
        if running:
            record.update(state="running", icon=":hourglass_flowing_sand:")
        elif key not in current:
            # A check with no value once the run is over did not execute -- an
            # optional dependency missing, or the check crashed. Saying
            # "running" there would be a lie that never resolves, and counting
            # it as clean would be worse.
            record.update(state="skipped", icon=":grey_question:")
        else:
            value = current[key]
            delta, icon, regression, fatal = _icon(value, base.get(key), kind)
            text = f"{value:.4g}" if kind == "growth" else str(int(value))
            improved = (kind != "growth" and key in base
                        and int(value) < int(base[key]))
            record.update(value=value, text=text, delta=delta, icon=icon,
                          regression=regression, fatal=fatal, improved=improved)
        records.append(record)
    return records


def _row_link(record, url_base):
    label = record["label"]
    return f"[{label}]({url_base}/README.md#{_slug(label)})" if url_base else label


def _row_value(record, detail_base):
    """A non-zero count links to its section in the findings file, which is
    addressed by the same slug as the README heading, so the link follows the
    label automatically."""
    text = record["text"]
    if detail_base and record["kind"] in ("count", "gate_count") and record["value"]:
        return f"[{text}]({detail_base}/qc_findings.md#{_slug(record['label'])})"
    return text


def _render_rows(records, url_base, detail_base):
    lines = []
    for record in records:
        name = _row_link(record, url_base)
        if record["state"] == "running":
            lines.append(f"| {name} | _running_ | | {record['icon']} |")
        elif record["state"] == "skipped":
            lines.append(f"| {name} | _not run_ | | {record['icon']} |")
        else:
            lines.append(f"| {name} | {_row_value(record, detail_base)} | "
                         f"{record['delta']} | {record['icon']} |")
    return lines


def read_status(directory: Path | None) -> dict[str, str]:
    """Read ``qc_status.tsv`` (``check<TAB>result``), if it exists."""
    path = directory / "qc_status.tsv" if directory else None
    if path is None or not path.is_file():
        return {}
    status = {}
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        key, _, value = line.partition("\t")
        if key.strip():
            status[key.strip()] = value.strip()
    return status


# essentialGenes.tsv classifications. The reference is the Stanford deletion collection:
# a true call agrees with it, a false non-essential gene is essential there but viable in
# the model, a false essential gene is viable there but not in the model.
_CALL_LABELS = {
    "TP": "true non-essential", "TN": "true essential",
    "FP": "false non-essential", "FN": "false essential",
}
_CORRECT_CALLS = {"TP", "TN"}
_MAX_LISTED = 20  # changed genes named in the comment; all are in the full report


def read_calls(directory: Path | None) -> dict[str, str]:
    """Read ``essentialGenes.tsv`` (``gene<TAB>classification``), if it exists."""
    path = directory / "essentialGenes.tsv" if directory else None
    if path is None or not path.is_file():
        return {}
    calls = {}
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        gene, _, call = line.partition("\t")
        if gene.strip():
            calls[gene.strip()] = call.strip()
    return calls


def read_gene_names(path: Path | None) -> dict[str, str]:
    """Systematic to standard gene names, from ``data/databases/SGDgeneNames.tsv``."""
    if path is None or not path.is_file():
        return {}
    names = {}
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        systematic, _, standard = line.partition("\t")
        if systematic and standard.strip():
            names[systematic] = standard.strip()
    return names


def essential_summary(calls: dict, base_calls: dict, names: dict) -> dict | None:
    """The genes that changed class since the target branch, and the wrong calls.

    ``changed`` rows are ``(gene, name, before, after, verdict)`` with the verdict
    ``better`` or ``worse`` (a call that became right or wrong), ``changed`` (wrong
    both times), ``new`` or ``removed`` (the gene is only in one model). ``None``
    when the branch has no calls."""
    if not calls:
        return None
    changed = []
    if base_calls:
        for gene in sorted(set(calls) | set(base_calls)):
            before, after = base_calls.get(gene), calls.get(gene)
            if before == after:
                continue
            if before is None:
                verdict = "new"
            elif after is None:
                verdict = "removed"
            elif (before in _CORRECT_CALLS) != (after in _CORRECT_CALLS):
                verdict = "better" if after in _CORRECT_CALLS else "worse"
            else:
                verdict = "changed"
            changed.append((gene, names.get(gene, ""), before, after, verdict))
    wrong = {call: [(gene, names.get(gene, "")) for gene, c in sorted(calls.items()) if c == call]
             for call in ("FP", "FN")}
    return {"changed": changed, "false_non_essential": wrong["FP"],
            "false_essential": wrong["FN"], "has_base": bool(base_calls)}


def _gene_label(gene: str, name: str) -> str:
    return f"`{gene}` ({name})" if name else f"`{gene}`"


def _call_text(call: str | None) -> str:
    return _CALL_LABELS.get(call, "not in the model") if call else "not in the model"


def _essential_cell(summary: dict | None, base_ref: str, running: bool) -> str:
    if running:
        return ":hourglass_flowing_sand: running"
    if summary is None:
        return ":grey_question: not run"
    wrong = (f"{len(summary['false_non_essential'])} false non-essential, "
             f"{len(summary['false_essential'])} false essential")
    if not summary["has_base"]:
        return f":grey_question: no baseline; {wrong}"
    count = {verdict: sum(1 for row in summary["changed"] if row[4] == verdict)
             for verdict in ("better", "worse", "changed", "new", "removed")}
    parts = [f"{n} {verdict}" for verdict, n in count.items() if n]
    if not parts:
        return f":white_check_mark: no gene changed class vs `{base_ref}`; {wrong}"
    icon = ":warning:" if count["worse"] else ":white_check_mark:"
    return f"{icon} {', '.join(parts)} vs `{base_ref}`; {wrong}"


def _change_line(row: tuple) -> str:
    gene, name, before, after, verdict = row
    return f"- {_gene_label(gene, name)}: {_call_text(before)} → {_call_text(after)} ({verdict})"


def _render_essential_full(summary: dict | None, base_ref: str) -> list[str]:
    lines = [
        "### Gene essentiality",
        "_Single-gene deletions against the Stanford deletion collection. A false "
        "non-essential gene is essential in the collection and viable in the model; a "
        "false essential gene is viable in the collection and not in the model. The "
        "per-gene calls are in [essentialGenes.tsv](essentialGenes.tsv)._",
        "",
    ]
    if summary is None:
        return [*lines, "_not run, the validation job reported no gene calls._", ""]
    if summary["has_base"]:
        lines += [f"#### Changed since `{base_ref}`", ""]
        if summary["changed"]:
            lines += ["| Gene | Name | Before | Now | Effect |", "| --- | --- | --- | --- | --- |"]
            lines += [f"| `{g}` | {n} | {_call_text(b)} | {_call_text(a)} | {v} |"
                      for g, n, b, a, v in summary["changed"]]
        else:
            lines += ["_No gene changed class._"]
        lines += [""]
    for title, key in (("False non-essential genes", "false_non_essential"),
                       ("False essential genes", "false_essential")):
        lines += [f"#### {title} ({len(summary[key])})", ""]
        lines += ["| Gene | Name |", "| --- | --- |"] if summary[key] else ["_None._"]
        lines += [f"| `{g}` | {n} |" for g, n in summary[key]] + [""]
    return lines


def _status_records(status, running):
    records = []
    for key, label in _STATUS_ROWS:
        record = {"key": key, "label": label, "kind": "status", "state": "measured",
                  "value": None, "text": "", "delta": "", "icon": "",
                  "regression": False, "fatal": False, "improved": False}
        if running:
            record.update(state="running", icon=":hourglass_flowing_sand:")
        elif key not in status:
            record.update(state="skipped", icon=":grey_question:")
        else:
            passed = status[key] == "pass"
            record.update(text=status[key], icon=":white_check_mark:" if passed else ":x:",
                          fatal=not passed)
        records.append(record)
    return records


def _fmt_delta(change: float) -> str:
    """4 significant digits, or a bare ``0`` for an exact (non-)change.

    Tolerance decides whether a move counts as a *regression* -- it must
    not also decide whether the reader gets to see the move at all. A
    change of 0.0008 sitting under a 0.005 tolerance is still a real
    change; showing it as literally ``0`` reads as "nothing moved" when
    something did, just not enough to matter for the verdict.
    """
    return f"{change:+.4g}" if abs(change) > 1e-9 else "0"


def _validation_records(current, base, running):
    """One record per validation metric; ``regression`` and ``improved`` follow
    the metric's direction, and a move within tolerance is neither."""
    records = []
    for key, label, direction, tol in _VALIDATION_ROWS:
        record = {"key": key, "label": label, "kind": "validation", "state": "measured",
                  "value": None, "text": "", "delta": "", "icon": "",
                  "regression": False, "fatal": False, "improved": False}
        if running:
            record.update(state="running", icon=":hourglass_flowing_sand:")
        elif key not in current:
            record.update(state="skipped", icon=":grey_question:")
        else:
            value = current[key]
            record.update(value=value, text=f"{value:.4g}")
            if key not in base:
                record.update(delta="new", icon=":grey_question:", state="new")
            else:
                change = value - base[key]
                record["delta"] = _fmt_delta(change)
                if abs(change) <= tol:
                    record["icon"] = (":warning:" if key in _ZERO_IS_IDEAL and value > 0
                                      else ":white_check_mark:")
                elif direction == "none":
                    record["icon"] = ":warning:"
                else:
                    improved = (change > 0) if direction == "higher" else (change < 0)
                    record["improved"] = improved
                    record["regression"] = not improved
                    if not improved:
                        record["icon"] = ":x:"
                    elif key in _ZERO_IS_IDEAL and value > 0:
                        record["icon"] = ":warning:"
                    else:
                        record["icon"] = ":white_check_mark:"
        records.append(record)
    return records


def _render_validation(records, url_base, detail_base):
    """The validation-metrics table lines."""
    lines = []
    for record in records:
        name = _row_link(record, url_base)
        if record["state"] == "running":
            lines.append(f"| {name} | _running_ | | {record['icon']} |")
        elif record["state"] == "skipped":
            lines.append(f"| {name} | _not run_ | | {record['icon']} |")
        else:
            shown = record["text"]
            if detail_base:
                # The value links to the measurements it summarises, the same
                # way a QC count links to the entries behind it.
                anchor = _slug(_DETAIL_ANCHOR.get(record["key"], record["label"]))
                shown = f"[{shown}]({detail_base}/validation_findings.md#{anchor})"
            lines.append(f"| {name} | {shown} | {record['delta']} | {record['icon']} |")
    return lines


def _parse_memote(directory: Path | None):
    """Parse ``memote_score.md`` -> ``(total_pct, {section: pct}, detailed)``.

    ``detailed`` is a list of ``(section, test, pct)`` triples, matched by
    its 3-column shape rather than by which heading it sits under -- the
    2-column section table above it cannot match a 3-column pattern, so
    this is unambiguous without needing to scope the search.

    None if ``directory`` is unset or the file is absent (the job did not
    run, or has not reached this directory yet) rather than [] or 0 -- a
    missing score and a score of zero must not read the same to the caller.
    """
    if directory is None:
        return None
    path = directory / "memote_score.md"
    if not path.is_file():
        return None
    text = path.read_text(encoding="utf-8")
    total = re.search(r"Total score:\s*([\d.]+)\s*%", text)
    if not total:
        return None
    sections = {m.group(1): float(m.group(2))
                for m in re.finditer(r"^\| (\w+) \| ([\d.]+)% \|$", text, re.M)}
    detailed = [(s, t, float(pct)) for s, t, pct in re.findall(
        r"^\| (.+?) \| (.+?) \| ([\d.]+)% \|$", text, re.M)]
    return float(total.group(1)), sections, detailed


def _memote_suite(directory: Path | None) -> str:
    """"full suite" when memote_score.md records a full run, else "fast subset"."""
    path = directory / "memote_score.md" if directory else None
    if path is not None and path.is_file() and "Full suite" in path.read_text(encoding="utf-8"):
        return "full suite"
    return "fast subset"


def _memote_delta(current: float, base: float | None) -> tuple[str, str]:
    """Percentage-point delta and icon, MEMOTE's own scale -- not the 4g
    convention used for validation metrics, since a score is already a
    rounded percent. Never :x:, unlike an ordinary count: this is compared
    against a stale baseline (see _render_memote), so a drop is shown, not
    scored as a build-blocking regression."""
    if base is None:
        return "new", ":grey_question:"
    change = current - base
    if abs(change) < 0.05:
        return "0", ":white_check_mark:"
    icon = ":white_check_mark:" if change > 0 else ":warning:"
    return f"{change:+.1f}", icon


def _render_memote(current: Path | None, base: Path | None, base_ref: str,
                    running: bool) -> list[str]:
    """The MEMOTE-score section, or a note if the job did not report.

    Compared against the target branch's own last-committed score rather
    than a fresh run on the base model, the way the other sections do it:
    MEMOTE (even the fast subset) is expensive enough that re-running it a
    second time in the same job would roughly double this job's cost for a
    score that mostly reflects annotation completeness, not something a
    typical pull request moves. See memote_snapshot.py's docstring.
    """
    header = [f"| Section | Score | &Delta; vs `{base_ref}` | |",
              "| --- | ---: | ---: | :---: |"]
    if running:
        return [*header, "| _running_ | | | :hourglass_flowing_sand: |", ""]
    parsed = _parse_memote(current)
    if parsed is None:
        return ["_not run — the MEMOTE job reported nothing._", ""]
    total, sections, detailed = parsed
    base_parsed = _parse_memote(base) if base is not None else None
    base_total, base_sections, _base_detailed = base_parsed or (None, {}, [])

    lines = list(header)
    total_delta, total_icon = _memote_delta(total, base_total)
    lines.append(f"| **Total** | {total:.1f}% | {total_delta} | {total_icon} |")
    for name, value in sections.items():
        delta, icon = _memote_delta(value, base_sections.get(name))
        lines.append(f"| {name} | {value:.1f}% | {delta} | {icon} |")
    lines.append("")
    # Collapsed by default: one row per MEMOTE test is too long to sit
    # in the open comment, but still worth having a click away rather
    # than only in the uploaded memote_result.json artifact.
    if detailed:
        lines += ["<details><summary>Per-test scores</summary>", "",
                  "| Section | Test | Score |", "| --- | --- | ---: |"]
        lines += [f"| {section} | {test} | {pct:.1f}% |"
                  for section, test, pct in detailed]
        lines += ["", "</details>", ""]
    return lines


# (metrics key, comment label). "wrong" splits into isolated vs. recurring
# the same way the report does; "drift"/"missing" are the two --fix can
# apply safely, kept for visibility rather than because a rise means
# anything is broken. Cross-annotation conflicts are a separate, occasional
# check (verify_annotation_consistency.py) that does not run in CI, so
# there is no annotation_conflict row here.
_ANNOTATION_SUMMARY_ROWS = [
    ("annotation_wrong_isolated", "`wrong`, isolated"),
    ("annotation_wrong_recurring", "`wrong`, recurring"),
    ("annotation_drift", "`drift` (safe to `--fix`)"),
    ("annotation_missing", "`missing` (safe to `--fix`)"),
]


def _annotation_delta(value: float, base: float | None) -> tuple[str, str]:
    """A count delta and icon, shown but never :x: -- see _render_annotation
    for why a rise here is not scored as a build-blocking regression the
    way an ordinary count's would be."""
    if base is None:
        return "new", (":warning:" if value else ":white_check_mark:")
    change = int(value) - int(base)
    if change == 0:
        return "0", ":white_check_mark:"
    text = f"+{change}" if change > 0 else str(change)
    icon = ":warning:" if change > 0 else ":white_check_mark:"
    return text, icon


def _render_annotation(current_dir: Path | None, base_dir: Path | None,
                        base_ref: str, running: bool, detail_base: str) -> list[str]:
    """The MetaNetX annotation-verification section, or a note if the job
    did not report.

    Like MEMOTE (see _render_memote), compared against the target branch's
    own last-committed summary rather than a fresh run on it -- the
    MetaNetX download makes a second run costly for little benefit -- and,
    unlike the other QC counts, a rise here is never treated as a
    regression: "wrong" needs a human decision, not a pass/fail rule (see
    the tool's own docstring), and MetaNetX's own reference data can move
    this count on its own, independent of anything this pull request did.
    """
    header = [f"| Finding | Count | &Delta; vs `{base_ref}` | |",
              "| --- | ---: | ---: | :---: |"]
    if running:
        return [*header, "| _running_ | | | :hourglass_flowing_sand: |", ""]
    current = read_metrics(current_dir, "annotation_metrics.tsv") if current_dir else {}
    if not current:
        return ["_not run — the annotation job reported nothing._", ""]
    base = read_metrics(base_dir, "annotation_metrics.tsv") if base_dir else {}
    lines = list(header)
    for key, label in _ANNOTATION_SUMMARY_ROWS:
        value = current.get(key)
        if value is None:
            lines.append(f"| {label} | _not run_ | | :grey_question: |")
            continue
        delta, icon = _annotation_delta(value, base.get(key))
        lines.append(f"| {label} | {int(value)} | {delta} | {icon} |")
    lines.append("")
    if detail_base:
        lines.append(f"[Full report]({detail_base}/annotation_report.md)")
        lines.append("")
    return lines


def collect(current: dict, base: dict, base_ref: str, running: bool,
            validation: dict | None = None, validation_base: dict | None = None,
            status: dict | None = None, essential: dict | None = None) -> dict:
    """Compute every row once; ``render_full`` and ``render_condensed`` format it."""
    groups = [
        ("Model checks",
         "_Growth, duplicate keys and reactions without metabolites are gates and "
         "block the merge; every other row is a non-blocking report._",
         _row_records(_MODEL_ROWS, current, base, running)),
        ("Mass/charge balance and MACAW", "",
         _row_records(_BALANCE_ROWS, current, base, running)),
        ("Annotations",
         "_Reaction/metabolite name agreement with the tsvs, and model/"
         "annotation-table consistency, are also gates; every other row "
         "is a non-blocking report._",
         _row_records(_ANNOTATION_ROWS, current, base, running)),
    ]
    model_files = _status_records(status or {}, running)
    if running or any(r["state"] == "measured" for r in model_files):
        groups.append(("Model files", "_All three are gates._", model_files))

    validation_records = _validation_records(validation or {}, validation_base or {}, running)
    records = [r for _, _, group in groups for r in group] + validation_records
    return {
        "groups": groups,
        "validation": validation_records,
        "essential": essential,
        "regressions": sum(r["regression"] for r in records),
        "warnings": sum(r["icon"] == ":warning:" for r in records
                        if r["kind"] != "validation" and r["state"] == "measured"),
        "pending": sum(r["state"] == "running" for r in records),
        "measured": sum(r["state"] == "measured" for r in records),
        "skipped": sum(r["state"] in ("skipped", "new") for r in records),
        "fatal_keys": [r["key"] for r in records if r["fatal"]],
    }


def _verdict(report: dict, base: dict, base_ref: str) -> str:
    regressions, pending = report["regressions"], report["pending"]
    if report["fatal_keys"]:
        reasons = " and ".join(
            _FATAL_MESSAGES.get(key, key) for key in dict.fromkeys(report["fatal_keys"]))
        verdict = f":x: **Merge blocked: {reasons}.** See the table(s) in the full report."
    elif regressions:
        extra = f" ({pending} check(s) still running)" if pending else ""
        verdict = (
            f":x: **{regressions} regression(s) vs `{base_ref}`** — this pull "
            f"request increased a finding count{extra}. Review the :x: rows."
        )
    elif pending:
        verdict = f":hourglass_flowing_sand: **{pending} check(s) still running.**"
        if report["measured"]:
            verdict += f" The rest are unchanged vs `{base_ref}`."
    elif not base:
        verdict = (
            ":information_source: First run for this comparison; no "
            "target-branch baseline yet."
        )
    elif report["warnings"]:
        verdict = (
            f":warning: **{report['warnings']} pre-existing finding(s), no regressions "
            f"vs `{base_ref}`.** Non-blocking."
        )
    else:
        verdict = (
            f":white_check_mark: **All checks clean, no regressions vs "
            f"`{base_ref}`.**"
        )
    if report["skipped"]:
        verdict += (
            f" {report['skipped']} check(s) did not run — their result is unknown, "
            "not clean."
        )
    return verdict


def render_full(report: dict, current: dict, base: dict, base_ref: str, url_base: str,
                running: bool, run_url: str, detail_base: str = "",
                memote_dir: Path | None = None, memote_base_dir: Path | None = None,
                annotation_dir: Path | None = None,
                annotation_base_dir: Path | None = None) -> str:
    """Every row of every check, with links to the entries behind each count."""
    body = []
    header = f"| Check | Result | &Delta; vs `{base_ref}` | |"
    separator = "| --- | ---: | ---: | :---: |"
    for title, note, records in report["groups"]:
        body += [f"### {title}"]
        if note:
            body += [note]
        body += ["", header, separator, *_render_rows(records, url_base, detail_base), ""]

    # Validation metrics share this comment: they are the same question
    # asked of the model's predictions rather than its structure, and
    # splitting them across two comments meant neither told the whole story.
    body += [
        "### Validation metrics",
        "_Predictions against measured data, to 4 significant digits. A "
        "change within tolerance still gets a checkmark; direction matters, "
        "so a rising R2 is a gain and a rising error is not._",
        "",
        f"| Metric | This branch | &Delta; vs `{base_ref}` | |",
        "| --- | ---: | ---: | :---: |",
        *_render_validation(report["validation"], url_base, detail_base),
        "",
        *_render_essential_full(report["essential"], base_ref),
    ]

    if _memote_suite(memote_dir) == "full suite":
        memote_note = (
            "_Full suite, run on this pull request. Compared against "
            f"`{base_ref}`'s own last-committed score._"
        )
    else:
        memote_note = (
            "_Fast subset -- flux-variability and matrix-rank tests are skipped, so "
            "this is not the same score `memote report snapshot` would give. Compared "
            f"against `{base_ref}`'s own last-committed score, not a fresh run on it. "
            "Comment `/run memote` on the pull request for the full suite._"
        )
    body += [
        "### MEMOTE score",
        memote_note,
        "",
        *_render_memote(memote_dir, memote_base_dir, base_ref, running),
    ]

    body += [
        "### Annotation verification (MetaNetX)",
        "_Metabolite/reaction identifiers cross-checked against MetaNetX's "
        "reference data -- see the full report below for detail._",
        "",
        *_render_annotation(annotation_dir, annotation_base_dir, base_ref,
                            running, detail_base),
    ]

    intro = (
        "_Each check name links to its explanation in the "
        f"[testResults README]({url_base}/README.md)._" if url_base else ""
    )
    footer = [_LEGEND, "", _FOOTNOTE]
    if run_url:
        footer += ["", f"[Full workflow run]({run_url})"]
    size_line = f"**Model size:** {_render_size(current, base, base_ref, running)}"
    return "\n".join(
        ["## Model quality report", "", _verdict(report, base, base_ref), "", size_line,
         "", intro, "", *body, *footer]
    ) + "\n"


_LEGEND = (
    ":x: = a count rose vs the target branch (regression) &middot; "
    ":warning: = a non-zero or worse-than-baseline finding (non-blocking) "
    "&middot; :hourglass_flowing_sand: = still running &middot; "
    ":grey_question: = the check did not run, or has no baseline yet."
)
_FOOTNOTE = (
    "_Most deltas are computed in one run -- the target branch is checked out and "
    "measured the same way -- so they reflect this pull request and nothing else. "
    "MEMOTE and the MetaNetX annotation check are the exception: compared against "
    "the target branch's own last-committed result instead._"
)


def _section_status(records: list[dict]) -> str:
    """One cell for a group of rows: the worst thing in it."""
    if any(r["state"] == "running" for r in records):
        return ":hourglass_flowing_sand: running"
    regressions = sum(r["regression"] for r in records)
    fatal = sum(r["fatal"] for r in records)
    skipped = sum(r["state"] == "skipped" for r in records)
    findings = sum(r["state"] == "measured" and r["icon"] == ":warning:" for r in records)
    if fatal:
        return f":x: **{fatal}** gate(s) failed"
    if regressions:
        return f":x: **{regressions}** regression(s)"
    parts = []
    if findings:
        parts.append(f":warning: **{findings}** pre-existing")
    if skipped:
        parts.append(f":grey_question: {skipped} not run")
    return ", ".join(parts) or ":white_check_mark: all clean"


def _list_item(record: dict, url_base: str, detail_base: str) -> str:
    delta = f" ({record['delta']})" if record["delta"] not in ("", "0") else ""
    return f"- {_row_link(record, url_base)}: {_row_value(record, detail_base)}{delta}"


def render_condensed(report: dict, current: dict, base: dict, base_ref: str, url_base: str,
                     running: bool, detail_base: str = "", summary_url: str = "",
                     memote_dir: Path | None = None,
                     memote_base_dir: Path | None = None,
                     annotation_dir: Path | None = None,
                     annotation_base_dir: Path | None = None) -> str:
    """One row per section and three lists, for the comment.

    The comment stays short enough to read in the pull request; the full tables
    are in the committed ``model_qc_summary.md``, linked from the last row.
    """
    table = ["| Section | Status |", "| --- | --- |"]
    for title, _, records in report["groups"]:
        table.append(f"| {title} ({len(records)}) | {_section_status(records)} |")
    table.append(f"| Validation metrics ({len(report['validation'])}) | "
                 f"{_section_status(report['validation'])} |")

    essential = report["essential"]
    table.append(f"| Gene essentiality | {_essential_cell(essential, base_ref, running)} |")

    parsed = _parse_memote(memote_dir)
    if running:
        memote_cell = ":hourglass_flowing_sand: running"
    elif parsed is None:
        memote_cell = ":grey_question: not run"
    else:
        base_parsed = _parse_memote(memote_base_dir) if memote_base_dir else None
        delta, icon = _memote_delta(parsed[0], base_parsed[0] if base_parsed else None)
        memote_cell = (f"{icon} **{parsed[0]:.1f}%** ({delta} vs `{base_ref}`, "
                       f"{_memote_suite(memote_dir)})")
    table.append(f"| MEMOTE | {memote_cell} |")

    annotation = read_metrics(annotation_dir, "annotation_metrics.tsv") if annotation_dir else {}
    if running:
        annotation_cell = ":hourglass_flowing_sand: running"
    elif not annotation:
        annotation_cell = ":grey_question: not run"
    else:
        wrong = int(annotation.get("annotation_wrong_isolated", 0)
                    + annotation.get("annotation_wrong_recurring", 0))
        annotation_cell = f"`wrong`: {wrong} (a report, never blocking)"
    table.append(f"| Annotation verification (MetaNetX) | {annotation_cell} |")
    if summary_url:
        table.append(f"| Full report | [model_qc_summary.md]({summary_url}) |")

    records = [r for _, _, group in report["groups"] for r in group] + report["validation"]
    measured = [r for r in records if r["state"] == "measured"]
    lists = []
    for heading, selected in (
        ("Regressions", [r for r in measured if r["regression"]]),
        ("Improved", [r for r in measured if r["improved"]]),
        ("Pre-existing", [r for r in measured if r["icon"] == ":warning:"
                          and not r["regression"] and not r["improved"]]),
    ):
        if selected:
            lists += [f"**{heading}**", "",
                      *[_list_item(r, url_base, detail_base) for r in selected], ""]

    if essential and essential["changed"] and not running:
        shown = essential["changed"][:_MAX_LISTED]
        more = len(essential["changed"]) - len(shown)
        lists += [f"**Gene essentiality, changed since `{base_ref}`**", "",
                  *[_change_line(row) for row in shown]]
        if more:
            lists += [f"- and {more} more, see the "
                      f"[full report]({summary_url}#gene-essentiality)" if summary_url
                      else f"- and {more} more in the full report"]
        lists += [""]

    size_line = f"**Model size:** {_render_size(current, base, base_ref, running)}"
    return "\n".join(
        ["## Model quality report", "", _verdict(report, base, base_ref), "", size_line, "",
         *table, "", *lists, _LEGEND]
    ) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--current", type=Path, required=True)
    parser.add_argument("--base", type=Path, default=None)
    parser.add_argument("--base-ref", default="the target branch")
    parser.add_argument("--url-base", default="",
                        help="base URL for the testResults README anchors")
    parser.add_argument("--detail-base", default="",
                        help="base URL for the findings files; when empty "
                             "the counts are plain text")
    parser.add_argument("--gene-names", type=Path,
                        default=Path(__file__).resolve().parents[4] / "data" / "databases"
                        / "SGDgeneNames.tsv",
                        help="systematic and standard gene names, for the essentiality lists")
    parser.add_argument("--run-url", default="")
    parser.add_argument("--running", action="store_true",
                        help="render every row as still running")
    parser.add_argument("--out", type=Path, required=True,
                        help="the condensed report, for the pull-request comment")
    parser.add_argument("--summary-out", type=Path, default=None,
                        help="the full report, committed next to the other results")
    args = parser.parse_args()

    current = read_metrics(args.current)
    base = read_metrics(args.base) if args.base else {}
    validation = read_metrics(args.current, "validation_metrics.tsv")
    validation_base = (
        read_metrics(args.base, "validation_metrics.tsv") if args.base else {}
    )
    status = read_status(args.current)
    essential = essential_summary(read_calls(args.current), read_calls(args.base),
                                  read_gene_names(args.gene_names))
    report = collect(current, base, args.base_ref, args.running, validation,
                     validation_base, status, essential)
    url_base, detail_base = args.url_base.rstrip("/"), args.detail_base.rstrip("/")
    directories = dict(memote_dir=args.current, memote_base_dir=args.base,
                       annotation_dir=args.current, annotation_base_dir=args.base)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    summary_url = f"{detail_base}/model_qc_summary.md" if detail_base else ""
    args.out.write_text(
        render_condensed(report, current, base, args.base_ref, url_base, args.running,
                         detail_base, summary_url, **directories)
        + (f"\n[Full workflow run]({args.run_url})\n" if args.run_url else ""),
        encoding="utf-8")
    print(f"Wrote {args.out}")
    if args.summary_out:
        # No run URL: the committed file must not change between two runs on the same model.
        args.summary_out.parent.mkdir(parents=True, exist_ok=True)
        args.summary_out.write_text(
            render_full(report, current, base, args.base_ref, url_base, args.running, "",
                        detail_base, **directories), encoding="utf-8")
        print(f"Wrote {args.summary_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
