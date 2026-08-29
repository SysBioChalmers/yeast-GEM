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
* every count -- :x: if it rose against the target branch (a regression),
  :warning: if it is non-zero but no worse (a pre-existing finding, not
  blocking), :white_check_mark: if it is zero.
* :hourglass_flowing_sand: -- the group is still being computed on this run.

Only the growth gate fails the build. Everything else is reported.

Usage:
    python build_qc_report.py --current data/testResults --base /tmp/base \\
        --base-ref develop --url-base https://github.com/.../data/testResults
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

# (metric key, comment label, kind)
_MODEL_ROWS = [
    ("growth", "Growth (biomass producible)", "growth"),
    ("duplicate_reactions", "Exact-duplicate reaction groups", "count"),
    ("orphan_metabolites", "Unused metabolites", "count"),
    ("unused_genes", "Unused genes", "count"),
    ("missing_formula", "Metabolites missing formula", "count"),
    ("missing_charge", "Metabolites missing charge", "count"),
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
]


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

    if base is None:
        return "new", (":warning:" if value else ":white_check_mark:"), False, False
    change = int(value) - int(base)
    if change > 0:
        return f"+{change}", ":x:", True, False
    text = "0" if change == 0 else str(change)
    if value > 0:
        return text, ":warning:", False, False
    return text, ":white_check_mark:", False, False


def _render_rows(rows, current, base, url_base, detail_base, running):
    lines, regressions, warnings, pending, skipped, fatal = [], 0, 0, 0, 0, False
    for key, label, kind in rows:
        name = (
            f"[{label}]({url_base}/README.md#{_slug(label)})"
            if url_base else label
        )
        if running:
            lines.append(f"| {name} | _running_ | | :hourglass_flowing_sand: |")
            pending += 1
            continue
        # A check with no value once the run is over did not execute -- an
        # optional dependency missing, or the check crashed. Saying "running"
        # there would be a lie that never resolves, and counting it as clean
        # would be worse.
        if key not in current:
            lines.append(f"| {name} | _not run_ | | :grey_question: |")
            skipped += 1
            continue
        value = current[key]
        delta, icon, regression, row_fatal = _icon(value, base.get(key), kind)
        text = f"{value:.4g}" if kind == "growth" else str(int(value))
        # A non-zero count links to its section in the findings file,
        # which is addressed by the same slug as the README heading, so
        # the link follows the label automatically.
        if detail_base and kind == "count" and value:
            text = f"[{text}]({detail_base}/qc_findings.md#{_slug(label)})"
        lines.append(f"| {name} | {text} | {delta} | {icon} |")
        regressions += regression
        warnings += icon == ":warning:"
        fatal = fatal or row_fatal
    return lines, regressions, warnings, pending, skipped, fatal


def _fmt_delta(change: float) -> str:
    """4 significant digits, or a bare ``0`` for an exact (non-)change.

    Tolerance decides whether a move counts as a *regression* -- it must
    not also decide whether the reader gets to see the move at all. A
    change of 0.0008 sitting under a 0.005 tolerance is still a real
    change; showing it as literally ``0`` reads as "nothing moved" when
    something did, just not enough to matter for the verdict.
    """
    return f"{change:+.4g}" if abs(change) > 1e-9 else "0"


def _render_validation(current, base, url_base, base_ref, detail_base=""):
    """The validation-metrics table, or a note if the job did not report."""
    if not current:
        return ["_not run — the validation metrics job reported nothing._"], 0, 1
    lines, regressions, skipped = [], 0, 0
    for key, label, direction, tol in _VALIDATION_ROWS:
        name = (
            f"[{label}]({url_base}/README.md#{_slug(label)})"
            if url_base else label
        )
        if key not in current:
            lines.append(f"| {name} | _not run_ | | | :grey_question: |")
            skipped += 1
            continue
        value = current[key]
        # The value links to the measurements it summarises, the same way a
        # QC count links to the entries behind it.
        shown = f"{value:.4g}"
        if detail_base:
            anchor = _slug(_DETAIL_ANCHOR.get(key, label))
            shown = (f"[{shown}]({detail_base}/validation_findings.md"
                     f"#{anchor})")
        if key not in base:
            lines.append(f"| {name} | {shown} | — | new | :grey_question: |")
            skipped += 1
            continue
        change = value - base[key]
        delta = _fmt_delta(change)
        if abs(change) <= tol:
            icon = (":warning:" if key in _ZERO_IS_IDEAL and value > 0
                    else ":white_check_mark:")
        elif direction == "none":
            icon = ":warning:"
        else:
            improved = (change > 0) if direction == "higher" else (change < 0)
            if not improved:
                icon = ":x:"
            elif key in _ZERO_IS_IDEAL and value > 0:
                icon = ":warning:"
            else:
                icon = ":white_check_mark:"
            regressions += not improved
        lines.append(
            f"| {name} | {shown} | {base[key]:.4g} | {delta} | {icon} |"
        )
    return lines, regressions, skipped


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


def _memote_delta(current: float, base: float | None) -> str:
    """Percentage-point delta, MEMOTE's own scale -- not the 4g convention
    used for validation metrics, since a score is already a rounded percent."""
    if base is None:
        return ""
    change = current - base
    if abs(change) < 0.05:
        return "0"
    return f"{change:+.1f} {':warning:' if change < 0 else ':white_check_mark:'}"


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
    if running:
        return ["_running_ &middot; :hourglass_flowing_sand:", ""]
    parsed = _parse_memote(current)
    if parsed is None:
        return ["_not run — the MEMOTE job reported nothing._", ""]
    total, sections, detailed = parsed
    base_parsed = _parse_memote(base) if base is not None else None
    base_total, base_sections, _base_detailed = base_parsed or (None, {}, [])

    lines = [
        f"**Total score: {total:.1f}%** &nbsp; "
        f"{_memote_delta(total, base_total)}".rstrip(),
        "",
    ]
    if sections:
        lines += [f"| Section | Score | &Delta; vs `{base_ref}` |",
                  "| --- | ---: | ---: |"]
        lines += [
            f"| {name} | {value:.1f}% | "
            f"{_memote_delta(value, base_sections.get(name))} |"
            for name, value in sections.items()
        ]
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


def build(current: dict, base: dict, base_ref: str, url_base: str,
          running: bool, run_url: str, detail_base: str = "",
          validation: dict | None = None,
          validation_base: dict | None = None,
          memote_dir: Path | None = None,
          memote_base_dir: Path | None = None) -> str:
    groups = [
        ("Model checks",
         "_Growth is a gate and blocks the merge; every other row is a "
         "non-blocking report._",
         _MODEL_ROWS),
        ("Mass/charge balance and MACAW", "", _BALANCE_ROWS),
        ("Annotations", "", _ANNOTATION_ROWS),
    ]

    body, regressions, warnings, pending, skipped, fatal = [], 0, 0, 0, 0, False
    header = f"| Check | Result | &Delta; vs `{base_ref}` | |"
    separator = "| --- | ---: | ---: | :---: |"
    for title, note, rows in groups:
        lines, reg, warn, pend, skip, fat = _render_rows(
            rows, current, base, url_base, detail_base, running
        )
        regressions += reg
        warnings += warn
        pending += pend
        skipped += skip
        fatal = fatal or fat
        body += [f"### {title}"]
        if note:
            body += [note]
        body += ["", header, separator, *lines, ""]

    # Validation metrics share this comment: they are the same question
    # asked of the model's predictions rather than its structure, and
    # splitting them across two comments meant neither told the whole story.
    if running:
        val_lines = ["| _running_ | | | | :hourglass_flowing_sand: |"]
        val_reg, val_skip = 0, 0
    else:
        val_lines, val_reg, val_skip = _render_validation(
            validation or {}, validation_base or {}, url_base, base_ref,
            detail_base,
        )
    regressions += val_reg
    skipped += val_skip
    body += [
        "### Validation metrics",
        "_Predictions against measured data, to 4 significant digits. A "
        "change within tolerance still gets a checkmark; direction matters, "
        "so a rising R2 is a gain and a rising error is not._",
        "",
        f"| Metric | This branch | `{base_ref}` | &Delta; | |",
        "| --- | ---: | ---: | ---: | :---: |",
        *val_lines,
        "",
    ]

    body += [
        "### MEMOTE score",
        "_Fast subset -- flux-variability and matrix-rank tests are "
        "skipped, so this is not the same score `memote report snapshot` "
        f"would give. Compared against `{base_ref}`'s own last-committed "
        "score, not a fresh run on it._",
        "",
        *_render_memote(memote_dir, memote_base_dir, base_ref, running),
    ]

    if fatal:
        verdict = (
            ":x: **Merge blocked: the model cannot grow.** "
            "See the Model checks table."
        )
    elif regressions:
        extra = f" ({pending} check(s) still running)" if pending else ""
        verdict = (
            f":x: **{regressions} regression(s) vs `{base_ref}`** — this pull "
            f"request increased a finding count{extra}. Review the :x: rows."
        )
    elif pending:
        verdict = (
            f":hourglass_flowing_sand: **{pending} check(s) still running.** "
            f"The rest are unchanged vs `{base_ref}`."
        )
    elif not base:
        verdict = (
            ":information_source: First run for this comparison; no "
            "target-branch baseline yet."
        )
    elif warnings:
        verdict = (
            f":warning: **{warnings} pre-existing finding(s), no regressions "
            f"vs `{base_ref}`.** Non-blocking."
        )
    else:
        verdict = (
            f":white_check_mark: **All checks clean, no regressions vs "
            f"`{base_ref}`.**"
        )
    if skipped:
        verdict += (
            f" {skipped} check(s) did not run — their result is unknown, "
            "not clean."
        )

    intro = (
        "_Each check name links to its explanation in the "
        f"[testResults README]({url_base}/README.md)._" if url_base else ""
    )
    footer = [
        ":x: = a count rose vs the target branch (regression) &middot; "
        ":warning: = a pre-existing non-zero finding (non-blocking) &middot; "
        ":hourglass_flowing_sand: = still running &middot; "
        ":grey_question: = the check did not run.",
        "",
        "_Both columns are computed in this run — the target branch is "
        "checked out and measured the same way — so the delta reflects this "
        "pull request and nothing else._",
    ]
    if run_url:
        footer += ["", f"[Full workflow run]({run_url})"]

    return "\n".join(
        ["## Model quality report", "", verdict, "", intro, "", *body, *footer]
    ) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--current", type=Path, required=True)
    parser.add_argument("--base", type=Path, default=None)
    parser.add_argument("--base-ref", default="the target branch")
    parser.add_argument("--url-base", default="",
                        help="base URL for the testResults README anchors")
    parser.add_argument("--detail-base", default="",
                        help="base URL for the per-check CSVs; when empty "
                             "the counts are plain text")
    parser.add_argument("--run-url", default="")
    parser.add_argument("--running", action="store_true",
                        help="render every row as still running")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    current = read_metrics(args.current)
    base = read_metrics(args.base) if args.base else {}
    validation = read_metrics(args.current, "validation_metrics.tsv")
    validation_base = (
        read_metrics(args.base, "validation_metrics.tsv") if args.base else {}
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        build(current, base, args.base_ref, args.url_base.rstrip("/"),
              args.running, args.run_url, args.detail_base.rstrip("/"),
              validation, validation_base,
              memote_dir=args.current, memote_base_dir=args.base),
        encoding="utf-8",
    )
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
