"""Python port of ``code/increaseVersion.m``, run on a ``release/X.Y.Z``
branch in CI.

Covers everything MATLAB does that has a Python equivalent: bump
validation, running the model tests, and refreshing every text artifact
those tests feed -- ``README.md``, ``data/testResults/README.md``,
``data/testResults/growth.md`` + ``growth.png``,
``data/testResults/essentialGenes.tsv``.

What this does **not** do: write ``model/yeast-GEM.xml``, ``.txt``,
``.xlsx`` or ``.mat``. Those four are RAVEN's ``exportForGit`` output,
derived *from* ``model/yeast-GEM.yml`` -- not from anything Python
writes -- and there is no Python writer for RAVEN's ``!!omap`` YAML, its
tab-separated ``.txt``, or its multi-sheet ``.xlsx`` layout. Building one
was explicitly out of scope here: getting a binary/structured export
format wrong is a much more expensive mistake than getting a markdown
table wrong, and there is no way to verify a hand-rolled writer's output
against RAVEN in this environment. This script stamps the three metadata
lines RAVEN's export reads (``id``, ``version``, ``date``) and stops; a
MATLAB step downstream runs ``exportForGit`` against the now-stamped
``.yml`` to produce the four derived files. See ``.github/workflows/release.yml``.

Two subcommands:

``validate``
    Check that ``--version`` is a legal bump from ``version.txt`` and
    that ``history.md`` already has a heading for it. No side effects;
    run before the release branch is cut, so a bad version number fails
    before anything is written.

``build``
    Run on the release branch. Stamps ``model/yeast-GEM.yml``, runs the
    tests, refreshes every file above, and guards that nothing else
    changed.

Run locally:
    python code/python/release/increase_version.py validate --version 9.1.1
    python code/python/release/increase_version.py build --version 9.1.1
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from datetime import UTC, datetime

import cobra
import matplotlib

matplotlib.use("Agg")

from yeastgem import conditions, model_tests
from yeastgem.io import REPO_PATH, _update_readme, read_yeast_model

_VERSION_TXT = REPO_PATH / "version.txt"
_HISTORY_MD = REPO_PATH / "history.md"
_YEAST_YML = REPO_PATH / "model" / "yeast-GEM.yml"
_ROOT_README = REPO_PATH / "README.md"
_RESULTS_README = REPO_PATH / "data" / "testResults" / "README.md"

# Files this step is allowed to change. Anything else in the working tree
# after it runs means develop and this branch have diverged in a way that
# was not accounted for -- mirrors the git diff --numstat guard at the end
# of increaseVersion.m, scoped to what this Python step is responsible
# for. model/yeast-GEM.xml/.txt/.xlsx/.mat are deliberately absent: they
# belong to the MATLAB export step that runs after this one, which has
# its own guard.
_EXPECTED_CHANGED = {
    "version.txt",
    "model/yeast-GEM.yml",
    "README.md",
    "data/testResults/README.md",
    "data/testResults/growth.md",
    "data/testResults/growth.png",
    "data/testResults/essentialGenes.tsv",
}


def _parse_version(text: str) -> tuple[int, int, int]:
    parts = text.strip().split(".")
    if len(parts) != 3 or not all(p.isdigit() for p in parts):
        raise SystemExit(f"'{text}' is not a X.Y.Z version")
    return tuple(int(p) for p in parts)  # type: ignore[return-value]


def _is_legal_bump(old: tuple[int, int, int], new: tuple[int, int, int]) -> bool:
    major, minor, patch = old
    return new in {
        (major + 1, 0, 0),
        (major, minor + 1, 0),
        (major, minor, patch + 1),
    }


def cmd_validate(args: argparse.Namespace) -> int:
    old = _parse_version(_VERSION_TXT.read_text(encoding="utf-8"))
    new = _parse_version(args.version)
    if not _is_legal_bump(old, new):
        raise SystemExit(
            f"{args.version} is not a legal bump from "
            f"{'.'.join(map(str, old))} (must be exactly one major, minor "
            "or patch step)"
        )

    heading = f"### yeast {args.version}:"
    if heading not in _HISTORY_MD.read_text(encoding="utf-8"):
        raise SystemExit(
            f"history.md has no '{heading}' heading. Update history.md "
            "by hand before cutting the release branch."
        )

    print(f"{args.version} is a legal bump from "
          f"{'.'.join(map(str, old))}, and history.md is ready.")
    return 0


def _stamp_yaml(version: str) -> None:
    """Rewrite the three metadata lines RAVEN's exportForGit reads.

    Everything else in the file -- every reaction, metabolite, gene --
    is untouched text, so this cannot introduce the kind of diff a full
    model rewrite would.

    The version line is a special case, handled separately below: RAVEN's
    writeYAMLmodel.m (writeMetadata) omits it entirely whenever
    model.version is empty, which is the normal state on develop -- so a
    release branch cut straight from develop, or from any commit since
    the last ordinary save, has no such line to match on. Confirmed
    directly: round-tripping the current develop yml through
    readYAMLmodel/writeYAMLmodel with no other change already drops a
    previously-committed `- version: ""` line. When that's the case,
    insert the line fresh right after `- name:`, matching writeMetadata's
    own field order (id, name, version, date).
    """
    text = _YEAST_YML.read_text(encoding="utf-8")
    today = datetime.now(UTC).strftime("%Y-%m-%d")
    substitutions = {
        r'^  - id: .*$': f'  - id: yeastGEM_v{version}',
        r'^  - date: ".*"$': f'  - date: "{today}"',
    }
    for pattern, replacement in substitutions.items():
        text, count = re.subn(pattern, replacement, text, count=1, flags=re.M)
        if count != 1:
            raise SystemExit(
                f"expected exactly one match for {pattern!r} in "
                f"{_YEAST_YML}, found {count}"
            )

    version_line = f'  - version: "{version}"'
    text, count = re.subn(r'^  - version: ".*"$', version_line, text, count=1, flags=re.M)
    if count == 0:
        text, count = re.subn(
            r'^(  - name: .*)$', rf'\1\n{version_line}', text, count=1, flags=re.M
        )
        if count != 1:
            raise SystemExit(
                f"expected exactly one '  - name: ...' line to insert the version "
                f"line after in {_YEAST_YML}, found {count}"
            )

    _YEAST_YML.write_text(text, encoding="utf-8")


def _update_root_readme_stats(accuracy, tp, tn, fp, fn, r2) -> None:
    """The five gene-essentiality/growth lines under '## Model overview'.

    The overview table row itself is handled separately by
    :func:`yeastgem.io._update_readme`, which already exists for the
    same purpose on every curation commit.
    """
    text = _ROOT_README.read_text(encoding="utf-8")
    substitutions = [
        (r"^(- Accuracy: )0\.\d+$", rf"\g<1>{accuracy:.3f}"),
        (r"^(- True non-essential genes: )\d+$", rf"\g<1>{len(tp)}"),
        (r"^(- True essential genes: )\d+$", rf"\g<1>{len(tn)}"),
        (r"^(- False non-essential genes: )\d+$", rf"\g<1>{len(fp)}"),
        (r"^(- False essential genes: )\d+$", rf"\g<1>{len(fn)}"),
        (r"^(- Correlation coefficient R<sup>2</sup>: )0\.\d+$",
         rf"\g<1>{r2:.3f}"),
    ]
    for pattern, replacement in substitutions:
        text, count = re.subn(pattern, replacement, text, count=1, flags=re.M)
        if count != 1:
            raise SystemExit(
                f"expected exactly one match for {pattern!r} in "
                f"{_ROOT_README}, found {count}"
            )
    _ROOT_README.write_text(text, encoding="utf-8")


def _update_results_readme(version: str, r2, flux, exchange, accuracy, tp,
                            tn, fp, fn) -> None:
    """Section 2 of data/testResults/README.md: the release summary table.

    Every row is rewritten in place by label, matching
    check_metrics.py's ``_METRICS`` labels exactly, so the two never
    silently disagree about what a row is called.
    """

    solver = f"Python {sys.version.split()[0]}, cobrapy {cobra.__version__}"
    text = _RESULTS_README.read_text(encoding="utf-8")

    metadata = [
        (r"^(- Model version: )[^\n]*$", rf"\g<1>{version}"),
        (r"^(- Software: )[^\n]*$", rf"\g<1>{solver}"),
    ]
    values = [
        ("Growth prediction R2", f"{r2:.16g}"),
        ("Anaerobic flux median fold error", f"{flux.median_fold_error:.16g}"),
        ("Anaerobic exchange mean relative error",
         f"{exchange.mean_relative_error:.16g}"),
        ("Ammonium per ATPase", f"{exchange.ammonium_per_atpase:.16g}"),
        ("Gene essentiality accuracy", f"{accuracy:.16g}"),
        ("True non-essential genes", str(len(tp))),
        ("True essential genes", str(len(tn))),
        ("False non-essential genes", str(len(fp))),
        ("False essential genes", str(len(fn))),
    ]
    for pattern, replacement in metadata:
        text, count = re.subn(pattern, replacement, text, count=1, flags=re.M)
        if count != 1:
            raise SystemExit(f"expected one match for {pattern!r}, found {count}")
    for label, value in values:
        pattern = rf"^(\| {re.escape(label)} \| )[\d.eE+-]+ \|$"
        text, count = re.subn(pattern, rf"\g<1>{value} |", text, count=1,
                               flags=re.M)
        if count != 1:
            raise SystemExit(
                f"expected one match for row {label!r}, found {count}. "
                "The label here and in check_metrics.py's _METRICS must "
                "agree."
            )

    # The per-product exchange table: the whole row is rewritten, not
    # just the number, so the row and the prediction cannot fall out of
    # step -- mirrors increaseVersion.m.
    for _i, row in exchange.results.iterrows():
        label = re.sub(r" (exchange|pseudoreaction)$", "", row["rxnName"])
        within = "yes" if row["withinError"] else "no"
        pattern = rf"^\| {re.escape(label)} \|[^\n]*$"
        replacement = (
            f"| {label} | {row['measured']:g} +/- {row['stdev']:g} | "
            f"{row['predicted']:.4f} | {within} |"
        )
        text, count = re.subn(pattern, replacement, text, count=1, flags=re.M)
        if count != 1:
            raise SystemExit(
                f"expected one match for the '{label}' exchange row, "
                f"found {count}"
            )
    _RESULTS_README.write_text(text, encoding="utf-8")


def _check_only_expected_changed() -> None:
    result = subprocess.run(
        ["git", "diff", "--numstat"], cwd=REPO_PATH,
        capture_output=True, text=True, check=True,
    )
    unexpected = []
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) != 3:
            continue
        path = parts[2].replace("\\", "/")
        if path not in _EXPECTED_CHANGED:
            unexpected.append(path)
    if unexpected:
        raise SystemExit(
            "Unexpected files changed: " + ", ".join(unexpected) + ". "
            "This usually means develop moved after the release branch "
            "was cut. Update develop, re-cut the branch, and try again."
        )
    print(f"Only expected files changed: {', '.join(sorted(_EXPECTED_CHANGED))}")


def cmd_build(args: argparse.Namespace) -> int:
    print("Loading model")
    model = read_yeast_model()
    model.id = f"yeastGEM_v{args.version}"

    print("Stamping model/yeast-GEM.yml")
    _stamp_yaml(args.version)

    print("Running gene essentiality analysis")
    essential = model_tests.essential_genes(model.copy(), write_output=True)

    print("Running growth analysis")
    r2 = model_tests.growth(model.copy(), write_output=True)

    print("Running anaerobic flux analysis")
    anaerobic = model.copy()
    conditions.apply(anaerobic, "anaerobic")
    flux = model_tests.anaerobic_flux_predictions(anaerobic)

    print("Running anaerobic exchange rate analysis")
    exchange = model_tests.plot_anaerobic(anaerobic, plot=False)

    print("Updating README.md")
    _update_readme(model)
    _update_root_readme_stats(
        essential.accuracy, essential.tp, essential.tn, essential.fp,
        essential.fn, r2,
    )

    print("Updating data/testResults/README.md")
    _update_results_readme(
        args.version, r2, flux, exchange, essential.accuracy, essential.tp,
        essential.tn, essential.fp, essential.fn,
    )

    _check_only_expected_changed()

    print(f"Bumping version.txt to {args.version}")
    _VERSION_TXT.write_text(args.version, encoding="utf-8")

    print(f"\nBuild complete for {args.version}. Remaining before this can "
          "be released: RAVEN's exportForGit against model/yeast-GEM.yml "
          "(writes .xml/.txt/.xlsx/.mat), then a memote snapshot.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("validate", help="check the version bump and history.md")
    p.add_argument("--version", required=True)
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("build", help="stamp, test, and refresh the release files")
    p.add_argument("--version", required=True)
    p.set_defaults(func=cmd_build)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
