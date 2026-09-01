"""Python port of ``code/increaseVersion.m``, run on a ``release/X.Y.Z``
branch in CI.

Covers everything MATLAB does that has a Python equivalent: bump
validation, running the model tests, and refreshing every text artifact
those tests feed -- ``README.md``, ``data/testResults/README.md``,
``data/testResults/growth.md`` + ``growth.png``,
``data/testResults/essentialGenes.tsv``.

What this does **not** do: write ``model/yeast-GEM.xml``, ``.txt``,
``.xlsx`` or ``.mat``. Those four are RAVEN's ``exportForGit`` output,
derived *from* ``model/yeast-GEM.yml``. ``raven_toolbox.io.write_yaml_model``
can write the ``.yml`` itself now (used by ``yeastgem.save_yeast_yaml``),
but ``_stamp_yaml`` below still edits the file's text directly rather than
loading the model and rewriting it whole: a full rewrite would touch
every line's formatting incidentally, while a release commit should only
ever change the three metadata lines it is actually stamping. There is
still no Python writer for RAVEN's tab-separated ``.txt`` or multi-sheet
``.xlsx`` layout -- those, and ``.xml``, stay a MATLAB step downstream,
which runs ``exportForGit`` against the now-stamped ``.yml`` to produce
the four derived files. See ``.github/workflows/release.yml``.

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

from yeastgem import conditions, load_yeast_yaml, model_tests
from yeastgem.io import REPO_PATH

_VERSION_TXT = REPO_PATH / "version.txt"
_HISTORY_MD = REPO_PATH / "history.md"
_YEAST_YML = REPO_PATH / "model" / "yeast-GEM.yml"
_ROOT_README = REPO_PATH / "README.md"
_RESULTS_README = REPO_PATH / "data" / "testResults" / "README.md"

# Matches root README.md's whole "# Model overview" section, from the
# heading up to (not including) the next top-level heading. On develop
# this is a static placeholder -- see the section itself -- so unlike
# the row-only rewrite this used to be, a release has to build the
# section from scratch rather than edit an existing table in place.
_MODEL_OVERVIEW_RE = re.compile(r"^# Model overview\n.*?(?=^# )", re.M | re.S)

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


def _current_version_text() -> str:
    """The currently-released version.

    version.txt lives on main only (develop does not track it -- see
    branch-hygiene.yml), but validate runs from a develop checkout when
    cutting a release. Read the local file if this happens to be run on
    a checkout that has one (main, or an already-cut release branch);
    otherwise fetch it from origin/main.
    """
    if _VERSION_TXT.is_file():
        return _VERSION_TXT.read_text(encoding="utf-8")
    subprocess.run(
        ["git", "fetch", "--depth=1", "origin", "main"], cwd=REPO_PATH,
        capture_output=True, text=True, check=True,
    )
    result = subprocess.run(
        ["git", "show", "origin/main:version.txt"], cwd=REPO_PATH,
        capture_output=True, text=True, check=True,
    )
    return result.stdout


def cmd_validate(args: argparse.Namespace) -> int:
    old = _parse_version(_current_version_text())
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


def _populate_root_readme_overview(version: str, n_rxns: int, n_mets: int,
                                    n_genes: int, accuracy, tp, tn, fp, fn,
                                    r2) -> None:
    """Replace root README.md's "# Model overview" placeholder with the
    full, stamped section: main is a specific released version, so
    unlike develop these numbers do not go stale between curations.
    """
    date = datetime.now(UTC).strftime("%d-%b-%Y")
    block = f"""# Model overview

| Taxonomy | Latest update | Version | Reactions | Metabolites | Genes |
|:-------|:--------------|:------|:------|:----------|:-----|
| _Saccharomyces cerevisiae_ | {date} | {version} | {n_rxns} | {n_mets} | {n_genes} |

Validated against experimental data on every pull request; see
[data/testResults/README.md](data/testResults/README.md) for the full
methodology and additional metrics, including anaerobic flux and
fermentation-product predictions.

### Gene essentiality prediction

Predicted gene essentiality vs. the Stanford yeast deletion collection
(1107 genes):

- Accuracy: {accuracy:.3f}
- True non-essential genes: {len(tp)}
- True essential genes: {len(tn)}
- False non-essential genes: {len(fp)}
- False essential genes: {len(fn)}

### Growth prediction

Predicted vs. measured growth rate across chemostat conditions from
[Österlund _et al._ (2013)](https://doi.org/10.1186/1752-0509-7-36):

- Correlation coefficient R<sup>2</sup>: {r2:.3f}

![Growth curve](data/testResults/growth.png)

"""
    text = _ROOT_README.read_text(encoding="utf-8")
    new_text, count = _MODEL_OVERVIEW_RE.subn(block, text, count=1)
    if count != 1:
        raise SystemExit(
            "expected exactly one '# Model overview' section in "
            f"{_ROOT_README}, found {count}"
        )
    _ROOT_README.write_text(new_text, encoding="utf-8")


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
    model = load_yeast_yaml()
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
    _populate_root_readme_overview(
        args.version, len(model.reactions), len(model.metabolites),
        len(model.genes), essential.accuracy, essential.tp, essential.tn,
        essential.fp, essential.fn, r2,
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
