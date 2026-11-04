"""The glue steps of the CI pipeline, in Python.

The workflow used to carry an awk one-liner for the gate, ``cp`` for
collecting results, a shell heredoc for the Gurobi licence and a block of
JavaScript for the comment. None of that could be run or tested outside
GitHub Actions, and each piece failed in its own way: the licence step
wrote a file nobody checked, the gate parsed a TSV with awk, and the
comment logic existed only inside a YAML string.

Every step is a subcommand here instead, so the whole pipeline is one
language, can be run locally, and fails loudly rather than silently.

Subcommands:
    licence    write the Gurobi WLS licence from the environment
    collect    copy the generated results into data/testResults
    comment    create or update the sticky pull-request comment
    stamp      record the pull request that last updated each group of results
    require    fail unless a job actually produced output
    gate       apply the merge gate to a metrics file
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

_API = "https://api.github.com"

# Files copied into the repository from a run's results. Named explicitly
# rather than globbed: a glob would sweep up whatever a future check
# happens to write, and the point of tracking these is that the set is
# predictable.
_TRACKED = (
    "qc_findings.md",
    "qc_metrics.tsv",
    "validation_findings.md",
    "validation_metrics.tsv",
    "essentialGenes.tsv",
    "memote_score.md",
    "annotation_report.md",
    "annotation_metrics.tsv",
    "qc_status.tsv",
)


def read_metrics(path: Path) -> dict[str, float]:
    """Read a ``key<TAB>value`` file. Missing file means no metrics."""
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


def cmd_licence(args) -> int:
    """Write the Gurobi WLS licence, if one is configured.

    Absence is not an error: secrets are unavailable to pull requests from
    forks, and cobra falls back to GLPK. Returns 0 either way, printing
    which case applied so the log says why a run used the solver it did.
    """
    licence = Path.home() / "gurobi.lic"
    file_form = os.environ.get("LICENCE_FILE", "").strip()
    fields = {
        name: os.environ.get(name, "").strip()
        for name in ("WLSACCESSID", "WLSSECRET", "LICENSEID")
    }

    if file_form:
        licence.write_text(file_form + "\n", encoding="utf-8")
    elif all(fields.values()):
        licence.write_text(
            "".join(f"{k}={v}\n" for k, v in fields.items()), encoding="utf-8"
        )
    else:
        print("No Gurobi licence configured; cobrapy will use GLPK.")
        return 0

    licence.chmod(0o600)
    # A named-user academic licence cannot authenticate on a hosted runner,
    # which otherwise looks like a solver bug rather than a licence one.
    if "WLSACCESSID" not in licence.read_text(encoding="utf-8"):
        print("::warning::the configured licence is not a WLS licence, so it "
              "will not authenticate on a hosted runner")
    else:
        print("Gurobi WLS licence installed")

    if args.github_env:
        with open(args.github_env, "a", encoding="utf-8") as handle:
            handle.write(f"GRB_LICENSE_FILE={licence}\n")
    return 0


def cmd_collect(args) -> int:
    """Copy the generated results into the repository.

    Reports what it copied and what was absent. A missing file is not
    fatal -- a job may have failed and the report says so -- but it is
    always named, because a result quietly not being updated is the thing
    that makes a tracked file untrustworthy.
    """
    args.dest.mkdir(parents=True, exist_ok=True)
    copied, missing = [], []
    for name in _TRACKED:
        source = args.source / name
        if source.is_file():
            shutil.copy2(source, args.dest / name)
            copied.append(name)
        else:
            missing.append(name)
    print(f"Copied: {', '.join(copied) or 'nothing'}")
    if missing:
        print(f"::warning::not produced by this run: {', '.join(missing)}")
    return 0


def _request(method: str, url: str, token: str, payload=None, attempts: int = 5):
    """One GitHub API call, retried with exponential backoff on a transient error.

    The API returns an occasional 5xx, or a 429 from the secondary rate limit.
    The run has finished its real work by the time it posts, so a blip there
    must not fail it.
    """
    data = json.dumps(payload).encode() if payload is not None else None
    for attempt in range(1, attempts + 1):
        request = urllib.request.Request(url, data=data, method=method)
        request.add_header("Authorization", f"Bearer {token}")
        request.add_header("Accept", "application/vnd.github+json")
        request.add_header("X-GitHub-Api-Version", "2022-11-28")
        if data is not None:
            request.add_header("Content-Type", "application/json")
        try:
            with urllib.request.urlopen(request) as response:
                return json.loads(response.read() or "null")
        except urllib.error.HTTPError as error:
            if (error.code != 429 and error.code < 500) or attempt == attempts:
                raise
            print(f"::warning::GitHub API returned {error.code}; retrying "
                  f"({attempt}/{attempts - 1})")
        except urllib.error.URLError as error:
            if attempt == attempts:
                raise
            print(f"::warning::GitHub API unreachable ({error.reason}); retrying "
                  f"({attempt}/{attempts - 1})")
        time.sleep(2 ** attempt)


def cmd_comment(args) -> int:
    """Create the sticky comment, or update it if it already exists.

    Identified by an HTML marker rather than by author or position, so it
    survives other comments being added and there is never more than one.
    """
    token = os.environ.get("GITHUB_TOKEN", "")
    if not token:
        print("::error::GITHUB_TOKEN is not set; cannot post the comment")
        return 1

    marker = f"<!-- {args.marker} -->"
    body = marker + "\n" + args.body.read_text(encoding="utf-8")
    base = f"{_API}/repos/{args.repo}/issues/{args.issue}/comments"

    try:
        existing = None
        page = 1
        while True:
            batch = _request("GET", f"{base}?per_page=100&page={page}", token)
            if not batch:
                break
            for comment in batch:
                if marker in (comment.get("body") or ""):
                    existing = comment["id"]
                    break
            if existing or len(batch) < 100:
                break
            page += 1

        if existing:
            _request(
                "PATCH",
                f"{_API}/repos/{args.repo}/issues/comments/{existing}",
                token, {"body": body},
            )
            print(f"Updated comment {existing}")
        else:
            created = _request("POST", base, token, {"body": body})
            print(f"Created comment {created['id']}")
    except urllib.error.HTTPError as error:
        print(f"::error::GitHub API returned {error.code}: {error.reason}")
        return 1
    except urllib.error.URLError as error:
        print(f"::error::GitHub API unreachable: {error.reason}")
        return 1
    return 0


def cmd_stamp(args) -> int:
    """Write the pull request number into the README's "last updated by" table.

    A table cell is marked ``<!-- last:GROUP -->**PR #N**``; the cells of the
    named groups get this pull request's number.
    """
    text = args.readme.read_text(encoding="utf-8")
    for group in args.group:
        pattern = re.compile(rf"(<!-- last:{re.escape(group)} -->\*\*PR #)\d+(\*\*)")
        text, count = pattern.subn(rf"\g<1>{args.pr}\g<2>", text)
        if not count:
            print(f"::warning::no 'last updated by' cell for group {group!r} in {args.readme}")
    args.readme.write_text(text, encoding="utf-8")
    return 0


def cmd_gate(args) -> int:
    """Fail the build only on the conditions that block a merge.

    Gates: growth; reaction/metabolite name agreement between
    model/yeast-GEM.yml and the annotation tsvs; id agreement between the model
    and the tsvs, with no deprecated id in active use (yeast-GEM#379); no
    duplicate keys in the model file; no reaction without metabolites; and, from
    ``--status``, the YAML and SBML round-trips and the YAML lint. Everything
    else in the report is a finding to review, not a reason to stop, so this
    deliberately reads specific values rather than looking at the report as a
    whole.
    """
    metrics = read_metrics(args.metrics)
    if not metrics:
        print(f"::error::{args.metrics} is missing or empty, so the checks "
              "did not produce a result")
        return 1

    failed = False

    if "growth" not in metrics:
        print("::error::the checks produced no growth value")
        failed = True
    else:
        growth = metrics["growth"]
        if growth <= 1e-6:
            print(f"::error::the model cannot produce biomass (growth = {growth:g})")
            failed = True
        else:
            print(f"Growth {growth:g}")

    # (metric, what a non-zero value means, where the detail is)
    zero_gates = (
        ("name_mismatches", "reaction/metabolite name(s) disagree between "
         "model/yeast-GEM.yml and the annotation tsvs"),
        ("annotation_consistency", "model/annotation-table inconsistenc(y/ies)"),
        ("duplicate_keys", "duplicate key(s) in model/yeast-GEM.yml"),
        ("empty_reactions", "reaction(s) without metabolites"),
    )
    for key, meaning in zero_gates:
        value = metrics.get(key, 0)
        if value > 0:
            print(f"::error::{value:g} {meaning} -- see qc_findings.md")
            failed = True
        else:
            print(f"{key}: 0")

    if args.status:
        status = {}
        if args.status.is_file():
            for line in args.status.read_text(encoding="utf-8").splitlines()[1:]:
                key, _, value = line.partition("\t")
                status[key.strip()] = value.strip()
        for key in ("roundtrip_yaml", "roundtrip_sbml", "yamllint"):
            if status.get(key) == "pass":
                print(f"{key}: pass")
            else:
                print(f"::error::{key}: {status.get(key, 'not run')}")
                failed = True

    return 1 if failed else 0


def cmd_require(args) -> int:
    """Fail unless at least one of the named files has content.

    The measuring steps are continue-on-error so that one failing side
    still uploads the other, which means the job's own status says nothing
    about whether it produced anything. This is the check that it did.
    """
    produced = [p for p in args.any if p.is_file() and p.stat().st_size > 0]
    if produced:
        print(f"Produced: {', '.join(str(p) for p in produced)}")
        return 0
    print("::error::none of these were produced: "
          + ", ".join(str(p) for p in args.any))
    return 1


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("licence", help="write the Gurobi WLS licence")
    p.add_argument("--github-env", default=os.environ.get("GITHUB_ENV"))
    p.set_defaults(func=cmd_licence)

    p = sub.add_parser("collect", help="copy results into the repository")
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--dest", type=Path, default=Path("data/testResults"))
    p.set_defaults(func=cmd_collect)

    p = sub.add_parser("comment", help="create or update the sticky comment")
    p.add_argument("--repo", required=True)
    p.add_argument("--issue", required=True)
    p.add_argument("--body", type=Path, required=True)
    p.add_argument("--marker", default="MODEL-QC-COMMENT")
    p.set_defaults(func=cmd_comment)

    p = sub.add_parser("require", help="fail unless a run produced output")
    p.add_argument("--any", type=Path, nargs="+", required=True)
    p.set_defaults(func=cmd_require)

    p = sub.add_parser("stamp", help="record the pull request in the results README")
    p.add_argument("--readme", type=Path, default=Path("data/testResults/README.md"))
    p.add_argument("--pr", required=True)
    p.add_argument("--group", nargs="+", required=True)
    p.set_defaults(func=cmd_stamp)

    p = sub.add_parser("gate", help="apply the merge gate")
    p.add_argument("--metrics", type=Path, required=True)
    p.add_argument("--status", type=Path, default=None,
                   help="qc_status.tsv, for the model-file gates")
    p.set_defaults(func=cmd_gate)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
