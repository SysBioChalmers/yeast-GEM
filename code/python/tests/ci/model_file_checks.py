"""Checks of the model file itself: YAML lint, YAML round-trip and SBML round-trip.

Each check records ``pass`` or ``fail`` under its key in ``qc_status.tsv`` in the
output directory (a ``check<TAB>result`` table, one row per key), which
``build_qc_report.py`` renders and ``pipeline.py gate`` applies:

* ``yamllint``: ``model/`` against the rules of ``yaml-validation.yml``;
* ``roundtrip_yaml``: ``model/yeast-GEM.yml`` read, written and read again gives a
  semantically equal model (``raven_toolbox.comparison.diff_models``);
* ``roundtrip_sbml``: the same through SBML, as in ``check_round_trip.py``.

A check that raises is recorded as ``fail`` with the error in the log, so the
status file always has all three keys.

Usage:
    python code/python/tests/ci/model_file_checks.py --out results/head
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
import traceback
from pathlib import Path

# Same configuration as .github/workflows/yaml-validation.yml: RAVEN and
# raven-toolbox write the model with sequences at the indentation of their key
# (yeast-GEM#397), which the default indent-sequences rule rejects.
YAMLLINT_CONFIG = ("{extends: default, rules: {line-length: disable, "
                   "indentation: {indent-sequences: false}}}")


def check_yamllint(repo: Path) -> bool:
    result = subprocess.run(["yamllint", "-d", YAMLLINT_CONFIG, str(repo / "model")],
                            capture_output=True, text=True, check=False)
    print(result.stdout[-2000:] + result.stderr[-2000:])
    return result.returncode == 0


def check_roundtrip_yaml(repo: Path) -> bool:
    from raven_toolbox.comparison import diff_models
    from raven_toolbox.io import read_yaml_model, write_yaml_model

    original = read_yaml_model(str(repo / "model" / "yeast-GEM.yml"))
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "yeast-GEM.yml"
        write_yaml_model(original, path)
        reloaded = read_yaml_model(str(path))
    report = diff_models(original, reloaded)
    print(report)
    return bool(report.equal)


def check_roundtrip_sbml(repo: Path) -> bool:
    from check_round_trip import main as round_trip

    return round_trip() == 0


CHECKS = (
    ("roundtrip_sbml", check_roundtrip_sbml),
    ("roundtrip_yaml", check_roundtrip_yaml),
    ("yamllint", check_yamllint),
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("."))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    status = {}
    for key, check in CHECKS:
        print(f"== {key}", flush=True)
        try:
            status[key] = "pass" if check(args.repo) else "fail"
        except Exception:  # a crashed check is a failed check
            traceback.print_exc()
            status[key] = "fail"
        print(f"{key}: {status[key]}", flush=True)

    args.out.mkdir(parents=True, exist_ok=True)
    lines = ["check\tresult"] + [f"{key}\t{status[key]}" for key in sorted(status)]
    (args.out / "qc_status.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0 if all(value == "pass" for value in status.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
