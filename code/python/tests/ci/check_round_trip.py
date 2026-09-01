"""Level-1 parity gate — Python SBML round-trip preserves the model.

Loads model/yeast-GEM.yml the same way curation does (load_yeast_yaml:
tsv cross-reference annotation merged in, deltaG restored), writes it out
as SBML to a temp file, reads that back, and asserts the two
``cobra.Model`` objects are semantically equal (delegated to
``raven_toolbox.comparison.diff_models``). Catches SBML library
regressions, annotation losses, and accidental ID rewrites -- for
whatever a curator's model actually contains, not just what already
happened to survive a previous xml round-trip.

Run locally:
    python code/python/tests/ci/check_round_trip.py
"""
from __future__ import annotations

import tempfile
from pathlib import Path

from cobra.io import read_sbml_model, write_sbml_model
from raven_toolbox.comparison import diff_models

from yeastgem import load_yeast_yaml


def main() -> int:
    model = load_yeast_yaml()
    with tempfile.TemporaryDirectory() as tmp:
        round_trip_path = Path(tmp) / "yeast-GEM.xml"
        write_sbml_model(model, str(round_trip_path))
        reloaded = read_sbml_model(str(round_trip_path))
    report = diff_models(model, reloaded)
    print(report)
    return 0 if report.equal else 1


if __name__ == "__main__":
    raise SystemExit(main())
