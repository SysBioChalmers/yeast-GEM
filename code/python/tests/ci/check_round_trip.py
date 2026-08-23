"""Level-1 parity gate — Python SBML round-trip preserves the model.

Loads the committed ``model/yeast-GEM.xml`` with cobrapy, writes it
back out to a temp file, reads that, and asserts the two
``cobra.Model`` objects are semantically equal (delegated to
``raven_toolbox.comparison.diff_models``). Catches SBML library
regressions, annotation losses, and accidental ID rewrites.

Run locally:
    python code/python/tests/ci/check_round_trip.py
"""
from __future__ import annotations

import tempfile
from pathlib import Path

from cobra.io import read_sbml_model, write_sbml_model
from raven_toolbox.comparison import diff_models

from yeastgem import MODEL_PATH


def main() -> int:
    model = read_sbml_model(str(MODEL_PATH))
    with tempfile.TemporaryDirectory() as tmp:
        round_trip_path = Path(tmp) / "yeast-GEM.xml"
        write_sbml_model(model, str(round_trip_path))
        reloaded = read_sbml_model(str(round_trip_path))
    report = diff_models(model, reloaded)
    print(report)
    return 0 if report.equal else 1


if __name__ == "__main__":
    raise SystemExit(main())
