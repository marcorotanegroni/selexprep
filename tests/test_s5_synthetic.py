"""Supplementary S5 is generated from the semi-synthetic run tables."""

from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "benchmarks" / "supplementary" / "s5_synthetic.py"


def _module():
    spec = importlib.util.spec_from_file_location("s5_synthetic", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_s5_reports_the_test_totals_and_every_configuration(tmp_path):
    s5 = _module()
    text = s5.build(legend=tmp_path / "absent.md")  # no legend: tables only
    assert (
        "1,500 runs, 1,480 with constants: 1,058 complete, 160 partial, 213 refused, 49 wrong"
        in text
    )
    assert "Negative controls: 20 of 20 refused." in text
    for config in ("BASE", "C80", "X3", "NC_donor", "NC_clone"):
        assert f"| {config} |" in text
    assert "| ID | Desirable | Expected from the current version |" in text
    assert "overrides equal the oracle in 30 of 30" in text


def test_s5_wrong_runs_are_only_in_the_two_known_conditions():
    import pandas as pd

    test = pd.read_csv(
        ROOT / "benchmarks" / "synthetic" / "runs" / "test" / "results.tsv", sep="\t"
    )
    assert set(test[test.outcome == "wrong"].config) == {"C80", "X3"}
