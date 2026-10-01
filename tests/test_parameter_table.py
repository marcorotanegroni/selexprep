"""The supplementary parameter table (benchmarks/parameters/parameter_table.py)."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "benchmarks" / "parameters" / "parameter_table.py"


@pytest.fixture(scope="module")
def table():
    spec = importlib.util.spec_from_file_location("parameter_table", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _public_constants(module) -> set[str]:
    return {n for n in vars(module) if n.isupper() and not n.startswith("_")}


def test_every_public_constant_is_listed_or_excluded_on_purpose(table):
    """A new constant in detect, adapters or extract cannot be left out unnoticed."""
    listed = {(m, n) for _, m, n, _ in table.PARAMETERS}
    for module_name, module in table.MODULES.items():
        for name in _public_constants(module):
            # Constants imported from another module are documented where they live.
            if (module_name, name) in listed or name in table.EXCLUDED:
                continue
            owner = next(
                (
                    m
                    for m, mod in table.MODULES.items()
                    if name in vars(mod) and (m, name) in listed
                ),
                None,
            )
            assert owner is not None, f"{module_name}.{name} is neither listed nor excluded"


def test_values_come_from_the_package(table):
    from selexprep.library import detect

    by_name = {r["parameter"]: r for r in table.rows()}
    assert by_name["BOUNDARY_HIGH_SUPPORT_POST_MAX"]["value"] == str(
        detect.BOUNDARY_HIGH_SUPPORT_POST_MAX
    )
    assert by_name["DEFAULT_TOP_N"]["value"] == "None (all)"
    assert "persistence 0.25" in by_name["COMPOSITE_WEIGHTS"]["value"]
    assert by_name["--error-rate"]["source"].startswith("cutadapt ")


def test_the_table_is_written(table, tmp_path):
    assert table.main(["--out", str(tmp_path)]) == 0
    md = (tmp_path / "parameters.md").read_text()
    assert md.startswith("# Parameters of `detect` and `extract` (selexprep ")
    assert (tmp_path / "parameters.tsv").read_text().count("\n") == len(table.rows()) + 1
