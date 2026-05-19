"""Smoke tests for the Phase 0 Typer CLI scaffold."""

from __future__ import annotations

from typer.testing import CliRunner

from selexprep.cli import app

runner = CliRunner()


def test_version_flag_emits_version_string() -> None:
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert "selexprep" in result.stdout


def test_help_lists_all_subcommands() -> None:
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    # Stub processing verbs + the Phase 1.5 catalog subapp
    for cmd in ("inspect", "fetch", "detect", "extract", "count", "qc", "run", "catalog"):
        assert cmd in result.stdout


def test_inspect_stub_exits_with_code_2() -> None:
    result = runner.invoke(app, ["inspect", "SRR000000"])
    assert result.exit_code == 2


def test_fetch_stub_exits_with_code_2() -> None:
    result = runner.invoke(app, ["fetch", "SRR000000", "--outdir", "/tmp/sx"])
    assert result.exit_code == 2


def test_no_args_shows_help() -> None:
    result = runner.invoke(app, [])
    assert "Usage:" in result.stdout
