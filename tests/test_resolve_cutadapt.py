"""Tests for ``selexprep._common.resolve_cutadapt``.

The locator lets selexprep find its declared cutadapt dependency even when the
environment isn't "activated" — ``pipx install`` (which exposes only
selexprep's own entry point on PATH), absolute-path invocation, or a workflow
runner with a sanitized PATH. It prefers ``$PATH`` and falls back to the
directory of the running interpreter.
"""

from __future__ import annotations

from pathlib import Path

from selexprep import _common


def test_resolve_cutadapt_prefers_path(monkeypatch) -> None:
    monkeypatch.setattr(_common.shutil, "which", lambda name: "/usr/bin/cutadapt")
    assert _common.resolve_cutadapt() == "/usr/bin/cutadapt"


def test_resolve_cutadapt_falls_back_to_interpreter_dir(monkeypatch, tmp_path: Path) -> None:
    # Nothing on PATH, but a cutadapt sits next to the interpreter (the pipx /
    # unactivated-venv case the fallback exists for).
    monkeypatch.setattr(_common.shutil, "which", lambda name: None)
    fake = tmp_path / "cutadapt"
    fake.write_text("#!/bin/sh\n", encoding="utf-8")
    monkeypatch.setattr(_common.sys, "executable", str(tmp_path / "python"))
    assert _common.resolve_cutadapt() == str(fake)


def test_resolve_cutadapt_absent_returns_none(monkeypatch, tmp_path: Path) -> None:
    # Not on PATH and not next to the interpreter -> callers decide (extract
    # raises; count soft-skips).
    monkeypatch.setattr(_common.shutil, "which", lambda name: None)
    monkeypatch.setattr(_common.sys, "executable", str(tmp_path / "python"))
    assert _common.resolve_cutadapt() is None
