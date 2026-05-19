"""Smoke tests for the Phase 0 Typer CLI scaffold."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

from typer.testing import CliRunner

from selexprep.cli import app

runner = CliRunner()


# Local copy of _synthetic_pool — keeps test_cli independent of test_detect.
def _synthetic_pool(
    primer_5p: str | None,
    primer_3p: str | None,
    n: int = 1000,
    random_len: int = 30,
) -> list[str]:
    bases = "ACGT"
    out: list[str] = []
    for i in range(n):
        rand = "".join(bases[(i * 7 + j * 13) % 4] for j in range(random_len))
        out.append((primer_5p or "") + rand + (primer_3p or ""))
    return out


def _write_fastq_gz(path: Path, seqs: list[str]) -> None:
    """Write a minimal FASTQ.gz: identical quality line per read."""
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for i, s in enumerate(seqs):
            fh.write(f"@read_{i}\n{s}\n+\n{'I' * len(s)}\n")


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


# ===========================================================================
# Phase 2 — `detect` CLI command (Wired)
# ===========================================================================


def test_detect_without_round_map_exits_with_code_2(tmp_path: Path) -> None:
    fq = tmp_path / "round0.fastq.gz"
    _write_fastq_gz(fq, _synthetic_pool("GGTAATACGACTCACTATAGGG", "CCATGCATGCATGCATGCAT"))
    result = runner.invoke(app, ["detect", str(fq), "--outdir", str(tmp_path / "out")])
    assert result.exit_code == 2
    assert "round-map" in result.stderr or "round-map" in result.output


def test_detect_with_round_map_emits_library_report(tmp_path: Path) -> None:
    primer_5p = "GGTAATACGACTCACTATAGGG"
    primer_3p = "CCATGCATGCATGCATGCAT"
    seqs = _synthetic_pool(primer_5p, primer_3p, n=1000)

    fq = tmp_path / "round0.fastq.gz"
    _write_fastq_gz(fq, seqs)

    round_map = tmp_path / "rounds.tsv"
    round_map.write_text("file\tround_number\nround0.fastq.gz\t0\n", encoding="utf-8")

    outdir = tmp_path / "out"
    result = runner.invoke(
        app, ["detect", str(fq), "--round-map", str(round_map), "--outdir", str(outdir)]
    )
    assert result.exit_code == 0, result.output

    out_path = outdir / "library_report.json"
    assert out_path.exists()
    payload = json.loads(out_path.read_text())
    assert payload["primer_5p"] is not None
    assert payload["primer_3p"] is not None
    assert payload["extraction_mode"] == "BOTH_PRIMERS_SINGLE_READ"


def test_detect_fastq_not_in_round_map_exits_with_code_2(tmp_path: Path) -> None:
    fq = tmp_path / "round0.fastq.gz"
    _write_fastq_gz(fq, _synthetic_pool("GGTAATACGACTCACTATAGGG", "CCATGCATGCATGCATGCAT"))

    round_map = tmp_path / "rounds.tsv"
    # Map references a different filename than the one we pass.
    round_map.write_text("file\tround_number\nelsewhere.fastq.gz\t0\n", encoding="utf-8")

    outdir = tmp_path / "out"
    result = runner.invoke(
        app, ["detect", str(fq), "--round-map", str(round_map), "--outdir", str(outdir)]
    )
    assert result.exit_code == 2


# ===========================================================================
# Phase 3 - `extract` CLI command
# ===========================================================================


def _write_library_report_json(path: Path, **overrides: object) -> None:
    """Hand-write a minimal LibraryReport JSON for CLI tests."""
    import json as _json

    base = {
        "primer_5p": "GGTAATACGACTCACTATAGGG",
        "primer_3p": "CCATGCATGCATGCATGCAT",
        "variants_5p": [],
        "variants_3p": [],
        "known_adapter_hits": {},
        "extraction_mode": "BOTH_PRIMERS_SINGLE_READ",
        "full_insert_recovered": True,
        "read_source": "R1",
        "required_action": "NONE",
        "orientation": "FORWARD",
        "n_length_mode": 30,
        "n_length_distribution": {"30": 100},
        "n_length_confidence": 1.0,
        "match_rate_5p": 0.95,
        "match_rate_3p": 0.95,
        "position_consistency_5p": 0.95,
        "position_consistency_3p": 0.95,
        "read_fraction_used_for_inference": 1.0,
        "sampling_seed": 42,
        "confidence": 0.85,
        "status": "HIGH",
        "failure_reason": None,
    }
    base.update(overrides)
    path.write_text(_json.dumps(base, indent=2), encoding="utf-8")


def test_extract_without_round_map_or_sample_sheet_errors(tmp_path: Path) -> None:
    fq = tmp_path / "r0.fastq.gz"
    _write_fastq_gz(fq, _synthetic_pool("GGTAATACGACTCACTATAGGG", "CCATGCATGCATGCATGCAT"))
    lr_path = tmp_path / "lr.json"
    _write_library_report_json(lr_path)

    result = runner.invoke(
        app,
        [
            "extract",
            str(fq),
            "--library-report",
            str(lr_path),
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code == 2


def test_extract_override_primer_flags_emit_phase4_error(tmp_path: Path) -> None:
    fq = tmp_path / "r0.fastq.gz"
    _write_fastq_gz(fq, _synthetic_pool("GGTAATACGACTCACTATAGGG", "CCATGCATGCATGCATGCAT"))
    lr_path = tmp_path / "lr.json"
    _write_library_report_json(lr_path)

    result = runner.invoke(
        app,
        [
            "extract",
            str(fq),
            "--library-report",
            str(lr_path),
            "--override-primer-5p",
            "AAAA",
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code == 2
    assert "Phase 4" in (result.stderr or result.output)


def test_extract_refuses_unable_library_report(tmp_path: Path) -> None:
    fq = tmp_path / "r0.fastq.gz"
    _write_fastq_gz(fq, _synthetic_pool("GGTAATACGACTCACTATAGGG", "CCATGCATGCATGCATGCAT"))
    lr_path = tmp_path / "lr.json"
    _write_library_report_json(
        lr_path,
        status="UNABLE_TO_INFER",
        extraction_mode="UNABLE_TO_EXTRACT",
        failure_reason="synthetic test failure",
    )
    rm = tmp_path / "rounds.tsv"
    rm.write_text(f"file\tround_number\n{fq.name}\t0\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "extract",
            str(fq),
            "--library-report",
            str(lr_path),
            "--round-map",
            str(rm),
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code == 2


def test_extract_end_to_end_both_primers(tmp_path: Path) -> None:
    """Smoke test: real cutadapt run via the CLI -> deterministic FASTA on disk."""
    import shutil as _shutil

    if _shutil.which("cutadapt") is None:
        return  # silently skip — covered by @pytest.mark.skipif in trim/runner tests

    fq = tmp_path / "r0.fastq.gz"
    _write_fastq_gz(fq, _synthetic_pool("GGTAATACGACTCACTATAGGG", "CCATGCATGCATGCATGCAT"))
    lr_path = tmp_path / "lr.json"
    _write_library_report_json(lr_path)
    rm = tmp_path / "rounds.tsv"
    rm.write_text(f"file\tround_number\n{fq.name}\t0\n", encoding="utf-8")
    outdir = tmp_path / "out"

    result = runner.invoke(
        app,
        [
            "extract",
            str(fq),
            "--library-report",
            str(lr_path),
            "--round-map",
            str(rm),
            "--outdir",
            str(outdir),
        ],
    )
    assert result.exit_code == 0, result.output
    assert (outdir / "round_00" / "extracted.fasta.gz").exists()
    assert (outdir / "trim_reports.json").exists()
