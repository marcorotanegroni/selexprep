"""The status note: printed wherever a status is shown next to inferred primers.

`status` rates how well the reads support the inferred primers, not a verified
boundary (benchmarks/synthetic test run: one sequence dominating every round is
taken for constant with status HIGH). The note must reach users of both
`detect` and `run`, and stay silent when nothing was inferred.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest
from typer.testing import CliRunner

from selexprep import cli
from selexprep.library.report import LibraryReport
from selexprep.run import RunReport, RunRowReport

runner = CliRunner()


def _report(status: str, primer_5p: str | None, primer_3p: str | None) -> LibraryReport:
    two_sided = bool(primer_5p and primer_3p)
    return LibraryReport(
        primer_5p=primer_5p,
        primer_3p=primer_3p,
        variants_5p=[],
        variants_3p=[],
        known_adapter_hits={},
        extraction_mode="BOTH_PRIMERS_SINGLE_READ" if two_sided else "UNABLE_TO_EXTRACT",
        full_insert_recovered=two_sided,
        read_source="R1",
        required_action="NONE" if two_sided else "MANUAL_PRIMERS_REQUIRED",
        orientation="FORWARD",
        n_length_mode=20 if two_sided else None,
        n_length_distribution={},
        n_length_confidence=1.0,
        match_rate_5p=1.0,
        match_rate_3p=1.0,
        position_consistency_5p=1.0,
        position_consistency_3p=1.0,
        read_fraction_used_for_inference=1.0,
        sampling_seed=42,
        confidence=0.95 if two_sided else 0.0,
        status=status,  # type: ignore[arg-type]
        failure_reason=None if two_sided else "no primer",
    )


def _inputs(tmp_path: Path) -> tuple[Path, Path]:
    fastq = tmp_path / "r1.fastq.gz"
    with gzip.open(fastq, "wt") as fh:
        fh.write("@r\nACGT\n+\nIIII\n")
    rounds = tmp_path / "rounds.tsv"
    rounds.write_text("file\tround_number\nr1.fastq.gz\t1\n")
    return fastq, rounds


@pytest.mark.parametrize(
    ("status", "p5", "p3", "shown"),
    [
        ("HIGH", "ACGTACGTACGTACGT", "TTGGCCAATTGGCCAA", True),
        ("MEDIUM", "ACGTACGTACGTACGT", "TTGGCCAATTGGCCAA", True),
        ("UNABLE_TO_INFER", None, None, False),
    ],
)
def test_detect_prints_the_note_when_primers_were_inferred(
    tmp_path, monkeypatch, status, p5, p3, shown
):
    fastq, rounds = _inputs(tmp_path)
    monkeypatch.setattr(cli, "compute_library_report", lambda *a, **k: _report(status, p5, p3))
    result = runner.invoke(
        cli.app,
        ["detect", str(fastq), "--round-map", str(rounds), "--outdir", str(tmp_path / "o")],
    )
    assert result.exit_code == 0, result.output
    assert (cli.STATUS_NOTE in result.output) is shown


@pytest.mark.parametrize(
    ("statuses", "shown"),
    [
        (["HIGH", "UNABLE_TO_INFER"], True),
        (["UNABLE_TO_INFER", None], False),
    ],
)
def test_run_prints_the_note_when_any_dataset_has_inferred_primers(
    tmp_path, monkeypatch, statuses, shown
):
    rows = [
        RunRowReport(
            accession=f"PRJ{i}",
            status="OK",
            last_stage_completed="qc",
            library_report_status=s,
        )
        for i, s in enumerate(statuses)
    ]
    report = RunReport(accessions_tsv=tmp_path / "a.tsv", outdir=tmp_path, rows=rows)
    monkeypatch.setattr("selexprep.run.run_batch", lambda *a, **k: report)
    accessions = tmp_path / "a.tsv"
    accessions.write_text("accession\nPRJ0\n")
    result = runner.invoke(cli.app, ["run", str(accessions), "--outdir", str(tmp_path / "o")])
    assert result.exit_code == 0, result.output
    assert (cli.STATUS_NOTE in result.output) is shown
