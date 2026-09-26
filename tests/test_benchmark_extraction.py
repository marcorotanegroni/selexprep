"""Tests for the per-run extraction-yield table (``selexprep.benchmark.extraction``)."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

from selexprep.benchmark.extraction import COLUMNS, extraction_rows, main


def _write_fasta_gz(path: Path, lengths: list[int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        for i, n in enumerate(lengths):
            fh.write(f">r{i}\n{'A' * n}\n")


def _deposit(results: Path, acc: str, *, mode: str = "BOTH_PRIMERS_SINGLE_READ") -> Path:
    acc_dir = results / acc
    (acc_dir / "extract").mkdir(parents=True)
    (acc_dir / "library_report.json").write_text(json.dumps({"extraction_mode": mode}))
    return acc_dir


def test_one_row_per_input_fastq_with_yield_and_round_length(tmp_path: Path) -> None:
    """Two runs pooled in one round: the second keeps nothing (a differing tag)."""
    acc_dir = _deposit(tmp_path, "PRJX")
    out = acc_dir / "extract" / "round_01" / "extracted.fasta.gz"
    _write_fasta_gz(out, [40] * 9 + [39])
    reports = [
        {
            "cutadapt_cmd": ["cutadapt", "-g", "A...C", "-o", "x.fa", f"/d/{srr}.fastq.gz"],
            "n_in": 100,
            "n_out": kept,
            "return_code": 0,
            "output_paths": [str(out)],
        }
        for srr, kept in (("SRR2", 0), ("SRR1", 10))
    ]
    (acc_dir / "extract" / "trim_reports.json").write_text(json.dumps(reports))
    (acc_dir / "extract" / "extract.rc").write_text("0\n")

    rows = extraction_rows(tmp_path, "PRJX", 40)

    assert [(r["fastq"], r["frac_kept"]) for r in rows] == [
        ("SRR1.fastq.gz", "0.1000"),
        ("SRR2.fastq.gz", "0.0000"),
    ]
    assert {r["round"] for r in rows} == {"01"}
    assert rows[0]["output_len_mode"] == 40
    assert rows[0]["output_frac_len_mode"] == "0.9000"
    assert rows[0]["extraction_mode"] == "BOTH_PRIMERS_SINGLE_READ"


def test_paired_split_call_reports_the_r1_input(tmp_path: Path) -> None:
    acc_dir = _deposit(tmp_path, "PRJP", mode="PAIRED_END_SPLIT_PRIMERS")
    r1_out = acc_dir / "extract" / "round_02" / "partial_5p_extracted_R1.fasta.gz"
    _write_fasta_gz(r1_out, [110, 110])
    cmd = [
        "cutadapt",
        "-g",
        "A",
        "-G",
        "C",
        "-o",
        "a",
        "-p",
        "b",
        "/d/S_1.fastq.gz",
        "/d/S_2.fastq.gz",
    ]
    entry = {
        "cutadapt_cmd": cmd,
        "n_in": 4,
        "n_out": 2,
        "return_code": 0,
        "output_paths": [str(r1_out), str(r1_out)],
    }
    (acc_dir / "extract" / "trim_reports.json").write_text(json.dumps([entry]))
    (acc_dir / "extract" / "extract.rc").write_text("0\n")

    (row,) = extraction_rows(tmp_path, "PRJP", 40)

    assert row["fastq"] == "S_1.fastq.gz"
    assert row["round"] == "02"
    assert row["frac_kept"] == "0.5000"


def test_refusal_and_missing_run_are_reported_not_raised(tmp_path: Path) -> None:
    refused = _deposit(tmp_path, "PRJR", mode="UNABLE_TO_EXTRACT")
    (refused / "extract" / "extract.rc").write_text("2\n")
    (refused / "extract" / "extract.log").write_text("extract: skipped - no primers\n\n")
    _deposit(tmp_path, "PRJN")

    (r_row,) = extraction_rows(tmp_path, "PRJR", 35)
    (n_row,) = extraction_rows(tmp_path, "PRJN", 35)

    assert r_row["note"] == "extract exit 2: extract: skipped - no primers"
    assert n_row["note"] == "extract not run"


def test_cli_covers_only_verified_recovery_rows(tmp_path: Path) -> None:
    gt = tmp_path / "ground_truth.tsv"
    header = "accession\tn_length_truth\tverified\tread_state"
    gt.write_text(
        f"{header}\n"
        "PRJA\t40\ttrue\traw_standard\n"
        "PRJB\t40\ttrue\tpre_trimmed\n"
        "PRJC\t40\tfalse\traw_standard\n"
    )
    _deposit(tmp_path / "results", "PRJA")
    out = tmp_path / "results" / "extraction_yield.tsv"

    argv = [
        "--ground-truth",
        str(gt),
        "--results-dir",
        str(tmp_path / "results"),
        "--out",
        str(out),
    ]
    assert main(argv) == 0

    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == COLUMNS
    assert [ln.split("\t")[0] for ln in lines[1:]] == ["PRJA"]
