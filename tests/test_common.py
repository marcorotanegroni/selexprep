"""Unit tests for selexprep._common."""

from __future__ import annotations

from pathlib import Path

from selexprep._common import iter_srr_files, load_csv, parse_round_number

# ----- iter_srr_files -----


def test_iter_srr_files_finds_single_end(tmp_path: Path) -> None:
    (tmp_path / "SRR123.fastq.gz").touch()
    assert iter_srr_files(tmp_path, "SRR123") == [tmp_path / "SRR123.fastq.gz"]


def test_iter_srr_files_finds_paired_end(tmp_path: Path) -> None:
    (tmp_path / "SRR456_1.fastq.gz").touch()
    (tmp_path / "SRR456_2.fastq.gz").touch()
    found = {p.name for p in iter_srr_files(tmp_path, "SRR456")}
    assert found == {"SRR456_1.fastq.gz", "SRR456_2.fastq.gz"}


def test_iter_srr_files_no_prefix_collision(tmp_path: Path) -> None:
    """SRR1234 must NOT match SRR12345.fastq.gz."""
    (tmp_path / "SRR1234.fastq.gz").touch()
    (tmp_path / "SRR12345.fastq.gz").touch()
    found = {p.name for p in iter_srr_files(tmp_path, "SRR1234")}
    assert found == {"SRR1234.fastq.gz"}


def test_iter_srr_files_missing_root_returns_empty(tmp_path: Path) -> None:
    assert iter_srr_files(tmp_path / "does_not_exist", "SRR1") == []


def test_iter_srr_files_recurses_subdirs(tmp_path: Path) -> None:
    sub = tmp_path / "round_01"
    sub.mkdir()
    (sub / "SRR789.fastq.gz").touch()
    assert iter_srr_files(tmp_path, "SRR789") == [sub / "SRR789.fastq.gz"]


# ----- load_csv -----


def test_load_csv_returns_list_of_dicts(tmp_path: Path) -> None:
    csv_path = tmp_path / "in.csv"
    csv_path.write_text("a,b\n1,2\n3,4\n", encoding="utf-8")
    rows = load_csv(csv_path)
    assert rows == [{"a": "1", "b": "2"}, {"a": "3", "b": "4"}]


def test_load_csv_missing_returns_empty(tmp_path: Path) -> None:
    assert load_csv(tmp_path / "missing.csv") == []


# ----- parse_round_number -----


def test_parse_round_number_extracts_int() -> None:
    assert parse_round_number(Path("round_03.counts.parquet")) == 3


def test_parse_round_number_strips_clusters_suffix() -> None:
    assert parse_round_number(Path("round_07.clusters.parquet")) == 7


def test_parse_round_number_falls_back_to_stem() -> None:
    assert parse_round_number(Path("enrich_final.parquet")) == "enrich_final"
