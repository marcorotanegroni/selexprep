"""Guards for a user running selexprep on their own FASTQ files.

Each test is a place where an outside user was stopped, or was given a
plausible but wrong result without a word: uncompressed FASTQ with a primer
override, a truncated gzip, an inferred "primer" that is a sequencing adapter,
a one-sided extraction mode with both primers supplied, and a reverse primer
typed as ordered.
"""

from __future__ import annotations

import gzip
import logging
import random
from pathlib import Path

import pytest
from typer.testing import CliRunner

from selexprep._common import resolve_cutadapt
from selexprep._io import InputReadError, open_fastx_text
from selexprep.cli import _warn_about_inferred_library, app
from selexprep.extract.runner import run_extract
from selexprep.library.adapters import KNOWN_ADAPTERS, known_adapter_within, reverse_complement
from selexprep.library.report import LibraryReport

runner = CliRunner()
needs_cutadapt = pytest.mark.skipif(resolve_cutadapt() is None, reason="cutadapt not available")

PRIMER_5P = "ATCCAGAGTGACGCAGCATTGC"
PRIMER_3P = "TGGACACGGTGGCTTAGTCAA"
TRUSEQ_TAIL = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"


def _inserts(n: int, length: int = 30, seed: int = 0) -> list[str]:
    rng = random.Random(seed)
    return ["".join(rng.choice("ACGT") for _ in range(length)) for _ in range(n)]


def _fastq_text(seqs: list[str]) -> str:
    return "".join(f"@read_{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in enumerate(seqs))


def _write_fastq(path: Path, seqs: list[str]) -> Path:
    if path.suffix == ".gz":
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            fh.write(_fastq_text(seqs))
    else:
        path.write_text(_fastq_text(seqs), encoding="utf-8")
    return path


def _report(
    *,
    primer_5p: str | None = PRIMER_5P,
    primer_3p: str | None = PRIMER_3P,
    extraction_mode: str = "BOTH_PRIMERS_SINGLE_READ",
) -> LibraryReport:
    return LibraryReport(
        primer_5p=primer_5p,
        primer_3p=primer_3p,
        variants_5p=[],
        variants_3p=[],
        known_adapter_hits={},
        extraction_mode=extraction_mode,  # type: ignore[arg-type]
        full_insert_recovered=extraction_mode == "BOTH_PRIMERS_SINGLE_READ",
        read_source="R1",
        required_action="NONE",
        orientation="FORWARD",
        n_length_mode=30,
        n_length_distribution={30: 20},
        n_length_confidence=1.0,
        match_rate_5p=0.95,
        match_rate_3p=0.95,
        position_consistency_5p=0.95,
        position_consistency_3p=0.95,
        read_fraction_used_for_inference=1.0,
        sampling_seed=42,
        confidence=0.85,
        status="HIGH",
        failure_reason=None,
    )


# ---------------------------------------------------------------------------
# Reading inputs: compressed or not, complete or not
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ["reads.fastq", "reads.fastq.gz", "plain_named.gz"])
def test_open_fastx_text_reads_by_content_not_by_name(tmp_path: Path, name: str) -> None:
    path = tmp_path / name
    if name == "plain_named.gz":
        path.write_text(_fastq_text(["ACGT", "TTTT"]), encoding="utf-8")
    else:
        _write_fastq(path, ["ACGT", "TTTT"])
    with open_fastx_text(path) as fh:
        lines = fh.read().splitlines()
    assert lines[1] == "ACGT" and lines[5] == "TTTT"


def _truncated_gz(tmp_path: Path) -> Path:
    whole = _write_fastq(tmp_path / "whole.fastq.gz", _inserts(2000, 60))
    cut = tmp_path / "cut.fastq.gz"
    data = whole.read_bytes()
    cut.write_bytes(data[: len(data) // 2])
    return cut


def test_truncated_gzip_raises_an_error_that_names_the_file(tmp_path: Path) -> None:
    cut = _truncated_gz(tmp_path)
    with pytest.raises(InputReadError) as err, open_fastx_text(cut) as fh:
        fh.read()
    assert str(cut) in str(err.value)
    assert "gzip -t" in str(err.value)


def test_detect_on_truncated_gzip_says_which_file(tmp_path: Path) -> None:
    """Was: a blank line and ``Aborted.`` with no file and no cause."""
    cut = _truncated_gz(tmp_path)
    round_map = tmp_path / "rounds.tsv"
    round_map.write_text(f"file\tround_number\n{cut.name}\t0\n", encoding="utf-8")
    result = runner.invoke(
        app, ["detect", str(cut), "--round-map", str(round_map), "--outdir", str(tmp_path / "out")]
    )
    assert result.exit_code == 2
    assert cut.name in result.output
    assert "truncated or corrupt" in result.output
    assert "Aborted" not in result.output


@needs_cutadapt
@pytest.mark.parametrize("suffix", [".fastq", ".fastq.gz"])
def test_override_works_on_uncompressed_and_compressed_fastq(tmp_path: Path, suffix: str) -> None:
    """Was: ``BadGzipFile`` on a plain FASTQ as soon as an override was given."""
    seqs = [PRIMER_5P + ins + PRIMER_3P for ins in _inserts(50)]
    fq = _write_fastq(tmp_path / f"round_00{suffix}", seqs)
    result = run_extract(
        _report(),
        [fq],
        tmp_path / "out",
        round_map={fq.name: 0},
        override_primer_5p=PRIMER_5P,
    )
    assert not result.skipped, result.skipped_reason
    assert [r.n_out for r in result.trim_reports] == [50]


@needs_cutadapt
def test_extract_on_truncated_gzip_exits_without_a_traceback(tmp_path: Path) -> None:
    cut = _truncated_gz(tmp_path)
    round_map = tmp_path / "rounds.tsv"
    round_map.write_text(f"file\tround_number\n{cut.name}\t0\n", encoding="utf-8")
    lr = tmp_path / "library_report.json"
    lr.write_text(_report().model_dump_json(), encoding="utf-8")
    result = runner.invoke(
        app,
        [
            "extract",
            str(cut),
            "--library-report",
            str(lr),
            "--round-map",
            str(round_map),
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code == 2
    assert cut.name in result.output
    assert "truncated or corrupt" in result.output
    assert "Traceback" not in result.output


# ---------------------------------------------------------------------------
# A called primer that is a sequencing adapter
# ---------------------------------------------------------------------------


def test_known_adapter_within_finds_the_adapter_anywhere() -> None:
    truseq = KNOWN_ADAPTERS["TRUSEQ_R1"]
    assert known_adapter_within(truseq + "ACAC") == "TRUSEQ_R1"
    assert known_adapter_within("AA" + truseq + "ACAC") == "TRUSEQ_R1"
    assert known_adapter_within(PRIMER_3P + truseq) == "TRUSEQ_R1"
    assert known_adapter_within(reverse_complement(truseq) + "GG") == "TRUSEQ_R1"
    assert known_adapter_within(PRIMER_3P) is None
    assert known_adapter_within(None) is None


def test_detect_warns_when_the_inferred_primer_contains_an_adapter(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The call the read-through datasets produced: two bases of constant, then TruSeq."""
    _warn_about_inferred_library(
        _report(primer_5p=None, primer_3p="AAAGATCGGAAGAGCACAC", extraction_mode="THREE_PRIME_ONLY")
    )
    err = capsys.readouterr().err
    assert "contains the TRUSEQ_R1 sequencing adapter" in err
    assert "--override-primer-3p" in err
    assert "only the 3' constant was found" in err


def test_detect_stays_quiet_for_a_clean_two_primer_library(
    capsys: pytest.CaptureFixture[str],
) -> None:
    _warn_about_inferred_library(_report())
    assert capsys.readouterr().err == ""


def test_detect_cli_warns_on_read_through_into_the_adapter(tmp_path: Path) -> None:
    """End to end: reads longer than the construct, running into TruSeq."""
    files = []
    for round_number in (0, 3):
        seqs = [
            PRIMER_5P + ins + PRIMER_3P + TRUSEQ_TAIL
            for ins in _inserts(1500, 30, seed=round_number)
        ]
        files.append(_write_fastq(tmp_path / f"pool_{round_number}.fastq.gz", seqs))
    round_map = tmp_path / "rounds.tsv"
    round_map.write_text(
        "file\tround_number\n"
        + "".join(f"{f.name}\t{r}\n" for f, r in zip(files, (0, 3), strict=True)),
        encoding="utf-8",
    )
    result = runner.invoke(
        app,
        [
            "detect",
            *map(str, files),
            "--round-map",
            str(round_map),
            "--outdir",
            str(tmp_path / "o"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "sequencing adapter" in result.output


# ---------------------------------------------------------------------------
# Supplied primers
# ---------------------------------------------------------------------------


def _warnings(caplog: pytest.LogCaptureFixture, needle: str) -> list[str]:
    return [
        r.getMessage()
        for r in caplog.records
        if r.levelname == "WARNING" and needle in r.getMessage()
    ]


@needs_cutadapt
def test_both_primers_with_one_sided_mode_warns_that_one_is_not_trimmed(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Was: the second primer was ignored for trimming, with no message."""
    seqs = [PRIMER_5P + ins + PRIMER_3P + TRUSEQ_TAIL for ins in _inserts(50)]
    fq = _write_fastq(tmp_path / "round_00.fastq.gz", seqs)
    with caplog.at_level(logging.WARNING):
        result = run_extract(
            _report(primer_3p=None, extraction_mode="FIVE_PRIME_ONLY"),
            [fq],
            tmp_path / "out",
            round_map={fq.name: 0},
            override_primer_5p=PRIMER_5P,
            override_primer_3p=PRIMER_3P,
        )
    assert not result.skipped, result.skipped_reason
    hits = _warnings(caplog, "both primers are set")
    assert len(hits) == 1
    assert "FIVE_PRIME_ONLY" in hits[0] and "BOTH_PRIMERS_SINGLE_READ" in hits[0]


@needs_cutadapt
def test_both_primers_with_two_sided_mode_does_not_warn(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    seqs = [PRIMER_5P + ins + PRIMER_3P for ins in _inserts(50)]
    fq = _write_fastq(tmp_path / "round_00.fastq.gz", seqs)
    with caplog.at_level(logging.WARNING):
        run_extract(
            _report(),
            [fq],
            tmp_path / "out",
            round_map={fq.name: 0},
            override_primer_5p=PRIMER_5P,
            override_primer_3p=PRIMER_3P,
        )
    assert _warnings(caplog, "both primers are set") == []
    assert _warnings(caplog, "of the earliest round's reads") == []


@needs_cutadapt
def test_reverse_primer_typed_as_ordered_is_named_in_the_warning(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Was: 0 reads kept, exit 0, and only the generic low-yield warning."""
    seqs = [PRIMER_5P + ins + PRIMER_3P for ins in _inserts(50)]
    fq = _write_fastq(tmp_path / "round_00.fastq.gz", seqs)
    as_ordered = reverse_complement(PRIMER_3P)
    with caplog.at_level(logging.WARNING):
        run_extract(
            _report(),
            [fq],
            tmp_path / "out",
            round_map={fq.name: 0},
            override_primer_3p=as_ordered,
        )
    hits = _warnings(caplog, "of the earliest round's reads")
    assert len(hits) == 1
    assert "supplied 3' primer" in hits[0]
    assert f"reverse complement, {PRIMER_3P}, is in 100.0%" in hits[0]


@needs_cutadapt
@pytest.mark.parametrize(
    "degenerate_5p",
    ["ATCCAGRGTGACGCAGCATTGC", "ATCCAGRGTGACGCAGCATYGC", "ATCCNGRGTGACGCAGCATYGC"],
    ids=["one_code", "two_codes", "three_codes"],
)
def test_degenerate_primer_that_fits_the_reads_does_not_warn(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, degenerate_5p: str
) -> None:
    """Was (0.4.8): two or more IUPAC codes gave "is in 0.0% of the reads"."""
    seqs = [PRIMER_5P + ins + PRIMER_3P for ins in _inserts(50)]
    fq = _write_fastq(tmp_path / "round_00.fastq.gz", seqs)
    with caplog.at_level(logging.WARNING):
        result = run_extract(
            _report(),
            [fq],
            tmp_path / "out",
            round_map={fq.name: 0},
            override_primer_5p=degenerate_5p,
        )
    assert [r.n_out for r in result.trim_reports] == [50]
    assert _warnings(caplog, "of the earliest round's reads") == []


# ---------------------------------------------------------------------------
# --extraction-mode
# ---------------------------------------------------------------------------


def _extracted_lengths(path: Path) -> set[int]:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return {len(line.strip()) for line in fh if not line.startswith(">")}


@needs_cutadapt
def test_extraction_mode_trims_both_primers_on_a_one_sided_report(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Read-through reads, a FIVE_PRIME_ONLY report, both primers known."""
    from selexprep.manifest import read_manifest_json

    seqs = [PRIMER_5P + ins + PRIMER_3P + TRUSEQ_TAIL for ins in _inserts(50)]
    fq = _write_fastq(tmp_path / "round_00.fastq.gz", seqs)
    with caplog.at_level(logging.WARNING):
        result = run_extract(
            _report(primer_3p=None, extraction_mode="FIVE_PRIME_ONLY"),
            [fq],
            tmp_path / "out",
            round_map={fq.name: 0},
            override_primer_5p=PRIMER_5P,
            override_primer_3p=PRIMER_3P,
            extraction_mode="BOTH_PRIMERS_SINGLE_READ",
        )
    assert not result.skipped, result.skipped_reason
    target = tmp_path / "out" / "overridden" / "round_00" / "extracted.fasta.gz"
    assert _extracted_lengths(target) == {30}
    report = read_manifest_json(result.manifest_path).library_report  # type: ignore[arg-type]
    assert report.extraction_mode == "BOTH_PRIMERS_SINGLE_READ"
    assert report.full_insert_recovered is True
    assert report.required_action == "NONE"
    assert report.n_length_mode == 30
    assert _warnings(caplog, "both primers are set") == []


@needs_cutadapt
def test_extraction_mode_alone_is_an_override(tmp_path: Path) -> None:
    """Without --rebuild it writes under ``overridden/``, like a primer override."""
    seqs = [PRIMER_5P + ins + PRIMER_3P for ins in _inserts(50)]
    fq = _write_fastq(tmp_path / "round_00.fastq.gz", seqs)
    result = run_extract(
        _report(),
        [fq],
        tmp_path / "out",
        round_map={fq.name: 0},
        extraction_mode="FIVE_PRIME_ONLY",
    )
    assert not result.skipped, result.skipped_reason
    target = tmp_path / "out" / "overridden" / "round_00" / "partial_5p_extracted.fasta.gz"
    assert _extracted_lengths(target) == {30 + len(PRIMER_3P)}


def test_extraction_mode_that_needs_a_missing_primer_is_refused(tmp_path: Path) -> None:
    fq = _write_fastq(tmp_path / "round_00.fastq.gz", [PRIMER_5P + "ACGT" * 8])
    result = run_extract(
        _report(primer_3p=None, extraction_mode="FIVE_PRIME_ONLY"),
        [fq],
        tmp_path / "out",
        round_map={fq.name: 0},
        extraction_mode="BOTH_PRIMERS_SINGLE_READ",
    )
    assert result.skipped
    assert "requires both primers" in (result.skipped_reason or "")


def test_extraction_mode_rejects_an_unknown_or_unextractable_value(tmp_path: Path) -> None:
    fq = _write_fastq(tmp_path / "round_00.fastq.gz", [PRIMER_5P + "ACGT" * 8 + PRIMER_3P])
    for bad in ("both", "UNABLE_TO_EXTRACT"):
        with pytest.raises(ValueError, match="is not one of"):
            run_extract(
                _report(), [fq], tmp_path / "out", round_map={fq.name: 0}, extraction_mode=bad
            )


def test_extract_cli_lists_the_modes_when_the_value_is_wrong(tmp_path: Path) -> None:
    fq = _write_fastq(tmp_path / "round_00.fastq.gz", [PRIMER_5P + "ACGT" * 8 + PRIMER_3P])
    round_map = tmp_path / "rounds.tsv"
    round_map.write_text(f"file\tround_number\n{fq.name}\t0\n", encoding="utf-8")
    lr = tmp_path / "library_report.json"
    lr.write_text(_report().model_dump_json(), encoding="utf-8")
    result = runner.invoke(
        app,
        [
            "extract",
            str(fq),
            "--library-report",
            str(lr),
            "--round-map",
            str(round_map),
            "--extraction-mode",
            "both",
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code == 2
    for mode in ("BOTH_PRIMERS_SINGLE_READ", "FIVE_PRIME_ONLY", "PAIRED_END_SPLIT_PRIMERS"):
        assert mode in result.output


@needs_cutadapt
def test_extract_cli_extraction_mode_end_to_end(tmp_path: Path) -> None:
    seqs = [PRIMER_5P + ins + PRIMER_3P + TRUSEQ_TAIL for ins in _inserts(50)]
    fq = _write_fastq(tmp_path / "round_00.fastq", seqs)
    round_map = tmp_path / "rounds.tsv"
    round_map.write_text(f"file\tround_number\n{fq.name}\t0\n", encoding="utf-8")
    lr = tmp_path / "library_report.json"
    lr.write_text(
        _report(primer_3p=None, extraction_mode="FIVE_PRIME_ONLY").model_dump_json(),
        encoding="utf-8",
    )
    result = runner.invoke(
        app,
        [
            "extract",
            str(fq),
            "--library-report",
            str(lr),
            "--round-map",
            str(round_map),
            "--override-primer-3p",
            PRIMER_3P,
            "--extraction-mode",
            "BOTH_PRIMERS_SINGLE_READ",
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code == 0, result.output
    target = tmp_path / "out" / "overridden" / "round_00" / "extracted.fasta.gz"
    assert _extracted_lengths(target) == {30}
