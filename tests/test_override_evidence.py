"""Five fixes ported from selexprepR (commits 820e5d9 and c640797).

1-2. Primer overrides carry the evidence of the supplied primers: match rates
     and random-region lengths measured on the earliest round, not the values
     detection reported for the primers it inferred.
3.   One IUPAC complement for reads, counting and primers.
4.   A side without a candidate flank is named as such, not reported as 0.00.
5.   The low-yield warning names reads in the opposite orientation.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import pytest

from selexprep._common import resolve_cutadapt, reverse_complement_iupac
from selexprep.count.counter import reverse_complement as count_revcomp
from selexprep.extract import runner
from selexprep.extract.strand import _revcomp_read
from selexprep.library.adapters import reverse_complement as primer_revcomp
from selexprep.library.report import LibraryReport, _classify, _low_match_reason

PRIMER_5P = "GGTAATACGACTCACTATAGGG"
PRIMER_3P = "CCATGCATGCATGCATGCAT"
N39 = "ACGTTGCAAGGCTTACGGATCCATGCAAGTCAGTCAGTA"
needs_cutadapt = pytest.mark.skipif(resolve_cutadapt() is None, reason="cutadapt not available")


def _fastq(path: Path, seqs: list[str]) -> Path:
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for i, s in enumerate(seqs):
            fh.write(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n")
    return path


def _report(**fields) -> LibraryReport:
    base = dict(
        primer_5p=PRIMER_5P,
        primer_3p=PRIMER_3P,
        variants_5p=[],
        variants_3p=[],
        known_adapter_hits={},
        extraction_mode="BOTH_PRIMERS_SINGLE_READ",
        full_insert_recovered=True,
        read_source="R1",
        required_action="NONE",
        orientation="FORWARD",
        n_length_mode=150,
        n_length_distribution={150: 100},
        n_length_confidence=1.0,
        match_rate_5p=0.0,
        match_rate_3p=0.0,
        position_consistency_5p=0.9,
        position_consistency_3p=0.9,
        read_fraction_used_for_inference=1.0,
        sampling_seed=42,
        confidence=0.7,
        status="MEDIUM",
        failure_reason=None,
    )
    base.update(fields)
    return LibraryReport(**base)


# ---------------------------------------------------------------------------
# 1. Match rates of the supplied primers
# ---------------------------------------------------------------------------


def test_override_match_rates_use_the_detection_definition(tmp_path):
    """Same case as the R test: reads shorter than a primer are left out of its rate."""
    one_substitution = "GGTAATACGACTCAGTATAGGG" + "ACGT" + PRIMER_3P
    reads = (
        [PRIMER_5P + "ACGT" + PRIMER_3P] * 60
        + [one_substitution] * 20
        + ["ACGT" + PRIMER_3P] * 20
        + [PRIMER_5P[1:]] * 20
    )
    rates = runner._override_match_rates(
        _report(),
        override_primer_5p=PRIMER_5P,
        override_primer_3p=PRIMER_3P,
        earliest_r1=[_fastq(tmp_path / "r0.fastq.gz", reads)],
        earliest_r2=None,
    )
    assert rates["match_rate_5p"] == pytest.approx(80 / 100)
    assert rates["match_rate_3p"] == pytest.approx(100 / 120)


def test_only_overridden_sides_are_measured(tmp_path):
    rates = runner._override_match_rates(
        _report(),
        override_primer_5p=PRIMER_5P,
        override_primer_3p=None,
        earliest_r1=[_fastq(tmp_path / "r0.fastq.gz", [PRIMER_5P + N39] * 10)],
        earliest_r2=None,
    )
    assert rates == {"match_rate_5p": 1.0}


def test_paired_split_measures_the_3p_primer_on_r2(tmp_path):
    r2 = [primer_revcomp(PRIMER_3P) + N39] * 30 + [N39 * 2] * 10
    rates = runner._override_match_rates(
        _report(extraction_mode="PAIRED_END_SPLIT_PRIMERS", full_insert_recovered=False),
        override_primer_5p=None,
        override_primer_3p=PRIMER_3P,
        earliest_r1=[_fastq(tmp_path / "r1.fastq.gz", [PRIMER_5P + N39] * 40)],
        earliest_r2=[_fastq(tmp_path / "r2.fastq.gz", r2)],
    )
    assert rates == {"match_rate_3p": pytest.approx(0.75)}


# ---------------------------------------------------------------------------
# 2. Random-region lengths after trimming the supplied primers
# ---------------------------------------------------------------------------


def test_insert_lengths_come_from_the_extracted_reads(tmp_path):
    fasta = tmp_path / "extracted.fasta.gz"
    with gzip.open(fasta, "wt") as fh:
        for i, s in enumerate([N39] * 9 + [N39[:-1]]):
            fh.write(f">r{i}\n{s}\n")
    assert runner._override_insert_lengths(fasta) == {
        "n_length_mode": 39,
        "n_length_distribution": {39: 9, 38: 1},
        "n_length_confidence": 0.9,
    }


@needs_cutadapt
def test_override_report_carries_the_supplied_primers_evidence(tmp_path):
    """A report inferred with a wrong boundary (150 nt, match rates 0) is
    overridden: the manifest reports the supplied primers' rates and the real
    39-nt random region, not the detection values."""
    from selexprep.manifest import read_manifest_json

    fq = _fastq(tmp_path / "round_00.fastq.gz", [PRIMER_5P + N39 + PRIMER_3P] * 50)
    result = runner.run_extract(
        _report(primer_5p="ACGTACGTACGTACGT", primer_3p="TTTTGGGGCCCCAAAA"),
        [fq],
        tmp_path / "out",
        round_map={fq.name: 0},
        override_primer_5p=PRIMER_5P,
        override_primer_3p=PRIMER_3P,
    )
    assert not result.skipped, result.skipped_reason
    lr = read_manifest_json(result.manifest_path).library_report
    assert (lr.match_rate_5p, lr.match_rate_3p) == (1.0, 1.0)
    assert lr.n_length_mode == 39
    assert lr.n_length_distribution == {39: 50}
    assert lr.n_length_confidence == 1.0


@needs_cutadapt
def test_override_rates_are_measured_in_the_extraction_orientation(tmp_path):
    from selexprep.manifest import read_manifest_json

    reversed_reads = [reverse_complement_iupac(PRIMER_5P + N39 + PRIMER_3P)] * 40
    fq = _fastq(tmp_path / "round_00.fastq.gz", reversed_reads)
    result = runner.run_extract(
        _report(orientation="REVERSE"),
        [fq],
        tmp_path / "out",
        round_map={fq.name: 0},
        override_primer_5p=PRIMER_5P,
    )
    assert not result.skipped, result.skipped_reason
    lr = read_manifest_json(result.manifest_path).library_report
    assert lr.match_rate_5p == 1.0
    assert lr.match_rate_3p == 0.0  # not overridden: detection's value stays
    assert lr.n_length_mode == 39


# ---------------------------------------------------------------------------
# 3. One IUPAC complement
# ---------------------------------------------------------------------------


def test_iupac_codes_are_complemented_keeping_case_and_other_symbols():
    assert count_revcomp("ACGTRYKMSWBVDHN") == "NDHBVWSKMRYACGT"
    assert count_revcomp("acgtrykm") == "kmryacgt"
    assert count_revcomp("AUu") == "aAT"
    assert count_revcomp("AC.-x") == "x-.GT"


def test_read_and_count_complements_agree():
    for seq in ("ACGTRYKMSWBVDHN", "acgtnry", "AC.GT-U"):
        assert _revcomp_read(seq) == count_revcomp(seq) == reverse_complement_iupac(seq)


def test_read_complement_keeps_unknown_symbols_instead_of_n():
    assert _revcomp_read("A.C") == "G.T"


def test_qc_kmer_complement_uses_the_same_table():
    from selexprep.qc.consistency import canonical_kmer

    assert canonical_kmer("TYAA") == canonical_kmer("TTRA") == "TTRA"
    assert canonical_kmer("ACGTN") == "ACGTN"


def test_primer_complement_agrees_on_uppercase_iupac():
    assert primer_revcomp("ACGTRYKMSWBVDHN") == count_revcomp("ACGTRYKMSWBVDHN")


# ---------------------------------------------------------------------------
# 4. Sides without a candidate flank
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("rates", "candidates", "reason"),
    [
        ((0, 0), (False, False), "No candidate 5' or 3' flank found in the earliest round"),
        ((0.3, 0), (True, False), "Both primer match rates below 0.40 (5'=0.30, 3'=no candidate)"),
        ((0, 0.1), (False, True), "Both primer match rates below 0.40 (5'=no candidate, 3'=0.10)"),
        ((0.3, 0.1), (True, True), "Both primer match rates below 0.40 (5'=0.30, 3'=0.10)"),
    ],
)
def test_low_match_reason_names_missing_candidates(rates, candidates, reason):
    assert (
        _low_match_reason(
            *rates, candidate_5p=candidates[0], candidate_3p=candidates[1], threshold=0.4
        )
        == reason
    )


def test_classify_passes_the_candidates_to_the_failure_reason():
    c = _classify(
        match_rate_5p=0.0,
        match_rate_3p=0.0,
        n_length_confidence=0.0,
        has_paired_split=False,
        paired_has_overlap=False,
        has_round_map=True,
        composite_confidence=0.0,
        primer_found_threshold=0.7,
        n_length_confident_threshold=0.8,
        unable_to_extract_threshold=0.4,
        status_high_cutoff=0.85,
        status_medium_cutoff=0.6,
        status_low_cutoff=0.3,
        candidate_5p=False,
        candidate_3p=False,
    )
    assert c.failure_reason == "No candidate 5' or 3' flank found in the earliest round"


def test_detect_on_random_reads_reports_no_candidate():
    import random

    from selexprep.library.detect import compute_library_report

    rng = random.Random(0)
    pool = ["".join(rng.choice("ACGT") for _ in range(60)) for _ in range(1000)]
    report = compute_library_report({0: pool, 1: pool[:]}, read_source="R1")
    assert report.status == "UNABLE_TO_INFER"
    assert report.failure_reason == "No candidate 5' or 3' flank found in the earliest round"


# ---------------------------------------------------------------------------
# 5. Low-yield warning
# ---------------------------------------------------------------------------


@needs_cutadapt
def test_low_yield_warning_names_the_opposite_orientation(tmp_path, caplog):
    forward = _fastq(tmp_path / "a.fastq.gz", [PRIMER_5P + N39 + PRIMER_3P] * 40)
    flipped = _fastq(
        tmp_path / "b.fastq.gz", [reverse_complement_iupac(PRIMER_5P + N39 + PRIMER_3P)] * 40
    )
    with caplog.at_level(logging.WARNING):
        result = runner.run_extract(
            _report(),
            [forward, flipped],
            tmp_path / "out",
            round_map={forward.name: 0, flipped.name: 1},
        )
    assert not result.skipped
    assert "reads in the opposite orientation" in caplog.text
    assert "b.fastq.gz" in caplog.text
