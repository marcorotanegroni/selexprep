"""Unit tests for selexprep.library.detect."""

from __future__ import annotations

import random
from pathlib import Path

import pandas as pd
import pytest

from selexprep.library.detect import (
    FlankResult,
    detect_flank,
    detect_from_parquet,
    detect_primers,
    earliest_round_parquet,
)


def _synthetic_pool(
    primer_5p: str | None,
    primer_3p: str | None,
    n: int = 1000,
    random_len: int = 30,
) -> list[str]:
    """Generate `n` sequences = primer_5p + random_region + primer_3p.

    The random region is deterministic per index so tests are reproducible;
    diversity is sufficient that no random k-mer dominates.
    """
    bases = "ACGT"
    seqs = []
    for i in range(n):
        # Deterministic but well-distributed random region
        random = "".join(bases[(i * 7 + j * 13) % 4] for j in range(random_len))
        seqs.append((primer_5p or "") + random + (primer_3p or ""))
    return seqs


# ----- detect_flank -----


def test_detect_flank_finds_strong_5p_consensus() -> None:
    seqs = _synthetic_pool(primer_5p="GGTAATACGACTCACTATAGGG", primer_3p=None, n=1000)
    result = detect_flank(seqs, is_prefix=True, min_seqs_for_detection=100)
    assert result.sequence is not None
    assert result.sequence.startswith("GGTAATACGACTCAC")
    assert result.confidence > 0.95
    assert result.length >= 15


def test_detect_flank_finds_strong_3p_consensus() -> None:
    seqs = _synthetic_pool(primer_5p=None, primer_3p="CCATGCATGCATGCATGCATGCATGCAT", n=1000)
    result = detect_flank(seqs, is_prefix=False, min_seqs_for_detection=100)
    assert result.sequence is not None
    assert "CCATGCATGCATGCA" in result.sequence
    assert result.confidence > 0.95


def test_detect_flank_returns_empty_when_no_consensus() -> None:
    # Truly random short reads — no shared flanks
    bases = "ACGT"
    seqs = ["".join(bases[(i * 7 + j) % 4] for j in range(40)) for i in range(1000)]
    result = detect_flank(seqs, is_prefix=True, min_seqs_for_detection=100)
    assert result.sequence is None
    assert result.confidence == 0.0


def test_detect_flank_tolerates_single_hamming_mismatch() -> None:
    """Hamming-≤1 variants still match the dominant consensus.

    Setup: 80% of reads carry the exact consensus primer; 20% carry it with
    one base flipped at a varying position. After Counter.most_common picks
    the exact form as the dominant L-mer, the variants are at Hamming = 1
    from it and therefore count as matches.
    """
    bases = "ACGT"
    seqs = []
    primer = "AAAAAAAAAAAAAAAA"  # 16 nt
    for i in range(1000):
        if i % 5 == 0:
            pos = (i // 5) % 16
            corrupted = primer[:pos] + "T" + primer[pos + 1 :]
        else:
            corrupted = primer
        random = "".join(bases[(i * 7 + j) % 4] for j in range(30))
        seqs.append(corrupted + random)
    result = detect_flank(seqs, is_prefix=True, min_seqs_for_detection=100, confidence=0.9)
    # 80% exact + 20% Hamming-1 from the consensus → 100% match
    assert result.sequence is not None
    assert result.confidence == pytest.approx(1.0)


def test_detect_primers_recovers_long_flanks_without_random_overextension() -> None:
    """Long T7-bearing constants must not be hard-capped at 30 nt.

    Regression for the Tier-1 PRJDB9110/PRJDB9111 finding: the old detector
    searched only up to 30 nt and, when max_len was raised naively, could
    absorb one adjacent random base because the whole candidate still matched
    within Hamming <= 1. The boundary should stop at the support cliff.
    """
    rng = random.Random(7)
    primer_5p = "TAATACGACTCACTATAGGGAGCAGGAGAGAGGTCAGATG"  # 40 nt
    primer_3p = "CCTATGCGTGCTAGTGTGA"  # 19 nt
    bases = "ACGT"
    seqs = [
        primer_5p + "".join(rng.choice(bases) for _ in range(30)) + primer_3p for _ in range(2_000)
    ]

    det = detect_primers(seqs, min_seqs_for_detection=100)

    assert det is not None
    assert det.primer_5p.sequence == primer_5p
    assert det.primer_5p.length == 40
    assert det.primer_3p.sequence == primer_3p
    assert det.primer_3p.length == 19
    assert det.estimated_random_region_len == 30


def test_detect_flank_stops_at_high_support_cliff_before_biased_tail() -> None:
    """A high-support primer core followed by biased random bases must not overextend.

    Regression for the PRJNA883192 R2 finding: the first 15 bases were
    a >95% consensus primer core, followed by a C-biased random-region
    shoulder. The boundary belongs at the support cliff after the core,
    not at the later point where the biased shoulder fully collapses.
    """
    primer = "TCGAGTTCAGAGCTC"  # read-level R2 prefix; revcomp = insert 3' flank
    bases = "ACGT"
    seqs: list[str] = []
    for i in range(2_000):
        biased_tail = [
            "C" if i % 5 != 0 else "A",  # 80% C
            "C" if i % 5 < 3 else "G",  # 60% C
            "C" if i % 5 < 3 else "T",  # 60% C
            "C" if i % 7 < 4 else "A",  # ~57% C
        ]
        random_tail = "".join(bases[(i * 7 + j * 13) % 4] for j in range(36))
        seqs.append(primer + "".join(biased_tail) + random_tail)

    result = detect_flank(seqs, is_prefix=True, min_seqs_for_detection=100)

    assert result.sequence == primer
    assert result.length == len(primer)
    assert result.confidence > 0.95


# ----- detect_primers -----


def test_detect_primers_both_flanks() -> None:
    seqs = _synthetic_pool(
        primer_5p="GGTAATACGACTCACTATAGGG",
        primer_3p="CCATGCATGCATGCATGCATGC",
        n=1000,
    )
    det = detect_primers(seqs, min_seqs_for_detection=100)
    assert det is not None
    assert det.primer_5p.sequence is not None
    assert det.primer_3p.sequence is not None
    assert det.estimated_random_region_len is not None
    assert det.estimated_random_region_len > 0


def test_detect_primers_only_5p() -> None:
    seqs = _synthetic_pool(primer_5p="GGTAATACGACTCACTATAGGG", primer_3p=None, n=1000)
    det = detect_primers(seqs, min_seqs_for_detection=100)
    assert det is not None
    assert det.primer_5p.sequence is not None
    assert det.primer_3p.sequence is None


def test_detect_primers_returns_none_when_too_few_sequences() -> None:
    seqs = _synthetic_pool(primer_5p="GGTAATACG", primer_3p=None, n=10)
    det = detect_primers(seqs, min_seqs_for_detection=500)
    assert det is None


# ----- detect_from_parquet -----


def test_detect_from_parquet(tmp_path: Path) -> None:
    seqs = _synthetic_pool(
        primer_5p="GGTAATACGACTCACTATAGGG",
        primer_3p="CCATGCATGCATGCATGCATGC",
        n=1000,
    )
    p = tmp_path / "round_00.counts.parquet"
    df = pd.DataFrame({"sequence": seqs, "reads": [1] * len(seqs), "rank": range(1, len(seqs) + 1)})
    df.to_parquet(p, index=False, compression="zstd")

    det = detect_from_parquet(p, top_n=1000, min_seqs_for_detection=100)
    assert det is not None
    assert det.source == p
    assert det.primer_5p.sequence is not None


def test_detect_from_parquet_default_uses_all_sequences(tmp_path: Path) -> None:
    """Default `top_n=None` must NOT subsample — every unique sequence is fed
    to the detection algorithm. This is the user's accuracy requirement: in
    a naive pool every read carries the primer, so the long tail of rare
    unique sequences also confirms the consensus."""
    seqs = _synthetic_pool(
        primer_5p="GGTAATACGACTCACTATAGGG",
        primer_3p="CCATGCATGCATGCATGCATGC",
        n=5_000,  # 5x the previous "top 1000" cap
    )
    p = tmp_path / "round_00.counts.parquet"
    df = pd.DataFrame({"sequence": seqs, "reads": [1] * len(seqs), "rank": range(1, len(seqs) + 1)})
    df.to_parquet(p, index=False, compression="zstd")

    det = detect_from_parquet(p, min_seqs_for_detection=100)  # no top_n
    assert det is not None
    assert det.n_sequences_analyzed == 5_000
    assert det.primer_5p.sequence is not None


def test_detect_from_parquet_explicit_top_n_subsamples(tmp_path: Path) -> None:
    """Passing an explicit positive `top_n` still caps the input — escape hatch
    for memory-constrained runs."""
    seqs = _synthetic_pool(primer_5p="GGTAATACG", primer_3p=None, n=2_000)
    p = tmp_path / "round_00.counts.parquet"
    df = pd.DataFrame(
        {"sequence": seqs, "reads": list(range(2_000, 0, -1)), "rank": range(1, 2_001)}
    )
    df.to_parquet(p, index=False, compression="zstd")

    det = detect_from_parquet(p, top_n=500, min_seqs_for_detection=100)
    assert det is not None
    assert det.n_sequences_analyzed == 500


# ----- earliest_round_parquet -----


def test_earliest_round_parquet_picks_lowest(tmp_path: Path) -> None:
    bp = tmp_path / "PRJ1"
    bp.mkdir()
    (bp / "round_03.counts.parquet").touch()
    (bp / "round_00.counts.parquet").touch()
    (bp / "round_07.counts.parquet").touch()
    chosen = earliest_round_parquet(bp)
    assert chosen is not None
    assert chosen.name == "round_00.counts.parquet"


def test_earliest_round_parquet_returns_none_when_empty(tmp_path: Path) -> None:
    bp = tmp_path / "PRJ1"
    bp.mkdir()
    assert earliest_round_parquet(bp) is None


# ----- FlankResult shape -----


def test_flank_result_dataclass_unpacking() -> None:
    fr = FlankResult(sequence="ACGT", length=4, confidence=0.95)
    assert fr.sequence == "ACGT"
    assert fr.length == 4
    assert fr.confidence == 0.95
