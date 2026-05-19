"""Unit tests for ``selexprep.library.adapters``."""

from __future__ import annotations

import pytest

from selexprep.library.adapters import (
    KNOWN_ADAPTERS,
    KNOWN_ADAPTERS_RC,
    count_adapter_hits,
    reverse_complement,
)

# ---------------------------------------------------------------------------
# reverse_complement
# ---------------------------------------------------------------------------


def test_reverse_complement_acgt() -> None:
    assert reverse_complement("ATCG") == "CGAT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement("GCAT") == "ATGC"


def test_reverse_complement_is_case_insensitive() -> None:
    assert reverse_complement("atcg") == "CGAT"
    assert reverse_complement("AtCg") == "CGAT"


def test_reverse_complement_treats_u_as_a_complement() -> None:
    # RNA primers are reported as DNA (locked plan line 296). U → A in the
    # complement (U pairs with A).
    assert reverse_complement("AUCG") == "CGAT"


def test_reverse_complement_rejects_ambiguous_bases() -> None:
    with pytest.raises(ValueError, match="non-ACGTU base 'N'"):
        reverse_complement("ANCG")
    with pytest.raises(ValueError, match="non-ACGTU base 'R'"):
        reverse_complement("ARCG")


def test_reverse_complement_is_involutive_on_acgt() -> None:
    for seq in ("ATCGATCG", "GCGCATAT", "TTAACCGG"):
        assert reverse_complement(reverse_complement(seq)) == seq


# ---------------------------------------------------------------------------
# KNOWN_ADAPTERS + RC table
# ---------------------------------------------------------------------------


def test_known_adapters_rc_table_matches_revcomp() -> None:
    for name, adapter in KNOWN_ADAPTERS.items():
        assert KNOWN_ADAPTERS_RC[name] == reverse_complement(adapter)


def test_known_adapters_includes_truseq_and_nextera() -> None:
    # v0.1 conservative set per the locked plan.
    assert "TRUSEQ_R1" in KNOWN_ADAPTERS
    assert "NEXTERA" in KNOWN_ADAPTERS
    assert KNOWN_ADAPTERS["TRUSEQ_R1"] == "AGATCGGAAGAGC"


# ---------------------------------------------------------------------------
# count_adapter_hits
# ---------------------------------------------------------------------------


def test_count_adapter_hits_finds_truseq_substring() -> None:
    truseq = KNOWN_ADAPTERS["TRUSEQ_R1"]
    seqs = [
        "ACGTACGT" + truseq + "TTTT",  # has TruSeq mid-read
        "ACGTACGTACGT",  # clean
        truseq + "GGGGCCCC",  # has TruSeq at start
    ]
    hits = count_adapter_hits(seqs)
    assert hits["TRUSEQ_R1"] == 2
    assert hits["NEXTERA"] == 0


def test_count_adapter_hits_finds_reverse_complement() -> None:
    truseq_rc = KNOWN_ADAPTERS_RC["TRUSEQ_R1"]
    seqs = [
        "ACGT" + truseq_rc + "GGGG",
        "ACGTACGTACGT",  # clean
    ]
    hits = count_adapter_hits(seqs)
    assert hits["TRUSEQ_R1"] == 1


def test_count_adapter_hits_does_not_double_count_per_sequence() -> None:
    # A sequence containing both the forward AND reverse complement still
    # counts as ONE hit for that adapter.
    truseq = KNOWN_ADAPTERS["TRUSEQ_R1"]
    truseq_rc = KNOWN_ADAPTERS_RC["TRUSEQ_R1"]
    seqs = [truseq + "AAAA" + truseq_rc]
    hits = count_adapter_hits(seqs)
    assert hits["TRUSEQ_R1"] == 1


def test_count_adapter_hits_returns_zero_for_clean_pool() -> None:
    seqs = ["ACGTACGTACGT", "TTTTAAAA", "GGGGCCCC"]
    hits = count_adapter_hits(seqs)
    assert hits == {name: 0 for name in KNOWN_ADAPTERS}


def test_count_adapter_hits_handles_empty_pool() -> None:
    hits = count_adapter_hits([])
    assert hits == {name: 0 for name in KNOWN_ADAPTERS}


def test_count_adapter_hits_uppercases_input() -> None:
    truseq = KNOWN_ADAPTERS["TRUSEQ_R1"]
    seqs = [("xx" + truseq + "yy").lower()]
    hits = count_adapter_hits(seqs)
    assert hits["TRUSEQ_R1"] == 1


def test_count_adapter_hits_rejects_nonpositive_k() -> None:
    with pytest.raises(ValueError, match="k must be positive"):
        count_adapter_hits(["ACGT"], k=0)
    with pytest.raises(ValueError, match="k must be positive"):
        count_adapter_hits(["ACGT"], k=-3)
