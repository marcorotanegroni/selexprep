"""Tests for ``selexprep.benchmark.equivalence``.

Coverage targets every kind in the ``EquivalenceKind`` Literal plus the
edge cases the metric aggregator depends on (None observed, IUPAC truth,
barcode stripping).
"""

from __future__ import annotations

from selexprep.benchmark.equivalence import EquivalenceResult, primer_equivalent


def test_exact_match_lowercase_input_returns_exact() -> None:
    r = primer_equivalent("acgtac", "ACGTAC")
    assert r == EquivalenceResult(matched=True, equivalence_kind="EXACT")


def test_revcomp_match_returns_revcomp_kind() -> None:
    # truth=ACGT → revcomp=ACGT (palindrome) → use a non-palindrome
    truth = "ACGTACG"
    rc = "CGTACGT"
    r = primer_equivalent(rc, truth)
    assert r.matched is True
    assert r.equivalence_kind == "REVCOMP"


def test_revcomp_disabled_drops_to_partial_or_mismatch() -> None:
    truth = "ACGTACG"
    rc = "CGTACGT"
    r = primer_equivalent(rc, truth, allow_revcomp=False)
    assert r.matched is False
    # ACGTACG and CGTACGT share no useful prefix/suffix → MISMATCH
    assert r.equivalence_kind == "MISMATCH"


def test_ut_normalized_match() -> None:
    # truth carries U (RNA); observed is T (DNA after sequencing)
    r = primer_equivalent("ACGTACGT", "ACGUACGU")
    assert r.matched is True
    assert r.equivalence_kind == "U_T_NORMALIZED"


def test_ut_disabled_with_u_in_truth_misses() -> None:
    # When U-T is disabled and revcomp is also disabled, an RNA-truth /
    # DNA-observed pair falls through to MISMATCH.
    r = primer_equivalent("ACGTACGT", "ACGUACGU", allow_ut=False, allow_revcomp=False)
    assert r.matched is False
    assert r.equivalence_kind == "MISMATCH"


def test_barcode_prefix_strip_match() -> None:
    r = primer_equivalent("AAAACGTACG", "ACGTACG", strip_barcodes=("AAA",))
    assert r.matched is True
    assert r.equivalence_kind == "BARCODE_STRIPPED"
    assert "AAA" in r.notes


def test_barcode_strip_picks_longest_matching_prefix() -> None:
    # Both "AAA" and "AAAACGT" prefix the observed; the longer barcode wins
    # so the stripped tail is the right "true primer".
    r = primer_equivalent("AAAACGTGCGT", "GCGT", strip_barcodes=("AAA", "AAAACGT"))
    assert r.matched is True
    assert r.equivalence_kind == "BARCODE_STRIPPED"
    assert "AAAACGT" in r.notes


def test_partial_5p_returns_unmatched_with_partial_kind() -> None:
    # Observed is a strict prefix of truth → PARTIAL_5P
    r = primer_equivalent("ACGTAC", "ACGTACGT")
    assert r.matched is False
    assert r.equivalence_kind == "PARTIAL_5P"


def test_partial_3p_returns_unmatched_with_partial_kind() -> None:
    # Observed is a suffix of truth — one captured only the 3' portion
    # of the other. ``_is_partial_3p`` requires one to be a strict suffix
    # of the other; shared-but-distinct suffixes alone fall through to
    # MISMATCH (which is also the correct accounting).
    r = primer_equivalent("ACGT", "GGGACGT")
    assert r.matched is False
    assert r.equivalence_kind == "PARTIAL_3P"


def test_mismatch_when_no_relationship() -> None:
    r = primer_equivalent("AAAAA", "GGGGG")
    assert r.matched is False
    assert r.equivalence_kind == "MISMATCH"


def test_iupac_in_truth_returns_unsupported() -> None:
    r = primer_equivalent("ACGTACGT", "ACGTACGN")
    assert r.matched is False
    assert r.equivalence_kind == "IUPAC_UNSUPPORTED"
    assert "IUPAC" in r.notes


def test_iupac_in_truth_lowercase_still_caught() -> None:
    # IUPAC detection is case-insensitive.
    r = primer_equivalent("ACGT", "acgrt")
    assert r.equivalence_kind == "IUPAC_UNSUPPORTED"


def test_observed_none_returns_mismatch() -> None:
    r = primer_equivalent(None, "ACGTACGT")
    assert r.matched is False
    assert r.equivalence_kind == "MISMATCH"


def test_empty_truth_returns_mismatch() -> None:
    r = primer_equivalent("ACGT", "")
    assert r.matched is False
    assert r.equivalence_kind == "MISMATCH"


def test_empty_observed_returns_mismatch() -> None:
    r = primer_equivalent("", "ACGT")
    assert r.equivalence_kind == "MISMATCH"


def test_case_folding_works_for_revcomp() -> None:
    # mixed-case observed and truth still match via revcomp
    r = primer_equivalent("CgTaCgT", "ACGTACG")
    assert r.equivalence_kind == "REVCOMP"
    assert r.matched is True
