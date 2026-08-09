"""Outer-edge core rescue (`_high_support_core`).

Regression cover for the failure mode found on PRJEB62495 and PRJNA315881: the
positional consensus walks inward from the *read edge*, so constant technical
sequence sitting outside the library constant (a truncated sequencing adapter,
an index remnant, the T7 start-G) is absorbed into the called flank. That
material is weakly conserved, the whole-primer ``match_rate`` collapses below
the primer-found threshold, and a two-sided library is downgraded to one-sided
extraction.

The rescue must trim only the *outer* edge, must never touch the inner
(random-region) boundary, and must leave well-supported flanks exactly as they
are.
"""

from selexprep.library.detect import (
    CORE_MIN_LEN,
    _high_support_core,
)


def _supports(*values: float) -> list[float]:
    return list(values)


class TestNoOpOnHealthyFlanks:
    """A uniformly well-supported flank must come back untouched (``None``)."""

    def test_healthy_prefix_is_not_trimmed(self):
        seq = "ACGTACGTACGTACGTACGT"
        assert _high_support_core(seq, [0.99] * len(seq), is_prefix=True) is None

    def test_healthy_suffix_is_not_trimmed(self):
        seq = "ACGTACGTACGTACGTACGT"
        assert _high_support_core(seq, [0.99] * len(seq), is_prefix=False) is None

    def test_support_exactly_at_floor_is_kept(self):
        # 0.90 is not < 0.90 — boundary case must not trim.
        seq = "ACGTACGTACGTACGT"
        assert _high_support_core(seq, [0.90] * len(seq), is_prefix=False) is None


# Per-position support actually measured on PRJEB62495 round_04 (368,381 reads),
# ordered 5'→3' like the called flank. The outer six positions carry a truncated
# TruSeq remnant; the 0.772 at index 3 is an isolated dip inside the genuine
# constant and must NOT drag the cut inward.
PRJEB62495_SUPPORTS = [
    0.921,
    0.932,
    0.931,
    0.772,
    0.932,
    0.978,
    0.923,
    0.920,
    0.911,
    0.906,
    0.903,
    0.902,
    0.898,
    0.899,
    0.884,
    0.875,
    0.976,
    0.909,
    0.901,
    0.612,
    0.483,
    0.886,
    0.896,
    0.897,
    0.749,
]
PRJEB62495_FLANK = "CGTGGTTACAGTCAGAGGACAGATT"


class TestTrimsOuterEdgeOnly:
    def test_real_prjeb62495_profile_yields_conserved_core(self):
        # Regression on the measured case: the 25 nt flank scores match_rate
        # 0.535 (below the 0.70 primer-found threshold) while the 19 nt core
        # scores 0.910, which is what unlocks two-sided extraction.
        assert (
            _high_support_core(PRJEB62495_FLANK, PRJEB62495_SUPPORTS, is_prefix=False)
            == "CGTGGTTACAGTCAGAGGA"
        )

    def test_interior_dip_alone_does_not_trim(self):
        # Same isolated 0.772 dip, but a clean outer edge: nothing to rescue.
        supports = list(PRJEB62495_SUPPORTS)
        supports[-6:] = [0.95] * 6
        assert _high_support_core(PRJEB62495_FLANK, supports, is_prefix=False) is None

    def test_prefix_drops_weak_head(self):
        # 5' flank: T7 start-G / index remnant ahead of the library constant.
        head = "GAC"
        core_seq = "TACACTGCACTGCGTTAGAG"
        seq = head + core_seq
        supports = _supports(0.55, 0.62, 0.71, *([0.99] * len(core_seq)))
        assert _high_support_core(seq, supports, is_prefix=True) == core_seq

    def test_inner_boundary_is_never_moved(self):
        # A weak position on the *inner* side must not cause a trim from the
        # outer side — the random-region boundary is off limits.
        seq = "ACGTACGTACGTACGTACGT"
        supports = _supports(*([0.99] * (len(seq) - 1)), 0.10)
        # Weak position is the inner edge for a prefix flank -> nothing to do.
        assert _high_support_core(seq, supports, is_prefix=True) is None

    def test_scan_stops_at_stable_run(self):
        # Weak outer pair, then a stable run, then another weak position: the
        # cut must stay on the outer pair and ignore what lies past the run.
        seq = "AA" + "CGTGGTTACAGTCAGAGGA"  # 21 nt
        supports = _supports(0.50, 0.60, 0.99, 0.99, 0.99, 0.40, *([0.99] * 15))
        assert _high_support_core(seq, supports, is_prefix=True) == "CGTGGTTACAGTCAGAGGA"


class TestGuardrails:
    """The specificity arm must keep making no primer call at all."""

    def test_core_below_min_len_is_refused(self):
        seq = "ACGTACGTACGTACGT"
        # Weak run leaves fewer than CORE_MIN_LEN well-supported bases.
        keep = CORE_MIN_LEN - 1
        cut = len(seq) - keep
        supports = _supports(*([0.20] * cut), *([0.99] * keep))
        assert _high_support_core(seq, supports, is_prefix=True) is None

    def test_no_sequence_is_refused(self):
        assert _high_support_core(None, [0.99, 0.99], is_prefix=True) is None

    def test_mismatched_supports_length_is_refused(self):
        assert _high_support_core("ACGTACGTACGTACGT", [0.99, 0.99], is_prefix=True) is None

    def test_all_weak_is_refused(self):
        # A pool with no conserved flank at all (pre-trimmed random region)
        # must not be rescued into a primer call.
        seq = "ACGTACGTACGTACGTACGT"
        assert _high_support_core(seq, [0.25] * len(seq), is_prefix=False) is None
