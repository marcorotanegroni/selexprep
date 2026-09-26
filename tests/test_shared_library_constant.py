"""Library constant = the part of a flank every round shares (`_shared_library_constant`).

The flank is called on the earliest round, from the read edge inward, so it
takes in whatever that round's reads carry outside the library constant. Found
on the Tier-1 benchmark with an all-read check:

- PRJNA809588: every run carries its own inline tag at both ends. The call took
  the earliest run's tags, and ``extract`` then kept under 1% of the reads of
  seven of the ten runs.
- PRJEB62495: one base of the 3' constant differs between runs.
- PRJNA1395820: a heterogeneous base sits ahead of the 5' constant.

The constant may lose sequence at its outer edge only: the inner core, and with
it the boundary with the random region, never moves.
"""

from __future__ import annotations

import random

from selexprep.library.detect import (
    CONSTANT_CORE_LEN,
    _shared_library_constant,
    compute_library_report,
)

C5 = "AGGACGATGCGGTACC"  # 16 nt library constants
C3 = "CAGACGACTCGCTGAGGATC"
N = 40


def _reads(n: int, *, tag5: str = "", tag3: str = "", c5: str = C5, c3: str = C3, seed: int = 0):
    rng = random.Random(seed)
    return [tag5 + c5 + "".join(rng.choice("ACGT") for _ in range(N)) + c3 + tag3 for _ in range(n)]


class TestUnchangedWithoutRunSpecificSequence:
    def test_flank_shared_by_every_round_is_kept_whole(self):
        pools = [_reads(1000, seed=s) for s in range(3)]
        assert _shared_library_constant(C5, pools, is_prefix=True) == (C5, [0, 0, 0])
        assert _shared_library_constant(C3, pools, is_prefix=False) == (C3, [0, 0, 0])

    def test_report_matches_a_clean_library(self):
        rounds = {r: _reads(1500, seed=r) for r in range(3)}
        report = compute_library_report(rounds, read_source="R1")
        assert (report.primer_5p, report.primer_3p) == (C5, C3)
        assert report.n_length_mode == N
        assert report.extraction_mode == "BOTH_PRIMERS_SINGLE_READ"


class TestRunSpecificTagsAreLeftOut:
    """PRJNA809588: each run has its own tags, of different lengths."""

    TAGS = (("ACTGAG", "GCATG"), ("TAAGCG", "CTGTC"), ("GTCAATG", "GTTAGATG"))

    def _rounds(self):
        return {r: _reads(1500, tag5=t5, tag3=t3, seed=r) for r, (t5, t3) in enumerate(self.TAGS)}

    def test_constants_are_the_shared_part(self):
        report = compute_library_report(self._rounds(), read_source="R1")
        # All three 5' tags end in G, so the shared 5' part is that G + C5.
        assert report.primer_5p == "G" + C5
        # The 3' tags share nothing, so the shared 3' part is C3 exactly.
        assert report.primer_3p == C3

    def test_boundary_and_signals_are_measured_where_the_constant_is(self):
        report = compute_library_report(self._rounds(), read_source="R1")
        assert report.n_length_mode == N
        assert report.position_consistency_5p > 0.95
        assert report.position_consistency_3p > 0.95
        assert report.extraction_mode == "BOTH_PRIMERS_SINGLE_READ"
        assert report.status == "HIGH"
        assert report.variants_5p[0][0] == "G" + C5

    def test_offsets_say_where_each_round_carries_the_constant(self):
        pools = list(self._rounds().values())
        constant, offsets = _shared_library_constant("ACTGAG" + C5, pools, is_prefix=True)
        assert constant == "G" + C5
        # Start of "G"+C5: one base before the end of each round's tag.
        assert offsets == [5, 5, 6]


def test_a_base_that_differs_between_rounds_is_left_out():
    """PRJEB62495: rounds disagree on the 3' constant's outermost base."""
    rounds = {
        0: _reads(1500, c3=C3[:-1] + "A", seed=0),
        1: _reads(1500, seed=1),
        2: _reads(1500, seed=2),
    }
    report = compute_library_report(rounds, read_source="R1")
    assert report.primer_3p == C3[:-1]
    assert report.n_length_mode == N


def test_heterogeneous_leading_base_is_left_out_even_in_one_round():
    """PRJNA1395820: a base ahead of the construct with no majority."""
    rng = random.Random(7)
    reads = [rng.choice("GAC") + r for r in _reads(2000, seed=7)]
    report = compute_library_report({9: reads}, read_source="R1")
    assert report.primer_5p == C5
    assert report.n_length_mode == N


def test_round_that_does_not_carry_the_core_does_not_vote():
    pools = [_reads(1000, tag5="ACTGAG", seed=0), _reads(1000, c5="T" * 16, seed=1)]
    constant, offsets = _shared_library_constant("ACTGAG" + C5, pools, is_prefix=True)
    assert constant == "ACTGAG" + C5
    assert offsets == [0, None]


def test_earliest_round_that_does_not_anchor_leaves_the_flank_alone():
    pools = [_reads(1000, c5="T" * 16, seed=0), _reads(1000, seed=1)]
    assert _shared_library_constant(C5, pools, is_prefix=True) == (C5, [0, 0])


def test_the_inner_core_is_never_trimmed():
    """Rounds that share nothing outside the core keep exactly the core."""
    core = C5[-CONSTANT_CORE_LEN:]
    pools = [_reads(1000, c5="AAAA" + core, seed=0), _reads(1000, c5="CCCC" + core, seed=1)]
    constant, _ = _shared_library_constant("AAAA" + core, pools, is_prefix=True)
    assert constant == core
