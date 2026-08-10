"""Refusal when the inferred random region has zero length.

Found by the adapter-control arm of the Tier-1 benchmark on PRJDB7022, a
non-SELEX small-RNA deposit in which 4.7M of 5.7M reads are one identical
51-mer. Positional consensus sees ~100% conservation at every position, so
there is no constant/random boundary to stop the inward walk: it consumes the
whole read from both ends, the two flanks meet, and the modal insert is 0 nt.

Before the guard, ``detect`` reported that as a healthy two-sided library
(``full_insert_recovered=True``, ``required_action="NONE"``, status MEDIUM),
which would have carried ``run`` through to ``count`` and written a table of
zero-length sequences without raising anything.

The guard must not touch the paired split-primer state, where the modal length
is ``None`` rather than 0 because no single read spans the insert.
"""

from selexprep.library.adapters import reverse_complement
from selexprep.library.detect import compute_library_report

PRIMER_5P = "GGTAATACGACTCACTATAGGG"
PRIMER_3P = "CCATGCATGCATGCATGCAT"

# The dominant read of PRJDB7022 (DRR129849 and mates): a single small-RNA
# insert followed by the TruSeq Small RNA 3' adapter, 51 nt total.
PRJDB7022_DOMINANT_READ = "TGCTTGGACTACATATGGTTGAGGGTTGTATGGAATTCTCGGGTGCCAAGG"


def _pool(primer_5p: str | None, primer_3p: str | None, *, n=600, random_len=30, offset=0):
    bases = "ACGT"
    return [
        (primer_5p or "")
        + "".join(bases[((i + offset) * 7 + j * 13) % 4] for j in range(random_len))
        + (primer_3p or "")
        for i in range(n)
    ]


def _rounds(primer_5p, primer_3p, *, random_len=30):
    return {
        r: _pool(primer_5p, primer_3p, random_len=random_len, offset=r * 1000) for r in range(3)
    }


class TestMonoclonalPoolIsRefused:
    def test_identical_reads_get_no_primer_call(self):
        pools = {r: [PRJDB7022_DOMINANT_READ] * 600 for r in range(3)}
        report = compute_library_report(pools, read_source="R1")

        assert report.status == "UNABLE_TO_INFER"
        assert report.primer_5p is None
        assert report.primer_3p is None
        assert report.extraction_mode == "UNABLE_TO_EXTRACT"
        assert report.required_action == "MANUAL_PRIMERS_REQUIRED"

    def test_refusal_states_the_reason(self):
        pools = {r: [PRJDB7022_DOMINANT_READ] * 600 for r in range(3)}
        report = compute_library_report(pools, read_source="R1")

        assert report.failure_reason
        assert "0 nt" in report.failure_reason
        # points the user at the escape hatch rather than just failing
        assert "--override-primer-5p" in report.failure_reason

    def test_converged_selex_pool_is_refused_the_same_way(self):
        """A late-round SELEX pool that has collapsed onto one winner produces
        the identical signal, and must be handled identically.

        The randomised region of a converged pool is as conserved as the
        primers flanking it, so no boundary exists to detect and the modal
        insert is 0 — exactly what a non-SELEX library looks like. `detect`
        cannot tell the two apart from reads alone, so it must not claim to:
        the refusal says a variable region is absent, never that the input is
        not SELEX. The user who knows better supplies the primers, which is
        what ``required_action`` asks for and why that advice is right here
        rather than merely harmless.
        """
        winner = PRIMER_5P + "ACGTACGTACGTACGTACGTACGTACGTAC" + PRIMER_3P
        pools = {r: [winner] * 600 for r in range(3)}

        report = compute_library_report(pools, read_source="R1")

        assert report.status == "UNABLE_TO_INFER"
        assert report.required_action == "MANUAL_PRIMERS_REQUIRED"
        assert report.primer_5p is None and report.primer_3p is None
        # the wording must stay agnostic about what the input actually is
        assert "not SELEX" not in report.failure_reason
        assert "If this really is a SELEX library" in report.failure_reason


class TestHealthyLibrariesAreUntouched:
    """The guard keys on a modal insert of exactly 0, which a real library
    cannot produce. Nothing with a randomised region may change."""

    def test_two_sided_library_still_recovered(self):
        report = compute_library_report(_rounds(PRIMER_5P, PRIMER_3P), read_source="R1")

        assert report.extraction_mode == "BOTH_PRIMERS_SINGLE_READ"
        assert report.required_action == "NONE"
        assert report.n_length_mode == 30
        assert report.primer_5p is not None and report.primer_3p is not None

    def test_one_sided_library_still_recovered(self):
        report = compute_library_report(_rounds(PRIMER_5P, None), read_source="R1")

        assert report.extraction_mode == "FIVE_PRIME_ONLY"
        assert report.primer_5p is not None

    def test_paired_split_none_length_is_not_caught(self):
        """``n_length_mode`` is None here, not 0 — a different state. A guard
        written as ``if not n_mode`` would swallow this one."""
        r1 = _rounds(PRIMER_5P, None, random_len=80)
        r2 = _rounds(reverse_complement(PRIMER_3P), None, random_len=80)

        report = compute_library_report(r1, read_source="R1_AND_R2", paired_mate_streams=r2)

        assert report.n_length_mode is None
        assert report.extraction_mode == "PAIRED_END_SPLIT_PRIMERS"
        assert report.required_action == "READ_MERGING_RECOMMENDED"
