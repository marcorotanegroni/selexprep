"""Unit tests for selexprep.qc.readiness."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from selexprep.qc.readiness import (
    FAIL,
    PASS,
    TRUSEQ_R1_PREFIX,
    WARN,
    BPReport,
    alphabet_violations,
    gc_content,
    length_distribution,
    positional_entropy,
    positional_top_base,
    review_bioproject,
    singleton_fraction,
    top_kmers,
    truseq_contamination,
)

# ----- pure primitives -----


def test_length_distribution_basic() -> None:
    mode_L, mode_frac, top5 = length_distribution(["AAAA", "AAAA", "AAAA", "CC", "GGG"])
    assert mode_L == 4
    assert mode_frac == pytest.approx(0.6)
    assert top5[0] == (4, pytest.approx(0.6))


def test_length_distribution_empty() -> None:
    mode_L, mode_frac, top5 = length_distribution([])
    assert mode_L == 0
    assert mode_frac == 0.0
    assert top5 == []


def test_alphabet_violations_counts_non_acgt() -> None:
    assert alphabet_violations(["ACGT", "ACGN", "RYKM"]) == 2


def test_gc_content_balanced() -> None:
    mean, _std = gc_content(["GCGC", "ATAT", "GCAT"])
    assert mean == pytest.approx(0.5)


def test_positional_top_base_constant_at_position_zero() -> None:
    out = positional_top_base(["AAAA", "AAGT", "AACC"], from_end=False, n_pos=4)
    # Position 0 is always 'A'
    assert out[0] == ("A", 1.0)


def test_positional_entropy_zero_when_constant() -> None:
    H = positional_entropy(["AAAA", "AAAA", "AAAA"], from_end=False, n_pos=4)
    assert all(h == 0.0 for h in H)


def test_positional_entropy_max_when_uniform() -> None:
    H = positional_entropy(["ACGT", "TGCA", "GTAC", "CATG"], from_end=False, n_pos=4)
    # Each position has equal A/C/G/T — H = log2(4) = 2
    assert all(h == pytest.approx(2.0) for h in H)


def test_top_kmers_returns_normalized_fractions() -> None:
    seqs = ["AAAAAA", "AAAAAA"]
    counts = [10, 5]
    out = top_kmers(seqs, counts, k=3, top_n=1)
    assert out[0][0] == "AAA"
    assert out[0][1] == pytest.approx(1.0)


def test_top_kmers_skips_short_sequences() -> None:
    out = top_kmers(["ACG"], [1], k=4)
    assert out == []


def test_truseq_contamination_detects_prefix() -> None:
    seqs = ["ACGT", "ACGT" + TRUSEQ_R1_PREFIX, "GGGG"]
    counts = [10, 100, 5]
    uniq_frac, reads_frac, n_reads = truseq_contamination(seqs, counts)
    assert n_reads == 100
    assert uniq_frac == pytest.approx(1 / 3)
    assert reads_frac == pytest.approx(100 / 115)


def test_singleton_fraction() -> None:
    assert singleton_fraction([1, 1, 1, 2, 5]) == pytest.approx(0.6)
    assert singleton_fraction([]) == 0.0


# ----- BPReport semantics -----


def test_bp_report_worst_is_fail_when_any_fail() -> None:
    r = BPReport(bp_id="PRJ1", tag="standard")
    r.add("pre", PASS)
    r.add("alphabet", WARN)
    r.add("trim_seq", FAIL)
    assert r.worst == FAIL


def test_bp_report_worst_is_warn_when_no_fail() -> None:
    r = BPReport(bp_id="PRJ1", tag="standard")
    r.add("pre", PASS)
    r.add("alphabet", WARN)
    assert r.worst == WARN


def test_bp_report_worst_pass_when_all_pass() -> None:
    r = BPReport(bp_id="PRJ1", tag="standard")
    r.add("pre", PASS)
    r.add("alphabet", PASS)
    assert r.worst == PASS


def test_bp_report_to_dict_roundtrips_json() -> None:
    r = BPReport(bp_id="PRJ1", tag="standard")
    r.add("pre", PASS, "all good")
    r.stats["n_reads"] = 1000
    d = r.to_dict()
    assert d["bp_id"] == "PRJ1"
    assert d["worst"] == PASS
    assert d["checks"][0]["section"] == "pre"
    assert d["stats"]["n_reads"] == 1000
    # Must be JSON-serialisable
    json.dumps(d)


# ----- review_bioproject integration -----


def _make_fixture_bp(
    bp_dir: Path,
    primer_5p: str = "GGTAATACGACTCACTATAGGG",
    primer_3p: str = "CCATGCATGCATGCATGCATGC",
    rrl: int = 30,
    n_rounds: int = 3,
    n_unique_per_round: int = 100,
) -> None:
    """Create a minimal processed BP directory with synthetic round parquets,
    a final enrich parquet, summary.json and cluster_stats.json."""
    bp_dir.mkdir(parents=True, exist_ok=True)
    bases = "ACGT"

    for r in range(n_rounds):
        sequences = []
        for i in range(n_unique_per_round):
            random_region = "".join(bases[(i * 7 + j * 13 + r) % 4] for j in range(rrl))
            sequences.append(random_region)
        df = pd.DataFrame(
            {
                "sequence": sequences,
                "reads": list(range(n_unique_per_round, 0, -1)),
                "rank": range(1, n_unique_per_round + 1),
            }
        )
        df.to_parquet(bp_dir / f"round_{r:02d}.counts.parquet", index=False, compression="zstd")
        # Cluster parquet — required by section_pre
        df.to_parquet(bp_dir / f"round_{r:02d}.clusters.parquet", index=False, compression="zstd")

    # Enrich: synthetic log2FC distribution with a clear winner
    n_seqs = 200
    enrich_df = pd.DataFrame(
        {
            "sequence": [f"SEQ{i}" for i in range(n_seqs)],
            "log2FC": [3.0 - i * 0.02 for i in range(n_seqs)],
        }
    )
    enrich_df.to_parquet(
        bp_dir / f"enrich_round_00_to_round_{n_rounds - 1:02d}.parquet",
        index=False,
        compression="zstd",
    )

    (bp_dir / "summary.json").write_text(
        json.dumps({"primer_5p": primer_5p, "primer_3p": primer_3p})
    )
    (bp_dir / "cluster_stats.json").write_text(json.dumps({}))


def test_review_bioproject_runs_all_sections(tmp_path: Path) -> None:
    bp = tmp_path / "PRJTEST1"
    _make_fixture_bp(bp)

    report = review_bioproject(
        bp,
        primer_5p="GGTAATACGACTCACTATAGGG",
        primer_3p="CCATGCATGCATGCATGCATGC",
        random_region_len=30,
        tag="standard",
        rrl_source="paper",
    )

    sections_run = {c.section for c in report.checks}
    # All eight diagnostic sections should fire (pre + 7 others)
    assert "pre" in sections_run
    assert "alphabet" in sections_run
    assert "lengths" in sections_run
    assert "trim_seq" in sections_run
    assert "composition" in sections_run
    assert "diversity" in sections_run
    assert "selection" in sections_run
    assert "consistency" in sections_run


def test_review_bioproject_missing_artifacts_fails_at_pre(tmp_path: Path) -> None:
    empty_bp = tmp_path / "EMPTY"
    empty_bp.mkdir()
    report = review_bioproject(empty_bp, primer_5p="X", primer_3p="X", random_region_len=30)
    # Only pre fires when artifacts are missing
    pre_check = next(c for c in report.checks if c.section == "pre")
    assert pre_check.status == FAIL
    assert "missing" in pre_check.detail


def test_review_bioproject_flags_short_primer(tmp_path: Path) -> None:
    bp = tmp_path / "PRJSHORT"
    _make_fixture_bp(bp)
    report = review_bioproject(
        bp,
        primer_5p="AC",  # < PRIMER_MIN_LEN
        primer_3p="CCATGCATGCATGCATGCATGC",
        random_region_len=30,
    )
    # Should have a second 'pre' check with FAIL ("incomplete seed")
    pre_checks = [c for c in report.checks if c.section == "pre"]
    assert any("incomplete seed" in c.detail for c in pre_checks)


def test_review_bioproject_no_primer_3p_for_multiplexed_origin(tmp_path: Path) -> None:
    """Multiplexed-origin tag accepts missing primer_3p (5p-only by design)."""
    bp = tmp_path / "PRJMULTI"
    _make_fixture_bp(bp, primer_3p="")
    report = review_bioproject(
        bp,
        primer_5p="GGTAATACGACTCACTATAGGG",
        primer_3p="",
        random_region_len=30,
        tag="multiplexed_origin",
    )
    pre_checks = [c for c in report.checks if c.section == "pre"]
    # Should NOT fail just because primer_3p is missing
    assert not any(c.status == FAIL and "incomplete seed" in c.detail for c in pre_checks)
