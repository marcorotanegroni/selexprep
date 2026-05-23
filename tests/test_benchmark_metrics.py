"""Tests for ``selexprep.benchmark.metrics``.

Covers the metric primitives (primer recovery, N-length recovery, count
correlation, extraction-mode distribution) plus three Codex-amendment
regression tests:

- ``test_aggregate_skips_unverified_rows_with_warning`` — amendment 1
- ``test_count_correlation_uses_union_with_zero_fill`` — amendment 5
- ``test_count_correlation_top_k_secondary_diagnostic`` — amendment 5
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import pytest

from selexprep.benchmark.metrics import (
    BenchmarkRow,
    aggregate_metrics,
    compute_count_correlation,
    compute_extraction_mode_distribution,
    compute_n_length_recovery,
    compute_pair_recovery_by_status,
    compute_primer_recovery,
    compute_required_action_distribution,
    compute_safe_failure_rate,
    load_ground_truth,
    write_metrics_json,
)
from selexprep.library.report import LibraryReport


def _make_lr(**overrides: Any) -> LibraryReport:
    """Minimal LibraryReport for tests; overrides apply to fields you care about."""
    defaults: dict[str, Any] = {
        "primer_5p": "ACGTACGT",
        "primer_3p": "TGCATGCA",
        "variants_5p": [],
        "variants_3p": [],
        "known_adapter_hits": {},
        "extraction_mode": "BOTH_PRIMERS_SINGLE_READ",
        "full_insert_recovered": True,
        "read_source": "R1",
        "required_action": "NONE",
        "orientation": "FORWARD",
        "n_length_mode": 30,
        "n_length_distribution": {},
        "n_length_confidence": 0.9,
        "match_rate_5p": 0.95,
        "match_rate_3p": 0.95,
        "position_consistency_5p": 0.95,
        "position_consistency_3p": 0.95,
        "read_fraction_used_for_inference": 1.0,
        "sampling_seed": 42,
        "confidence": 0.9,
        "status": "HIGH",
        "failure_reason": None,
    }
    defaults.update(overrides)
    return LibraryReport(**defaults)


def _row(accession: str = "PRJ_TEST", *, verified: bool = True, **overrides: Any) -> BenchmarkRow:
    defaults: dict[str, Any] = {
        "accession": accession,
        "library_kind": "DNA",
        "target_kind": "protein",
        "primer_5p_truth": "ACGTACGT",
        "primer_3p_truth": "TGCATGCA",
        "n_length_truth": 30,
        "paper_doi": "10.1000/test",
        "paper_pmid": "0",
        "round_map_source": "auto",
        "round_map_path": "",
        "verified": verified,
        "notes": "synthetic test row",
        "library_report": _make_lr(),
        "observed_counts": None,
        "reference_counts": None,
    }
    defaults.update(overrides)
    return BenchmarkRow(**defaults)


# ---------------------------------------------------------------------------
# Primer recovery
# ---------------------------------------------------------------------------


def test_primer_recovery_exact_match_on_both_sides() -> None:
    rows = [_row()]
    report = compute_primer_recovery(rows)
    assert report.n_evaluated == 1
    assert report.counts_5p == {"EXACT": 1}
    assert report.counts_3p == {"EXACT": 1}
    assert report.counts_by_status_5p == {"HIGH": {"EXACT": 1}}


def test_primer_recovery_mismatch_when_library_report_missing() -> None:
    rows = [_row(library_report=None)]
    report = compute_primer_recovery(rows)
    assert report.counts_5p == {"MISMATCH": 1}
    assert report.counts_3p == {"MISMATCH": 1}
    assert "NO_REPORT" in report.counts_by_status_5p


def test_primer_recovery_revcomp_kind_recorded() -> None:
    # Truth = ACGTACG, observed = CGTACGT (revcomp)
    rows = [
        _row(
            primer_5p_truth="ACGTACG",
            library_report=_make_lr(primer_5p="CGTACGT", primer_3p="TGCATGCA"),
        )
    ]
    report = compute_primer_recovery(rows)
    assert report.counts_5p == {"REVCOMP": 1}


# ---------------------------------------------------------------------------
# N-length recovery
# ---------------------------------------------------------------------------


def test_n_length_recovery_within_tolerance() -> None:
    rows = [
        _row(n_length_truth=30, library_report=_make_lr(n_length_mode=29)),
        _row(n_length_truth=30, library_report=_make_lr(n_length_mode=33)),  # out
        _row(n_length_truth=30, library_report=_make_lr(n_length_mode=None)),  # unmeasurable
    ]
    report = compute_n_length_recovery(rows, tolerance=2)
    assert report.n_in_tolerance == 1
    assert report.n_out_of_tolerance == 1
    assert report.n_unmeasurable == 1


def test_n_length_recovery_zero_truth_treated_unmeasurable() -> None:
    rows = [_row(n_length_truth=0, library_report=_make_lr(n_length_mode=30))]
    report = compute_n_length_recovery(rows)
    assert report.n_unmeasurable == 1


# ---------------------------------------------------------------------------
# Extraction mode distribution
# ---------------------------------------------------------------------------


def test_extraction_mode_distribution_counts() -> None:
    rows = [
        _row(library_report=_make_lr(extraction_mode="BOTH_PRIMERS_SINGLE_READ")),
        _row(library_report=_make_lr(extraction_mode="BOTH_PRIMERS_SINGLE_READ")),
        _row(library_report=_make_lr(extraction_mode="PAIRED_END_SPLIT_PRIMERS")),
        _row(library_report=None),
    ]
    dist = compute_extraction_mode_distribution(rows)
    assert dist.counts == {
        "BOTH_PRIMERS_SINGLE_READ": 2,
        "PAIRED_END_SPLIT_PRIMERS": 1,
        "NO_REPORT": 1,
    }


# ---------------------------------------------------------------------------
# Count correlation — Codex amendment 5 regressions
# ---------------------------------------------------------------------------


def test_count_correlation_uses_union_with_zero_fill() -> None:
    """Codex amendment 5: correlation must be computed on the union of sequences.

    Demonstration: a shared subset with perfectly-agreeing counts +
    disjoint unique-to-each entries.

    - Intersection-only would consider just the shared subset:
      Pearson(obs_shared, ref_shared) = 1.0 (perfect agreement).
    - Union + zero-fill considers the disjoint entries too:
      Pearson(obs_union, ref_union) < 1.0 — honest accounting of the
      fact that each tool saw sequences the other didn't.

    The function must return the union (< 1.0) Pearson, not the
    intersection (= 1.0) Pearson.
    """
    rows = [
        _row(
            accession="PRJ_X",
            # Two sequences with identical counts in both, plus one
            # unique to each side.
            observed_counts={"shared_A": 100, "shared_B": 50, "obs_only": 30},
            reference_counts={"shared_A": 100, "shared_B": 50, "ref_only": 30},
        )
    ]
    report = compute_count_correlation(rows, top_k=2)
    assert report.n_union == 4
    assert report.n_observed_only == 1
    assert report.n_reference_only == 1
    assert report.n_in_both == 2
    # If the implementation had used intersection-only, pearson would
    # be exactly 1.0 here (perfect agreement on shared_A + shared_B).
    # The union+zero-fill correlation is strictly less than 1.0.
    assert report.pearson is not None
    assert 0.0 < report.pearson < 1.0
    # And materially less — disjoint mass should drop pearson well
    # below 0.9 on this fixture (analytical: 4400 / 5300 ≈ 0.83).
    assert report.pearson < 0.9


def test_count_correlation_top_k_secondary_diagnostic() -> None:
    """Top-K Pearson is reported under ``top_k_pearson``, NOT the headline ``pearson``."""
    rows = [
        _row(
            accession="PRJ_Y",
            observed_counts={f"seq{i}": (10 - i) * 10 for i in range(10)},
            reference_counts={f"seq{i}": (10 - i) * 10 + (1 if i == 0 else 0) for i in range(10)},
        )
    ]
    report = compute_count_correlation(rows, top_k=5)
    assert report.pearson is not None
    assert report.top_k_pearson is not None
    # The top_k_pearson field exists and is labeled correctly; the
    # primary pearson is the union-on-everything correlation.
    assert report.top_k == 5


def test_count_correlation_skips_rows_without_counts_cleanly() -> None:
    rows = [
        _row(accession="PRJ_A", observed_counts=None, reference_counts=None),
        _row(accession="PRJ_B", observed_counts={"x": 1}, reference_counts=None),
        _row(accession="PRJ_C", observed_counts={"x": 1}, reference_counts={"x": 1}),
    ]
    report = compute_count_correlation(rows)
    skip_reasons = {r["accession"]: r for r in report.per_row}
    assert skip_reasons["PRJ_A"]["skipped"] == "no_counts"
    assert skip_reasons["PRJ_B"]["skipped"] == "no_reference_counts"
    assert "skipped" not in skip_reasons["PRJ_C"]


# ---------------------------------------------------------------------------
# Aggregator — Codex amendment 1 regression
# ---------------------------------------------------------------------------


def test_aggregate_skips_unverified_rows_with_warning(
    caplog: pytest.LogCaptureFixture, capsys: pytest.CaptureFixture[str]
) -> None:
    """``verified=false`` rows are excluded from numeric metrics AND warned.

    The warning surfaces via both `logger.warning` and stderr so the
    Snakefile + CI logs and the running terminal both see it.
    """
    rows = [
        _row(accession="PRJ_VERIFIED", verified=True),
        _row(accession="PRJ_TBD_1", verified=False),
        _row(accession="PRJ_TBD_2", verified=False),
    ]
    with caplog.at_level(logging.WARNING):
        report = aggregate_metrics(rows)

    assert report.n_verified == 1
    assert report.n_unverified == 2
    assert report.n_total == 3
    assert report.skipped_unverified_accessions == ["PRJ_TBD_1", "PRJ_TBD_2"]
    # Each unverified accession appears in the logged warnings.
    log_text = "\n".join(rec.getMessage() for rec in caplog.records)
    assert "PRJ_TBD_1" in log_text
    assert "PRJ_TBD_2" in log_text
    # And on stderr (so the Snakefile rule prints the skip).
    err = capsys.readouterr().err
    assert "PRJ_TBD_1" in err
    assert "PRJ_TBD_2" in err
    # The single verified row IS counted in primer recovery.
    assert report.primer_recovery.n_evaluated == 1


# ---------------------------------------------------------------------------
# JSON serialization
# ---------------------------------------------------------------------------


def test_write_metrics_json_is_deterministic(tmp_path: Path) -> None:
    rows = [_row(accession="PRJ_A"), _row(accession="PRJ_B")]
    report = aggregate_metrics(rows)
    out1 = tmp_path / "a.json"
    out2 = tmp_path / "b.json"
    write_metrics_json(report, out1)
    write_metrics_json(report, out2)
    assert out1.read_bytes() == out2.read_bytes()
    # And it's sorted-keys at the top level.
    payload = json.loads(out1.read_text())
    keys = list(payload.keys())
    assert keys == sorted(keys)


def test_load_ground_truth_parses_verified_column(tmp_path: Path) -> None:
    tsv = tmp_path / "gt.tsv"
    tsv.write_text(
        "accession\tlibrary_kind\ttarget_kind\tprimer_5p_truth\tprimer_3p_truth"
        "\tn_length_truth\tpaper_doi\tpaper_pmid\tround_map_source\tround_map_path"
        "\tverified\tnotes\n"
        "PRJ_OK\tDNA\tprotein\tACGT\tTGCA\t30\t10.1/abc\t1\tauto\t\ttrue\tok\n"
        "PRJ_TBD\tDNA\tprotein\t\t\t0\t10.1/xyz\t2\tauto\t\tfalse\tnot verified\n",
        encoding="utf-8",
    )
    rows = load_ground_truth(tsv)
    assert {r.accession for r in rows} == {"PRJ_OK", "PRJ_TBD"}
    by_acc = {r.accession: r for r in rows}
    assert by_acc["PRJ_OK"].verified is True
    assert by_acc["PRJ_TBD"].verified is False


# ===========================================================================
# Scope pivot (2026-05-22) — primer-inference-first benchmark
# ===========================================================================


def test_pair_recovery_by_status_classifies_pair_exact() -> None:
    rows = [_row(library_report=_make_lr(status="HIGH"))]
    report = compute_pair_recovery_by_status(rows)
    assert report.n_evaluated == 1
    assert report.counts == {"HIGH": {"pair_exact": 1}}


def test_pair_recovery_by_status_classifies_pair_equivalent() -> None:
    """Both sides match via U-T normalization — pair_equivalent, not pair_exact."""
    rows = [
        _row(
            library_report=_make_lr(
                primer_5p="ACGTACGT",
                primer_3p="TGCATGCA",
                status="HIGH",
            ),
            primer_5p_truth="ACGUACGU",
            primer_3p_truth="UGCAUGCA",
        )
    ]
    report = compute_pair_recovery_by_status(rows)
    assert report.counts == {"HIGH": {"pair_equivalent": 1}}


def test_pair_recovery_by_status_pair_partial_when_one_side_fails() -> None:
    rows = [
        _row(
            library_report=_make_lr(primer_5p="ACGTACGT", primer_3p="GGGGGGGG", status="MEDIUM"),
            primer_5p_truth="ACGTACGT",
            primer_3p_truth="TGCATGCA",
        )
    ]
    report = compute_pair_recovery_by_status(rows)
    assert report.counts == {"MEDIUM": {"pair_partial": 1}}


def test_pair_recovery_by_status_pair_failed_when_no_report() -> None:
    rows = [_row(library_report=None)]
    report = compute_pair_recovery_by_status(rows)
    assert report.counts == {"NO_REPORT": {"pair_failed": 1}}


def test_safe_failure_rate_counts_status_unable_to_infer() -> None:
    rows = [
        _row(accession="PRJ_OK", library_report=_make_lr(status="HIGH")),
        _row(
            accession="PRJ_REFUSE",
            library_report=_make_lr(
                status="UNABLE_TO_INFER",
                primer_5p=None,
                primer_3p=None,
                extraction_mode="UNABLE_TO_EXTRACT",
                required_action="MANUAL_PRIMERS_REQUIRED",
                failure_reason="primers below confidence floor",
            ),
        ),
    ]
    report = compute_safe_failure_rate(rows)
    assert report.n_evaluated == 2
    assert report.n_safe_failures == 1
    assert report.safe_failure_accessions == ["PRJ_REFUSE"]
    assert report.rate == 0.5
    # All three refusal signals triggered for the single PRJ_REFUSE row,
    # so each appears in the by_reason breakdown.
    assert report.by_reason == {
        "extraction_mode_UNABLE_TO_EXTRACT": 1,
        "required_action_MANUAL_PRIMERS_REQUIRED": 1,
        "status_UNABLE_TO_INFER": 1,
    }


def test_safe_failure_rate_treats_missing_library_report_as_safe_failure() -> None:
    """No library_report = the inference run itself failed; counts as a safe failure."""
    rows = [_row(accession="PRJ_NOREPORT", library_report=None)]
    report = compute_safe_failure_rate(rows)
    assert report.n_safe_failures == 1
    assert report.by_reason == {"no_library_report": 1}
    assert report.rate == 1.0


def test_safe_failure_rate_zero_when_all_high_confidence() -> None:
    rows = [_row() for _ in range(3)]
    report = compute_safe_failure_rate(rows)
    assert report.n_safe_failures == 0
    assert report.rate == 0.0
    assert report.by_reason == {}


def test_required_action_distribution_counts() -> None:
    rows = [
        _row(library_report=_make_lr(required_action="NONE")),
        _row(library_report=_make_lr(required_action="NONE")),
        _row(library_report=_make_lr(required_action="READ_MERGING_RECOMMENDED")),
        _row(library_report=None),
    ]
    dist = compute_required_action_distribution(rows)
    assert dist.counts == {
        "NONE": 2,
        "READ_MERGING_RECOMMENDED": 1,
        "NO_REPORT": 1,
    }


def test_aggregate_metrics_does_not_populate_count_correlation() -> None:
    """Scope pivot: count_correlation field stays at its default (empty) values.

    Phase 6c populates it via the self-consistency checker; 6b.1's
    aggregate_metrics deliberately skips it.
    """
    rows = [_row(observed_counts={"x": 5}, reference_counts={"x": 5})]
    report = aggregate_metrics(rows)
    # Default-constructed: pearson and spearman are None, counts are 0.
    assert report.count_correlation.pearson is None
    assert report.count_correlation.spearman is None
    assert report.count_correlation.n_union == 0


def test_aggregate_metrics_populates_pair_recovery_safe_failure_required_action() -> None:
    """End-to-end: the new scope-pivot fields are wired through aggregate_metrics."""
    rows = [
        _row(accession="PRJ_OK", library_report=_make_lr(status="HIGH", required_action="NONE")),
        _row(
            accession="PRJ_REFUSE",
            library_report=_make_lr(
                status="UNABLE_TO_INFER",
                primer_5p=None,
                primer_3p=None,
                extraction_mode="UNABLE_TO_EXTRACT",
                required_action="MANUAL_PRIMERS_REQUIRED",
            ),
        ),
    ]
    report = aggregate_metrics(rows)
    # Pair recovery cross-tab is non-empty.
    assert report.pair_recovery_by_status.n_evaluated == 2
    assert "HIGH" in report.pair_recovery_by_status.counts
    # Safe failure rate is populated.
    assert report.safe_failure_rate.n_safe_failures == 1
    assert report.safe_failure_rate.rate == 0.5
    # Required action distribution is non-empty.
    assert "NONE" in report.required_action_distribution.counts
    assert "MANUAL_PRIMERS_REQUIRED" in report.required_action_distribution.counts


def test_write_metrics_json_includes_new_pivot_fields(tmp_path: Path) -> None:
    rows = [_row(accession="PRJ_A")]
    report = aggregate_metrics(rows)
    out = tmp_path / "metrics.json"
    write_metrics_json(report, out)
    payload = json.loads(out.read_text())
    # Schema keys for the scope pivot are present.
    assert "pair_recovery_by_status" in payload
    assert "safe_failure_rate" in payload
    assert "required_action_distribution" in payload
    # count_correlation remains in the schema (as the Phase 6c entry).
    assert "count_correlation" in payload
    # safe_failure_rate has the documented shape.
    sf = payload["safe_failure_rate"]
    assert set(sf.keys()) == {
        "n_evaluated",
        "n_safe_failures",
        "rate",
        "by_reason",
        "safe_failure_accessions",
    }
