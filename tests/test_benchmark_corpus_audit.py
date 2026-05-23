"""Tests for the Tier 2 corpus-audit module (Phase 6b.3a)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from selexprep.benchmark.corpus_audit import (
    CorpusAuditReport,
    accessions_sha256,
    aggregate_audit_from_run_outputs,
    main,
    sample_corpus,
    write_accessions_tsv,
    write_audit_json,
)

# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------


def test_sample_corpus_deterministic_for_same_seed() -> None:
    """Same (n, seed) on the bundled catalog → identical accession list."""
    a = sample_corpus(n=5, seed=42)
    b = sample_corpus(n=5, seed=42)
    assert a == b
    assert len(a) <= 5
    assert all(isinstance(x, str) and x for x in a)


def test_sample_corpus_different_seed_different_sample() -> None:
    """Different seed → different draw on a non-trivial catalog."""
    a = sample_corpus(n=10, seed=1)
    b = sample_corpus(n=10, seed=2)
    # On a 200+-row INSDC pool a 10-element sample with two distinct
    # seeds should almost certainly differ; assert weakly to keep this
    # robust against future catalog size changes.
    assert a != b


def test_sample_corpus_excludes_ground_truth() -> None:
    """Excluded accessions never appear in the sample."""
    candidates = sample_corpus(n=20, seed=42)
    assert candidates  # sanity — there are INSDC rows
    exclude = (candidates[0],)
    re_sampled = sample_corpus(n=20, seed=42, exclude=exclude)
    assert candidates[0] not in re_sampled


def test_sample_corpus_insdc_only() -> None:
    """INSDC filter drops figshare:* / zenodo:* / utexas:* prefixes."""
    sample = sample_corpus(n=50, seed=42)
    for accession in sample:
        assert not accession.startswith(("figshare:", "zenodo:", "utexas:")), accession


def test_sample_corpus_source_substring_filter() -> None:
    """The ``sources`` substring filter narrows the pool to that adapter."""
    # ena substring should keep at least one row but never more than the
    # unfiltered total. We assert the pool shrinks (or stays equal).
    full = sample_corpus(n=10, seed=42)
    ena_only = sample_corpus(n=10, seed=42, sources="ena")
    assert len(ena_only) <= len(full)


def test_sample_corpus_empty_pool_returns_empty() -> None:
    """No matching INSDC accessions → empty list, never an error."""
    sample = sample_corpus(n=5, seed=42, sources="zzzz_no_such_source")
    assert sample == []


# ---------------------------------------------------------------------------
# accessions_sha256 + write_accessions_tsv
# ---------------------------------------------------------------------------


def test_accessions_sha256_independent_of_order() -> None:
    """The hash is over sorted order, so input order does not matter."""
    a = ["PRJNA1", "PRJNA2", "PRJNA3"]
    b = ["PRJNA3", "PRJNA1", "PRJNA2"]
    assert accessions_sha256(a) == accessions_sha256(b)


def test_write_accessions_tsv_emits_two_columns_sorted(tmp_path: Path) -> None:
    """The TSV format matches what ``selexprep run`` consumes (locked plan)."""
    out = tmp_path / "accessions.tsv"
    write_accessions_tsv(["PRJNA42", "PRJNA1"], out)
    text = out.read_text(encoding="utf-8")
    lines = text.splitlines()
    assert lines[0] == "accession\tnotes"
    assert lines[1].startswith("PRJNA1\t")
    assert lines[2].startswith("PRJNA42\t")


# ---------------------------------------------------------------------------
# Aggregator
# ---------------------------------------------------------------------------


def _write_run_summary(path: Path, rows: list[dict[str, str]]) -> None:
    """Helper: emit a synthetic run_summary.tsv with the canonical column set."""
    columns = (
        "accession",
        "status",
        "last_stage_completed",
        "library_report_status",
        "extraction_mode",
        "required_action",
        "confidence",
        "flags_raised",
        "notes",
    )
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(r.get(c, "") for c in columns) + "\n")


def test_aggregator_inference_safe_failure_rate_excludes_fetch_failures(
    tmp_path: Path,
) -> None:
    """Methodological correction (Codex + user): denominator is ONLY rows with a LibraryReport.

    Three sampled accessions; one fetch-fails outright. The safe-failure
    rate is computed over the remaining two rows that produced a
    LibraryReport — NOT 1/3 (which would inflate the metric).
    """
    summary = tmp_path / "run_summary.tsv"
    _write_run_summary(
        summary,
        [
            {
                "accession": "PRJA_OK",
                "status": "OK",
                "library_report_status": "HIGH",
                "extraction_mode": "BOTH_PRIMERS_SINGLE_READ",
                "required_action": "NONE",
                "flags_raised": "0",
            },
            {
                "accession": "PRJA_REFUSE",
                "status": "EXTRACT_REFUSED",
                "library_report_status": "UNABLE_TO_INFER",
                "extraction_mode": "UNABLE_TO_EXTRACT",
                "required_action": "MANUAL_PRIMERS_REQUIRED",
            },
            {
                "accession": "PRJA_FETCH_FAIL",
                "status": "FETCH_FAILED",
                # No library_report fields populated.
            },
        ],
    )
    report = aggregate_audit_from_run_outputs(
        run_summary_tsv=summary,
        ground_truth_tsv=None,
        catalog_version="test-v1",
        sample_seed=42,
        sample_accessions_sha256="dummy",
    )
    assert report.n_sampled == 3
    assert report.n_with_library_report == 2
    assert report.n_inference_safe_failures == 1
    # 1 / 2, not 1 / 3.
    assert report.inference_safe_failure_rate == pytest.approx(0.5)
    assert report.fetch_outcome_distribution["FETCH_FAILED"] == 1
    assert report.fetch_outcome_distribution["OK"] == 1
    assert report.fetch_outcome_distribution["EXTRACT_REFUSED"] == 1


def test_aggregator_counts_fetchable_correctly(tmp_path: Path) -> None:
    """Anything past the fetch stage counts as fetchable (locked plan partition)."""
    summary = tmp_path / "run_summary.tsv"
    _write_run_summary(
        summary,
        [
            {
                "accession": "A",
                "status": "OK",
                "library_report_status": "HIGH",
                "extraction_mode": "BOTH_PRIMERS_SINGLE_READ",
                "required_action": "NONE",
                "flags_raised": "0",
            },
            {"accession": "B", "status": "DETECT_FAILED"},
            {"accession": "C", "status": "FETCH_REFUSED"},
            {"accession": "D", "status": "FETCH_FAILED"},
            {
                "accession": "E",
                "status": "EXTRACT_FAILED",
                "library_report_status": "HIGH",
                "extraction_mode": "BOTH_PRIMERS_SINGLE_READ",
                "required_action": "NONE",
            },
        ],
    )
    report = aggregate_audit_from_run_outputs(
        run_summary_tsv=summary,
        ground_truth_tsv=None,
        catalog_version=None,
        sample_seed=0,
        sample_accessions_sha256="",
    )
    # OK + DETECT_FAILED + EXTRACT_FAILED = 3 fetchable; FETCH_REFUSED +
    # FETCH_FAILED = 2 not fetchable.
    assert report.n_fetchable == 3
    assert report.fetch_outcome_distribution["DETECT_FAILED"] == 1
    assert report.fetch_outcome_distribution["FETCH_REFUSED"] == 1
    assert report.fetch_outcome_distribution["FETCH_FAILED"] == 1


def test_aggregator_qc_flags_histogram(tmp_path: Path) -> None:
    """Only OK rows with a flags_raised value contribute to the histogram."""
    summary = tmp_path / "run_summary.tsv"
    _write_run_summary(
        summary,
        [
            {
                "accession": "A",
                "status": "OK",
                "library_report_status": "HIGH",
                "extraction_mode": "BOTH_PRIMERS_SINGLE_READ",
                "required_action": "NONE",
                "flags_raised": "2",
            },
            {
                "accession": "B",
                "status": "OK",
                "library_report_status": "HIGH",
                "extraction_mode": "BOTH_PRIMERS_SINGLE_READ",
                "required_action": "NONE",
                "flags_raised": "0",
            },
            {
                "accession": "C",
                "status": "OK",
                "library_report_status": "HIGH",
                "extraction_mode": "BOTH_PRIMERS_SINGLE_READ",
                "required_action": "NONE",
                "flags_raised": "2",
            },
            # QC didn't run — should NOT enter the histogram.
            {"accession": "D", "status": "EXTRACT_FAILED"},
        ],
    )
    report = aggregate_audit_from_run_outputs(
        run_summary_tsv=summary,
        ground_truth_tsv=None,
        catalog_version=None,
        sample_seed=0,
        sample_accessions_sha256="",
    )
    assert report.n_with_qc_run == 3
    assert report.flags_raised_histogram[0] == 1
    assert report.flags_raised_histogram[2] == 2


def test_aggregator_ground_truth_overlap(tmp_path: Path) -> None:
    """``is_in_ground_truth`` annotation + overlap count."""
    summary = tmp_path / "run_summary.tsv"
    _write_run_summary(
        summary,
        [
            {
                "accession": "PRJ_IN_GT",
                "status": "OK",
                "library_report_status": "HIGH",
                "extraction_mode": "BOTH_PRIMERS_SINGLE_READ",
                "required_action": "NONE",
                "flags_raised": "0",
            },
            {"accession": "PRJ_NOT_IN_GT", "status": "FETCH_FAILED"},
        ],
    )
    gt = tmp_path / "ground_truth.tsv"
    gt.write_text(
        "accession\tlibrary_kind\nPRJ_IN_GT\tDNA\nPRJ_OTHER\tDNA\n",
        encoding="utf-8",
    )
    report = aggregate_audit_from_run_outputs(
        run_summary_tsv=summary,
        ground_truth_tsv=gt,
        catalog_version="test",
        sample_seed=42,
        sample_accessions_sha256="x",
    )
    assert report.n_in_ground_truth_overlap == 1
    accessions_in_gt = {r["accession"] for r in report.per_accession if r["is_in_ground_truth"]}
    assert accessions_in_gt == {"PRJ_IN_GT"}


# ---------------------------------------------------------------------------
# JSON writer determinism
# ---------------------------------------------------------------------------


def test_write_audit_json_deterministic(tmp_path: Path) -> None:
    """sorted keys + stable per_accession order → bit-identical output across runs.

    Also covers the Codex peer-review defensive sort: the writer re-sorts
    ``per_accession`` even if the caller already sorted upstream, so direct
    callers (tests, hand-rolled scripts) can't accidentally write
    non-deterministic ordering. We pass intentionally unsorted input here.
    """
    report = CorpusAuditReport(
        catalog_version="x",
        sample_seed=7,
        sample_accessions_sha256="abc",
        n_sampled=2,
        fetch_outcome_distribution={"OK": 1, "FETCH_FAILED": 1},
        n_fetchable=1,
        n_with_library_report=1,
        library_report_status_distribution={"HIGH": 1},
        extraction_mode_distribution={"BOTH_PRIMERS_SINGLE_READ": 1},
        required_action_distribution={"NONE": 1},
        inference_safe_failure_rate=0.0,
        n_inference_safe_failures=0,
        flags_raised_histogram={0: 1},
        n_with_qc_run=1,
        per_accession=[
            {"accession": "B", "status": "OK", "is_in_ground_truth": False},
            {"accession": "A", "status": "FETCH_FAILED", "is_in_ground_truth": False},
        ],
    )
    out1 = tmp_path / "audit1.json"
    out2 = tmp_path / "audit2.json"
    write_audit_json(report, out1)
    write_audit_json(report, out2)
    assert out1.read_bytes() == out2.read_bytes()
    parsed = json.loads(out1.read_text(encoding="utf-8"))
    assert list(parsed.keys()) == sorted(parsed.keys())
    # flags_raised_histogram keys must be stringified ints.
    assert "0" in parsed["flags_raised_histogram"]
    # Codex peer-review defensive sort: per_accession ordering is enforced
    # by write_audit_json itself, NOT trusted from the caller. The input
    # report's per_accession was B, A (unsorted); the JSON must show A, B.
    assert [r["accession"] for r in parsed["per_accession"]] == ["A", "B"]


# ---------------------------------------------------------------------------
# CLI: sample + aggregate
# ---------------------------------------------------------------------------


def test_main_sample_writes_tsv_and_manifest(tmp_path: Path) -> None:
    """``corpus_audit sample`` writes accessions TSV + sidecar manifest with envelope."""
    out = tmp_path / "accs.tsv"
    rc = main(
        [
            "sample",
            "--n",
            "5",
            "--seed",
            "42",
            "--out",
            str(out),
        ]
    )
    assert rc == 0
    assert out.exists()
    sidecar = out.with_suffix(".manifest.json")
    assert sidecar.exists()
    manifest = json.loads(sidecar.read_text(encoding="utf-8"))
    assert manifest["sample_seed"] == 42
    assert manifest["catalog_version"]
    assert len(manifest["sample_accessions_sha256"]) == 64


def test_main_aggregate_picks_up_sidecar(tmp_path: Path) -> None:
    """The aggregator reads catalog_version + seed + sha from the sidecar manifest."""
    summary = tmp_path / "run_summary.tsv"
    _write_run_summary(
        summary,
        [
            {
                "accession": "A",
                "status": "OK",
                "library_report_status": "HIGH",
                "extraction_mode": "BOTH_PRIMERS_SINGLE_READ",
                "required_action": "NONE",
                "flags_raised": "0",
            },
        ],
    )
    sidecar = tmp_path / "audit_accessions.manifest.json"
    sidecar.write_text(
        json.dumps(
            {
                "catalog_version": "v0.1.5-test",
                "sample_seed": 99,
                "sample_accessions_sha256": "deadbeef" * 8,
                "n_sampled": 1,
            }
        ),
        encoding="utf-8",
    )
    out = tmp_path / "audit_metrics.json"
    rc = main(
        [
            "aggregate",
            "--run-summary",
            str(summary),
            "--sample-manifest",
            str(sidecar),
            "--out",
            str(out),
        ]
    )
    assert rc == 0
    parsed = json.loads(out.read_text(encoding="utf-8"))
    assert parsed["catalog_version"] == "v0.1.5-test"
    assert parsed["sample_seed"] == 99
    assert parsed["sample_accessions_sha256"] == "deadbeef" * 8


def test_main_aggregate_refuses_when_no_sha_available(tmp_path: Path) -> None:
    """Codex peer-review fix: empty reproducibility envelope is a hard refusal.

    Neither a sidecar manifest nor an explicit ``--sample-sha`` is provided.
    The CLI must refuse to emit an audit JSON with empty provenance — the
    whole point of the envelope is that the audit can be reproduced
    across catalog refreshes by fingerprinting the sampled accession list.
    """
    summary = tmp_path / "run_summary.tsv"
    _write_run_summary(
        summary,
        [
            {
                "accession": "A",
                "status": "OK",
                "library_report_status": "HIGH",
                "extraction_mode": "BOTH_PRIMERS_SINGLE_READ",
                "required_action": "NONE",
                "flags_raised": "0",
            },
        ],
    )
    out = tmp_path / "audit_metrics.json"
    with pytest.raises(SystemExit) as exc_info:
        main(
            [
                "aggregate",
                "--run-summary",
                str(summary),
                "--out",
                str(out),
            ]
        )
    # SystemExit message carries the explanation; the audit JSON must NOT
    # have been written.
    assert "sample sha256" in str(exc_info.value)
    assert not out.exists()


def test_main_aggregate_warns_on_n_sampled_mismatch(
    tmp_path: Path, capsys: pytest.CaptureFixture
) -> None:
    """Codex peer-review fix: sidecar n_sampled != run summary rows → loud warning.

    The audit JSON's ``n_sampled`` field semantically means "what was
    actually processed" — len(run_summary.tsv). If the sidecar says
    something else, that's either a selexprep run bug (rows lost) or
    operator error (re-ran with a different TSV). Surface it; don't fail
    silently.
    """
    summary = tmp_path / "run_summary.tsv"
    _write_run_summary(
        summary,
        [
            {
                "accession": "A",
                "status": "OK",
                "library_report_status": "HIGH",
                "extraction_mode": "BOTH_PRIMERS_SINGLE_READ",
                "required_action": "NONE",
                "flags_raised": "0",
            },
        ],
    )
    sidecar = tmp_path / "audit_accessions.manifest.json"
    sidecar.write_text(
        json.dumps(
            {
                "catalog_version": "v0.1.5-test",
                "sample_seed": 42,
                "sample_accessions_sha256": "deadbeef" * 8,
                "n_sampled": 30,  # sidecar says 30, summary has 1 — mismatch
            }
        ),
        encoding="utf-8",
    )
    out = tmp_path / "audit_metrics.json"
    rc = main(
        [
            "aggregate",
            "--run-summary",
            str(summary),
            "--sample-manifest",
            str(sidecar),
            "--out",
            str(out),
        ]
    )
    assert rc == 0
    # WARNING printed to stdout, audit JSON still emitted with the
    # actually-processed row count.
    captured = capsys.readouterr()
    assert "WARNING" in captured.out
    assert "30" in captured.out  # sidecar's n_sampled
    parsed = json.loads(out.read_text(encoding="utf-8"))
    assert parsed["n_sampled"] == 1  # summary's row count is the source of truth
