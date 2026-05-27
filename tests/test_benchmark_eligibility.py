"""Unit tests for selexprep.benchmark.eligibility (Phase 6b.5b)."""

from __future__ import annotations

import csv
from pathlib import Path
from unittest.mock import patch

import requests

from selexprep.benchmark.eligibility import (
    AuditEligibility,
    EligibilityReport,
    _normalize_library_name,
    classification_distribution,
    classify_accession,
    classify_plan,
    eligible_accessions,
    main,
    read_eligibility_tsv,
    write_eligibility_tsv,
)
from selexprep.fetch.metadata import RoundRecord
from selexprep.fetch.plan import FetchPlan, FetchRun


def _make_run(
    srr: str,
    library_strategy: str = "OTHER",
    library_name: str = "",
    sample_title: str = "",
    round_number: int | None = None,
) -> FetchRun:
    """Factory: synthesize a FetchRun without going through ENA."""
    confidence = "HIGH" if round_number is not None else "NONE"
    return FetchRun(
        srr=srr,
        sample_accession=f"SAM_{srr}",
        sample_title=sample_title,
        library_name=library_name,
        experiment_title="",
        read_count=1000,
        base_count=75000,
        fastq_urls=[f"ftp.ena.example/{srr}.fastq.gz"],
        fastq_md5s=["deadbeef"],
        fastq_bytes=[1000],
        paired_end=False,
        library_strategy=library_strategy,
        round_record=RoundRecord(
            srr=srr,
            round_number=round_number,
            confidence=confidence,
            source_field="sample_title" if round_number is not None else "none",
            matched_pattern="round_word_digit" if round_number is not None else "none",
            round_candidates=[round_number] if round_number is not None else [],
        ),
    )


def _plan(accession: str, runs: list[FetchRun]) -> FetchPlan:
    return FetchPlan(
        accession=accession,
        bioproject_id=accession,
        study_title="test",
        library_strategy="OTHER",
        library_source="SYNTHETIC",
        runs=runs,
    )


# ----- _normalize_library_name -----


def test_normalize_library_name_digits_replaced() -> None:
    assert _normalize_library_name("RAPT26-1R") == "RAPT{N}-{N}R"
    assert _normalize_library_name("SELEX_round26_N40") == "SELEX_round{N}_N{N}"
    assert _normalize_library_name("DNAFOXR00") == "DNAFOXR{N}"
    # Empty / whitespace → empty (groups collapse to one)
    assert _normalize_library_name("") == ""
    assert _normalize_library_name("   ") == ""


def test_normalize_library_name_collapses_sub_libraries() -> None:
    """PRJNA1246497 in the pilot: 4 runs across 2 sub-libraries x 2 rounds."""
    names = [
        "SELEX_round26_N40",
        "SELEX_round2_N40",
        "SELEX_round26_Stem2",
        "SELEX_round2_Stem2",
    ]
    normalized = {_normalize_library_name(n) for n in names}
    assert normalized == {"SELEX_round{N}_N{N}", "SELEX_round{N}_Stem{N}"}


# ----- classify_plan: all 5 buckets -----


def test_classify_eligible_ht_selex_rounds_clean() -> None:
    """Standard SELEX deposit: ≥2 rounds, single trajectory, all compatible."""
    plan = _plan(
        "PRJ_HTSELEX",
        runs=[
            _make_run("SRR1", library_name="RAPT26-1R", round_number=1),
            _make_run("SRR2", library_name="RAPT26-2R", round_number=2),
            _make_run("SRR3", library_name="RAPT26-3R", round_number=3),
            _make_run("SRR4", library_name="RAPT26-4R", round_number=4),
        ],
    )
    r = classify_plan(plan)
    assert r.classification == AuditEligibility.ELIGIBLE_HT_SELEX_ROUNDS
    assert r.n_runs == 4
    assert r.n_runs_with_round == 4
    assert r.distinct_rounds_parsed == [1, 2, 3, 4]
    assert r.mixed_group_count == 1  # single trajectory
    assert r.n_runs_strategy_compatible == 4
    assert r.n_runs_strategy_blocklisted == 0


def test_classify_non_selex_assay_all_blocklisted() -> None:
    """All runs RNA-Seq → NON_SELEX_ASSAY."""
    plan = _plan(
        "PRJ_RNASEQ",
        runs=[
            _make_run("SRR1", library_strategy="RNA-Seq"),
            _make_run("SRR2", library_strategy="RNA-Seq"),
        ],
    )
    r = classify_plan(plan)
    assert r.classification == AuditEligibility.NON_SELEX_ASSAY
    assert r.n_runs_strategy_blocklisted == 2
    assert r.n_runs_strategy_compatible == 0
    assert "RNA-Seq" in r.reason


def test_classify_no_round_structure_too_few_runs() -> None:
    """Only 1 run with parseable round → NO_ROUND_STRUCTURE."""
    plan = _plan(
        "PRJ_ONE_ROUND",
        runs=[
            _make_run("SRR1", library_name="round1", round_number=1),
            _make_run("SRR2", library_name="unrelated"),  # no round
        ],
    )
    r = classify_plan(plan)
    assert r.classification == AuditEligibility.NO_ROUND_STRUCTURE
    assert r.n_runs_with_round == 1
    assert "≥ 2" in r.reason or "needs" in r.reason


def test_classify_no_round_structure_only_one_distinct_round() -> None:
    """Multiple runs all parsing to the same round → no trajectory."""
    plan = _plan(
        "PRJ_SAME_ROUND",
        runs=[
            _make_run("SRR1", round_number=1),
            _make_run("SRR2", round_number=1),
            _make_run("SRR3", round_number=1),
        ],
    )
    r = classify_plan(plan)
    assert r.classification == AuditEligibility.NO_ROUND_STRUCTURE
    assert r.distinct_rounds_parsed == [1]


def test_classify_mixed_project_via_sub_library_groups() -> None:
    """PRJNA1246497-like case: ≥2 normalized library_name groups → MIXED."""
    plan = _plan(
        "PRJ_MIXED_LIBS",
        runs=[
            _make_run("SRR1", library_name="SELEX_round26_N40", round_number=26),
            _make_run("SRR2", library_name="SELEX_round2_N40", round_number=2),
            _make_run("SRR3", library_name="SELEX_round26_Stem2", round_number=26),
            _make_run("SRR4", library_name="SELEX_round2_Stem2", round_number=2),
        ],
    )
    r = classify_plan(plan)
    assert r.classification == AuditEligibility.MIXED_PROJECT_NEEDS_GROUPING
    assert r.mixed_group_count == 2
    assert r.distinct_rounds_parsed == [2, 26]


def test_classify_mixed_project_via_mixed_strategy() -> None:
    """Mixed library_strategy (SELEX + RNA-Seq control) → MIXED, not NON_SELEX."""
    plan = _plan(
        "PRJ_MIXED_STRATEGY",
        runs=[
            _make_run("SRR1", library_strategy="OTHER", library_name="round1", round_number=1),
            _make_run("SRR2", library_strategy="OTHER", library_name="round2", round_number=2),
            _make_run("SRR3", library_strategy="OTHER", library_name="round3", round_number=3),
            _make_run("SRR4", library_strategy="RNA-Seq", library_name="control"),
        ],
    )
    r = classify_plan(plan)
    assert r.classification == AuditEligibility.MIXED_PROJECT_NEEDS_GROUPING
    assert r.n_runs_strategy_compatible == 3
    assert r.n_runs_strategy_blocklisted == 1
    assert "mixed library_strategy" in r.reason


def test_classify_fetch_dead_empty_plan() -> None:
    """Empty plan → FETCH_DEAD (defensive — build_fetch_plan raises ValueError
    on empty ENA response, but a directly-constructed empty plan is valid
    Python)."""
    plan = _plan("PRJ_EMPTY", runs=[])
    r = classify_plan(plan)
    assert r.classification == AuditEligibility.FETCH_DEAD


# ----- classify_accession: network mocking -----


def test_classify_accession_fetch_failure_yields_fetch_dead() -> None:
    """build_fetch_plan raising HTTPError → FETCH_DEAD with error in reason."""
    with patch(
        "selexprep.benchmark.eligibility.build_fetch_plan",
        side_effect=requests.HTTPError("HTTP 500"),
    ):
        r = classify_accession("PRJ_DEAD")
    assert r.classification == AuditEligibility.FETCH_DEAD
    assert "HTTP 500" in r.reason or "HTTPError" in r.reason


def test_classify_accession_no_records_yields_fetch_dead() -> None:
    """build_fetch_plan raising ValueError ('ENA returned no records') → FETCH_DEAD."""
    with patch(
        "selexprep.benchmark.eligibility.build_fetch_plan",
        side_effect=ValueError("ENA returned no records for accession 'PRJ_X'"),
    ):
        r = classify_accession("PRJ_X")
    assert r.classification == AuditEligibility.FETCH_DEAD
    assert "ENA returned no records" in r.reason


def test_classify_accession_success_delegates_to_classify_plan() -> None:
    """Network call succeeds → classification matches classify_plan output."""
    fake_plan = _plan(
        "PRJ_OK",
        runs=[
            _make_run("SRR1", library_name="R1", round_number=1),
            _make_run("SRR2", library_name="R2", round_number=2),
        ],
    )
    with patch(
        "selexprep.benchmark.eligibility.build_fetch_plan",
        return_value=fake_plan,
    ):
        r = classify_accession("PRJ_OK")
    assert r.classification == AuditEligibility.ELIGIBLE_HT_SELEX_ROUNDS
    assert r.n_runs_with_round == 2


# ----- TSV round-trip + helpers -----


def test_write_and_read_eligibility_tsv_roundtrip(tmp_path: Path) -> None:
    """Inverse property: write → read returns the same EligibilityReports."""
    reports = [
        EligibilityReport(
            accession="PRJ_A",
            classification=AuditEligibility.ELIGIBLE_HT_SELEX_ROUNDS,
            n_runs=4,
            n_runs_with_round=4,
            distinct_rounds_parsed=[0, 1, 2, 3],
            library_strategies=["OTHER"],
            n_runs_strategy_compatible=4,
            n_runs_strategy_blocklisted=0,
            mixed_group_count=1,
            reason="ok",
        ),
        EligibilityReport(
            accession="PRJ_B",
            classification=AuditEligibility.FETCH_DEAD,
            reason="ENA fetch failed: HTTPError",
        ),
    ]
    out = tmp_path / "eligibility.tsv"
    write_eligibility_tsv(reports, out)
    roundtrip = read_eligibility_tsv(out)
    assert len(roundtrip) == 2
    by_acc = {r.accession: r for r in roundtrip}
    assert by_acc["PRJ_A"].classification == AuditEligibility.ELIGIBLE_HT_SELEX_ROUNDS
    assert by_acc["PRJ_A"].distinct_rounds_parsed == [0, 1, 2, 3]
    assert by_acc["PRJ_A"].library_strategies == ["OTHER"]
    assert by_acc["PRJ_B"].classification == AuditEligibility.FETCH_DEAD


def test_write_eligibility_tsv_sorted_by_accession(tmp_path: Path) -> None:
    """TSV rows are sorted by accession for deterministic diffs."""
    reports = [
        EligibilityReport(accession="PRJ_Z", classification=AuditEligibility.FETCH_DEAD),
        EligibilityReport(accession="PRJ_A", classification=AuditEligibility.FETCH_DEAD),
        EligibilityReport(accession="PRJ_M", classification=AuditEligibility.FETCH_DEAD),
    ]
    out = tmp_path / "eligibility.tsv"
    write_eligibility_tsv(reports, out)
    with out.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert [r["accession"] for r in rows] == ["PRJ_A", "PRJ_M", "PRJ_Z"]


def test_eligible_accessions_returns_only_eligible_sorted() -> None:
    reports = [
        EligibilityReport(
            accession="PRJ_C", classification=AuditEligibility.ELIGIBLE_HT_SELEX_ROUNDS
        ),
        EligibilityReport(
            accession="PRJ_A", classification=AuditEligibility.ELIGIBLE_HT_SELEX_ROUNDS
        ),
        EligibilityReport(accession="PRJ_B", classification=AuditEligibility.NO_ROUND_STRUCTURE),
        EligibilityReport(accession="PRJ_D", classification=AuditEligibility.NON_SELEX_ASSAY),
        EligibilityReport(
            accession="PRJ_E", classification=AuditEligibility.MIXED_PROJECT_NEEDS_GROUPING
        ),
        EligibilityReport(accession="PRJ_F", classification=AuditEligibility.FETCH_DEAD),
    ]
    assert eligible_accessions(reports) == ["PRJ_A", "PRJ_C"]


def test_classification_distribution_sorted_alpha() -> None:
    reports = [
        EligibilityReport(accession="A", classification=AuditEligibility.NO_ROUND_STRUCTURE),
        EligibilityReport(accession="B", classification=AuditEligibility.NO_ROUND_STRUCTURE),
        EligibilityReport(accession="C", classification=AuditEligibility.ELIGIBLE_HT_SELEX_ROUNDS),
        EligibilityReport(accession="D", classification=AuditEligibility.FETCH_DEAD),
    ]
    d = classification_distribution(reports)
    # Sorted alphabetically (stable JSON shape).
    assert list(d.keys()) == ["ELIGIBLE_HT_SELEX_ROUNDS", "FETCH_DEAD", "NO_ROUND_STRUCTURE"]
    assert d["NO_ROUND_STRUCTURE"] == 2


# ----- CLI: classify-catalog -----


def test_cli_classify_catalog_writes_tsv_and_distribution(tmp_path: Path) -> None:
    """End-to-end CLI smoke: classify a tiny synthetic catalog with mocked
    ENA, verify the TSV is written + the stderr distribution is printed."""
    catalog = tmp_path / "bioprojects.csv"
    catalog.write_text(
        "bioproject_id,source,study_title\n"
        "PRJNA_ELIG,ena,Real SELEX\n"
        "PRJNA_RNASEQ,ena,RNA-seq false positive\n"
        "zenodo:99,zenodo,Non-INSDC (filtered out by classifier)\n",
        encoding="utf-8",
    )
    elig_plan = _plan(
        "PRJNA_ELIG",
        runs=[
            _make_run("SRR1", library_name="R1", round_number=1),
            _make_run("SRR2", library_name="R2", round_number=2),
        ],
    )
    rnaseq_plan = _plan(
        "PRJNA_RNASEQ",
        runs=[_make_run("SRR3", library_strategy="RNA-Seq")],
    )
    plans_by_acc = {"PRJNA_ELIG": elig_plan, "PRJNA_RNASEQ": rnaseq_plan}

    def fake_fetch(acc: str, *, timeout_s: int = 30) -> FetchPlan:
        return plans_by_acc[acc]

    out = tmp_path / "eligibility.tsv"
    with patch("selexprep.benchmark.eligibility.build_fetch_plan", side_effect=fake_fetch):
        rc = main(["classify-catalog", "--catalog", str(catalog), "--out", str(out)])
    assert rc == 0
    reports = read_eligibility_tsv(out)
    # Only the 2 INSDC accessions get classified; zenodo:99 is skipped.
    by_acc = {r.accession: r for r in reports}
    assert set(by_acc) == {"PRJNA_ELIG", "PRJNA_RNASEQ"}
    assert by_acc["PRJNA_ELIG"].classification == AuditEligibility.ELIGIBLE_HT_SELEX_ROUNDS
    assert by_acc["PRJNA_RNASEQ"].classification == AuditEligibility.NON_SELEX_ASSAY


# ----- Empirical pilot regression tests -----
#
# The Phase 6b.4 audit pilot's 30 accessions, post-Phase 6b.5a refresh
# of the catalog. These cases pin classify_plan against the empirical
# data Codex surfaced in the per-row fetch_metadata.json inspection.


def test_pilot_rapt26_is_eligible() -> None:
    """PRJDB19138 (RAPT26-1R..5R, all MEDIUM single-match) → ELIGIBLE."""
    plan = _plan(
        "PRJDB19138",
        runs=[
            _make_run("DRR618526", library_name="RAPT26-1R", round_number=1),
            _make_run("DRR618527", library_name="RAPT26-2R", round_number=2),
            _make_run("DRR618528", library_name="RAPT26-3R", round_number=3),
            _make_run("DRR618529", library_name="RAPT26-4R", round_number=4),
            _make_run("DRR618530", library_name="RAPT26-5R", round_number=5),
        ],
    )
    r = classify_plan(plan)
    assert r.classification == AuditEligibility.ELIGIBLE_HT_SELEX_ROUNDS
    assert r.mixed_group_count == 1


def test_pilot_pyoverdine_single_run_is_no_round_structure() -> None:
    """PRJNA932049 (1 run only, no round info) → NO_ROUND_STRUCTURE."""
    plan = _plan(
        "PRJNA932049",
        runs=[
            _make_run(
                "SRR23354534",
                library_name="GSM7029482",
                sample_title="Aptamers, 2' FY-RNA, pyoverdine",
            )
        ],
    )
    r = classify_plan(plan)
    assert r.classification == AuditEligibility.NO_ROUND_STRUCTURE
    assert r.n_runs_with_round == 0
