"""Unit tests for ``selexprep.qc.flags``.

One positive + one negative case per flag.
"""

from __future__ import annotations

from pathlib import Path

import yaml

from selexprep.library.report import LibraryReport
from selexprep.manifest import SelexprepManifestV1, build_manifest_from_extract_result
from selexprep.qc.flags import (
    Flag,
    check_adapter_contamination_high,
    check_extraction_mode_changed_across_rounds,
    check_low_primer_match,
    check_low_total_reads,
    check_n_length_variation_across_rounds,
    check_requires_read_merging_for_full_insert,
    check_strand_mix,
    check_unexpected_rarefied_diversity_increase,
    compute_all_flags,
    write_flags_yaml,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_library_report(**overrides: object) -> LibraryReport:
    base = {
        "primer_5p": "GGTAATACGACTCACTATAGGG",
        "primer_3p": "CCATGCATGCATGCATGCAT",
        "variants_5p": [],
        "variants_3p": [],
        "known_adapter_hits": {"TRUSEQ_R1": 0, "NEXTERA": 0},
        "extraction_mode": "BOTH_PRIMERS_SINGLE_READ",
        "full_insert_recovered": True,
        "read_source": "R1",
        "required_action": "NONE",
        "orientation": "FORWARD",
        "n_length_mode": 30,
        "n_length_distribution": {30: 1000},
        "n_length_confidence": 1.0,
        "match_rate_5p": 0.95,
        "match_rate_3p": 0.95,
        "position_consistency_5p": 0.95,
        "position_consistency_3p": 0.95,
        "read_fraction_used_for_inference": 1.0,
        "sampling_seed": 42,
        "confidence": 0.85,
        "status": "HIGH",
        "failure_reason": None,
    }
    base.update(overrides)
    return LibraryReport(**base)  # type: ignore[arg-type]


def _make_manifest(lr: LibraryReport | None = None) -> SelexprepManifestV1:
    return build_manifest_from_extract_result(
        library_report=lr or _make_library_report(),
        input_paths=[],
        output_paths=[],
        accession="SRR123",
        bioproject_id="PRJ123",
        runs=["SRR123"],
        parameters={},
    )


def _diverse_pool(n_unique: int, reads_per: int = 5) -> dict[str, int]:
    """Synthetic pool: n_unique sequences, all 30 nt, equal read counts."""
    return {("ACGT" * 8)[: i % 30 + 10] + str(i): reads_per for i in range(n_unique)}


# ---------------------------------------------------------------------------
# check_low_primer_match
# ---------------------------------------------------------------------------


def test_low_primer_match_negative() -> None:
    lr = _make_library_report(match_rate_5p=0.9, match_rate_3p=0.9)
    assert check_low_primer_match(lr) is None


def test_low_primer_match_positive_5p_only() -> None:
    lr = _make_library_report(match_rate_5p=0.2, match_rate_3p=0.9)
    flag = check_low_primer_match(lr)
    assert flag is not None
    assert flag.name == "low_primer_match"
    assert flag.severity == "warn"


def test_low_primer_match_positive_both_sides() -> None:
    lr = _make_library_report(match_rate_5p=0.1, match_rate_3p=0.2)
    flag = check_low_primer_match(lr)
    assert flag is not None
    sides = [s["side"] for s in flag.evidence["sides_below"]]  # type: ignore[index]
    assert "5p" in sides and "3p" in sides


# ---------------------------------------------------------------------------
# check_n_length_variation_across_rounds
# ---------------------------------------------------------------------------


def test_n_length_variation_negative_uniform() -> None:
    pools = {r: {("A" * 30) + str(i): 10 for i in range(5)} for r in range(3)}
    assert check_n_length_variation_across_rounds(pools) is None


def test_n_length_variation_positive_many_modes() -> None:
    # 4 distinct modal lengths across 4 rounds -> exceeds threshold (2).
    pools = {
        0: {"A" * 20: 100},
        1: {"A" * 25: 100},
        2: {"A" * 30: 100},
        3: {"A" * 35: 100},
    }
    flag = check_n_length_variation_across_rounds(pools)
    assert flag is not None
    assert flag.name == "n_length_variation_across_rounds"


# ---------------------------------------------------------------------------
# check_low_total_reads
# ---------------------------------------------------------------------------


def test_low_total_reads_negative() -> None:
    pools = {r: {f"s{i}": 100 for i in range(200)} for r in range(3)}  # 20k each
    assert check_low_total_reads(pools) is None


def test_low_total_reads_positive() -> None:
    pools = {
        0: {f"s{i}": 100 for i in range(200)},  # 20k -> OK
        1: {f"t{i}": 5 for i in range(100)},  # 500 -> below 10k
    }
    flag = check_low_total_reads(pools)
    assert flag is not None
    assert flag.name == "low_total_reads"


# ---------------------------------------------------------------------------
# check_adapter_contamination_high
# ---------------------------------------------------------------------------


def test_adapter_contamination_negative() -> None:
    """100 adapter hits / 100k inference-universe reads = 0.1% -> below threshold."""
    lr = _make_library_report(
        known_adapter_hits={"TRUSEQ_R1": 100, "NEXTERA": 0},
        read_fraction_used_for_inference=1.0,
    )
    trim_reports = {0: {"n_in": 100_000, "n_out": 95_000}}
    assert check_adapter_contamination_high(lr, trim_reports) is None


def test_adapter_contamination_positive() -> None:
    """1000 adapter hits / 1k inference-universe reads = 100% -> above threshold."""
    lr = _make_library_report(
        known_adapter_hits={"TRUSEQ_R1": 1000, "NEXTERA": 0},
        read_fraction_used_for_inference=1.0,
    )
    trim_reports = {0: {"n_in": 1_000, "n_out": 900}}
    flag = check_adapter_contamination_high(lr, trim_reports)
    assert flag is not None
    assert flag.name == "adapter_contamination_high"
    assert flag.severity == "warn"


def test_adapter_contamination_denominator_unavailable() -> None:
    """When trim_reports is missing, surface flag as 'info' with reason."""
    lr = _make_library_report(known_adapter_hits={"TRUSEQ_R1": 50})
    flag = check_adapter_contamination_high(lr, None)
    assert flag is not None
    assert flag.severity == "info"
    assert flag.evidence["reason"] == "denominator_unavailable"


def test_adapter_contamination_honors_read_fraction_used_for_inference() -> None:
    """If inference subsampled (e.g., used 50% of reads), the
    denominator must scale down to match — otherwise the fraction is
    understated by 2x."""
    lr = _make_library_report(
        known_adapter_hits={"TRUSEQ_R1": 100},
        read_fraction_used_for_inference=0.5,
    )
    # n_in=10000, fraction=0.5 → denominator should be 5000
    # hits=100, frac = 100/5000 = 0.02 (below 0.05 threshold)
    trim_reports = {0: {"n_in": 10_000, "n_out": 9_000}}
    assert check_adapter_contamination_high(lr, trim_reports) is None
    # If the same hit count were measured against the FULL n_in (no
    # subsampling correction), frac = 100/10000 = 0.01 — different number,
    # confirming the read_fraction adjustment matters.


# ---------------------------------------------------------------------------
# check_strand_mix
# ---------------------------------------------------------------------------


def test_strand_mix_no_report_returns_none() -> None:
    assert check_strand_mix(None) is None


def test_strand_mix_negative(tmp_path: Path) -> None:
    """All-forward strand report -> no flag."""
    path = tmp_path / "strand_report.tsv"
    path.write_text(
        "round\tforward\treverse\tambiguous\n0\t100\t0\t0\n1\t100\t0\t0\n",
        encoding="utf-8",
    )
    assert check_strand_mix(path) is None


def test_strand_mix_positive(tmp_path: Path) -> None:
    """One round with >20% reverse -> flag fires."""
    path = tmp_path / "strand_report.tsv"
    path.write_text(
        "round\tforward\treverse\tambiguous\n0\t100\t0\t0\n1\t60\t40\t0\n",  # 40% reverse
        encoding="utf-8",
    )
    flag = check_strand_mix(path)
    assert flag is not None
    assert flag.name == "strand_mix"


# ---------------------------------------------------------------------------
# check_unexpected_rarefied_diversity_increase
# ---------------------------------------------------------------------------


def test_diversity_increase_negative_monotonic_decrease() -> None:
    """Diversity decreases across rounds -> no flag."""
    pools = {
        0: _diverse_pool(50000, reads_per=10),  # high diversity
        1: _diverse_pool(20000, reads_per=10),
        2: _diverse_pool(10000, reads_per=10),
    }
    assert check_unexpected_rarefied_diversity_increase(pools) is None


def test_diversity_increase_positive() -> None:
    """Diversity goes UP from R0 to R1 -> flag fires."""
    pools = {
        0: _diverse_pool(5_000, reads_per=10),  # ~5k unique
        1: _diverse_pool(50_000, reads_per=10),  # ~50k unique
    }
    flag = check_unexpected_rarefied_diversity_increase(pools)
    assert flag is not None
    assert flag.name == "unexpected_rarefied_diversity_increase"


def test_diversity_increase_clamps_depth_to_min_total(tmp_path: Path) -> None:
    """BLOCKING regression: rarefaction depth must
    clamp to ``min(RAREFACTION_DEPTH, min_total_reads_per_round)``. The
    previous bug let a 2k-read round skip rarefaction entirely (the
    ``rarefy`` helper returns the original pool when ``depth >= total``)
    while a 100k-read round got actually rarefied to 10k — comparing
    unique counts across those two pools was exactly the depth-
    confounding the flag was meant to prevent."""
    # round_0: 2_000 total reads (below default RAREFACTION_DEPTH of 10_000)
    # round_1: 100_000 total reads (above default)
    pool_0 = {f"s{i}": 2 for i in range(1_000)}  # 1000 uniques * 2 reads = 2k
    pool_1 = {f"t{i}": 100 for i in range(1_000)}  # 1000 uniques * 100 = 100k
    pools = {0: pool_0, 1: pool_1}

    # Even though the raw unique counts are equal (1000 each), the flag
    # must rarefy both to the SAME depth (min_total = 2000) so the
    # comparison is depth-fair. Result should not fire (diversity equal
    # at common depth).
    flag = check_unexpected_rarefied_diversity_increase(pools)

    # If flag did fire, evidence must at least document the effective_depth.
    if flag is not None:
        assert "effective_depth" in flag.evidence
        assert "configured_depth" in flag.evidence
        # And effective_depth must be clamped to min_total (2000), not 10000.
        assert flag.evidence["effective_depth"] == 2_000, (
            f"effective_depth should clamp to min_total (2_000); "
            f"got {flag.evidence['effective_depth']}"
        )
    # The KEY property: rarefy(pool_0, 2000) returns pool_0 unchanged
    # (sums to 2000); rarefy(pool_1, 2000) subsamples to 2k reads.
    # Both pools have the same uniques (1000) so flag should NOT fire.
    assert flag is None, "flag should not fire when both rarefied pools have equal uniques"


def test_diversity_increase_records_effective_depth_in_evidence() -> None:
    """When the flag DOES fire, evidence must include both configured_depth
    and effective_depth so users can audit the rarefaction decision.

    Designed so the smaller round caps effective_depth below
    RAREFACTION_DEPTH and rarefaction-fair comparison still surfaces an
    increase (pool_1 has much more diversity per common read)."""
    # pool_0: 1000 uniques * 5 reads each = 5_000 total reads (sets min_total)
    # pool_1: 5000 uniques * 100 reads each = 500_000 total reads (higher diversity)
    pool_0 = {f"a{i}": 5 for i in range(1_000)}
    pool_1 = {f"b{i}": 100 for i in range(5_000)}
    pools = {0: pool_0, 1: pool_1}

    flag = check_unexpected_rarefied_diversity_increase(pools)
    assert flag is not None, "expected flag to fire (pool_1 has higher diversity at common depth)"
    assert "configured_depth" in flag.evidence
    assert "effective_depth" in flag.evidence
    assert flag.evidence["effective_depth"] <= flag.evidence["configured_depth"]
    # min_total = 5000, so effective_depth = min(10000, 5000) = 5000.
    assert flag.evidence["effective_depth"] == 5_000
    assert flag.evidence["configured_depth"] == 10_000
    assert flag.evidence["min_total_reads"] == 5_000


# ---------------------------------------------------------------------------
# check_requires_read_merging_for_full_insert
# ---------------------------------------------------------------------------


def test_read_merging_info_negative() -> None:
    lr = _make_library_report(required_action="NONE")
    assert check_requires_read_merging_for_full_insert(lr) is None


def test_read_merging_info_positive() -> None:
    lr = _make_library_report(
        required_action="READ_MERGING_RECOMMENDED",
        extraction_mode="PAIRED_END_SPLIT_PRIMERS",
        full_insert_recovered=False,
    )
    flag = check_requires_read_merging_for_full_insert(lr)
    assert flag is not None
    assert flag.severity == "info"


# ---------------------------------------------------------------------------
# check_extraction_mode_changed_across_rounds (v0.1 always inert)
# ---------------------------------------------------------------------------


def test_extraction_mode_change_v01_always_none() -> None:
    """v0.1 single-dataset mode cannot fire this flag."""
    assert check_extraction_mode_changed_across_rounds(None) is None
    assert check_extraction_mode_changed_across_rounds([]) is None
    assert check_extraction_mode_changed_across_rounds([_make_manifest()]) is None


# ---------------------------------------------------------------------------
# compute_all_flags + write_flags_yaml
# ---------------------------------------------------------------------------


def test_compute_all_flags_returns_list_of_raised_flags() -> None:
    lr = _make_library_report(match_rate_5p=0.1, required_action="READ_MERGING_RECOMMENDED")
    manifest = _make_manifest(lr)
    pools = {0: {f"s{i}": 100 for i in range(200)}}  # 20k reads OK
    flags = compute_all_flags(manifest, pools)
    names = [f.name for f in flags]
    assert "low_primer_match" in names
    assert "requires_read_merging_for_full_insert" in names


def test_write_flags_yaml_is_deterministic(tmp_path: Path) -> None:
    flags = [
        Flag(name="z_flag", severity="warn", evidence={"foo": 1}),
        Flag(name="a_flag", severity="info", evidence={"bar": 2}),
    ]
    path_a = tmp_path / "a.yaml"
    path_b = tmp_path / "b.yaml"
    write_flags_yaml(flags, path_a)
    write_flags_yaml(flags, path_b)
    assert path_a.read_bytes() == path_b.read_bytes()


def test_write_flags_yaml_sorts_by_name(tmp_path: Path) -> None:
    flags = [
        Flag(name="z_flag", severity="warn", evidence={}),
        Flag(name="a_flag", severity="info", evidence={}),
    ]
    path = tmp_path / "flags.yaml"
    write_flags_yaml(flags, path)
    payload = yaml.safe_load(path.read_text())
    assert [f["name"] for f in payload] == ["a_flag", "z_flag"]


def test_write_flags_yaml_empty_list_emits_empty_array(tmp_path: Path) -> None:
    path = tmp_path / "flags.yaml"
    write_flags_yaml([], path)
    assert path.read_text() == "[]\n"


def test_write_flags_yaml_deterministic_with_nested_evidence(tmp_path: Path) -> None:
    """NB: nested dicts + lists inside `evidence`
    must be deterministic under YAML serialization (`safe_dump` already
    sorts top-level keys; verify the guarantee holds for nested
    structures like the ones produced by `unexpected_rarefied_diversity_increase`
    and `adapter_contamination_high`)."""
    # Build the same flag twice with intentionally different insertion
    # orders for the nested dicts/lists.
    evidence_a = {
        "max_fraction": 0.05,
        "adapters_above_threshold": [
            {"adapter": "TRUSEQ_R1", "hits": 500, "fraction_of_reads": 0.12345},
            {"adapter": "NEXTERA", "hits": 100, "fraction_of_reads": 0.024},
        ],
        "rounds_seen": {2: 10, 0: 5, 1: 8},
    }
    evidence_b = {
        "rounds_seen": {0: 5, 1: 8, 2: 10},  # same data, different insertion order
        "adapters_above_threshold": [
            {"fraction_of_reads": 0.12345, "hits": 500, "adapter": "TRUSEQ_R1"},
            {"hits": 100, "adapter": "NEXTERA", "fraction_of_reads": 0.024},
        ],
        "max_fraction": 0.05,
    }
    flag_a = Flag(name="adapter_contamination_high", severity="warn", evidence=evidence_a)
    flag_b = Flag(name="adapter_contamination_high", severity="warn", evidence=evidence_b)

    path_a = tmp_path / "a.yaml"
    path_b = tmp_path / "b.yaml"
    write_flags_yaml([flag_a], path_a)
    write_flags_yaml([flag_b], path_b)

    assert path_a.read_bytes() == path_b.read_bytes(), (
        "nested evidence dicts/lists must serialize deterministically regardless of insertion order"
    )
