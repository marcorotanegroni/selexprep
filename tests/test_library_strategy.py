"""Unit tests for selexprep.fetch.library_strategy."""

from __future__ import annotations

import pytest

from selexprep.fetch.library_strategy import (
    LIBRARY_STRATEGY_BLOCKLIST,
    classify_study_by_library_strategies,
    is_library_strategy_compatible_with_selex,
)

# ----- per-run helper -----


def test_empty_string_is_compatible() -> None:
    """Missing annotation → treat as compatible (don't punish blank fields)."""
    assert is_library_strategy_compatible_with_selex("") is True
    assert is_library_strategy_compatible_with_selex("   ") is True


def test_other_is_compatible() -> None:
    """``OTHER`` is the standard CV code used by SELEX deposits."""
    assert is_library_strategy_compatible_with_selex("OTHER") is True


def test_selex_is_compatible() -> None:
    """``SELEX`` (rare, not officially registered) is also compatible."""
    assert is_library_strategy_compatible_with_selex("SELEX") is True


def test_targeted_capture_is_compatible() -> None:
    assert is_library_strategy_compatible_with_selex("Targeted-Capture") is True


@pytest.mark.parametrize(
    "strategy",
    [
        "RNA-Seq",
        "ChIP-Seq",
        "ChIP",
        "miRNA-Seq",
        "ATAC-seq",
        "Bisulfite-Seq",
        "WGS",
        "WGA",
        "WXS",
        "Hi-C",
        "RIP-Seq",
        "Tn-Seq",
        "ncRNA-Seq",
        "ssRNA-seq",
    ],
)
def test_known_blocklisted_strategies(strategy: str) -> None:
    assert is_library_strategy_compatible_with_selex(strategy) is False


def test_amplicon_is_compatible_empirical() -> None:
    """calibration: ENA's live API showed 176/10000 real
    HT-SELEX runs tagged ``AMPLICON`` (e.g. PRJDB40017 NRd2 series).
    SELEX rounds ARE PCR-amplified random libraries — AMPLICON is a
    legitimate tag. The initial blocklist assumption that AMPLICON
    implied a non-SELEX assay was empirically wrong; this regression
    pins the corrected behavior."""
    assert is_library_strategy_compatible_with_selex("AMPLICON") is True


def test_blocklist_is_frozenset() -> None:
    """Blocklist is immutable to prevent accidental mutation during runtime."""
    assert isinstance(LIBRARY_STRATEGY_BLOCKLIST, frozenset)


# ----- per-study aggregator -----


def test_all_compatible_runs_keeps_study_not_mixed() -> None:
    """ALL compatible → study kept, not flagged as mixed."""
    c = classify_study_by_library_strategies("PRJ_OK", ["OTHER", "OTHER", ""])
    assert c.should_exclude is False
    assert c.is_mixed_strategy is False
    assert c.n_runs_total == 3
    assert c.n_runs_compatible == 3
    assert c.n_runs_blocklisted == 0
    assert c.blocklisted_strategies == {}


def test_all_blocklisted_runs_excludes_study() -> None:
    """ALL blocklisted → study dropped with exclusion reason."""
    c = classify_study_by_library_strategies(
        "PRJ_RNASEQ",
        ["RNA-Seq", "RNA-Seq", "RNA-Seq"],
    )
    assert c.should_exclude is True
    assert c.is_mixed_strategy is False
    assert c.n_runs_total == 3
    assert c.n_runs_compatible == 0
    assert c.n_runs_blocklisted == 3
    assert c.blocklisted_strategies == {"RNA-Seq": 3}
    assert "RNA-Seq" in c.exclusion_reason
    assert "3 runs" in c.exclusion_reason


def test_mixed_runs_keep_study_flagged_mixed() -> None:
    """SOME compatible + SOME blocklisted → study kept, flagged as mixed
    (audit eligibility layer in will classify as
    MIXED_PROJECT_NEEDS_GROUPING)."""
    c = classify_study_by_library_strategies(
        "PRJ_MIXED",
        ["OTHER", "RNA-Seq", "OTHER"],
    )
    assert c.should_exclude is False
    assert c.is_mixed_strategy is True
    assert c.n_runs_compatible == 2
    assert c.n_runs_blocklisted == 1
    assert c.blocklisted_strategies == {"RNA-Seq": 1}


def test_mixed_multi_strategy_breakdown_sorted_alpha() -> None:
    """When multiple blocklisted strategies present, breakdown is alpha-sorted."""
    c = classify_study_by_library_strategies(
        "PRJ_MULTIBLOCK",
        ["OTHER", "RNA-Seq", "ChIP-Seq", "ChIP-Seq"],
    )
    assert c.is_mixed_strategy is True
    # Insertion-ordered dict; alphabetical: ChIP-Seq before RNA-Seq.
    assert list(c.blocklisted_strategies.keys()) == ["ChIP-Seq", "RNA-Seq"]
    assert c.blocklisted_strategies == {"ChIP-Seq": 2, "RNA-Seq": 1}


def test_empty_run_list_does_not_exclude() -> None:
    """Degenerate case: no runs → keep (defer to downstream classification)."""
    c = classify_study_by_library_strategies("PRJ_EMPTY", [])
    assert c.should_exclude is False
    assert c.is_mixed_strategy is False
    assert c.n_runs_total == 0


def test_singleton_blocklisted_singular_grammar() -> None:
    """Exclusion reason uses 'run' (singular) when n=1."""
    c = classify_study_by_library_strategies("PRJ_X", ["RNA-Seq"])
    assert c.should_exclude is True
    assert "1 run" in c.exclusion_reason


# ----- Empirical pilot regression tests -----
#
# The audit pilot identified five INSDC studies that were
# pulled into the catalog despite their runs being unambiguously not
# SELEX. Pin those exact cases here.


def test_empirical_chip_seq_excluded() -> None:
    """PRJNA1244400 in the pilot was 12-run ChIP-seq input controls."""
    c = classify_study_by_library_strategies(
        "PRJNA1244400",
        ["ChIP-Seq"] * 12,
    )
    assert c.should_exclude is True


def test_empirical_rna_seq_bladder_cancer_excluded() -> None:
    """PRJNA998371 in the pilot was 63 runs of bladder cancer RNA-Seq."""
    c = classify_study_by_library_strategies("PRJNA998371", ["RNA-Seq"] * 63)
    assert c.should_exclude is True


def test_empirical_real_selex_with_other_kept() -> None:
    """A real SELEX deposit (typically library_strategy=OTHER) must be kept."""
    c = classify_study_by_library_strategies(
        "PRJDB19138",
        ["OTHER"] * 5,  # RAPT26-1R..5R were all OTHER in the pilot
    )
    assert c.should_exclude is False
    assert c.is_mixed_strategy is False


def test_empirical_selex_with_one_rna_control_kept_as_mixed() -> None:
    """Hypothetical: a SELEX deposit that included one RNA-Seq control run
    must NOT be dropped — the SELEX runs are real data. Mark mixed
    instead and let the audit-eligibility layer handle it.

    This is the key rule: "do not treat ANY run blocklisted
    as NON_SELEX_ASSAY for the whole BioProject."
    """
    c = classify_study_by_library_strategies(
        "PRJ_REAL_SELEX_WITH_CONTROL",
        ["OTHER", "OTHER", "OTHER", "OTHER", "RNA-Seq"],
    )
    assert c.should_exclude is False
    assert c.is_mixed_strategy is True
    assert c.n_runs_compatible == 4
    assert c.n_runs_blocklisted == 1
