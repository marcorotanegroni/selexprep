"""Per-run + per-BioProject library_strategy classification for catalog hygiene.

Phase 6b.5a — empirically motivated. The N=30 Tier 2 audit pilot found
that ENA studies whose runs are unambiguously RNA-Seq / ChIP-Seq /
miRNA-Seq / cell-treatment timecourses were ending up in the bundled
discovery catalog because the broad-recall text search matched their
abstracts. This module classifies them out at refresh time.

**Per-run first, then per-BioProject** (Phase 6b.5a user amendment).
Mixed BioProjects exist (SELEX runs alongside controls or adjacent
assays); treating "ANY run blocklisted → exclude the whole study"
would throw away real SELEX runs. The decision rule:

- ALL runs blocklisted              → ``should_exclude=True`` (drop +
                                         record exclusion reason)
- SOME blocklisted + SOME compatible → ``should_exclude=False,
                                         is_mixed_strategy=True``
                                         (keep; the audit eligibility
                                         layer in 6b.5b will classify
                                         as MIXED_PROJECT_NEEDS_GROUPING)
- ALL compatible                     → ``should_exclude=False,
                                         is_mixed_strategy=False``

Used by:

- :mod:`selexprep.catalog.rebuild` — the ``selexprep catalog refresh``
  code path that regenerates the bundled ``bioprojects.csv``.
- :class:`selexprep.fetch.discover.ENAAdapter` — the initial-bootstrap
  discovery pipeline (multi-source).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

#: ENA's `library_strategy` controlled-vocabulary codes that cannot
#: plausibly be SELEX data.
#:
#: **Empirical calibration (Phase 6b.5a — ENA Portal API live query on
#: 2026-05-24).** Real HT-SELEX deposits use the following CV codes:
#:
#: - ``SELEX``           (most common — 9824/10000 in the HT-SELEX query)
#: - ``OTHER``           (typical fallback when no formal code matches)
#: - ``AMPLICON``        (legit — SELEX rounds ARE PCR-amplified random
#:                        libraries; e.g. PRJDB40017 NRd2 series is tagged
#:                        AMPLICON; 176/10000 in HT-SELEX, 30/327 in
#:                        SELEX-seq, 261/810 in aptamer)
#: - ``Targeted-Capture`` (rare; e.g. 20/810 in aptamer query)
#: - ``""``              (field not annotated)
#:
#: AMPLICON was initially in the blocklist (wrong assumption that
#: amplicon-sequencing implied a different assay); empirical evidence
#: removed it.
#:
#: The blocklist enumerates CV codes that ARE formal labels for *other*
#: assay types. A run carrying one of these labels cannot be a SELEX
#: round by construction (the depositor marked it as something else).
#: Adding here is conservative; removing requires evidence that real
#: SELEX deposits use the code.
LIBRARY_STRATEGY_BLOCKLIST: frozenset[str] = frozenset(
    {
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
        "RAD-Seq",
        "CLONE",
        "CLONEEND",
        "EST",
        "FINISHING",
        "FL-cDNA",
        "MeDIP-Seq",
        "MNase-Seq",
        "MRE-Seq",
        "MBD-Seq",
        "POOLCLONE",
        "RIP-Seq",
        "Tn-Seq",
        "VALIDATION",
        "ncRNA-Seq",
        "ssRNA-seq",
        "DNase-Hypersensitivity",
    }
)


def is_library_strategy_compatible_with_selex(strategy: str) -> bool:
    """Per-run check: could this ``library_strategy`` plausibly be SELEX?

    Returns True for:

    - Empty/whitespace strings (field not annotated).
    - Any value NOT in :data:`LIBRARY_STRATEGY_BLOCKLIST`.

    Returns False for:

    - Any value in :data:`LIBRARY_STRATEGY_BLOCKLIST` (unambiguously
      not SELEX).
    """
    s = strategy.strip()
    return not s or s not in LIBRARY_STRATEGY_BLOCKLIST


@dataclass
class StudyStrategyClassification:
    """Per-BioProject summary of library_strategy across all its runs.

    See module docstring for the decision rule. ``bioprojects_excluded.csv``
    rows are produced from instances with ``should_exclude=True``; the
    ``exclusion_reason`` field is the human-readable summary written
    there.
    """

    bioproject_id: str
    n_runs_total: int = 0
    n_runs_compatible: int = 0
    n_runs_blocklisted: int = 0
    #: ``{strategy_value: run_count}`` for blocklisted strategies only,
    #: sorted alphabetically. Empty for all-compatible studies.
    blocklisted_strategies: dict[str, int] = field(default_factory=dict)
    #: True iff every run was blocklisted — the study is dropped from
    #: the catalog and recorded in ``bioprojects_excluded.csv``.
    should_exclude: bool = False
    #: True iff at least one run was compatible AND at least one was
    #: blocklisted. The study is KEPT in the catalog (its compatible
    #: runs are real SELEX data), but flagged for the audit's
    #: MIXED_PROJECT_NEEDS_GROUPING bucket in Phase 6b.5b.
    is_mixed_strategy: bool = False
    #: Populated iff ``should_exclude`` — written to the sidecar CSV.
    exclusion_reason: str = ""


def classify_study_by_library_strategies(
    bioproject_id: str,
    library_strategies: list[str],
) -> StudyStrategyClassification:
    """Aggregate per-run library_strategy values into a study-level decision.

    Parameters
    ----------
    bioproject_id
        The study accession (e.g. ``"PRJNA1244400"``).
    library_strategies
        One entry per run in the study, in any order. Empty / missing
        values are treated as compatible (see
        :func:`is_library_strategy_compatible_with_selex`).

    Returns
    -------
    :class:`StudyStrategyClassification`
    """
    result = StudyStrategyClassification(bioproject_id=bioproject_id)
    result.n_runs_total = len(library_strategies)
    blocklisted_counter: Counter[str] = Counter()

    for s in library_strategies:
        if is_library_strategy_compatible_with_selex(s):
            result.n_runs_compatible += 1
        else:
            result.n_runs_blocklisted += 1
            blocklisted_counter[s.strip()] += 1

    result.blocklisted_strategies = dict(sorted(blocklisted_counter.items()))

    # Degenerate case: ENA returned the study but no runs (malformed
    # response or upstream filter dropped them). Keep — the audit's
    # downstream classification will mark it as NO_ROUND_STRUCTURE or
    # FETCH_DEAD when it tries to actually use the study.
    if result.n_runs_total == 0:
        return result

    if result.n_runs_compatible == 0:
        # 100% blocklisted — drop the whole study.
        result.should_exclude = True
        strategies_summary = ", ".join(
            f"{s} ({n} run{'s' if n != 1 else ''})"
            for s, n in result.blocklisted_strategies.items()
        )
        result.exclusion_reason = (
            f"all {result.n_runs_total} runs use blocklisted "
            f"library_strategy values: {strategies_summary}"
        )
        return result

    # At least one compatible run — keep. Flag mixed if blocklisted
    # runs are also present.
    result.is_mixed_strategy = result.n_runs_blocklisted > 0
    return result


__all__ = [
    "LIBRARY_STRATEGY_BLOCKLIST",
    "StudyStrategyClassification",
    "classify_study_by_library_strategies",
    "is_library_strategy_compatible_with_selex",
]
