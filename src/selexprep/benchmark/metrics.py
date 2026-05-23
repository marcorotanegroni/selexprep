"""Benchmark metrics for the primer-inference benchmark (locked plan Phase 6, revised).

Scope pivot (2026-05-22): the v0.1 benchmark headline is **paper-reported
primer recovery from accession-derived reads**, not comparator-tool
agreement. AptaPLEX / EasyDIVER+ both require known primers as input;
they cannot benchmark the unique selexprep claim (primer inference from
accession reads). The "do downstream counts agree" question can be
answered intrinsically by comparing inferred-primer runs to
``--override-primer``-driven runs of selexprep itself — see
:func:`compute_count_correlation`, kept here as a Phase 6c entry point
but NOT called from :func:`aggregate_metrics`.

Inputs are :class:`BenchmarkRow` instances that pair a ground-truth entry
(from ``benchmarks/ground_truth.tsv``) with a ``LibraryReport`` (from
``selexprep detect`` output). The top-level :func:`aggregate_metrics`
returns a :class:`BenchmarkMetricsReport` that the Figure A plot consumes
and that ``metrics.json`` serializes.

Headline metrics:

- **Primer recovery**: per-side equivalence outcomes
  (EXACT / REVCOMP / U_T_NORMALIZED / BARCODE_STRIPPED / partials /
  MISMATCH / IUPAC_UNSUPPORTED) — :func:`compute_primer_recovery`.
- **Pair recovery by status**: cross-tab of pair-level recovery
  (pair_exact / pair_equivalent / pair_partial / pair_failed) against
  ``LibraryReport.status`` — :func:`compute_pair_recovery_by_status`.
- **Safe-failure rate**: count of rows where selexprep honestly refused
  (status=UNABLE_TO_INFER or extraction_mode=UNABLE_TO_EXTRACT or
  required_action=MANUAL_PRIMERS_REQUIRED) — :func:`compute_safe_failure_rate`.
  THIS is selexprep's unique distinguishing capability: comparator tools
  cannot refuse because they require primers up-front.
- **N-length recovery**: ±tolerance match vs paper-reported length —
  :func:`compute_n_length_recovery`.
- **Distributions**: ``extraction_mode`` + ``required_action`` —
  :func:`compute_extraction_mode_distribution`,
  :func:`compute_required_action_distribution` (honest-accounting per
  locked plan line 369).

Codex pre-implementation amendments still folded in:

- **Amendment 1**: ``verified=false`` rows are filtered out for all
  numerical metrics; each skip emits a stderr warning; the report
  records ``n_verified`` / ``n_unverified`` counts.
- **Amendment 4** (count_correlation, Phase 6c entry only): Pearson via
  ``pandas.Series.corr``, Spearman via Pearson-of-ranks (no SciPy).
- **Amendment 5** (count_correlation, Phase 6c entry only): union of
  observed + reference sequences with zero-fill; top-K-by-reference
  Pearson labeled as secondary diagnostic only.
"""

from __future__ import annotations

import json
import logging
import sys
from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from selexprep.benchmark.equivalence import EquivalenceResult, primer_equivalent
from selexprep.library.report import LibraryReport, read_library_report_json

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BenchmarkRow:
    """One accession's ground-truth + inferred-LibraryReport pair.

    ``observed_counts`` and ``reference_counts`` are optional per-sequence
    count maps used by :func:`compute_count_correlation`. They are absent
    in Phase 6b.1 (real comparator runs land in 6b.2); the data class
    accepts ``None`` so the scaffolding can be exercised without them.
    """

    accession: str
    library_kind: str
    target_kind: str
    primer_5p_truth: str
    primer_3p_truth: str
    n_length_truth: int
    paper_doi: str
    paper_pmid: str
    round_map_source: str  # "auto" | "curated"
    round_map_path: str
    verified: bool
    notes: str
    library_report: LibraryReport | None
    observed_counts: dict[str, int] | None = None
    reference_counts: dict[str, int] | None = None


@dataclass(frozen=True)
class PrimerSidePair:
    """The equivalence outcome for one row's 5' and 3' primer comparisons."""

    accession: str
    status_5p: EquivalenceResult
    status_3p: EquivalenceResult


@dataclass
class PrimerRecoveryReport:
    """Aggregate primer-recovery counts.

    Each ``EquivalenceKind`` gets its own count for 5' and 3', plus a
    breakdown by ``LibraryReport.status`` so the Figure A panels can
    label "recovery rate by inference confidence". The exact-match count
    on both sides is the headline number; partials and IUPAC-unsupported
    are reported separately (locked plan line 365).
    """

    counts_5p: dict[str, int] = field(default_factory=dict)
    counts_3p: dict[str, int] = field(default_factory=dict)
    counts_by_status_5p: dict[str, dict[str, int]] = field(default_factory=dict)
    counts_by_status_3p: dict[str, dict[str, int]] = field(default_factory=dict)
    n_evaluated: int = 0
    pairs: list[PrimerSidePair] = field(default_factory=list)


@dataclass
class NLengthRecoveryReport:
    """N-length recovery vs paper-reported truth.

    ``tolerance`` is symmetric: ``observed_in_tolerance`` counts rows
    whose ``|observed - truth| <= tolerance``. Rows where the
    LibraryReport's ``n_length_mode`` is ``None`` (e.g.,
    ``PAIRED_END_SPLIT_PRIMERS`` — the full insert spans both reads and
    cannot be measured from R1 alone) count toward ``n_unmeasurable``.
    """

    n_in_tolerance: int = 0
    n_out_of_tolerance: int = 0
    n_unmeasurable: int = 0
    tolerance: int = 2
    per_row: list[dict[str, Any]] = field(default_factory=list)


@dataclass
class CountCorrelationReport:
    """Union+zero-fill correlation between selexprep counts and a reference.

    Phase 6c entry point. The "reference" is selexprep run with
    ``--override-primer-{5p,3p}`` set to the paper-reported truth — the
    self-consistency check (inferred-primer counts vs override-primer
    counts) replaces the original comparator-tool oracle. This dataclass
    keeps the data plumbing so the function is correct by the time the
    6c self-consistency runner lands.

    ``pearson`` / ``spearman`` are the headline numbers on the **union**
    of sequences observed in either source (missing side zero-filled).
    ``top_k_pearson`` is a secondary diagnostic computed on the top-K
    sequences by reference count — clearly labeled, never the headline.
    """

    pearson: float | None = None
    spearman: float | None = None
    top_k_pearson: float | None = None
    top_k: int = 100
    n_union: int = 0
    n_observed_only: int = 0
    n_reference_only: int = 0
    n_in_both: int = 0
    per_row: list[dict[str, Any]] = field(default_factory=list)


@dataclass
class ExtractionModeDistribution:
    """Honest-accounting distribution of ``LibraryReport.extraction_mode``."""

    counts: dict[str, int] = field(default_factory=dict)


@dataclass
class RequiredActionDistribution:
    """Honest-accounting distribution of ``LibraryReport.required_action``.

    Companion to ``ExtractionModeDistribution`` — extraction_mode is
    biology, required_action is workflow guidance (locked plan line 19).
    Both belong in the benchmark accounting because they answer different
    questions.
    """

    counts: dict[str, int] = field(default_factory=dict)


@dataclass
class PairRecoveryByStatus:
    """Cross-tab of pair-level recovery vs ``LibraryReport.status``.

    For each row, both primer comparisons are aggregated into one
    pair-level outcome:

    - ``pair_exact``      — both sides ``EXACT``.
    - ``pair_equivalent`` — both sides matched via any equivalence
      rule (``EXACT`` / ``REVCOMP`` / ``U_T_NORMALIZED`` /
      ``BARCODE_STRIPPED``).
    - ``pair_partial``    — exactly one side matched; the other is
      ``PARTIAL_5P`` / ``PARTIAL_3P`` / ``MISMATCH`` /
      ``IUPAC_UNSUPPORTED``.
    - ``pair_failed``     — neither side matched.

    Bucketed by ``LibraryReport.status`` so the Figure A headline panel
    shows "of HIGH-confidence calls, what fraction recovered the pair
    exactly?" — the Codex-honest selling point.
    """

    counts: dict[str, dict[str, int]] = field(default_factory=dict)
    n_evaluated: int = 0


@dataclass
class SafeFailureRate:
    """Counts of rows where selexprep honestly refused vs miscalled.

    A "safe failure" is any of:

    - ``status == "UNABLE_TO_INFER"``           — inference confidence
      below the threshold.
    - ``extraction_mode == "UNABLE_TO_EXTRACT"`` — both primer match
      rates below ``UNABLE_TO_EXTRACT_MATCH_RATE``.
    - ``required_action == "MANUAL_PRIMERS_REQUIRED"`` — workflow
      flagged as needing curator input.

    These are the distinguishing capability vs known-primer pipelines
    (AptaPLEX / EasyDIVER+ cannot refuse — they require primers up
    front). The benchmark reports the count + a per-reason breakdown
    so the paper can claim "on N of M ambiguous datasets, selexprep
    correctly refused to call rather than silently miscalling".

    ``rate`` is ``n_safe_failures / n_evaluated`` (zero when
    ``n_evaluated == 0``).
    """

    n_evaluated: int = 0
    n_safe_failures: int = 0
    rate: float = 0.0
    by_reason: dict[str, int] = field(default_factory=dict)
    safe_failure_accessions: list[str] = field(default_factory=list)


@dataclass
class BenchmarkMetricsReport:
    """Top-level metrics report. Serializes to deterministic JSON.

    Scope pivot (2026-05-22): ``count_correlation`` is no longer
    populated by :func:`aggregate_metrics`; the field stays in the
    schema as a Phase 6c entry point so the serialized JSON keeps a
    stable shape across versions, but it always has the default
    (empty) values. Phase 6c will populate it via a separate
    self-consistency check (inferred-primer vs paper-override-primer
    selexprep runs).
    """

    n_verified: int = 0
    n_unverified: int = 0
    n_total: int = 0
    skipped_unverified_accessions: list[str] = field(default_factory=list)
    primer_recovery: PrimerRecoveryReport = field(default_factory=PrimerRecoveryReport)
    pair_recovery_by_status: PairRecoveryByStatus = field(default_factory=PairRecoveryByStatus)
    safe_failure_rate: SafeFailureRate = field(default_factory=SafeFailureRate)
    n_length_recovery: NLengthRecoveryReport = field(default_factory=NLengthRecoveryReport)
    extraction_mode_distribution: ExtractionModeDistribution = field(
        default_factory=ExtractionModeDistribution
    )
    required_action_distribution: RequiredActionDistribution = field(
        default_factory=RequiredActionDistribution
    )
    # Phase 6c entry point — populated by an out-of-band self-consistency
    # checker, NOT by aggregate_metrics in 6b.1.
    count_correlation: CountCorrelationReport = field(default_factory=CountCorrelationReport)


# ---------------------------------------------------------------------------
# Loading + row filtering
# ---------------------------------------------------------------------------


def load_ground_truth(path: Path) -> list[BenchmarkRow]:
    """Read ``ground_truth.tsv`` into :class:`BenchmarkRow` instances.

    The ``library_report`` field is ``None`` until paired with the
    selexprep output via :func:`attach_library_reports`. The aggregator
    refuses to score a row whose ``library_report`` is still ``None``
    (matches "report not produced" — surface as MISMATCH at the
    primer-comparison level).
    """
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    rows: list[BenchmarkRow] = []
    for _, r in df.iterrows():
        rows.append(
            BenchmarkRow(
                accession=r["accession"],
                library_kind=r.get("library_kind", ""),
                target_kind=r.get("target_kind", ""),
                primer_5p_truth=r.get("primer_5p_truth", ""),
                primer_3p_truth=r.get("primer_3p_truth", ""),
                n_length_truth=int(r["n_length_truth"]) if r.get("n_length_truth") else 0,
                paper_doi=r.get("paper_doi", ""),
                paper_pmid=r.get("paper_pmid", ""),
                round_map_source=r.get("round_map_source", "auto") or "auto",
                round_map_path=r.get("round_map_path", ""),
                verified=str(r.get("verified", "")).strip().lower() == "true",
                notes=r.get("notes", ""),
                library_report=None,
            )
        )
    return rows


def attach_library_reports(rows: list[BenchmarkRow], reports_dir: Path) -> list[BenchmarkRow]:
    """For each row, attach the matching ``library_report.json`` if present.

    Looks for ``reports_dir/<accession>/library_report.json``. Missing
    reports leave ``library_report=None``; the aggregator handles that
    case by counting the row's primer comparison as a MISMATCH on both
    sides (the inference simply didn't produce a report).
    """
    out: list[BenchmarkRow] = []
    for row in rows:
        lr_path = reports_dir / row.accession / "library_report.json"
        lr = read_library_report_json(lr_path) if lr_path.exists() else None
        out.append(
            BenchmarkRow(
                accession=row.accession,
                library_kind=row.library_kind,
                target_kind=row.target_kind,
                primer_5p_truth=row.primer_5p_truth,
                primer_3p_truth=row.primer_3p_truth,
                n_length_truth=row.n_length_truth,
                paper_doi=row.paper_doi,
                paper_pmid=row.paper_pmid,
                round_map_source=row.round_map_source,
                round_map_path=row.round_map_path,
                verified=row.verified,
                notes=row.notes,
                library_report=lr,
                observed_counts=row.observed_counts,
                reference_counts=row.reference_counts,
            )
        )
    return out


def _filter_verified(rows: list[BenchmarkRow]) -> tuple[list[BenchmarkRow], list[str]]:
    """Drop unverified rows, log a warning per skip, return the skipped accessions."""
    verified: list[BenchmarkRow] = []
    skipped: list[str] = []
    for r in rows:
        if not r.verified:
            logger.warning("benchmark: skipping unverified row accession=%s", r.accession)
            print(
                f"benchmark: skipping unverified row accession={r.accession}",
                file=sys.stderr,
            )
            skipped.append(r.accession)
        else:
            verified.append(r)
    return verified, skipped


# ---------------------------------------------------------------------------
# Metric primitives
# ---------------------------------------------------------------------------


def compute_primer_recovery(rows: list[BenchmarkRow]) -> PrimerRecoveryReport:
    """Aggregate per-row 5' + 3' equivalence outcomes.

    Rows with ``library_report is None`` (no inference output) are
    counted as ``MISMATCH`` on both sides — selexprep failed to recover
    them at all.
    """
    report = PrimerRecoveryReport()
    for r in rows:
        report.n_evaluated += 1
        lr = r.library_report
        primer_5p = lr.primer_5p if lr is not None else None
        primer_3p = lr.primer_3p if lr is not None else None
        status_label = lr.status if lr is not None else "NO_REPORT"

        eq_5p = primer_equivalent(primer_5p, r.primer_5p_truth)
        eq_3p = primer_equivalent(primer_3p, r.primer_3p_truth)

        report.counts_5p[eq_5p.equivalence_kind] = (
            report.counts_5p.get(eq_5p.equivalence_kind, 0) + 1
        )
        report.counts_3p[eq_3p.equivalence_kind] = (
            report.counts_3p.get(eq_3p.equivalence_kind, 0) + 1
        )
        bucket_5p = report.counts_by_status_5p.setdefault(status_label, {})
        bucket_5p[eq_5p.equivalence_kind] = bucket_5p.get(eq_5p.equivalence_kind, 0) + 1
        bucket_3p = report.counts_by_status_3p.setdefault(status_label, {})
        bucket_3p[eq_3p.equivalence_kind] = bucket_3p.get(eq_3p.equivalence_kind, 0) + 1

        report.pairs.append(PrimerSidePair(accession=r.accession, status_5p=eq_5p, status_3p=eq_3p))
    return report


# Equivalence kinds that count as "matched" for pair-level aggregation.
_MATCHING_KINDS: frozenset[str] = frozenset(
    {"EXACT", "REVCOMP", "U_T_NORMALIZED", "BARCODE_STRIPPED"}
)


def compute_pair_recovery_by_status(rows: list[BenchmarkRow]) -> PairRecoveryByStatus:
    """Cross-tab pair-level recovery against ``LibraryReport.status``.

    Per-row classification:

    - both sides ``EXACT``                    → ``pair_exact``
    - both sides matched (any equivalence)    → ``pair_equivalent``
    - exactly one side matched                → ``pair_partial``
    - neither side matched                    → ``pair_failed``

    ``pair_equivalent`` is the superset bucket: it includes
    ``pair_exact`` if you sum them, so the Figure A panel should display
    them as stacked bars (exact at the bottom, equivalent-but-not-exact
    on top) to avoid double-counting.

    Rows with no ``library_report`` are bucketed under
    ``status == "NO_REPORT"`` and counted as ``pair_failed`` —
    selexprep produced no inference at all.
    """
    report = PairRecoveryByStatus()
    for r in rows:
        report.n_evaluated += 1
        lr = r.library_report
        status = lr.status if lr is not None else "NO_REPORT"
        primer_5p = lr.primer_5p if lr is not None else None
        primer_3p = lr.primer_3p if lr is not None else None

        eq_5p = primer_equivalent(primer_5p, r.primer_5p_truth)
        eq_3p = primer_equivalent(primer_3p, r.primer_3p_truth)

        both_exact = eq_5p.equivalence_kind == "EXACT" and eq_3p.equivalence_kind == "EXACT"
        both_matched = (
            eq_5p.equivalence_kind in _MATCHING_KINDS and eq_3p.equivalence_kind in _MATCHING_KINDS
        )
        one_matched = (
            eq_5p.equivalence_kind in _MATCHING_KINDS or eq_3p.equivalence_kind in _MATCHING_KINDS
        )

        if both_exact:
            bucket = "pair_exact"
        elif both_matched:
            bucket = "pair_equivalent"
        elif one_matched:
            bucket = "pair_partial"
        else:
            bucket = "pair_failed"

        status_bucket = report.counts.setdefault(status, {})
        status_bucket[bucket] = status_bucket.get(bucket, 0) + 1

    return report


# Reasons that count as "safe failure" — selexprep honestly refused
# rather than silently miscalling. This is the unique distinguishing
# capability vs known-primer pipelines.
_SAFE_FAILURE_REASONS: tuple[tuple[str, str], ...] = (
    ("status_UNABLE_TO_INFER", "status"),
    ("extraction_mode_UNABLE_TO_EXTRACT", "extraction_mode"),
    ("required_action_MANUAL_PRIMERS_REQUIRED", "required_action"),
)


def compute_safe_failure_rate(rows: list[BenchmarkRow]) -> SafeFailureRate:
    """Count rows where selexprep honestly refused vs silently miscalled.

    A row counts as a safe failure if ANY of the following hold:

    - ``LibraryReport.status == "UNABLE_TO_INFER"``
    - ``LibraryReport.extraction_mode == "UNABLE_TO_EXTRACT"``
    - ``LibraryReport.required_action == "MANUAL_PRIMERS_REQUIRED"``

    A row with no ``library_report`` (the inference run failed
    entirely) is also counted as a safe failure — the downstream
    pipeline refuses to act on a missing report.

    The per-reason breakdown records each individual refusal signal
    (a single row can trigger multiple if e.g. status=UNABLE_TO_INFER
    AND extraction_mode=UNABLE_TO_EXTRACT). The headline
    ``n_safe_failures`` counts rows where AT LEAST ONE reason fired.
    """
    report = SafeFailureRate()
    for r in rows:
        report.n_evaluated += 1
        lr = r.library_report

        # Missing library_report — the inference run itself failed.
        if lr is None:
            report.n_safe_failures += 1
            report.by_reason["no_library_report"] = report.by_reason.get("no_library_report", 0) + 1
            report.safe_failure_accessions.append(r.accession)
            continue

        triggered: list[str] = []
        if lr.status == "UNABLE_TO_INFER":
            triggered.append("status_UNABLE_TO_INFER")
        if lr.extraction_mode == "UNABLE_TO_EXTRACT":
            triggered.append("extraction_mode_UNABLE_TO_EXTRACT")
        if lr.required_action == "MANUAL_PRIMERS_REQUIRED":
            triggered.append("required_action_MANUAL_PRIMERS_REQUIRED")

        if triggered:
            report.n_safe_failures += 1
            report.safe_failure_accessions.append(r.accession)
            for reason in triggered:
                report.by_reason[reason] = report.by_reason.get(reason, 0) + 1

    report.safe_failure_accessions = sorted(report.safe_failure_accessions)
    report.rate = report.n_safe_failures / report.n_evaluated if report.n_evaluated > 0 else 0.0
    return report


def compute_required_action_distribution(rows: list[BenchmarkRow]) -> RequiredActionDistribution:
    """Tally ``LibraryReport.required_action`` values across rows.

    Companion to :func:`compute_extraction_mode_distribution` —
    extraction_mode is biology, required_action is workflow guidance
    (locked plan line 19). Both belong in the honest-accounting
    section because they answer different questions a reviewer would
    ask.
    """
    dist = RequiredActionDistribution()
    for r in rows:
        lr = r.library_report
        key = lr.required_action if lr is not None else "NO_REPORT"
        dist.counts[key] = dist.counts.get(key, 0) + 1
    return dist


def compute_n_length_recovery(
    rows: list[BenchmarkRow], *, tolerance: int = 2
) -> NLengthRecoveryReport:
    """Compare inferred ``n_length_mode`` to ``n_length_truth`` within ±tolerance."""
    report = NLengthRecoveryReport(tolerance=tolerance)
    for r in rows:
        lr = r.library_report
        observed = lr.n_length_mode if lr is not None else None
        per_row_entry: dict[str, Any] = {
            "accession": r.accession,
            "truth": r.n_length_truth,
            "observed": observed,
        }
        if observed is None or r.n_length_truth <= 0:
            report.n_unmeasurable += 1
            per_row_entry["bucket"] = "unmeasurable"
        elif abs(observed - r.n_length_truth) <= tolerance:
            report.n_in_tolerance += 1
            per_row_entry["bucket"] = "in_tolerance"
        else:
            report.n_out_of_tolerance += 1
            per_row_entry["bucket"] = "out_of_tolerance"
        report.per_row.append(per_row_entry)
    return report


def compute_count_correlation(
    rows: list[BenchmarkRow], *, top_k: int = 100
) -> CountCorrelationReport:
    """Pearson + Spearman on the union of observed + reference sequences (zero-filled).

    Codex amendment 5: NEVER use the intersection — it hides the
    sequences only one tool emitted and biases agreement upward. The
    top-K Pearson is a secondary diagnostic only, computed on the
    top-K sequences by reference count.

    For 6b.1 scaffolding most rows lack one or both count maps; rows
    without both are excluded from the headline correlations but still
    enumerated in ``per_row`` so the metrics JSON reflects coverage
    honestly.
    """
    report = CountCorrelationReport(top_k=top_k)

    observed_union: dict[str, int] = Counter()
    reference_union: dict[str, int] = Counter()

    for r in rows:
        oc = r.observed_counts or {}
        rc = r.reference_counts or {}
        if not oc and not rc:
            report.per_row.append({"accession": r.accession, "skipped": "no_counts"})
            continue
        if not oc:
            report.per_row.append({"accession": r.accession, "skipped": "no_observed_counts"})
            continue
        if not rc:
            report.per_row.append({"accession": r.accession, "skipped": "no_reference_counts"})
            continue
        # Per-row also accumulate into the corpus-level union (for the
        # headline correlation across all accessions).
        observed_union.update(oc)
        reference_union.update(rc)
        report.per_row.append(
            {
                "accession": r.accession,
                "n_observed": len(oc),
                "n_reference": len(rc),
                "in_both": len(set(oc) & set(rc)),
            }
        )

    union_keys = set(observed_union) | set(reference_union)
    if len(union_keys) < 2:
        # Single-point (or empty) union has no defined correlation —
        # return the union counts but leave correlation fields as None.
        # Also skips a noisy numpy RuntimeWarning about zero variance.
        report.n_union = len(union_keys)
        for s in union_keys:
            if s in observed_union and s not in reference_union:
                report.n_observed_only += 1
            elif s in reference_union and s not in observed_union:
                report.n_reference_only += 1
            else:
                report.n_in_both += 1
        return report

    obs_series = pd.Series(
        [observed_union.get(seq, 0) for seq in union_keys],
        index=sorted(union_keys),
        dtype=float,
    )
    ref_series = pd.Series(
        [reference_union.get(seq, 0) for seq in union_keys],
        index=sorted(union_keys),
        dtype=float,
    )

    report.n_union = len(union_keys)
    report.n_observed_only = sum(
        1 for s in union_keys if s in observed_union and s not in reference_union
    )
    report.n_reference_only = sum(
        1 for s in union_keys if s in reference_union and s not in observed_union
    )
    report.n_in_both = sum(1 for s in union_keys if s in observed_union and s in reference_union)

    # pandas.Series.corr returns NaN when variance is zero (e.g., all-zero
    # column); convert to None for clean JSON.
    pearson = obs_series.corr(ref_series, method="pearson")
    # Codex amendment 4: pandas.Series.corr(method="spearman") actually
    # imports scipy.stats.spearmanr lazily. To honor "no SciPy" we
    # compute Spearman ourselves as Pearson-of-ranks (the mathematical
    # definition; pandas .rank() handles ties by averaging by default).
    spearman = obs_series.rank().corr(ref_series.rank(), method="pearson")
    report.pearson = float(pearson) if pd.notna(pearson) else None
    report.spearman = float(spearman) if pd.notna(spearman) else None

    # Secondary: top-K by reference count.
    top_keys = ref_series.sort_values(ascending=False).head(top_k).index
    if len(top_keys) >= 2:
        top_obs = obs_series.loc[top_keys]
        top_ref = ref_series.loc[top_keys]
        top_corr = top_obs.corr(top_ref, method="pearson")
        report.top_k_pearson = float(top_corr) if pd.notna(top_corr) else None

    return report


def compute_extraction_mode_distribution(
    rows: list[BenchmarkRow],
) -> ExtractionModeDistribution:
    """Tally ``LibraryReport.extraction_mode`` values across rows (locked plan line 369)."""
    dist = ExtractionModeDistribution()
    for r in rows:
        lr = r.library_report
        key = lr.extraction_mode if lr is not None else "NO_REPORT"
        dist.counts[key] = dist.counts.get(key, 0) + 1
    return dist


# ---------------------------------------------------------------------------
# Top-level aggregator + JSON writer
# ---------------------------------------------------------------------------


def aggregate_metrics(
    rows: list[BenchmarkRow], *, n_length_tolerance: int = 2
) -> BenchmarkMetricsReport:
    """Top-level combiner. Filters ``verified=true`` rows and aggregates the primer-inference metrics.

    Scope pivot (2026-05-22): does NOT compute ``count_correlation`` —
    that's a Phase 6c entry, populated by a separate self-consistency
    checker comparing inferred-primer vs ``--override-primer`` runs.
    The field remains in :class:`BenchmarkMetricsReport` for schema
    stability but stays at default (empty) values.

    The Figure A plot consumes the returned report (via the
    ``metrics.json`` serialization).
    """
    verified, skipped = _filter_verified(rows)
    report = BenchmarkMetricsReport(
        n_verified=len(verified),
        n_unverified=len(skipped),
        n_total=len(rows),
        skipped_unverified_accessions=sorted(skipped),
    )
    report.primer_recovery = compute_primer_recovery(verified)
    report.pair_recovery_by_status = compute_pair_recovery_by_status(verified)
    report.safe_failure_rate = compute_safe_failure_rate(verified)
    report.n_length_recovery = compute_n_length_recovery(verified, tolerance=n_length_tolerance)
    report.extraction_mode_distribution = compute_extraction_mode_distribution(verified)
    report.required_action_distribution = compute_required_action_distribution(verified)
    # count_correlation deliberately not populated — Phase 6c entry point.
    return report


def _equivalence_result_to_dict(eq: EquivalenceResult) -> dict[str, Any]:
    return {
        "matched": eq.matched,
        "equivalence_kind": eq.equivalence_kind,
        "notes": eq.notes,
    }


def _primer_recovery_to_dict(pr: PrimerRecoveryReport) -> dict[str, Any]:
    return {
        "n_evaluated": pr.n_evaluated,
        "counts_5p": dict(sorted(pr.counts_5p.items())),
        "counts_3p": dict(sorted(pr.counts_3p.items())),
        "counts_by_status_5p": {
            k: dict(sorted(v.items())) for k, v in sorted(pr.counts_by_status_5p.items())
        },
        "counts_by_status_3p": {
            k: dict(sorted(v.items())) for k, v in sorted(pr.counts_by_status_3p.items())
        },
        "pairs": [
            {
                "accession": p.accession,
                "status_5p": _equivalence_result_to_dict(p.status_5p),
                "status_3p": _equivalence_result_to_dict(p.status_3p),
            }
            for p in sorted(pr.pairs, key=lambda x: x.accession)
        ],
    }


def _pair_recovery_to_dict(pr: PairRecoveryByStatus) -> dict[str, Any]:
    return {
        "n_evaluated": pr.n_evaluated,
        "counts": {
            status: dict(sorted(bucket.items())) for status, bucket in sorted(pr.counts.items())
        },
    }


def _safe_failure_to_dict(sf: SafeFailureRate) -> dict[str, Any]:
    return {
        "n_evaluated": sf.n_evaluated,
        "n_safe_failures": sf.n_safe_failures,
        "rate": round(sf.rate, 6),
        "by_reason": dict(sorted(sf.by_reason.items())),
        "safe_failure_accessions": sf.safe_failure_accessions,
    }


def write_metrics_json(report: BenchmarkMetricsReport, path: Path) -> None:
    """Serialize ``report`` to deterministic JSON (sorted keys + stable ordering)."""
    payload: dict[str, Any] = {
        "n_verified": report.n_verified,
        "n_unverified": report.n_unverified,
        "n_total": report.n_total,
        "skipped_unverified_accessions": report.skipped_unverified_accessions,
        "primer_recovery": _primer_recovery_to_dict(report.primer_recovery),
        "pair_recovery_by_status": _pair_recovery_to_dict(report.pair_recovery_by_status),
        "safe_failure_rate": _safe_failure_to_dict(report.safe_failure_rate),
        "n_length_recovery": asdict(report.n_length_recovery),
        "extraction_mode_distribution": {
            "counts": dict(sorted(report.extraction_mode_distribution.counts.items()))
        },
        "required_action_distribution": {
            "counts": dict(sorted(report.required_action_distribution.counts.items()))
        },
        # Phase 6c entry point — present in schema for stability, populated
        # by a separate self-consistency checker.
        "count_correlation": asdict(report.count_correlation),
    }
    text = json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    path.write_text(text, encoding="utf-8")


# ---------------------------------------------------------------------------
# CLI entry (called by the Snakefile's compute_metrics rule)
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    import argparse

    p = argparse.ArgumentParser(description="Aggregate Phase 6 benchmark metrics.")
    p.add_argument("--ground-truth", required=True, type=Path)
    p.add_argument("--reports-dir", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--n-length-tolerance", type=int, default=2)
    args = p.parse_args(argv)

    rows = load_ground_truth(args.ground_truth)
    rows = attach_library_reports(rows, args.reports_dir)
    report = aggregate_metrics(rows, n_length_tolerance=args.n_length_tolerance)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_metrics_json(report, args.out)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
