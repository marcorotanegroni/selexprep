"""Benchmark metrics for the primer-inference benchmark.

Scope pivot (2026-05-22): the v0.1 benchmark headline is **paper-reported
primer recovery from accession-derived reads**, not comparator-tool
agreement. AptaPLEX / EasyDIVER+ both require known primers as input;
they cannot benchmark the unique selexprep claim (primer inference from
accession reads). The "do downstream counts agree" question can be
answered intrinsically by comparing inferred-primer runs to
``--override-primer``-driven runs of selexprep itself — see
:func:`compute_count_correlation`, kept here as an entry point
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
  the design).

Additional behaviors:

- ``verified=false`` rows are filtered out for all
  numerical metrics; each skip emits a stderr warning; the report
  records ``n_verified`` / ``n_unverified`` counts.
- ``count_correlation`` (entry only): Pearson via
  ``pandas.Series.corr``, Spearman via Pearson-of-ranks (no SciPy).
- ``count_correlation`` (entry only): union of
  observed + reference sequences with zero-fill; top-K-by-reference
  Pearson labeled as secondary diagnostic only.
"""

from __future__ import annotations

import json
import logging
import sys
from collections import Counter
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any

import pandas as pd

from selexprep.benchmark.equivalence import EquivalenceResult, primer_equivalent
from selexprep.fetch.plan import FetchPlan, fastq_filenames_for_run
from selexprep.fetch.runner import read_fetch_metadata_json
from selexprep.library.report import LibraryReport, read_library_report_json

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FetchStats:
    """Per-accession fetch accounting, honestly inferable post-fetch.

    Derived from the :class:`~selexprep.fetch.plan.FetchPlan` audit trail
    (``fetch_metadata.json``: the PLAN, i.e. what ENA said *should* exist)
    cross-referenced with the on-disk primary-read ``fastqs.manifest``
    (R1/single-end FASTQs; R2 mates live in ``fastqs.r2.manifest`` for
    paired-aware detection). the plan is a
    PLAN, not an outcome log, so the field names stay honest — we report
    expected / available / missing / no-fastq-url, and NEVER claim
    ``failed_runs`` (a run with no on-disk FASTQ might be embargoed,
    no-URL, or a genuine download failure — the plan can't tell us which).

    - ``fetch_expected_runs`` — runs in the plan.
    - ``fetch_available_runs`` — expected runs whose R1 FASTQ is on disk.
    - ``fetch_missing_runs`` — expected minus available (includes no-URL runs).
    - ``runs_with_no_fastq_url`` — runs whose ENA record had no
      ``fastq_ftp`` URL (an upstream-availability explanation for some
      missing runs; cf. PRJEB70964's 10 no-URL runs).
    - ``unassigned_missing_or_skipped`` — missing runs that DID have a URL
      (``missing minus runs_with_no_fastq_url``). honestly
      named — from the plan + manifest alone we cannot prove a run was
      *skipped* for round-unassignment vs *failed* mid-download; true
      failed-vs-skipped would need the Snakefile to persist a
      ``fetch_result.json`` (a deferred option). cf. PRJEB22637 / PRJNA315881
      where the missing run is an unassigned-round skip.
    - ``partial_fetch`` — some but not all expected runs landed
      (``0 < available < expected``). A total miss (available == 0) is
      NOT a partial fetch; a complete fetch (available == expected) is not
      either.
    """

    accession: str
    fetch_expected_runs: int = 0
    fetch_available_runs: int = 0
    fetch_missing_runs: int = 0
    runs_with_no_fastq_url: int = 0
    unassigned_missing_or_skipped: int = 0
    partial_fetch: bool = False


@dataclass(frozen=True)
class BenchmarkRow:
    """One accession's ground-truth + inferred-LibraryReport pair.

    Adds the read-state arm label + orthogonal limitation
    flags: ``read_state`` is ``raw_standard`` (recovery
    arm) or ``pre_trimmed`` (specificity arm), labeled from independent
    read-length + architecture evidence — NEVER from detect output. The
    flags (``mono_round`` / ``partial_fetch`` / ``paired_end_r1_only`` /
    ``demultiplexed``) carry limitations orthogonally so they aren't
    conflated with read-state. ``fetch_stats`` is attached post-fetch by
    :func:`attach_fetch_stats`.

    ``observed_counts`` and ``reference_counts`` are optional per-sequence
    count maps used by :func:`compute_count_correlation`. They are absent
    in the standard flow (only the optional comparator path populates them); the data class
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
    # read-state arm + orthogonal limitation flags.
    read_state: str = ""  # "raw_standard" | "pre_trimmed" | ""
    mono_round: bool = False
    partial_fetch: bool = False
    paired_end_r1_only: bool = False
    demultiplexed: bool = False
    # score the 3' side only when its truth is
    # paper-grounded. False ⇒ the 3' truth is read-resolved/diagnostic
    # (e.g. PRJNA883192, whose 3' truth was derived from the deposited
    # reads), so scoring it against detect's read-derived 3' would be
    # circular — the 5' still scores normally.
    score_3p: bool = True
    fetch_stats: FetchStats | None = None


@dataclass(frozen=True)
class PrimerSidePair:
    """The equivalence outcome for one row's 5' and 3' primer comparisons."""

    accession: str
    status_5p: EquivalenceResult
    status_3p: EquivalenceResult
    # False ⇒ the 3' truth is read-resolved/diagnostic and was NOT scored
    # toward paper-grounded recovery (the pair was classified on the 5').
    score_3p: bool = True


@dataclass
class PrimerRecoveryReport:
    """Aggregate primer-recovery counts.

    Each ``EquivalenceKind`` gets its own count for 5' and 3', plus a
    breakdown by ``LibraryReport.status`` so the Figure A panels can
    label "recovery rate by inference confidence". The exact-match count
    on both sides is the headline number; partials and IUPAC-unsupported
    are reported separately.
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

    entry point. The "reference" is selexprep run with
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
    biology, required_action is workflow guidance.
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
    - ``pair_partial``    — at least one side carried informative recovery
      (matched OR a boundary ``PARTIAL_5P`` / ``PARTIAL_3P``) but the pair
      did not both fully match. Includes both-sides-PARTIAL (correct
      regions, fuzzy boundaries).
    - ``pair_failed``     — NO informative recovery on either side (both
      ``MISMATCH`` / ``IUPAC_UNSUPPORTED`` / no report).

    Bucketed by ``LibraryReport.status`` so the Figure A headline panel
    shows "of HIGH-confidence calls, what fraction recovered the pair
    exactly?" — the honest selling point.
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
class SpecificityReport:
    """Specificity arm — no-false-call on pre-trimmed deposits.

    . Evaluated ONLY over ``read_state ==
    "pre_trimmed"`` rows: deposits whose reads are the N-region only
    (flanking constants trimmed off before deposit), independently
    classified from read-length + architecture evidence. The expected
    behavior is that selexprep emits NO primer call.

    - ``n_no_false_call`` — ``primer_5p is None and primer_3p is None``
      (or no report at all): selexprep correctly refused to fabricate a
      primer that isn't in the reads. This is the specificity "gold".
    - ``n_false_positive`` — selexprep emitted any primer on a
      pre-trimmed deposit (a fabricated call).

    Headline: *"On independently classified
    pre-trimmed deposits, selexprep made N_false/N false-positive primer
    calls."* This does NOT claim detect *detects* pre-trimming — the
    benchmark layer knows the read-state; detect only reports that no
    primers were inferred.
    """

    n_evaluated: int = 0  # rows WITH a report (no_false_call + false_positive)
    n_no_false_call: int = 0
    n_false_positive: int = 0
    # a pre_trimmed row with no library_report is NOT a
    # no-call success — it's not_evaluable (the inference run produced
    # nothing to judge). It is excluded from the n_evaluated denominator.
    n_not_evaluable: int = 0
    false_positive_accessions: list[str] = field(default_factory=list)
    not_evaluable_accessions: list[str] = field(default_factory=list)
    per_row: list[dict[str, Any]] = field(default_factory=list)


@dataclass
class BenchmarkMetricsReport:
    """Top-level metrics report. Serializes to deterministic JSON.

    Scope pivot (2026-05-22): ``count_correlation`` is no longer
    populated by :func:`aggregate_metrics`; the field stays in the
    schema as an entry point so the serialized JSON keeps a
    stable shape across versions, but it always has the default
    (empty) values. will populate it via a separate
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
    # two-arm sensitivity/specificity reframe.
    # recovery_denominator = n raw_standard rows (the recovery-arm
    # denominator); recovery metrics above run ONLY on this arm.
    recovery_denominator: int = 0
    # specificity arm — no-false-call on pre_trimmed deposits.
    specificity: SpecificityReport = field(default_factory=SpecificityReport)
    # adapter-control panel — no-false-call on adapter_control deposits
    # (constant collides with a known adapter; out of the recovery denominator).
    adapter_control: SpecificityReport = field(default_factory=SpecificityReport)
    # multi-round-only sensitivity sub-report (recovery arm minus
    # mono_round rows) — so confidence-limited single-round rows don't
    # silently deflate the headline (don't drop hard rows
    # post-hoc; report them in a sub-arm).
    multi_round_sensitivity: PairRecoveryByStatus = field(default_factory=PairRecoveryByStatus)
    # per-accession fetch accounting (honest field names; no failed_runs).
    fetch_stats: dict[str, FetchStats] = field(default_factory=dict)
    # entry point — populated by an out-of-band self-consistency
    # checker, NOT by aggregate_metrics in .
    count_correlation: CountCorrelationReport = field(default_factory=CountCorrelationReport)


# ---------------------------------------------------------------------------
# Loading + row filtering
# ---------------------------------------------------------------------------


def _truthy(value: Any) -> bool:
    """Parse a TSV cell as a boolean (``"true"`` / ``"false"`` / blank)."""
    return str(value).strip().lower() == "true"


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
                verified=_truthy(r.get("verified", "")),
                notes=r.get("notes", ""),
                library_report=None,
                # read-state arm + flags (default to ""/False
                # so older ground_truth.tsv files without these columns still
                # parse; the arm split treats "" as neither arm).
                read_state=str(r.get("read_state", "") or "").strip(),
                mono_round=_truthy(r.get("mono_round", "")),
                partial_fetch=_truthy(r.get("partial_fetch", "")),
                paired_end_r1_only=_truthy(r.get("paired_end_r1_only", "")),
                demultiplexed=_truthy(r.get("demultiplexed", "")),
                # score_3p defaults True; only an explicit "false" disables
                # 3' scoring (read-resolved/diagnostic 3' truth).
                score_3p=str(r.get("score_3p", "") or "").strip().lower() != "false",
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
        # ``replace`` copies every field (incl. read_state + flags +
        # fetch_stats) so this can't silently drop new fields as the
        # schema grows.
        out.append(replace(row, library_report=lr))
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
# Fetch accounting (per-accession; honest field names, no failed_runs)
# ---------------------------------------------------------------------------


def compute_fetch_stats_from_plan(
    accession: str, plan: FetchPlan, manifest_basenames: set[str]
) -> FetchStats:
    """Cross-reference a FetchPlan (what ENA said exists) against primary FASTQs.

    ``manifest_basenames`` is the set of FASTQ *basenames* present in the
    primary ``fastqs.manifest`` (R1 for paired-end, the only FASTQ for
    single-end). R2 mates are tracked separately for detect via
    ``fastqs.r2.manifest`` and are intentionally not needed to decide whether
    a run landed. A run counts as *available* iff at least one of its non-R2
    FASTQ filenames is in that set.

    Honest semantics: the plan is a PLAN, not an outcome
    log, so we report expected / available / missing / no-URL only —
    never ``failed_runs`` (a missing run might be embargoed, no-URL, or a
    genuine failure; the plan cannot disambiguate).
    """
    expected = len(plan.runs)
    runs_with_no_url = sum(1 for run in plan.runs if not run.fastq_urls)
    available = 0
    for run in plan.runs:
        r1_names = [n for n in fastq_filenames_for_run(run) if not n.endswith("_2.fastq.gz")]
        if any(name in manifest_basenames for name in r1_names):
            available += 1
    missing = max(0, expected - available)
    # Of the missing runs, the ones that HAD a URL (so not an upstream
    # no-URL case) — skipped for round-unassignment or failed mid-download;
    # the plan can't disambiguate, hence the honest name.
    unassigned_missing_or_skipped = max(0, missing - runs_with_no_url)
    partial = 0 < available < expected
    return FetchStats(
        accession=accession,
        fetch_expected_runs=expected,
        fetch_available_runs=available,
        fetch_missing_runs=missing,
        runs_with_no_fastq_url=runs_with_no_url,
        unassigned_missing_or_skipped=unassigned_missing_or_skipped,
        partial_fetch=partial,
    )


def _read_manifest_basenames(manifest_path: Path) -> set[str]:
    """Read a ``fastqs.manifest`` (one FASTQ path per line) into basenames."""
    basenames: set[str] = set()
    for line in manifest_path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped:
            basenames.add(Path(stripped).name)
    return basenames


def load_fetch_stats(reports_dir: Path, accession: str) -> FetchStats | None:
    """Load ``fetch_metadata.json`` + ``fastqs.manifest`` for one accession.

    Returns ``None`` when no ``fetch_metadata.json`` exists — fetch_stats
    is then simply absent for that accession rather than fabricated.
    """
    meta_path = reports_dir / accession / "fetch_metadata.json"
    if not meta_path.exists():
        return None
    plan = read_fetch_metadata_json(meta_path)
    manifest_path = reports_dir / accession / "fastqs.manifest"
    manifest_basenames = (
        _read_manifest_basenames(manifest_path) if manifest_path.exists() else set()
    )
    return compute_fetch_stats_from_plan(accession, plan, manifest_basenames)


def attach_fetch_stats(rows: list[BenchmarkRow], reports_dir: Path) -> list[BenchmarkRow]:
    """Attach per-accession :class:`FetchStats` parsed from the fetch audit trail.

    Mirrors :func:`attach_library_reports`. Rows with no
    ``fetch_metadata.json`` keep ``fetch_stats=None``.
    """
    return [replace(row, fetch_stats=load_fetch_stats(reports_dir, row.accession)) for row in rows]


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
        bucket_5p = report.counts_by_status_5p.setdefault(status_label, {})
        bucket_5p[eq_5p.equivalence_kind] = bucket_5p.get(eq_5p.equivalence_kind, 0) + 1
        # The 3' side is tallied only when it's a paper-grounded truth. A
        # read-resolved/diagnostic 3' (score_3p=false, e.g. PRJNA883192) would
        # be circular to score against detect's read-derived 3', so it is
        # excluded from the 3' recovery counts (the 5' still scores normally).
        if r.score_3p:
            report.counts_3p[eq_3p.equivalence_kind] = (
                report.counts_3p.get(eq_3p.equivalence_kind, 0) + 1
            )
            bucket_3p = report.counts_by_status_3p.setdefault(status_label, {})
            bucket_3p[eq_3p.equivalence_kind] = bucket_3p.get(eq_3p.equivalence_kind, 0) + 1

        report.pairs.append(
            PrimerSidePair(
                accession=r.accession, status_5p=eq_5p, status_3p=eq_3p, score_3p=r.score_3p
            )
        )
    return report


# Equivalence kinds that count as "matched" for pair-level aggregation.
_MATCHING_KINDS: frozenset[str] = frozenset(
    {"EXACT", "REVCOMP", "U_T_NORMALIZED", "BARCODE_STRIPPED"}
)
# Boundary partials: the correct constant region was localized but its extent
# is off (e.g. detect under-extends a long T7-containing 5' flank). Not a full
# match, but informative recovery — NOT pair_failed.
_PARTIAL_KINDS: frozenset[str] = frozenset({"PARTIAL_5P", "PARTIAL_3P"})
# A side carries informative recovery if it matched OR was a boundary partial.
_INFORMATIVE_KINDS: frozenset[str] = _MATCHING_KINDS | _PARTIAL_KINDS


def compute_pair_recovery_by_status(rows: list[BenchmarkRow]) -> PairRecoveryByStatus:
    """Cross-tab pair-level recovery against ``LibraryReport.status``.

    Per-row classification:

    - both sides ``EXACT``                    → ``pair_exact``
    - both sides matched (any equivalence)    → ``pair_equivalent``
    - ≥1 side informative (matched/partial)   → ``pair_partial``
    - no informative recovery either side     → ``pair_failed``

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

        if r.score_3p:
            both_exact = eq_5p.equivalence_kind == "EXACT" and eq_3p.equivalence_kind == "EXACT"
            both_matched = (
                eq_5p.equivalence_kind in _MATCHING_KINDS
                and eq_3p.equivalence_kind in _MATCHING_KINDS
            )
            # pair_failed is reserved for NO informative recovery
            # on either side (both MISS: MISMATCH / IUPAC_UNSUPPORTED / no
            # report). A side that matched OR was a boundary partial (correct
            # region, fuzzy extent) is informative → the pair is at least
            # pair_partial. This fixes the both-sides-PARTIAL case, which
            # previously fell through to pair_failed.
            any_informative = (
                eq_5p.equivalence_kind in _INFORMATIVE_KINDS
                or eq_3p.equivalence_kind in _INFORMATIVE_KINDS
            )
        else:
            # Read-resolved/diagnostic 3' (score_3p=false): only the 5' is
            # paper-grounded, so never claim a both-sides pair — classify on
            # the 5' alone (informative 5' → pair_partial; else pair_failed).
            both_exact = False
            both_matched = False
            any_informative = eq_5p.equivalence_kind in _INFORMATIVE_KINDS

        if both_exact:
            bucket = "pair_exact"
        elif both_matched:
            bucket = "pair_equivalent"
        elif any_informative:
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
    . Both belong in the honest-accounting
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
        extraction_mode = lr.extraction_mode if lr is not None else None
        per_row_entry: dict[str, Any] = {
            "accession": r.accession,
            "truth": r.n_length_truth,
            "observed": observed,
        }
        # N-length recovery is only meaningful when BOTH
        # flanks were recovered in a single read (the random region is
        # bounded on both sides). When primers weren't recovered —
        # UNABLE_TO_EXTRACT (null primers; n_length_mode degenerates to the
        # raw read length) or FIVE/THREE_PRIME_ONLY (one boundary unknown)
        # — the row is not_evaluable for N, NOT in/out of tolerance.
        # Otherwise a null-primer row could be spuriously credited as
        # "in tolerance" (the PRJEB70964 N=35 artifact).
        if extraction_mode != "BOTH_PRIMERS_SINGLE_READ":
            report.n_unmeasurable += 1
            per_row_entry["bucket"] = "unmeasurable"
            per_row_entry["reason"] = (
                "no_report" if lr is None else f"extraction_mode={extraction_mode}"
            )
        elif observed is None or r.n_length_truth <= 0:
            report.n_unmeasurable += 1
            per_row_entry["bucket"] = "unmeasurable"
            per_row_entry["reason"] = "missing_observed_or_truth"
        elif abs(observed - r.n_length_truth) <= tolerance:
            report.n_in_tolerance += 1
            per_row_entry["bucket"] = "in_tolerance"
        else:
            report.n_out_of_tolerance += 1
            per_row_entry["bucket"] = "out_of_tolerance"
        report.per_row.append(per_row_entry)
    return report


def _compute_no_false_call(rows: list[BenchmarkRow], read_state: str) -> SpecificityReport:
    """No-false-call report over rows of a given ``read_state``.

    Correct (no false positive) iff selexprep emitted no primer —
    ``primer_5p is None and primer_3p is None``. Any emitted primer is a
    false positive (selexprep fabricated a constant that isn't recoverable).
    A row with no ``library_report`` is ``not_evaluable`` and
    is excluded from the ``n_evaluated`` denominator.

    Shared by the specificity arm (``pre_trimmed``) and the adapter-control
    panel (``adapter_control``): both expect a refusal, differing only in
    WHY (pre-trimming vs adapter-collision), so the measurement is identical.
    """
    report = SpecificityReport()
    for r in rows:
        if r.read_state != read_state:
            continue
        lr = r.library_report
        if lr is None:
            # No report ⇒ not_evaluable (NOT a no-call success).
            report.n_not_evaluable += 1
            report.not_evaluable_accessions.append(r.accession)
            report.per_row.append(
                {
                    "accession": r.accession,
                    "bucket": "not_evaluable",
                    "primer_5p": None,
                    "primer_3p": None,
                }
            )
            continue
        report.n_evaluated += 1
        primer_5p = lr.primer_5p
        primer_3p = lr.primer_3p
        emitted = primer_5p is not None or primer_3p is not None
        if emitted:
            report.n_false_positive += 1
            report.false_positive_accessions.append(r.accession)
            bucket = "false_positive"
        else:
            report.n_no_false_call += 1
            bucket = "no_false_call"
        report.per_row.append(
            {
                "accession": r.accession,
                "bucket": bucket,
                "primer_5p": primer_5p,
                "primer_3p": primer_3p,
            }
        )
    report.false_positive_accessions = sorted(report.false_positive_accessions)
    report.not_evaluable_accessions = sorted(report.not_evaluable_accessions)
    return report


def compute_specificity(rows: list[BenchmarkRow]) -> SpecificityReport:
    """Specificity arm — no-false-call on ``read_state == "pre_trimmed"`` deposits.

    Pre-trimmed reads are the N-region only, so the correct behavior is to
    emit NO primer. Delegates to :func:`_compute_no_false_call`.
    """
    return _compute_no_false_call(rows, "pre_trimmed")


def compute_adapter_control(rows: list[BenchmarkRow]) -> SpecificityReport:
    """Adapter-control panel — no-false-call on ``read_state == "adapter_control"``.

    Deposits whose constant collides with a known sequencing adapter (e.g.
    PRJEB70964: 5' constant = revcomp(TruSeq R1)). Excluded from the recovery
    denominator — scoring "refused to call an adapter a primer" as a recovery
    miss would be unfair — and reported here as a NEGATIVE CONTROL where the
    correct behavior is no primer call. Delegates to :func:`_compute_no_false_call`.
    """
    return _compute_no_false_call(rows, "adapter_control")


def compute_count_correlation(
    rows: list[BenchmarkRow], *, top_k: int = 100
) -> CountCorrelationReport:
    """Pearson + Spearman on the union of observed + reference sequences (zero-filled).

    NEVER use the intersection — it hides the
    sequences only one tool emitted and biases agreement upward. The
    top-K Pearson is a secondary diagnostic only, computed on the
    top-K sequences by reference count.

    For scaffolding most rows lack one or both count maps; rows
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
    # pandas.Series.corr(method="spearman") actually
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
    """Tally ``LibraryReport.extraction_mode`` values across rows."""
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
    that's an entry, populated by a separate self-consistency
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

    # split by read_state into benchmark roles. Recovery
    # metrics (sensitivity) run ONLY on the raw_standard arm, where primers
    # are physically present and recovery is even applicable. The
    # specificity arm (pre_trimmed) and adapter-control panel
    # (adapter_control) self-filter inside their no-false-call computations.
    # Rows with an empty/unknown read_state belong to no role and are
    # excluded from all of them (but still count in the honest-accounting
    # distributions below).
    raw_standard = [r for r in verified if r.read_state == "raw_standard"]

    report.recovery_denominator = len(raw_standard)
    report.primer_recovery = compute_primer_recovery(raw_standard)
    report.pair_recovery_by_status = compute_pair_recovery_by_status(raw_standard)
    report.n_length_recovery = compute_n_length_recovery(raw_standard, tolerance=n_length_tolerance)

    # Multi-round-only sensitivity sub-report: recovery arm minus
    # mono_round rows: don't drop confidence-limited
    # single-round rows post-hoc; report them in a sub-arm.
    multi_round = [r for r in raw_standard if not r.mono_round]
    report.multi_round_sensitivity = compute_pair_recovery_by_status(multi_round)

    # Specificity arm — no-false-call on pre_trimmed deposits.
    report.specificity = compute_specificity(verified)
    # Adapter-control panel — no-false-call on adapter_control deposits
    # (e.g. PRJEB70964, 5' constant = revcomp(TruSeq R1)); out of recovery.
    report.adapter_control = compute_adapter_control(verified)

    # Honest accounting runs across ALL verified rows (both arms +
    # any unlabeled), so the distributions + safe-failure tally reflect
    # the whole benchmark set, not just one arm.
    report.safe_failure_rate = compute_safe_failure_rate(verified)
    report.extraction_mode_distribution = compute_extraction_mode_distribution(verified)
    report.required_action_distribution = compute_required_action_distribution(verified)

    # Per-accession fetch accounting (attached by attach_fetch_stats).
    fetch_stats: dict[str, FetchStats] = {}
    for r in verified:
        if r.fetch_stats is not None:
            fetch_stats[r.accession] = r.fetch_stats
    report.fetch_stats = fetch_stats

    # count_correlation deliberately not populated — entry point.
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
                "score_3p": p.score_3p,
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


def _specificity_to_dict(sp: SpecificityReport) -> dict[str, Any]:
    return {
        "n_evaluated": sp.n_evaluated,
        "n_no_false_call": sp.n_no_false_call,
        "n_false_positive": sp.n_false_positive,
        "n_not_evaluable": sp.n_not_evaluable,
        "false_positive_accessions": sp.false_positive_accessions,
        "not_evaluable_accessions": sp.not_evaluable_accessions,
        "per_row": sorted(sp.per_row, key=lambda x: x["accession"]),
    }


def _fetch_stats_to_dict(fs: FetchStats) -> dict[str, Any]:
    return {
        "accession": fs.accession,
        "fetch_expected_runs": fs.fetch_expected_runs,
        "fetch_available_runs": fs.fetch_available_runs,
        "fetch_missing_runs": fs.fetch_missing_runs,
        "runs_with_no_fastq_url": fs.runs_with_no_fastq_url,
        "unassigned_missing_or_skipped": fs.unassigned_missing_or_skipped,
        "partial_fetch": fs.partial_fetch,
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
        # two-arm sensitivity/specificity fields.
        "recovery_denominator": report.recovery_denominator,
        "specificity": _specificity_to_dict(report.specificity),
        "adapter_control": _specificity_to_dict(report.adapter_control),
        "multi_round_sensitivity": _pair_recovery_to_dict(report.multi_round_sensitivity),
        "fetch_stats": {
            acc: _fetch_stats_to_dict(fs) for acc, fs in sorted(report.fetch_stats.items())
        },
        # entry point — present in schema for stability, populated
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

    p = argparse.ArgumentParser(description="Aggregate benchmark metrics.")
    p.add_argument("--ground-truth", required=True, type=Path)
    p.add_argument("--reports-dir", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--n-length-tolerance", type=int, default=2)
    args = p.parse_args(argv)

    rows = load_ground_truth(args.ground_truth)
    rows = attach_library_reports(rows, args.reports_dir)
    rows = attach_fetch_stats(rows, args.reports_dir)
    report = aggregate_metrics(rows, n_length_tolerance=args.n_length_tolerance)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_metrics_json(report, args.out)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
