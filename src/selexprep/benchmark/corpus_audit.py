"""Tier 2 corpus-audit scaffolding (Phase 6b.3a — pipeline only, no results).

The audit answers a different question than Tier 1's curated primer-recovery
benchmark (``benchmark.metrics`` + Figure A): given the messy public
HT-SELEX corpus surfaced by selexprep's bundled discovery catalog, what
fraction of accessions can be fetched at all, what fraction reach
HIGH-confidence primer inference, what fraction safely refuse, and how
do ``extraction_mode`` / ``required_action`` distribute?

These are **distributional** metrics — there is no per-row ground truth.
The methods text says so explicitly. Comparator tools (AptaPLEX,
EasyDIVER+) cannot produce this kind of audit by construction: they need
the primer sequences as input, so they cannot be pointed at a random
catalog accession and asked "what does the LibraryReport say?".

Scope (locked plan + Phase 6b.3a):

- ``sample_corpus(n, sources=..., exclude=..., seed=...)`` — deterministic
  uniform-random sample from the bundled catalog, INSDC-only by
  construction (``filter_catalog(insdc_only=True)``).
- ``write_accessions_tsv`` — emits the 2-column TSV that
  ``selexprep run`` consumes.
- ``aggregate_audit_from_run_outputs`` — parses ``run_summary.tsv``
  (Phase 6a batch driver) into a :class:`CorpusAuditReport` with
  EXPLICIT denominators per metric family. The methodological
  correction (Codex peer-review + user) is that
  ``inference_safe_failure_rate`` is computed only among rows with a
  LibraryReport — NEVER mixed with fetch failures, which would inflate
  the metric with ENA/network problems.
- ``write_audit_json`` — deterministic sorted-keys JSON, same
  discipline as ``benchmark.metrics.write_metrics_json``.
- ``main(argv)`` — CLI with ``sample`` and ``aggregate`` subcommands so
  the Snakefile ``rule sample_corpus`` and ``rule aggregate_audit`` can
  invoke this module directly.

6b.3a ships **scaffolding only**. The actual HPC run that produces
``audit_metrics.json`` + ``audit_accessions.tsv`` + ``figure_b.{pdf,png}``
is the 6b.4 follow-up commit (no code changes). CI does NOT execute the
Snakefile.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from selexprep.catalog.filter import filter_catalog
from selexprep.catalog.reader import catalog_version as _catalog_version
from selexprep.catalog.reader import load_catalog

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# RunStatus partition for fetchability accounting
# ---------------------------------------------------------------------------

# Mirrors ``selexprep.run.runner.RunStatus``. Kept local so the audit module
# doesn't pull the runner's heavy dependency chain just to compute a few
# distributions. Fetch is considered to have succeeded for any row whose
# RunStatus indicates the pipeline got past the fetch stage — which means
# DETECT_FAILED + EXTRACT_REFUSED + EXTRACT_FAILED + COUNT_FAILED +
# QC_FAILED + SKIPPED_READ_MERGING_RECOMMENDED + OK all count as fetchable.
_FETCH_SUCCEEDED_STATUSES: frozenset[str] = frozenset(
    {
        "OK",
        "SKIPPED_READ_MERGING_RECOMMENDED",
        "DETECT_FAILED",
        "EXTRACT_REFUSED",
        "EXTRACT_FAILED",
        "COUNT_FAILED",
        "QC_FAILED",
    }
)


# ---------------------------------------------------------------------------
# Dataclass
# ---------------------------------------------------------------------------


@dataclass
class CorpusAuditReport:
    """Tier 2 audit report — distributional metrics with explicit denominators.

    The denominator partition is the methodological correction from the
    Codex peer-review + user pass: a single rate over all metrics would
    inflate ``inference_safe_failure_rate`` with fetch failures (network /
    ENA / regional restrictions are not safe failures, they're operational
    failures). Each metric family is therefore bucketed under the
    denominator it answers a question about:

    - **fetchability** (``n_sampled``): can public data be obtained?
    - **inference**   (``n_with_library_report``): given selexprep ran
      ``detect``, what did the report say?
    - **qc**          (``n_with_qc_run``): for rows that reached qc, how
      many flags were raised?

    Each panel in Figure B labels its denominator in the subtitle so a
    reviewer never has to guess the normalization.
    """

    # Reproducibility envelope.
    catalog_version: str | None = None
    sample_seed: int = 42
    sample_accessions_sha256: str = ""
    n_sampled: int = 0
    n_in_ground_truth_overlap: int = 0

    # Fetchability metrics (denominator = ``n_sampled``).
    fetch_outcome_distribution: dict[str, int] = field(default_factory=dict)
    n_fetchable: int = 0

    # Inference metrics (denominator = ``n_with_library_report``).
    n_with_library_report: int = 0
    library_report_status_distribution: dict[str, int] = field(default_factory=dict)
    extraction_mode_distribution: dict[str, int] = field(default_factory=dict)
    required_action_distribution: dict[str, int] = field(default_factory=dict)
    inference_safe_failure_rate: float = 0.0
    n_inference_safe_failures: int = 0
    inference_safe_failure_by_reason: dict[str, int] = field(default_factory=dict)

    # QC metrics (denominator = ``n_with_qc_run``).
    n_with_qc_run: int = 0
    flags_raised_histogram: dict[int, int] = field(default_factory=dict)

    # Phase 6b.5b: catalog-level eligibility breakdown (denominator =
    # the full classified catalog, NOT n_sampled). Populated when an
    # eligibility TSV is passed to the aggregator. Provides the
    # "audit-eligible vs not" layer-1 story alongside the layer-2
    # in-sample distributions above.
    catalog_classification_distribution: dict[str, int] = field(default_factory=dict)
    n_catalog_classified: int = 0
    n_catalog_eligible: int = 0

    # Phase 6b.5d: full-catalog denominator (INSDC-classified + passthrough
    # rows). The eligibility classifier only sees INSDC rows because
    # figshare/zenodo passthrough deposits don't have ENA filereport
    # endpoints or per-run library_strategy metadata. Surfacing the total
    # here keeps Figure B's title honest — without it, a reviewer reads
    # "X of N audit-eligible" as "X of all catalog rows" rather than
    # "X of the N rows that the classifier can act on".
    n_catalog_total: int = 0
    n_catalog_non_insdc_passthrough: int = 0

    # Phase 6b.5d: free-form caveats block surfaced verbatim in the audit
    # JSON. Currently carries the multiplex-detection caveat (NO_ROUND_STRUCTURE
    # may include single-FASTQ inline-barcoded SELEX deposits that v0.1
    # cannot detect without a user-supplied sample sheet). The aggregator
    # populates this; consumers (Figure B, Application Note) read it.
    caveats: dict[str, str] = field(default_factory=dict)

    # Traceability — the parsed per-accession rows with ground-truth
    # overlap annotation. Stable sort by accession.
    per_accession: list[dict[str, Any]] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------


def sample_corpus(
    n: int,
    *,
    sources: str | None = None,
    exclude: tuple[str, ...] = (),
    seed: int = 42,
    eligible_only: tuple[str, ...] | None = None,
) -> list[str]:
    """Sample ``n`` INSDC accessions from the bundled catalog deterministically.

    Parameters
    ----------
    n
        Sample size. If fewer INSDC accessions remain after filtering,
        the entire remaining pool is returned (no padding).
    sources
        Substring (case-insensitive) matched against the catalog's
        ``source`` column. ``None`` keeps all INSDC sources — ENA, SRA,
        DDBJ. Useful for restricting the audit to a single discovery
        adapter.
    exclude
        Accessions to drop before sampling. Used by the Snakefile to
        pass ``ground_truth.tsv``'s accession list so Tier 1 rows don't
        double-count in Tier 2. Also useful for chaining multiple
        sweeps without repeats.
    seed
        Seed for ``random.Random``. The same ``(catalog, n, sources,
        exclude, seed)`` quadruple always returns the same accession
        list (bit-deterministic).
    eligible_only
        Phase 6b.5b: restrict the sampling pool to this set of
        accessions (intersected with the INSDC catalog). The audit
        Snakefile passes the ``ELIGIBLE_HT_SELEX_ROUNDS`` rows from
        the eligibility TSV here so non-eligible rows (NON_SELEX_ASSAY,
        NO_ROUND_STRUCTURE, MIXED, FETCH_DEAD) can't be drawn. ``None``
        means no restriction (backward compat with pre-6b.5b callers).

    Returns
    -------
    A sorted list of accession identifiers. Sorted output (rather than
    the random.sample call order) is the bit-deterministic contract —
    downstream consumers should not rely on draw order.

    Notes
    -----
    The catalog version is NOT pinned by this function; callers that
    need cross-snapshot reproducibility should record
    :func:`selexprep.catalog.reader.catalog_version` and the returned
    sample's sha256 (see :func:`aggregate_audit_from_run_outputs`).
    """
    df = load_catalog()
    df = filter_catalog(df, insdc_only=True, source_contains=sources)
    excluded = set(exclude)
    eligible_set = set(eligible_only) if eligible_only is not None else None
    pool = sorted(
        a
        for a in df["bioproject_id"].tolist()
        if a and a not in excluded and (eligible_set is None or a in eligible_set)
    )
    if not pool:
        return []
    take = min(n, len(pool))
    rng = random.Random(seed)
    sampled = rng.sample(pool, take)
    return sorted(sampled)


def write_accessions_tsv(accessions: list[str], path: Path, *, notes: str = "") -> None:
    """Emit the 2-column ``accession\\tnotes`` TSV that ``selexprep run`` consumes.

    The TSV is sorted by accession (stable across runs).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("accession\tnotes\n")
        for acc in sorted(accessions):
            fh.write(f"{acc}\t{notes}\n")


def accessions_sha256(accessions: list[str]) -> str:
    """sha256 of the sorted accession list, joined with ``\\n`` and trailing newline.

    Used to fingerprint a sampled audit population so the audit JSON +
    figure are reproducible across future catalog refreshes (which
    shift row indices and break seed-based reproducibility on the
    catalog DataFrame alone). Records the fingerprint of the resulting
    accession list itself.
    """
    payload = "\n".join(sorted(accessions)) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _read_ground_truth_accessions(path: Path) -> list[str]:
    """Read the ``accession`` column from ``ground_truth.tsv`` (Tier 1).

    Used by the Snakefile to compute ``--exclude-ground-truth`` and by
    the aggregator to annotate ``is_in_ground_truth`` per row.
    """
    if not path.exists():
        return []
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return sorted([a for a in df["accession"].tolist() if a])


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------


_MULTIPLEX_CAVEAT_NO_ROUND_STRUCTURE = (
    "The NO_ROUND_STRUCTURE bucket may include single-FASTQ inline-barcoded "
    "multiplexed SELEX deposits that selexprep v0.1 cannot detect without a "
    "user-supplied sample sheet (multiplex auto-detection is a v0.2 "
    "deferral per the locked plan). Visible v0.1 signals are: few runs + "
    "library_strategy=SELEX + 0 parseable rounds. Single-run deposits with "
    "n_runs_with_round=1 are genuinely single-round and correctly classified."
)


def aggregate_audit_from_run_outputs(
    run_summary_tsv: Path,
    ground_truth_tsv: Path | None,
    *,
    catalog_version: str | None,
    sample_seed: int,
    sample_accessions_sha256: str,
    eligibility_tsv: Path | None = None,
    catalog_csv: Path | None = None,
) -> CorpusAuditReport:
    """Parse ``run_summary.tsv`` (Phase 6a batch driver) into the audit report.

    The methodological correction folded in here (Codex peer-review +
    user pass): ``inference_safe_failure_rate`` is computed ONLY among
    rows whose ``library_report_status`` is non-empty — i.e., detect
    produced a report at all. Mixing fetch failures into this rate
    would inflate the metric with ENA/network problems and conflate
    "selexprep refused" (a feature) with "the dataset was
    unreachable" (an external problem).

    Parameters
    ----------
    run_summary_tsv
        The ``run_summary.tsv`` emitted by ``selexprep run``. Required
        columns: ``accession``, ``status``, ``library_report_status``,
        ``extraction_mode``, ``required_action``, ``flags_raised``.
    ground_truth_tsv
        Optional Tier 1 ground-truth TSV. When present, each
        per-accession row is annotated with ``is_in_ground_truth: bool``
        and ``n_in_ground_truth_overlap`` is filled in. ``None`` is
        accepted so the aggregator can run without a Tier 1 set
        available.
    catalog_version, sample_seed, sample_accessions_sha256
        Reproducibility envelope values stamped by the Snakefile
        ``rule aggregate_audit`` from the values the sampler recorded.
        Required (not defaulted) so callers cannot forget to thread
        them through.
    eligibility_tsv
        Phase 6b.5b: optional path to the eligibility TSV emitted by
        ``selexprep.benchmark.eligibility.classify-catalog``. When
        present, the audit JSON gains
        ``catalog_classification_distribution`` (counts per
        AuditEligibility value across the full classified catalog) +
        ``n_catalog_classified`` + ``n_catalog_eligible`` for the
        layer-1 "what fraction of the catalog is audit-eligible" story.
    catalog_csv
        Phase 6b.5d: optional path to ``bioprojects.csv``. When given,
        the audit JSON gains ``n_catalog_total`` +
        ``n_catalog_non_insdc_passthrough`` so a reviewer can read
        Figure B's "X of N audit-eligible" segment against the full
        catalog denominator rather than the implicit INSDC-only one.
        ``None`` leaves both fields at 0 (backward compat).
    """
    df = pd.read_csv(run_summary_tsv, sep="\t", dtype=str).fillna("")
    report = CorpusAuditReport(
        catalog_version=catalog_version,
        sample_seed=sample_seed,
        sample_accessions_sha256=sample_accessions_sha256,
    )

    # Phase 6b.5d: populate full-catalog denominators (when a catalog
    # path is supplied) + the multiplex caveat under NO_ROUND_STRUCTURE.
    # Both feed Figure B's title + the paper's "what the public corpus
    # looks like" narrative.
    if catalog_csv is not None and catalog_csv.exists():
        from selexprep.catalog.filter import is_insdc_accession

        catalog_df = pd.read_csv(catalog_csv)
        ids = catalog_df["bioproject_id"].fillna("").astype(str).tolist()
        report.n_catalog_total = len(ids)
        report.n_catalog_non_insdc_passthrough = sum(1 for a in ids if not is_insdc_accession(a))

    report.caveats["NO_ROUND_STRUCTURE"] = _MULTIPLEX_CAVEAT_NO_ROUND_STRUCTURE

    if eligibility_tsv is not None and eligibility_tsv.exists():
        # Local import to avoid circular: eligibility imports from
        # fetch.plan which is fine, but corpus_audit is imported by
        # the Snakefile too. Lazy keeps module-level imports tight.
        from selexprep.benchmark.eligibility import (
            classification_distribution,
            read_eligibility_tsv,
        )

        elig_reports = read_eligibility_tsv(eligibility_tsv)
        report.catalog_classification_distribution = classification_distribution(elig_reports)
        report.n_catalog_classified = len(elig_reports)
        report.n_catalog_eligible = report.catalog_classification_distribution.get(
            "ELIGIBLE_HT_SELEX_ROUNDS", 0
        )

    gt_set = set(_read_ground_truth_accessions(ground_truth_tsv)) if ground_truth_tsv else set()
    report.n_sampled = len(df)

    per_accession: list[dict[str, Any]] = []
    for _, row in df.iterrows():
        status = row.get("status", "")
        lr_status = row.get("library_report_status", "")
        extraction_mode = row.get("extraction_mode", "")
        required_action = row.get("required_action", "")
        flags_raised_raw = row.get("flags_raised", "")
        accession = row.get("accession", "")

        # Fetch-outcome distribution: every status counted.
        report.fetch_outcome_distribution[status] = (
            report.fetch_outcome_distribution.get(status, 0) + 1
        )
        if status in _FETCH_SUCCEEDED_STATUSES:
            report.n_fetchable += 1

        # Inference metrics: only rows with a LibraryReport (denominator).
        has_library_report = bool(lr_status)
        if has_library_report:
            report.n_with_library_report += 1
            report.library_report_status_distribution[lr_status] = (
                report.library_report_status_distribution.get(lr_status, 0) + 1
            )
            if extraction_mode:
                report.extraction_mode_distribution[extraction_mode] = (
                    report.extraction_mode_distribution.get(extraction_mode, 0) + 1
                )
            if required_action:
                report.required_action_distribution[required_action] = (
                    report.required_action_distribution.get(required_action, 0) + 1
                )
            # Safe-failure trigger list (a row can satisfy multiple
            # reasons; n_inference_safe_failures counts the row once,
            # by_reason tallies each reason independently).
            triggered: list[str] = []
            if lr_status == "UNABLE_TO_INFER":
                triggered.append("status_UNABLE_TO_INFER")
            if extraction_mode == "UNABLE_TO_EXTRACT":
                triggered.append("extraction_mode_UNABLE_TO_EXTRACT")
            if required_action == "MANUAL_PRIMERS_REQUIRED":
                triggered.append("required_action_MANUAL_PRIMERS_REQUIRED")
            if triggered:
                report.n_inference_safe_failures += 1
                for reason in triggered:
                    report.inference_safe_failure_by_reason[reason] = (
                        report.inference_safe_failure_by_reason.get(reason, 0) + 1
                    )

        # QC metrics: only rows that reached qc (status == "OK" since
        # qc is the last stage). flags_raised is a stringified int (or
        # empty if no qc).
        if status == "OK" and flags_raised_raw != "":
            try:
                flags_raised = int(flags_raised_raw)
            except ValueError:
                logger.warning(
                    "audit: non-integer flags_raised %r for accession %s — skipping",
                    flags_raised_raw,
                    accession,
                )
            else:
                report.n_with_qc_run += 1
                report.flags_raised_histogram[flags_raised] = (
                    report.flags_raised_histogram.get(flags_raised, 0) + 1
                )

        per_accession.append(
            {
                "accession": accession,
                "status": status,
                "library_report_status": lr_status,
                "extraction_mode": extraction_mode,
                "required_action": required_action,
                "flags_raised": flags_raised_raw,
                "is_in_ground_truth": accession in gt_set,
            }
        )

    report.per_accession = sorted(per_accession, key=lambda r: r["accession"])
    report.n_in_ground_truth_overlap = sum(
        1 for r in report.per_accession if r["is_in_ground_truth"]
    )
    report.inference_safe_failure_rate = (
        report.n_inference_safe_failures / report.n_with_library_report
        if report.n_with_library_report > 0
        else 0.0
    )
    return report


# ---------------------------------------------------------------------------
# JSON writer
# ---------------------------------------------------------------------------


def write_audit_json(report: CorpusAuditReport, path: Path) -> None:
    """Serialize the audit report to deterministic JSON (sorted keys + stable ordering).

    The JSON is the source of truth for Figure B (matplotlib byte-output is
    not guaranteed deterministic; the audit JSON is, by the same discipline
    used in ``benchmark.metrics.write_metrics_json``).

    Defensive sort (Codex peer-review): ``per_accession`` is sorted by
    accession here even if the caller already sorted upstream, so direct
    callers (tests, future hand-rolled scripts) can't accidentally write a
    non-deterministic ordering.
    """
    per_accession_sorted = sorted(report.per_accession, key=lambda r: r.get("accession", ""))
    payload: dict[str, Any] = {
        "catalog_version": report.catalog_version,
        "sample_seed": report.sample_seed,
        "sample_accessions_sha256": report.sample_accessions_sha256,
        "n_sampled": report.n_sampled,
        "n_in_ground_truth_overlap": report.n_in_ground_truth_overlap,
        # Phase 6b.5d: full-catalog denominator (zero unless the
        # aggregator was given --catalog). Lets Figure B's title quote
        # the honest "X of N INSDC rows · M non-INSDC passthrough · K
        # catalog total" breakdown.
        "n_catalog_total": report.n_catalog_total,
        "n_catalog_non_insdc_passthrough": report.n_catalog_non_insdc_passthrough,
        # Phase 6b.5b layer-1 view: catalog-level eligibility breakdown.
        # Empty when the aggregator wasn't given an eligibility TSV
        # (backward compat with pre-6b.5b audit runs).
        "n_catalog_classified": report.n_catalog_classified,
        "n_catalog_eligible": report.n_catalog_eligible,
        "catalog_classification_distribution": dict(
            sorted(report.catalog_classification_distribution.items())
        ),
        # Phase 6b.5d: free-form caveats keyed by bucket / topic.
        "caveats": dict(sorted(report.caveats.items())),
        "fetch_outcome_distribution": dict(sorted(report.fetch_outcome_distribution.items())),
        "n_fetchable": report.n_fetchable,
        "n_with_library_report": report.n_with_library_report,
        "library_report_status_distribution": dict(
            sorted(report.library_report_status_distribution.items())
        ),
        "extraction_mode_distribution": dict(sorted(report.extraction_mode_distribution.items())),
        "required_action_distribution": dict(sorted(report.required_action_distribution.items())),
        "inference_safe_failure_rate": round(report.inference_safe_failure_rate, 6),
        "n_inference_safe_failures": report.n_inference_safe_failures,
        "inference_safe_failure_by_reason": dict(
            sorted(report.inference_safe_failure_by_reason.items())
        ),
        "n_with_qc_run": report.n_with_qc_run,
        # Histogram keys are ints; JSON keys are strings — stringify
        # explicitly with a numeric sort so the on-disk shape is stable.
        "flags_raised_histogram": {
            str(k): report.flags_raised_histogram[k]
            for k in sorted(report.flags_raised_histogram.keys())
        },
        "per_accession": per_accession_sorted,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    path.write_text(text, encoding="utf-8")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _cmd_sample(args: argparse.Namespace) -> int:
    exclude_tuple: tuple[str, ...] = ()
    if args.exclude_ground_truth is not None:
        exclude_tuple = tuple(_read_ground_truth_accessions(args.exclude_ground_truth))
    # Phase 6b.5b: when an eligibility TSV is provided, restrict the
    # sampling pool to ``ELIGIBLE_HT_SELEX_ROUNDS`` accessions only.
    eligible_tuple: tuple[str, ...] | None = None
    if args.eligibility is not None:
        from selexprep.benchmark.eligibility import (
            eligible_accessions,
            read_eligibility_tsv,
        )

        elig_reports = read_eligibility_tsv(args.eligibility)
        eligible_tuple = tuple(eligible_accessions(elig_reports))
    accessions = sample_corpus(
        n=args.n,
        sources=args.sources,
        exclude=exclude_tuple,
        seed=args.seed,
        eligible_only=eligible_tuple,
    )
    write_accessions_tsv(accessions, args.out)
    # Also emit a sidecar manifest with the reproducibility envelope so the
    # aggregator's --catalog-version / --sample-seed / --sample-sha can
    # source from a file rather than depending on the runner passing them
    # through the shell. The Snakefile reads this sidecar when invoking
    # ``aggregate``.
    sidecar = args.out.with_suffix(".manifest.json")
    payload = {
        "catalog_version": _catalog_version(),
        "sample_seed": args.seed,
        "sample_accessions_sha256": accessions_sha256(accessions),
        "n_sampled": len(accessions),
    }
    sidecar.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.out} ({len(accessions)} accessions)")
    print(f"wrote {sidecar}")
    return 0


def _cmd_aggregate(args: argparse.Namespace) -> int:
    # Reproducibility envelope: prefer the sidecar emitted by ``sample``
    # but allow CLI overrides (useful for re-aggregating from a
    # pre-existing run_summary.tsv when the sidecar isn't on disk).
    catalog_version_value: str | None = args.catalog_version
    sample_seed_value: int = args.sample_seed
    sample_sha_value: str = args.sample_sha
    sidecar_n_sampled: int | None = None
    if args.sample_manifest is not None and args.sample_manifest.exists():
        manifest = json.loads(args.sample_manifest.read_text(encoding="utf-8"))
        catalog_version_value = catalog_version_value or manifest.get("catalog_version")
        if args.sample_seed == 42 and "sample_seed" in manifest:
            sample_seed_value = int(manifest["sample_seed"])
        if not sample_sha_value:
            sample_sha_value = manifest.get("sample_accessions_sha256", "")
        if "n_sampled" in manifest:
            sidecar_n_sampled = int(manifest["n_sampled"])

    # Codex peer-review fix: the function signature requires
    # ``sample_accessions_sha256`` (it's a positional kwarg, not a default)
    # so the audit JSON's reproducibility envelope cannot silently land
    # empty. We enforce the same at the CLI boundary: if neither the
    # sidecar nor an explicit ``--sample-sha`` provides a value, refuse to
    # emit a JSON with empty provenance. That's the whole point of the
    # envelope — a paper-defensible audit needs to fingerprint its sample.
    if not sample_sha_value:
        raise SystemExit(
            "audit aggregate: no sample sha256 available. Pass either "
            "--sample-manifest pointing at the sidecar emitted by "
            "'corpus_audit sample', or --sample-sha <hex> explicitly. "
            "Refusing to emit an audit JSON with empty provenance."
        )

    report = aggregate_audit_from_run_outputs(
        run_summary_tsv=args.run_summary,
        ground_truth_tsv=args.ground_truth,
        catalog_version=catalog_version_value,
        sample_seed=sample_seed_value,
        sample_accessions_sha256=sample_sha_value,
        eligibility_tsv=args.eligibility,
        catalog_csv=args.catalog,
    )

    # Codex peer-review fix: the run summary's row count should match the
    # sidecar's ``n_sampled``. If selexprep run lost a row (a runner bug)
    # or someone re-ran with a different TSV (operator error), warn loudly
    # — the audit JSON's ``n_sampled`` represents what was actually
    # processed (len(run_summary)), but the discrepancy needs to be
    # visible. Don't fail: the per-accession breakdown is still usable.
    if sidecar_n_sampled is not None and sidecar_n_sampled != report.n_sampled:
        logger.warning(
            "audit aggregate: run_summary has %d rows but sidecar manifest "
            "says n_sampled=%d. Either selexprep run lost rows or the run "
            "summary was generated from a different accessions TSV.",
            report.n_sampled,
            sidecar_n_sampled,
        )
        print(
            f"WARNING: run_summary rows ({report.n_sampled}) != sidecar "
            f"n_sampled ({sidecar_n_sampled}); audit JSON reports "
            f"{report.n_sampled} (rows actually processed)."
        )

    write_audit_json(report, args.out)
    print(f"wrote {args.out}")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Tier 2 corpus-audit pipeline (Phase 6b.3a).")
    sub = p.add_subparsers(dest="cmd", required=True)

    sample = sub.add_parser("sample", help="Sample N INSDC accessions from the catalog.")
    sample.add_argument("--n", type=int, default=30)
    sample.add_argument("--seed", type=int, default=42)
    sample.add_argument(
        "--sources",
        type=str,
        default=None,
        help="Substring match on catalog 'source' column (case-insensitive).",
    )
    sample.add_argument(
        "--exclude-ground-truth",
        type=Path,
        default=None,
        help="Path to ground_truth.tsv; its accessions are excluded from the sample.",
    )
    sample.add_argument(
        "--eligibility",
        type=Path,
        default=None,
        help=(
            "Phase 6b.5b: path to eligibility.tsv emitted by "
            "``selexprep.benchmark.eligibility classify-catalog``. "
            "When provided, only accessions classified as "
            "ELIGIBLE_HT_SELEX_ROUNDS are eligible for sampling."
        ),
    )
    sample.add_argument("--out", type=Path, required=True)
    sample.set_defaults(func=_cmd_sample)

    aggregate = sub.add_parser(
        "aggregate",
        help="Aggregate run_summary.tsv into a Tier 2 audit_metrics.json.",
    )
    aggregate.add_argument("--run-summary", type=Path, required=True)
    aggregate.add_argument("--ground-truth", type=Path, default=None)
    aggregate.add_argument("--out", type=Path, required=True)
    aggregate.add_argument(
        "--sample-manifest",
        type=Path,
        default=None,
        help="Sidecar manifest emitted by 'sample' subcommand.",
    )
    aggregate.add_argument("--catalog-version", type=str, default=None)
    aggregate.add_argument("--sample-seed", type=int, default=42)
    aggregate.add_argument("--sample-sha", type=str, default="")
    aggregate.add_argument(
        "--eligibility",
        type=Path,
        default=None,
        help=(
            "Phase 6b.5b: path to eligibility.tsv. When provided, the "
            "audit JSON gains ``catalog_classification_distribution`` "
            "(layer-1 view) alongside the in-sample metrics."
        ),
    )
    aggregate.add_argument(
        "--catalog",
        type=Path,
        default=None,
        help=(
            "Phase 6b.5d: path to bioprojects.csv. When provided, the "
            "audit JSON gains ``n_catalog_total`` + "
            "``n_catalog_non_insdc_passthrough`` so Figure B's title can "
            "report the honest full-catalog denominator."
        ),
    )
    aggregate.set_defaults(func=_cmd_aggregate)

    args = p.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())


__all__ = [
    "CorpusAuditReport",
    "accessions_sha256",
    "aggregate_audit_from_run_outputs",
    "main",
    "sample_corpus",
    "write_accessions_tsv",
    "write_audit_json",
]
