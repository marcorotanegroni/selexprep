"""Phase 6b.5b — audit eligibility classifier for Tier 2 corpus audit.

Catalog rows go through a per-accession ENA fetch + classification BEFORE
the audit samples. The audit only samples from
``ELIGIBLE_HT_SELEX_ROUNDS``; other buckets are counted and reported in
``audit_metrics.json``'s ``catalog_classification_distribution``.

Decision tree (Phase 6b.5b user amendment — per-run first, then
per-BioProject; mixed projects exist):

1. ENA fetch fails (HTTPError / ValueError / no records)
   → ``FETCH_DEAD``
2. ALL runs use blocklisted ``library_strategy`` (defense in depth — the
   Phase 6b.5a catalog filter should already have caught these)
   → ``NON_SELEX_ASSAY``
3. Among the compatible runs: < 2 distinct round numbers parsed
   → ``NO_ROUND_STRUCTURE`` (cross-round persistence needs ≥ 2)
4. Compatible runs partition into ≥ 2 sub-trajectory groups by
   normalized ``library_name`` OR the catalog/strategy classifier
   flagged the BioProject as mixed-strategy
   → ``MIXED_PROJECT_NEEDS_GROUPING``
5. Otherwise
   → ``ELIGIBLE_HT_SELEX_ROUNDS``

Sub-trajectory detection (#4): each compatible run's ``library_name`` has
its digit substrings normalized to a literal ``{N}`` placeholder; runs
with the same normalized string are one trajectory. Example:

- ``RAPT26-1R``, ``RAPT26-2R``, ``RAPT26-3R`` all normalize to
  ``RAPT{N}-{N}R`` → 1 group → not MIXED.
- ``SELEX_round26_N40`` + ``SELEX_round2_N40`` + ``SELEX_round26_Stem2``
  + ``SELEX_round2_Stem2`` normalize to two distinct strings
  (``SELEX_round{N}_N{N}`` + ``SELEX_round{N}_Stem{N}``) → 2 groups →
  MIXED.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import sys
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

import requests

from selexprep.fetch.library_strategy import (
    classify_study_by_library_strategies,
    is_library_strategy_compatible_with_selex,
)
from selexprep.fetch.plan import FetchPlan, build_fetch_plan

logger = logging.getLogger(__name__)


class AuditEligibility(str, Enum):
    """Per-accession eligibility classification for Tier 2 audit sampling.

    ``(str, Enum)`` rather than ``StrEnum`` because StrEnum is Python
    3.11+; the project supports 3.10. The behavior is equivalent —
    instances are ``str`` subclasses and compare equal to their values.
    """

    ELIGIBLE_HT_SELEX_ROUNDS = "ELIGIBLE_HT_SELEX_ROUNDS"
    MIXED_PROJECT_NEEDS_GROUPING = "MIXED_PROJECT_NEEDS_GROUPING"
    NO_ROUND_STRUCTURE = "NO_ROUND_STRUCTURE"
    NON_SELEX_ASSAY = "NON_SELEX_ASSAY"
    FETCH_DEAD = "FETCH_DEAD"

    def __str__(self) -> str:  # match StrEnum semantics
        return self.value


@dataclass
class EligibilityReport:
    """Classification result for one catalog accession.

    Carries enough diagnostic detail to populate ``eligibility.tsv``
    (the per-accession audit-time view) and ``audit_metrics.json``'s
    ``catalog_classification_distribution`` block (the aggregate view).
    """

    accession: str
    classification: AuditEligibility
    n_runs: int = 0
    n_runs_with_round: int = 0
    distinct_rounds_parsed: list[int] = field(default_factory=list)
    library_strategies: list[str] = field(default_factory=list)
    n_runs_strategy_compatible: int = 0
    n_runs_strategy_blocklisted: int = 0
    #: Number of distinct sub-trajectory groups detected by normalizing
    #: each compatible run's ``library_name`` (digit sequences → ``{N}``).
    #: 0 when no runs; 1 for a single trajectory; ≥ 2 → MIXED.
    mixed_group_count: int = 0
    reason: str = ""


_DIGIT_NORMALIZER = re.compile(r"\d+")


def _normalize_library_name(name: str) -> str:
    """Replace digit sequences with the literal ``{N}`` placeholder.

    Used by sub-trajectory grouping (#4 in the decision tree). Empty /
    whitespace strings normalize to an empty string and collapse to a
    single group, which is the desired behavior — runs without
    library_name annotation are treated as one trajectory by default.
    """
    return _DIGIT_NORMALIZER.sub("{N}", name.strip())


def _count_distinct_library_name_groups(runs: list) -> int:
    """Count distinct sub-trajectory groups among a list of FetchRuns."""
    if not runs:
        return 0
    normalized = {_normalize_library_name(r.library_name) for r in runs}
    return max(1, len(normalized))


def classify_plan(plan: FetchPlan) -> EligibilityReport:
    """Pure classification logic from a populated :class:`FetchPlan`.

    Network-free; the caller is responsible for fetching the plan.
    Used by :func:`classify_accession` and directly in unit tests.
    """
    runs = plan.runs
    n_runs = len(runs)
    library_strategies = [r.library_strategy for r in runs]

    # Defense in depth: empty plan → FETCH_DEAD (build_fetch_plan
    # normally raises ValueError on empty ENA response; this branch
    # is reachable only if a caller constructs an empty FetchPlan
    # directly).
    if n_runs == 0:
        return EligibilityReport(
            accession=plan.accession,
            classification=AuditEligibility.FETCH_DEAD,
            reason="empty FetchPlan (no runs); ENA returned no records",
        )

    # Per-run + per-BioProject library_strategy classification (reuses
    # the Phase 6b.5a discovery filter).
    strategy_class = classify_study_by_library_strategies(plan.accession, library_strategies)

    if strategy_class.should_exclude:
        return EligibilityReport(
            accession=plan.accession,
            classification=AuditEligibility.NON_SELEX_ASSAY,
            n_runs=n_runs,
            library_strategies=library_strategies,
            n_runs_strategy_compatible=0,
            n_runs_strategy_blocklisted=n_runs,
            reason=strategy_class.exclusion_reason,
        )

    # Partition runs by strategy compatibility. Only compatible runs
    # contribute to round-structure + sub-trajectory checks.
    compatible_runs = [
        r for r in runs if is_library_strategy_compatible_with_selex(r.library_strategy)
    ]
    n_runs_compat = len(compatible_runs)
    n_runs_blocked = n_runs - n_runs_compat

    # NO_ROUND_STRUCTURE: < 2 compatible runs with parseable round, OR
    # only 1 distinct round across them. detect's cross-round
    # persistence inference needs ≥ 2 distinct rounds.
    runs_with_round = [r for r in compatible_runs if r.round_record.round_number is not None]
    distinct_rounds = sorted(
        {
            r.round_record.round_number
            for r in runs_with_round
            if r.round_record.round_number is not None
        }
    )

    if len(runs_with_round) < 2:
        return EligibilityReport(
            accession=plan.accession,
            classification=AuditEligibility.NO_ROUND_STRUCTURE,
            n_runs=n_runs,
            n_runs_with_round=len(runs_with_round),
            distinct_rounds_parsed=distinct_rounds,
            library_strategies=library_strategies,
            n_runs_strategy_compatible=n_runs_compat,
            n_runs_strategy_blocklisted=n_runs_blocked,
            reason=(
                f"only {len(runs_with_round)} compatible run(s) have a "
                "parseable round number; cross-round persistence needs "
                "≥ 2 to infer primers safely"
            ),
        )

    if len(distinct_rounds) < 2:
        return EligibilityReport(
            accession=plan.accession,
            classification=AuditEligibility.NO_ROUND_STRUCTURE,
            n_runs=n_runs,
            n_runs_with_round=len(runs_with_round),
            distinct_rounds_parsed=distinct_rounds,
            library_strategies=library_strategies,
            n_runs_strategy_compatible=n_runs_compat,
            n_runs_strategy_blocklisted=n_runs_blocked,
            reason=(
                f"only 1 distinct round parsed across {len(runs_with_round)} "
                f"runs (round={distinct_rounds[0]}); cross-round persistence "
                "needs ≥ 2 distinct rounds"
            ),
        )

    # MIXED detection: either the strategy classifier flagged mixed
    # (some compatible + some blocklisted runs) OR the compatible runs
    # partition into ≥ 2 normalized-library_name groups.
    group_count = _count_distinct_library_name_groups(compatible_runs)

    if strategy_class.is_mixed_strategy or group_count >= 2:
        reason_parts: list[str] = []
        if strategy_class.is_mixed_strategy:
            reason_parts.append(
                f"mixed library_strategy: {n_runs_compat} compatible + "
                f"{n_runs_blocked} blocklisted runs"
            )
        if group_count >= 2:
            reason_parts.append(
                f"{group_count} distinct sub-trajectories detected from "
                "normalized library_name (e.g. multiple sub-libraries x "
                "rounds — needs per-group extraction logic beyond v0.1 scope)"
            )
        return EligibilityReport(
            accession=plan.accession,
            classification=AuditEligibility.MIXED_PROJECT_NEEDS_GROUPING,
            n_runs=n_runs,
            n_runs_with_round=len(runs_with_round),
            distinct_rounds_parsed=distinct_rounds,
            library_strategies=library_strategies,
            n_runs_strategy_compatible=n_runs_compat,
            n_runs_strategy_blocklisted=n_runs_blocked,
            mixed_group_count=group_count,
            reason="; ".join(reason_parts),
        )

    return EligibilityReport(
        accession=plan.accession,
        classification=AuditEligibility.ELIGIBLE_HT_SELEX_ROUNDS,
        n_runs=n_runs,
        n_runs_with_round=len(runs_with_round),
        distinct_rounds_parsed=distinct_rounds,
        library_strategies=library_strategies,
        n_runs_strategy_compatible=n_runs_compat,
        n_runs_strategy_blocklisted=n_runs_blocked,
        mixed_group_count=group_count,
        reason=(
            f"{len(distinct_rounds)} distinct rounds across "
            f"{len(runs_with_round)} runs, single trajectory, all "
            "library_strategy values SELEX-compatible"
        ),
    )


def classify_accession(accession: str, *, timeout_s: int = 30) -> EligibilityReport:
    """Fetch the accession's :class:`FetchPlan` from ENA and classify it.

    Catches the network/parse errors that :func:`build_fetch_plan` can
    raise (``requests.RequestException`` / ``ValueError``) and maps
    them to ``FETCH_DEAD`` with the underlying error message recorded
    in ``reason``. This is the lazy-evaluation aspect the user's
    eligibility plan calls for — fetch failures aren't decided
    pre-classification; they surface as the FETCH_DEAD bucket.
    """
    try:
        plan = build_fetch_plan(accession, timeout_s=timeout_s)
    except (requests.RequestException, ValueError) as e:
        return EligibilityReport(
            accession=accession,
            classification=AuditEligibility.FETCH_DEAD,
            reason=f"ENA fetch failed: {type(e).__name__}: {e}",
        )
    return classify_plan(plan)


# ---------------------------------------------------------------------------
# TSV writer (per-accession view consumed by the Snakefile)
# ---------------------------------------------------------------------------


_TSV_COLUMNS = (
    "accession",
    "classification",
    "n_runs",
    "n_runs_with_round",
    "distinct_rounds_parsed",
    "library_strategies",
    "n_runs_strategy_compatible",
    "n_runs_strategy_blocklisted",
    "mixed_group_count",
    "reason",
)


def write_eligibility_tsv(reports: list[EligibilityReport], out_path: Path) -> None:
    """Emit per-accession eligibility classifications to a TSV.

    Sorted by accession for deterministic diffs. ``distinct_rounds_parsed``
    and ``library_strategies`` are JSON-encoded lists so downstream
    consumers can parse them losslessly.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sorted_reports = sorted(reports, key=lambda r: r.accession)

    with out_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(_TSV_COLUMNS), delimiter="\t")
        writer.writeheader()
        for r in sorted_reports:
            writer.writerow(
                {
                    "accession": r.accession,
                    "classification": str(r.classification),
                    "n_runs": r.n_runs,
                    "n_runs_with_round": r.n_runs_with_round,
                    "distinct_rounds_parsed": json.dumps(r.distinct_rounds_parsed),
                    "library_strategies": json.dumps(sorted(set(r.library_strategies))),
                    "n_runs_strategy_compatible": r.n_runs_strategy_compatible,
                    "n_runs_strategy_blocklisted": r.n_runs_strategy_blocklisted,
                    "mixed_group_count": r.mixed_group_count,
                    "reason": r.reason.replace("\t", " ").replace("\n", " "),
                }
            )


def read_eligibility_tsv(path: Path) -> list[EligibilityReport]:
    """Inverse of :func:`write_eligibility_tsv` — round-trip for tests."""
    out: list[EligibilityReport] = []
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            out.append(
                EligibilityReport(
                    accession=row["accession"],
                    classification=AuditEligibility(row["classification"]),
                    n_runs=int(row["n_runs"]),
                    n_runs_with_round=int(row["n_runs_with_round"]),
                    distinct_rounds_parsed=json.loads(row.get("distinct_rounds_parsed") or "[]"),
                    library_strategies=json.loads(row.get("library_strategies") or "[]"),
                    n_runs_strategy_compatible=int(row["n_runs_strategy_compatible"]),
                    n_runs_strategy_blocklisted=int(row["n_runs_strategy_blocklisted"]),
                    mixed_group_count=int(row["mixed_group_count"]),
                    reason=row.get("reason", ""),
                )
            )
    return out


def eligible_accessions(reports: list[EligibilityReport]) -> list[str]:
    """Return only the accessions classified as ELIGIBLE_HT_SELEX_ROUNDS.

    Sorted for determinism. Consumed by the Tier 2 audit sampler so
    only fully-eligible rows enter ``sample_corpus``.
    """
    return sorted(
        r.accession
        for r in reports
        if r.classification == AuditEligibility.ELIGIBLE_HT_SELEX_ROUNDS
    )


def classification_distribution(
    reports: list[EligibilityReport],
) -> dict[str, int]:
    """Aggregate distribution {AuditEligibility value: count}.

    Used by the audit aggregator to populate
    ``catalog_classification_distribution`` in ``audit_metrics.json``.
    Sorted alphabetically for stable JSON shape.
    """
    counts: dict[str, int] = {}
    for r in reports:
        key = str(r.classification)
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


# ---------------------------------------------------------------------------
# CLI: classify-catalog subcommand
# ---------------------------------------------------------------------------


def _cmd_classify_catalog(args: argparse.Namespace) -> int:
    import pandas as pd  # local import — pandas is core dep, kept lazy

    catalog_df = pd.read_csv(args.catalog)
    insdc_only = catalog_df[
        catalog_df["bioproject_id"].str.match(r"^(PRJ[NDE][ABM]|SRP|ERP|DRP)", na=False)
    ]
    accessions = sorted(set(insdc_only["bioproject_id"].tolist()))
    if args.limit is not None and args.limit > 0:
        accessions = accessions[: args.limit]

    print(f"Classifying {len(accessions)} INSDC accessions...", file=sys.stderr)
    reports: list[EligibilityReport] = []
    for i, acc in enumerate(accessions):
        report = classify_accession(acc, timeout_s=args.timeout_s)
        reports.append(report)
        if (i + 1) % 25 == 0:
            print(f"  classified {i + 1}/{len(accessions)}", file=sys.stderr)

    write_eligibility_tsv(reports, args.out)
    dist = classification_distribution(reports)
    print(f"wrote {args.out}", file=sys.stderr)
    for k, n in dist.items():
        print(f"  {k}: {n}", file=sys.stderr)
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Tier 2 audit eligibility classifier (Phase 6b.5b).")
    sub = p.add_subparsers(dest="cmd", required=True)

    classify = sub.add_parser(
        "classify-catalog",
        help=(
            "Read the bundled bioprojects.csv, classify each INSDC "
            "accession via an ENA fetch, write a per-accession TSV."
        ),
    )
    classify.add_argument(
        "--catalog",
        type=Path,
        required=True,
        help="Path to bioprojects.csv (the catalog snapshot to classify).",
    )
    classify.add_argument("--out", type=Path, required=True, help="Output eligibility.tsv path.")
    classify.add_argument(
        "--timeout-s",
        type=int,
        default=30,
        help="HTTP timeout per ENA query (seconds).",
    )
    classify.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional cap on accessions to classify (for smoke tests).",
    )
    classify.set_defaults(func=_cmd_classify_catalog)

    args = p.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())


__all__ = [
    "AuditEligibility",
    "EligibilityReport",
    "classification_distribution",
    "classify_accession",
    "classify_plan",
    "eligible_accessions",
    "main",
    "read_eligibility_tsv",
    "write_eligibility_tsv",
]
