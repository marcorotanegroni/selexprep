"""Per-BioProject round-coverage classification.

Classifies whether a BioProject exposes all SELEX rounds as separate public
runs, only a subset, or bundled/multiplexed runs that cannot be cleanly
separated. The classification feeds downstream filtering (a multiplexed-
unrecoverable project can't be used for per-round counts) and the QC report.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Callable
from datetime import datetime

STATUS_ALL = "all_rounds_public"
STATUS_PARTIAL = "partial_rounds_public"
STATUS_MULTIPLEXED = "rounds_multiplexed_unrecoverable"
STATUS_UNKNOWN = "unknown_round_coverage"


def _parse_int(value: str | None) -> int | None:
    raw = (value or "").strip()
    if raw.lstrip("-").isdigit():
        return int(raw)
    return None


def _has_multiplex_hint(*texts: str | None) -> bool:
    text = " ".join((t or "").strip().lower() for t in texts if t)
    if not text:
        return False

    hints = (
        "multiplex",
        "contains rounds",
        "contained rounds",
        "multiple rounds in one",
        "single run contains",
        "single srr contains",
        "bundled rounds",
        "combined rounds",
        "pooled rounds",
    )
    return any(h in text for h in hints)


def classify_bioproject_round_coverage(
    bioproject: dict,
    samples: list[dict],
    rounds_by_srr: dict[str, dict],
) -> dict:
    """Classify round-coverage status for a single BioProject.

    Returns a dict with `round_coverage_status` (one of `STATUS_*`),
    `reason`, observed-vs-declared round counts, multiplex-hint flag, and
    per-confidence histogram. Routing logic (in order):

    1. No public SRRs → `STATUS_UNKNOWN` ("no_public_srrs").
    2. No declared round count + multiplex hint → `STATUS_MULTIPLEXED`.
    3. No declared round count → `STATUS_UNKNOWN`.
    4. Observed unique rounds ≥ declared → `STATUS_ALL`.
    5. Multiplex hint → `STATUS_MULTIPLEXED`.
    6. Any unresolved SRRs → `STATUS_UNKNOWN`.
    7. Otherwise → `STATUS_PARTIAL`.
    """
    bp_id = bioproject.get("bioproject_id", "")
    declared = _parse_int(bioproject.get("n_rounds_declared"))
    sample_srrs = [(s.get("srr") or "").strip() for s in samples if (s.get("srr") or "").strip()]

    assigned_rounds: list[int] = []
    unknown_srrs: list[str] = []
    confidence_counts: Counter = Counter()

    for srr in sample_srrs:
        round_rec = rounds_by_srr.get(srr)
        if not round_rec:
            unknown_srrs.append(srr)
            confidence_counts["missing"] += 1
            continue

        confidence = round_rec.get("confidence") or "missing"
        confidence_counts[confidence] += 1

        rn = _parse_int(round_rec.get("round_number"))
        if rn is None:
            unknown_srrs.append(srr)
        else:
            assigned_rounds.append(rn)

    assigned_unique = sorted(set(assigned_rounds))
    multiplex_hint = _has_multiplex_hint(
        bioproject.get("manual_curation_notes"),
        bioproject.get("study_title"),
    )

    if not sample_srrs:
        status = STATUS_UNKNOWN
        reason = "no_public_srrs"
    elif declared is None:
        if multiplex_hint:
            status = STATUS_MULTIPLEXED
            reason = "multiplex_hint_without_declared_round_count"
        else:
            status = STATUS_UNKNOWN
            reason = "declared_round_count_missing"
    elif len(assigned_unique) >= declared:
        status = STATUS_ALL
        reason = "declared_rounds_observed"
    elif multiplex_hint:
        status = STATUS_MULTIPLEXED
        reason = "manual_notes_indicate_multiplexed_rounds"
    elif unknown_srrs:
        status = STATUS_UNKNOWN
        reason = "unresolved_round_assignments_present"
    else:
        status = STATUS_PARTIAL
        reason = "declared_rounds_exceed_observed_rounds"

    return {
        "bioproject_id": bp_id,
        "round_coverage_status": status,
        "reason": reason,
        "n_rounds_declared": declared,
        "n_srr": len(sample_srrs),
        "n_rounds_assigned": len(assigned_unique),
        "assigned_rounds": assigned_unique,
        "unknown_round_srrs": unknown_srrs,
        "unknown_round_count": len(unknown_srrs),
        "multiplex_hint": multiplex_hint,
        "round_confidence_counts": dict(confidence_counts),
    }


def _default_include_filter(bp: dict) -> bool:
    return (bp.get("include", "") or "").strip().lower() == "y"


def build_round_coverage_report(
    bioprojects: list[dict],
    samples: list[dict],
    rounds: list[dict],
    include_filter: Callable[[dict], bool] | None = None,
) -> dict:
    """Aggregate per-BioProject classifications into a single report.

    Parameters
    ----------
    bioprojects, samples, rounds
        The three relational tables as lists of row dicts.
    include_filter
        Predicate determining which BioProjects feed the `included_*` summary.
        Default: BioProjects whose `include` column is `y`. Pass your own
        callable for thesis-specific filtering (e.g. RNA-only, target-specific).
    """
    samples_by_bp: dict[str, list[dict]] = defaultdict(list)
    for sample in samples:
        samples_by_bp[sample.get("bioproject_id", "")].append(sample)

    rounds_by_srr = {r["srr"]: r for r in rounds if r.get("srr")}

    entries = [
        classify_bioproject_round_coverage(
            bioproject=bp,
            samples=samples_by_bp.get(bp.get("bioproject_id", ""), []),
            rounds_by_srr=rounds_by_srr,
        )
        for bp in bioprojects
    ]

    if include_filter is None:
        include_filter = _default_include_filter

    status_counts = Counter(entry["round_coverage_status"] for entry in entries)
    included_entries = [
        entry for entry, bp in zip(entries, bioprojects, strict=True) if include_filter(bp)
    ]
    included_status_counts = Counter(entry["round_coverage_status"] for entry in included_entries)

    return {
        "generated_at": datetime.now().isoformat(),
        "summary": {
            "n_bioprojects": len(entries),
            "status_counts": dict(status_counts),
            "included_status_counts": dict(included_status_counts),
        },
        "included_bioprojects": included_entries,
        "entries": entries,
    }
