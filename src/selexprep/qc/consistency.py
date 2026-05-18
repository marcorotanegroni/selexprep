"""K-mer Jaccard consistency check across SELEX rounds (diagnostic only).

Computes k-mer Jaccard distance between all per-round count parquets of a
BioProject, then verifies that the assigned round ordering is monotonic with
respect to pairwise distances (neighbouring rounds should be more similar
than distant rounds).

**Epistemological rule: this module NEVER infers or corrects round numbers
from enrichment signals.** Doing so would be CIRCULAR — the signal we want
to measure is being used as the label, contaminating downstream training.
The module is STRICTLY diagnostic: it flags inconsistencies for manual
review and writes them to a report. Callers decide how to resolve them.

Algorithm:
  1. For each BioProject with ≥3 rounds assigned:
     - Load each `round_XX.counts.parquet`, sample top-N (default 10k) sequences
     - Extract canonical k-mers (k=6) from each round's top sequences
     - Compute pairwise Jaccard distance → distance matrix
  2. Check: dist(R_i, R_{i+1}) < dist(R_i, R_{i+k}) for all k > 1 (10% tolerance)
  3. Flag BioProjects where this monotonicity fails.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


DEFAULT_K = 6
DEFAULT_TOP_N = 10_000
_MONOTONICITY_TOLERANCE = 1.1


# ---------------------------------------------------------------------------
# K-mer extraction
# ---------------------------------------------------------------------------


def canonical_kmer(kmer: str) -> str:
    """Return the lexicographically smaller of `kmer` and its reverse complement."""
    complement = str.maketrans("ACGTN", "TGCAN")
    rc = kmer.translate(complement)[::-1]
    return kmer if kmer < rc else rc


def extract_kmers(sequence: str, k: int) -> set[str]:
    """Extract canonical k-mers from a single sequence."""
    if len(sequence) < k:
        return set()
    return {canonical_kmer(sequence[i : i + k]) for i in range(len(sequence) - k + 1)}


def _kmer_set_for_round(parquet: Path, k: int, top_n: int) -> set[str]:
    """Take top-N sequences from a counts parquet; return union of canonical k-mers."""
    df = pd.read_parquet(parquet, columns=["sequence", "reads"])
    df = df.nlargest(top_n, "reads")
    kmers: set[str] = set()
    for seq in df["sequence"]:
        kmers.update(extract_kmers(seq, k))
    return kmers


# ---------------------------------------------------------------------------
# Jaccard distance
# ---------------------------------------------------------------------------


def jaccard_distance(a: set[str], b: set[str]) -> float:
    """Jaccard distance between two sets. Returns 0.0 if both empty, 1.0 if disjoint."""
    if not a and not b:
        return 0.0
    union = len(a | b)
    if union == 0:
        return 1.0
    return 1.0 - (len(a & b) / union)


# ---------------------------------------------------------------------------
# Monotonicity check
# ---------------------------------------------------------------------------


def check_monotonicity(
    ordered_rounds: list[int],
    distance_matrix: dict[tuple[int, int], float],
) -> dict:
    """Check if consecutive rounds have smaller distances than distant rounds.

    Returns a dict with ``monotonic: bool`` and ``violations: list[dict]``.
    10% tolerance on the consecutive < distant inequality.
    """
    violations = []
    is_monotonic = True

    for i, ri in enumerate(ordered_rounds):
        for j in range(i + 2, len(ordered_rounds)):
            rj = ordered_rounds[j]
            r_next = ordered_rounds[i + 1]
            d_consecutive = distance_matrix.get((min(ri, r_next), max(ri, r_next)))
            d_distant = distance_matrix.get((min(ri, rj), max(ri, rj)))
            if d_consecutive is None or d_distant is None:
                continue
            if d_consecutive > d_distant * _MONOTONICITY_TOLERANCE:
                is_monotonic = False
                violations.append(
                    {
                        "from_round": ri,
                        "neighbour_round": r_next,
                        "distant_round": rj,
                        "dist_consecutive": round(d_consecutive, 4),
                        "dist_distant": round(d_distant, 4),
                        "reason": "consecutive distance > distant distance",
                    }
                )

    return {"monotonic": is_monotonic, "violations": violations}


# ---------------------------------------------------------------------------
# Primer-trim provenance
# ---------------------------------------------------------------------------


def _primer_trim_status(bp_dir: Path) -> dict:
    """Inspect `summary.json` for primer-trim provenance flags.

    Untrimmed primer flanks inflate k-mer Jaccard similarity between rounds
    (shared universal primers), which can mask convergence or create false
    monotonicity. This surfaces a diagnostic flag — never a correction.
    """
    summary_path = bp_dir / "summary.json"
    if not summary_path.exists():
        return {"status": "no_summary"}
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {"status": "summary_unreadable"}

    primer_known = not summary.get("primer_unknown", True)
    per_round = {r["round"]: r.get("primer_trim_applied") for r in summary.get("rounds", [])}
    states = set(per_round.values())

    flags = []
    if primer_known and summary.get("primer_trim_applied") is False:
        flags.append("primers_known_but_not_trimmed")
    if primer_known and len(states - {None}) > 1:
        flags.append("primer_trim_mixed_across_rounds")
    if primer_known and None in states and (True in states or False in states):
        flags.append("primer_trim_partially_cached")

    return {
        "primer_known": primer_known,
        "aggregate_trim_applied": summary.get("primer_trim_applied"),
        "per_round_trim_applied": per_round,
        "flags": flags,
    }


# ---------------------------------------------------------------------------
# Per-BioProject check
# ---------------------------------------------------------------------------


def check_bioproject(
    bp_dir: Path,
    bioproject_id: str | None = None,
    k: int = DEFAULT_K,
    top_n: int = DEFAULT_TOP_N,
) -> dict:
    """Diagnostic consistency check for a single processed BioProject directory.

    `bp_dir` must contain `round_*.counts.parquet` files. Returns a dict with
    distance matrix, monotonicity verdict, violations, and primer-trim
    provenance flags. `bioproject_id` defaults to `bp_dir.name`.
    """
    if bioproject_id is None:
        bioproject_id = bp_dir.name

    parquets = sorted(bp_dir.glob("round_*.counts.parquet"))
    if len(parquets) < 3:
        return {
            "bioproject_id": bioproject_id,
            "status": "insufficient_rounds",
            "n_rounds": len(parquets),
        }

    round_kmers: dict[int, set[str]] = {}
    for p in parquets:
        stem = p.stem.replace(".counts", "")
        try:
            rn = int(stem.split("_")[-1])
        except ValueError:
            logger.warning("[%s] cannot parse round number from %s", bioproject_id, p.name)
            continue
        round_kmers[rn] = _kmer_set_for_round(p, k, top_n)

    if len(round_kmers) < 3:
        return {
            "bioproject_id": bioproject_id,
            "status": "insufficient_parseable_rounds",
            "n_rounds": len(round_kmers),
        }

    distance_matrix: dict[tuple[int, int], float] = {}
    round_nums = sorted(round_kmers.keys())
    for i, r1 in enumerate(round_nums):
        for r2 in round_nums[i + 1 :]:
            distance_matrix[(r1, r2)] = jaccard_distance(round_kmers[r1], round_kmers[r2])

    monotonicity = check_monotonicity(round_nums, distance_matrix)
    trim_status = _primer_trim_status(bp_dir)

    result = {
        "bioproject_id": bioproject_id,
        "status": "ok" if monotonicity["monotonic"] else "suspicious",
        "n_rounds": len(round_nums),
        "round_numbers": round_nums,
        "k": k,
        "top_n_per_round": top_n,
        "distance_matrix": {f"R{r1}-R{r2}": round(d, 4) for (r1, r2), d in distance_matrix.items()},
        "monotonic": monotonicity["monotonic"],
        "violations": monotonicity["violations"],
        "primer_trim_status": trim_status,
    }

    flags: list[str] = []
    if not monotonicity["monotonic"]:
        flags.append("round_ordering_suspicious")
    flags.extend(trim_status.get("flags", []))
    if flags:
        result["flags"] = flags

    return result


# ---------------------------------------------------------------------------
# Multi-BP runner
# ---------------------------------------------------------------------------


def run_consistency_check(
    processed_root: Path,
    bioproject_id: str | None = None,
    k: int = DEFAULT_K,
    top_n: int = DEFAULT_TOP_N,
    output_path: Path | None = None,
) -> dict:
    """Run consistency check on every (or one) BioProject under `processed_root`.

    Each BP directory must contain a `summary.json` to be processed. Result is
    optionally written to `output_path` as JSON.
    """
    if bioproject_id:
        bp_ids = [bioproject_id]
    else:
        bp_ids = sorted(
            p.name for p in processed_root.iterdir() if p.is_dir() and (p / "summary.json").exists()
        )

    if not bp_ids:
        logger.warning("No processed BioProjects found in %s", processed_root)
        return {}

    logger.info("Checking consistency for %d BioProject(s)", len(bp_ids))
    reports: dict[str, dict] = {}
    for bp_id in bp_ids:
        logger.info("[%s] k-mer Jaccard consistency check", bp_id)
        try:
            reports[bp_id] = check_bioproject(processed_root / bp_id, bp_id, k=k, top_n=top_n)
        except Exception as e:
            logger.error("[%s] failed: %s", bp_id, e)
            reports[bp_id] = {"bioproject_id": bp_id, "status": "error", "error": str(e)}

    n_ok = sum(1 for r in reports.values() if r.get("status") == "ok")
    n_suspicious = sum(1 for r in reports.values() if r.get("status") == "suspicious")
    n_insufficient = sum(1 for r in reports.values() if "insufficient" in r.get("status", ""))
    n_primer_trim_flagged = sum(
        1 for r in reports.values() if r.get("primer_trim_status", {}).get("flags")
    )

    output = {
        "generated_at": datetime.now().isoformat(),
        "k": k,
        "top_n_per_round": top_n,
        "summary": {
            "n_bioprojects": len(reports),
            "n_ok": n_ok,
            "n_suspicious": n_suspicious,
            "n_insufficient_data": n_insufficient,
            "n_primer_trim_flagged": n_primer_trim_flagged,
        },
        "reports": reports,
    }

    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(output, indent=2), encoding="utf-8")
        logger.info(
            "Consistency report written → %s (%d ok, %d suspicious, %d insufficient, %d primer-trim flagged)",
            output_path,
            n_ok,
            n_suspicious,
            n_insufficient,
            n_primer_trim_flagged,
        )

    if n_suspicious > 0:
        logger.warning(
            "%d BioProjects have round ordering inconsistencies — REVIEW MANUALLY. "
            "Do NOT auto-reassign rounds based on this report.",
            n_suspicious,
        )

    return output
