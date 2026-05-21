"""Depth-aware suspicion flags for a single SELEX dataset.

Locked plan lines 350-358 spec eight flags. Each is computed from the
manifest's ``LibraryReport`` + the per-round ``counts.parquet`` outputs
emitted by ``selexprep count``. The output is ``flags.yaml`` (sorted by
flag name; deterministic).

The match-rate threshold is **imported** from ``library.detect`` to keep
the two thresholds in sync — when Codex calibration tunes
``UNABLE_TO_EXTRACT_MATCH_RATE`` it automatically tightens the QC flag
as well.

**v0.1 scope**: ``extraction_mode_changed_across_rounds`` cannot fire
in single-dataset mode (only one ``LibraryReport`` per manifest). It is
a placeholder for the ``selexprep run`` batch driver (Phase 6).
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from typing import Literal

import yaml

from selexprep.library.detect import UNABLE_TO_EXTRACT_MATCH_RATE
from selexprep.library.report import LibraryReport
from selexprep.manifest import SelexprepManifestV1
from selexprep.qc.diversity import rarefy, unique_count

logger = logging.getLogger(__name__)


Severity = Literal["warn", "info"]


# ---------------------------------------------------------------------------
# CALIBRATION-TODO constants
# ---------------------------------------------------------------------------

# CALIBRATION-TODO: locked plan line 351 ("rarefied R_n > R_{n-1}").
# Codex confirms rarefaction depth after Phase 6 benchmark numbers.
RAREFACTION_DEPTH = 10_000

# Match-rate threshold for "low primer match" is the SAME number as the
# UNABLE_TO_EXTRACT classifier threshold - imported, not redeclared.

# CALIBRATION-TODO: locked plan line 353 ("> 2 modal lengths").
N_LENGTH_MAX_MODAL_LENGTHS = 2

# CALIBRATION-TODO: locked plan line 354 ("> 20% wrong orientation").
STRAND_MIX_MAX_REVERSE_FRACTION = 0.20

# CALIBRATION-TODO: locked plan line 355 ("< 10k any round").
LOW_TOTAL_READS_MIN = 10_000

# CALIBRATION-TODO: locked plan line 356 ("> 5% of reads").
ADAPTER_CONTAMINATION_MAX_FRACTION = 0.05

# Diversity rarefaction RNG seed - deterministic across reruns.
_RAREFY_SEED = 42


# ---------------------------------------------------------------------------
# Flag dataclass
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Flag:
    """One suspicion flag raised by the QC layer."""

    name: str
    severity: Severity
    evidence: dict[str, object]


# ---------------------------------------------------------------------------
# Individual checks
# ---------------------------------------------------------------------------


def check_unexpected_rarefied_diversity_increase(
    counts_by_round: dict[int, dict[str, int]],
) -> Flag | None:
    """Rarefied unique counts should monotonically decrease across rounds.

    Locked plan line 351: an increase is suspicious - it suggests
    contamination or that a "round" was actually a re-sequencing of the
    same pool. Uses rarefaction (not raw counts) so depth differences
    do not confound the comparison.

    Phase 5 Codex pass 1 fix: rarefy to ``effective_depth = min(
    RAREFACTION_DEPTH, min_total_reads_per_round)``. The previous version
    used ``RAREFACTION_DEPTH`` as a fixed target, but ``rarefy`` returns
    the original pool unchanged when ``depth >= total`` - so a 2k-read
    round was being compared against a 10k-rarefied round, reintroducing
    exactly the depth-confounding the flag is supposed to avoid.
    ``low_total_reads`` covers the absolute-depth concern separately.
    """
    if len(counts_by_round) < 2:
        return None

    total_per_round = {r: sum(counts_by_round[r].values()) for r in counts_by_round}
    min_total = min(total_per_round.values())
    if min_total <= 0:
        # An empty round is a different problem (the low_total_reads flag
        # will surface it); rarefaction comparison is undefined.
        return None
    effective_depth = min(RAREFACTION_DEPTH, min_total)

    rarefied_uniques: dict[int, int] = {}
    for r in sorted(counts_by_round):
        rarefied = rarefy(counts_by_round[r], depth=effective_depth, seed=_RAREFY_SEED)
        rarefied_uniques[r] = unique_count(rarefied)

    rounds = sorted(rarefied_uniques)
    increases: list[dict[str, object]] = []
    for prev, curr in pairwise(rounds):
        u_prev = rarefied_uniques[prev]
        u_curr = rarefied_uniques[curr]
        if u_curr > u_prev:
            increases.append(
                {
                    "from_round": prev,
                    "to_round": curr,
                    "rarefied_unique_prev": u_prev,
                    "rarefied_unique_curr": u_curr,
                }
            )
    if not increases:
        return None
    return Flag(
        name="unexpected_rarefied_diversity_increase",
        severity="warn",
        evidence={
            "configured_depth": RAREFACTION_DEPTH,
            "effective_depth": effective_depth,
            "min_total_reads": min_total,
            "rarefied_uniques_per_round": rarefied_uniques,
            "increases": increases,
        },
    )


def check_low_primer_match(lr: LibraryReport) -> Flag | None:
    """Either 5' or 3' match rate below UNABLE_TO_EXTRACT_MATCH_RATE."""
    rate_5p = lr.match_rate_5p
    rate_3p = lr.match_rate_3p
    below = []
    if rate_5p < UNABLE_TO_EXTRACT_MATCH_RATE:
        below.append({"side": "5p", "match_rate": rate_5p})
    if rate_3p < UNABLE_TO_EXTRACT_MATCH_RATE:
        below.append({"side": "3p", "match_rate": rate_3p})
    if not below:
        return None
    return Flag(
        name="low_primer_match",
        severity="warn",
        evidence={"threshold": UNABLE_TO_EXTRACT_MATCH_RATE, "sides_below": below},
    )


def check_n_length_variation_across_rounds(
    counts_by_round: dict[int, dict[str, int]],
) -> Flag | None:
    """More than ``N_LENGTH_MAX_MODAL_LENGTHS`` distinct modal lengths across rounds."""
    if not counts_by_round:
        return None
    modal_per_round: dict[int, int] = {}
    for r, counts in counts_by_round.items():
        if not counts:
            continue
        # Modal length = most common sequence length weighted by reads.
        length_counter: Counter[int] = Counter()
        for seq, n in counts.items():
            length_counter[len(seq)] += n
        if not length_counter:
            continue
        modal_per_round[r] = length_counter.most_common(1)[0][0]

    distinct_modes = set(modal_per_round.values())
    if len(distinct_modes) <= N_LENGTH_MAX_MODAL_LENGTHS:
        return None
    return Flag(
        name="n_length_variation_across_rounds",
        severity="warn",
        evidence={
            "max_modal_lengths": N_LENGTH_MAX_MODAL_LENGTHS,
            "modal_length_per_round": modal_per_round,
            "distinct_modes": sorted(distinct_modes),
        },
    )


def check_strand_mix(strand_report_path: Path | None) -> Flag | None:
    """> STRAND_MIX_MAX_REVERSE_FRACTION reverse reads in any round.

    The strand report is the TSV emitted by Phase 3's extract pipeline
    (only present when ``orientation in {MIXED, REVERSE}``). Missing
    file -> no flag (FORWARD orientation never emits a strand report).
    """
    if strand_report_path is None or not strand_report_path.exists():
        return None
    high_rounds: list[dict[str, object]] = []
    with strand_report_path.open(encoding="utf-8") as fh:
        header = fh.readline().strip().split("\t")
        for raw in fh:
            row = raw.rstrip("\n").split("\t")
            if len(row) != len(header):
                continue
            d = dict(zip(header, row, strict=False))
            try:
                fwd = int(d["forward"])
                rev = int(d["reverse"])
                amb = int(d["ambiguous"])
            except (KeyError, ValueError):
                continue
            total = fwd + rev + amb
            if total == 0:
                continue
            frac_reverse = rev / total
            if frac_reverse > STRAND_MIX_MAX_REVERSE_FRACTION:
                high_rounds.append(
                    {
                        "round": int(d["round"]),
                        "reverse_fraction": frac_reverse,
                        "forward": fwd,
                        "reverse": rev,
                        "ambiguous": amb,
                    }
                )
    if not high_rounds:
        return None
    return Flag(
        name="strand_mix",
        severity="warn",
        evidence={
            "max_reverse_fraction": STRAND_MIX_MAX_REVERSE_FRACTION,
            "rounds_above_threshold": high_rounds,
        },
    )


def check_low_total_reads(counts_by_round: dict[int, dict[str, int]]) -> Flag | None:
    """Any round below ``LOW_TOTAL_READS_MIN`` total reads."""
    if not counts_by_round:
        return None
    low_rounds: list[dict[str, object]] = []
    for r in sorted(counts_by_round):
        total = sum(counts_by_round[r].values())
        if total < LOW_TOTAL_READS_MIN:
            low_rounds.append({"round": r, "total_reads": total})
    if not low_rounds:
        return None
    return Flag(
        name="low_total_reads",
        severity="warn",
        evidence={"min_threshold": LOW_TOTAL_READS_MIN, "low_rounds": low_rounds},
    )


def check_adapter_contamination_high(
    lr: LibraryReport,
    trim_reports_by_round: dict[int, dict[str, int]] | None = None,
) -> Flag | None:
    """Any adapter's hits > ADAPTER_CONTAMINATION_MAX_FRACTION of reads.

    Phase 5 Codex pass 1 fix: the numerator (``lr.known_adapter_hits``)
    comes from Phase 2's inference pass on the EARLIEST round's
    subsampled pre-extraction reads. The denominator must match that
    measurement universe — using post-extraction ``counts.parquet``
    totals (the previous behavior) compared apples to oranges and
    could inflate or deflate the fraction arbitrarily.

    Correct denominator: ``trim_reports_by_round[earliest_round]["n_in"]
    * lr.read_fraction_used_for_inference``. When ``trim_reports`` is
    not available (e.g. manifest hand-edited, qc run on a partial
    artifact set), surface the flag as ``severity="info"`` with
    ``reason="denominator_unavailable"`` rather than silently guessing.
    """
    if not lr.known_adapter_hits:
        return None

    if not trim_reports_by_round:
        return Flag(
            name="adapter_contamination_high",
            severity="info",
            evidence={
                "reason": "denominator_unavailable",
                "note": (
                    "trim_reports.json is missing; cannot reconstruct the "
                    "Phase 2 inference-universe denominator. Hits are "
                    "preserved as raw counts."
                ),
                "known_adapter_hits": {a: int(h) for a, h in lr.known_adapter_hits.items()},
            },
        )

    # ASSUMPTION (v0.1): the lowest round number in trim_reports.json is
    # the same round Phase 2's primer inference used. This holds in the
    # documented v0.1 CLI flow because the user passes the same FASTQ set
    # to both `detect` and `extract`, so detect's "earliest round" equals
    # min(extract's trim_reports.rounds). The assumption can BREAK in:
    #   - partial extract runs (user processed only a subset of detect's
    #     rounds)
    #   - manual artifact stitching (hand-assembled trim_reports.json)
    #   - v0.2 batch driver (`selexprep run`) where multiple extract
    #     invocations may produce a Frankenstein trim_reports.json
    # Bulletproof fix (deferred to v0.2): add ``earliest_inference_round``
    # to the LibraryReport schema in Phase 2's compute_library_report; use
    # ``lr.earliest_inference_round`` here instead of min(trim_reports).
    # A schema bump is disproportionate for an edge case unreachable
    # through the v0.1 single-dataset CLI flow.
    earliest_round = min(trim_reports_by_round)
    n_in_earliest = int(trim_reports_by_round[earliest_round].get("n_in", 0))
    denominator = int(n_in_earliest * lr.read_fraction_used_for_inference)
    if denominator <= 0:
        return Flag(
            name="adapter_contamination_high",
            severity="info",
            evidence={
                "reason": "denominator_zero",
                "note": (
                    "Inference-universe denominator computed as 0 "
                    f"(n_in_earliest={n_in_earliest}, "
                    f"read_fraction_used_for_inference={lr.read_fraction_used_for_inference}). "
                    "Cannot evaluate adapter fraction."
                ),
                "known_adapter_hits": {a: int(h) for a, h in lr.known_adapter_hits.items()},
            },
        )

    over: list[dict[str, object]] = []
    for adapter, hits in lr.known_adapter_hits.items():
        frac = hits / denominator
        if frac > ADAPTER_CONTAMINATION_MAX_FRACTION:
            over.append(
                {
                    "adapter": adapter,
                    "hits": int(hits),
                    "fraction_of_reads": frac,
                }
            )
    if not over:
        return None
    return Flag(
        name="adapter_contamination_high",
        severity="warn",
        evidence={
            "max_fraction": ADAPTER_CONTAMINATION_MAX_FRACTION,
            "denominator_basis": "trim_reports.json n_in[earliest_round] * read_fraction_used_for_inference",
            "earliest_round": earliest_round,
            "denominator": denominator,
            "adapters_above_threshold": over,
        },
    )


def check_extraction_mode_changed_across_rounds(
    manifests: list[SelexprepManifestV1] | None = None,
) -> Flag | None:
    """v0.1 single-dataset scope: cannot fire (one LR per manifest).

    Placeholder for the ``selexprep run`` batch driver (Phase 6) which
    will compare ``extraction_mode`` across multiple manifests from the
    same SRA project.
    """
    if not manifests or len(manifests) < 2:
        return None
    modes = {m.extraction_mode for m in manifests}
    if len(modes) <= 1:
        return None
    return Flag(
        name="extraction_mode_changed_across_rounds",
        severity="warn",
        evidence={"distinct_extraction_modes": sorted(modes)},
    )


def check_requires_read_merging_for_full_insert(lr: LibraryReport) -> Flag | None:
    """Informational: paired-end split-primer datasets need read merging
    to recover the full insert (locked plan line 358)."""
    if lr.required_action != "READ_MERGING_RECOMMENDED":
        return None
    return Flag(
        name="requires_read_merging_for_full_insert",
        severity="info",
        evidence={
            "extraction_mode": lr.extraction_mode,
            "required_action": lr.required_action,
            "note": (
                "Paired-end split-primer libraries need merging (v0.2) "
                "to recover the full insert; v0.1 emits R1/R2 partials only."
            ),
        },
    )


# ---------------------------------------------------------------------------
# Aggregator + YAML emitter
# ---------------------------------------------------------------------------


def compute_all_flags(
    manifest: SelexprepManifestV1,
    counts_by_round: dict[int, dict[str, int]],
    *,
    strand_report_path: Path | None = None,
    trim_reports_by_round: dict[int, dict[str, int]] | None = None,
) -> list[Flag]:
    """Run every flag check; return the list of flags that fired.

    ``trim_reports_by_round`` is the aggregated per-round ``{n_in, n_out}``
    derived from ``trim_reports.json`` (qc/runner does the parsing). Used
    by ``check_adapter_contamination_high`` to derive the correct
    denominator (Phase 2's inference universe, not the post-extraction
    universe). When absent, that flag degrades gracefully.
    """
    lr = manifest.library_report
    candidates: list[Flag | None] = [
        check_unexpected_rarefied_diversity_increase(counts_by_round),
        check_low_primer_match(lr),
        check_n_length_variation_across_rounds(counts_by_round),
        check_strand_mix(strand_report_path),
        check_low_total_reads(counts_by_round),
        check_adapter_contamination_high(lr, trim_reports_by_round),
        check_extraction_mode_changed_across_rounds(None),  # v0.1 always None
        check_requires_read_merging_for_full_insert(lr),
    ]
    return [f for f in candidates if f is not None]


def write_flags_yaml(flags: list[Flag], path: Path) -> None:
    """Write `flags` as deterministic YAML (sorted by flag name)."""
    payload = [
        {
            "name": f.name,
            "severity": f.severity,
            "evidence": _normalize_for_yaml(f.evidence),
        }
        for f in sorted(flags, key=lambda x: x.name)
    ]
    text = yaml.safe_dump(payload, sort_keys=True, default_flow_style=False)
    if not text or text == "[]\n":
        text = "[]\n"
    if not text.endswith("\n"):
        text += "\n"
    path.write_text(text, encoding="utf-8")


def _normalize_for_yaml(obj: object) -> object:
    """Recursively convert dict-of-int-keys to dict-of-str-keys for stable YAML."""
    if isinstance(obj, dict):
        return {str(k): _normalize_for_yaml(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_normalize_for_yaml(x) for x in obj]
    if isinstance(obj, float):
        # Round floats to 6 decimal places so YAML output is stable
        # across float-repr platform quirks.
        return round(obj, 6)
    return obj
