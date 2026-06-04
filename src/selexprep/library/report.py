"""Public `LibraryReport` contract + locked classification table.

This module is the **public contract** between primer inference and every
downstream stage (extract, count, qc, manifest). It is declared
**strict-mypy** in ``pyproject.toml`` because schema drift here breaks
the JSON-on-disk format the manifest hashes.

Three things live here and only here:

1. **Literal aliases** for the five categorical fields. Importing these
   aliases (rather than re-typing each ``Literal[...]``) keeps the schema
   and the classifier in sync.
2. **``LibraryReport``** — the pydantic v2 ``BaseModel``. Fields are
   verbatim from the design.
3. **``_classify``** — the pure-function decision table from the design,
   plus the status-cap rule (no round map ⇒
   status ≤ MEDIUM). Takes thresholds as kwargs so the algorithm module
   (``library/detect.py``) holds the single source of truth for
   calibration constants and tests can drive the table directly.
4. **Deterministic JSON I/O** — ``write_library_report_json`` /
   ``read_library_report_json``. Bit-identical output across reruns; the
   manifest's ``output_sha256`` depends on it.

**Out of scope for this module:** the inference algorithm itself
(per-signal helpers, persistence score, composite confidence) — those
live in ``library/detect.py``. This module is declarative.
"""

from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path
from typing import Literal, NamedTuple

from pydantic import BaseModel, ConfigDict

# ---------------------------------------------------------------------------
# Literal aliases — one source of truth for the five categorical fields.
# ---------------------------------------------------------------------------

ExtractionMode = Literal[
    "BOTH_PRIMERS_SINGLE_READ",
    "FIVE_PRIME_ONLY",
    "THREE_PRIME_ONLY",
    "PAIRED_END_SPLIT_PRIMERS",
    "UNABLE_TO_EXTRACT",
]
"""Biology-only descriptor of what the reads contain."""

ReadSource = Literal["R1", "R2", "R1_AND_R2", "INTERLEAVED", "UNKNOWN"]
"""Which physical reads carry the random region."""

RequiredAction = Literal["NONE", "MANUAL_PRIMERS_REQUIRED", "READ_MERGING_RECOMMENDED"]
"""Workflow guidance separate from biology (extraction_mode)."""

Orientation = Literal["FORWARD", "REVERSE", "MIXED"]
"""Strand orientation summary (diagnostic; acts on this)."""

Status = Literal["HIGH", "MEDIUM", "LOW", "UNABLE_TO_INFER"]
"""Composite-confidence status. Capped at MEDIUM when no round map is
available."""


# ---------------------------------------------------------------------------
# LibraryReport — public schema. Fields verbatim from the design 233-285.
# ---------------------------------------------------------------------------


class LibraryReport(BaseModel):
    """Inferred primer + library-structure report for one SELEX dataset.

    Locked schema. Field order and types match the plan exactly; changes
    here are a schema-version bump (out of scope for v0.1).
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    # ---- Primer detection ------------------------------------------------
    primer_5p: str | None
    primer_3p: str | None
    variants_5p: list[tuple[str, int]]
    variants_3p: list[tuple[str, int]]

    # ---- Adapter blacklist (records hits; does NOT filter reads) --------
    known_adapter_hits: dict[str, int]

    # ---- Biology ---------------------------------------------------------
    extraction_mode: ExtractionMode
    full_insert_recovered: bool

    # ---- File layout -----------------------------------------------------
    read_source: ReadSource

    # ---- Workflow guidance ----------------------------------------------
    required_action: RequiredAction

    # ---- Orientation ----------------------------------------------------
    orientation: Orientation

    # ---- N-region -------------------------------------------------------
    n_length_mode: int | None
    n_length_distribution: dict[int, int]
    n_length_confidence: float

    # ---- Match quality --------------------------------------------------
    match_rate_5p: float
    match_rate_3p: float
    position_consistency_5p: float
    position_consistency_3p: float

    # ---- Reproducibility ------------------------------------------------
    read_fraction_used_for_inference: float
    sampling_seed: int

    # ---- Composite ------------------------------------------------------
    confidence: float
    status: Status
    failure_reason: str | None


# ---------------------------------------------------------------------------
# Classifier — the design + status cap.
# ---------------------------------------------------------------------------


class Classification(NamedTuple):
    """Result of ``_classify``: the four output fields driven by the decision table."""

    extraction_mode: ExtractionMode
    full_insert_recovered: bool
    required_action: RequiredAction
    status: Status
    failure_reason: str | None


def _classify(
    *,
    match_rate_5p: float,
    match_rate_3p: float,
    n_length_confidence: float,
    has_paired_split: bool,
    paired_has_overlap: bool,
    has_round_map: bool,
    composite_confidence: float,
    primer_found_threshold: float,
    n_length_confident_threshold: float,
    unable_to_extract_threshold: float,
    status_high_cutoff: float,
    status_medium_cutoff: float,
    status_low_cutoff: float,
) -> Classification:
    """Map detection signals to (extraction_mode, full_insert_recovered, required_action, status).

    Pure function. Thresholds are passed in (not read from globals) so the
    algorithm module owns the single source of truth and tests can drive
    the decision table directly.

    Implements the decision table plus the cap "no round map ⇒ status ≤ MEDIUM".
    """
    # Row 8 (catch-all failure): both rates below the 0.4 floor. Wins over
    # every other case — even a paired-split signal is meaningless if
    # neither side has primer evidence.
    if match_rate_5p < unable_to_extract_threshold and match_rate_3p < unable_to_extract_threshold:
        return Classification(
            extraction_mode="UNABLE_TO_EXTRACT",
            full_insert_recovered=False,
            required_action="MANUAL_PRIMERS_REQUIRED",
            status="UNABLE_TO_INFER",
            failure_reason=(
                f"Both primer match rates below {unable_to_extract_threshold:.2f} "
                f"(5'={match_rate_5p:.2f}, 3'={match_rate_3p:.2f})"
            ),
        )

    # Rows 6 + 7: paired-end split (R1 carries 5', R2 carries 3').
    if has_paired_split:
        if paired_has_overlap:
            # Row 7 — v0.2 territory; v0.1 cannot reach this branch
            # (paired_has_overlap is hard-coded False in the orchestrator
            # until read merging ships).
            return Classification(
                extraction_mode="PAIRED_END_SPLIT_PRIMERS",
                full_insert_recovered=True,
                required_action="NONE",
                status=_assign_status(
                    extraction_mode="PAIRED_END_SPLIT_PRIMERS",
                    composite_confidence=composite_confidence,
                    has_round_map=has_round_map,
                    status_high_cutoff=status_high_cutoff,
                    status_medium_cutoff=status_medium_cutoff,
                    status_low_cutoff=status_low_cutoff,
                ),
                failure_reason=None,
            )
        return Classification(
            extraction_mode="PAIRED_END_SPLIT_PRIMERS",
            full_insert_recovered=False,
            required_action="READ_MERGING_RECOMMENDED",
            status=_assign_status(
                extraction_mode="PAIRED_END_SPLIT_PRIMERS",
                composite_confidence=composite_confidence,
                has_round_map=has_round_map,
                status_high_cutoff=status_high_cutoff,
                status_medium_cutoff=status_medium_cutoff,
                status_low_cutoff=status_low_cutoff,
            ),
            failure_reason=None,
        )

    # Row 1: both primers strong on the same read.
    if match_rate_5p > primer_found_threshold and match_rate_3p > primer_found_threshold:
        return Classification(
            extraction_mode="BOTH_PRIMERS_SINGLE_READ",
            full_insert_recovered=True,
            required_action="NONE",
            status=_assign_status(
                extraction_mode="BOTH_PRIMERS_SINGLE_READ",
                composite_confidence=composite_confidence,
                has_round_map=has_round_map,
                status_high_cutoff=status_high_cutoff,
                status_medium_cutoff=status_medium_cutoff,
                status_low_cutoff=status_low_cutoff,
            ),
            failure_reason=None,
        )

    # Rows 2 + 3: only 5' strong.
    if match_rate_5p > primer_found_threshold:
        if n_length_confidence > n_length_confident_threshold:
            return Classification(
                extraction_mode="FIVE_PRIME_ONLY",
                full_insert_recovered=False,
                required_action="NONE",
                status=_assign_status(
                    extraction_mode="FIVE_PRIME_ONLY",
                    composite_confidence=composite_confidence,
                    has_round_map=has_round_map,
                    status_high_cutoff=status_high_cutoff,
                    status_medium_cutoff=status_medium_cutoff,
                    status_low_cutoff=status_low_cutoff,
                ),
                failure_reason=None,
            )
        return Classification(
            extraction_mode="UNABLE_TO_EXTRACT",
            full_insert_recovered=False,
            required_action="MANUAL_PRIMERS_REQUIRED",
            status="UNABLE_TO_INFER",
            failure_reason=(
                f"Only 5' primer detected but N-length confidence "
                f"({n_length_confidence:.2f}) ≤ threshold "
                f"({n_length_confident_threshold:.2f})"
            ),
        )

    # Rows 4 + 5: only 3' strong.
    if match_rate_3p > primer_found_threshold:
        if n_length_confidence > n_length_confident_threshold:
            return Classification(
                extraction_mode="THREE_PRIME_ONLY",
                full_insert_recovered=False,
                required_action="NONE",
                status=_assign_status(
                    extraction_mode="THREE_PRIME_ONLY",
                    composite_confidence=composite_confidence,
                    has_round_map=has_round_map,
                    status_high_cutoff=status_high_cutoff,
                    status_medium_cutoff=status_medium_cutoff,
                    status_low_cutoff=status_low_cutoff,
                ),
                failure_reason=None,
            )
        return Classification(
            extraction_mode="UNABLE_TO_EXTRACT",
            full_insert_recovered=False,
            required_action="MANUAL_PRIMERS_REQUIRED",
            status="UNABLE_TO_INFER",
            failure_reason=(
                f"Only 3' primer detected but N-length confidence "
                f"({n_length_confidence:.2f}) ≤ threshold "
                f"({n_length_confident_threshold:.2f})"
            ),
        )

    # Ambiguous middle ground: both rates between unable_to_extract (0.4)
    # and primer_found (0.7). Conservative: refuse to extract.
    return Classification(
        extraction_mode="UNABLE_TO_EXTRACT",
        full_insert_recovered=False,
        required_action="MANUAL_PRIMERS_REQUIRED",
        status="UNABLE_TO_INFER",
        failure_reason=(
            f"Ambiguous primer evidence: 5'={match_rate_5p:.2f}, "
            f"3'={match_rate_3p:.2f}, neither side passes the "
            f"{primer_found_threshold:.2f} 'primer found' threshold"
        ),
    )


def _assign_status(
    *,
    extraction_mode: ExtractionMode,
    composite_confidence: float,
    has_round_map: bool,
    status_high_cutoff: float,
    status_medium_cutoff: float,
    status_low_cutoff: float,
) -> Status:
    """Map composite confidence to a Status, applying the no-round-map cap.

    the design: "If no round map is available, cross-round
    persistence cannot run → ``status`` is capped at ``MEDIUM``."
    """
    if extraction_mode == "UNABLE_TO_EXTRACT":
        return "UNABLE_TO_INFER"

    raw: Status
    if composite_confidence >= status_high_cutoff:
        raw = "HIGH"
    elif composite_confidence >= status_medium_cutoff:
        raw = "MEDIUM"
    elif composite_confidence >= status_low_cutoff:
        raw = "LOW"
    else:
        raw = "UNABLE_TO_INFER"

    if not has_round_map and raw == "HIGH":
        return "MEDIUM"
    return raw


# ---------------------------------------------------------------------------
# Deterministic JSON I/O — bit-identical output across reruns.
# ---------------------------------------------------------------------------


def _to_deterministic_payload(report: LibraryReport) -> dict[str, object]:
    """Convert a LibraryReport to a dict ready for deterministic JSON dump.

    - Top-level field order follows pydantic's declaration order (stable
      across runs because pydantic v2 preserves declaration order).
    - ``n_length_distribution`` keys are pre-sorted numerically (otherwise
      ``json.dumps(sort_keys=True)`` would sort them lexically: "10",
      "100", "20", … — deterministic but unreadable).
    - ``variants_5p`` / ``variants_3p`` retain their input order (caller's
      responsibility to pre-rank).
    """
    payload: dict[str, object] = report.model_dump(mode="json")

    # Numeric sort for the int-keyed dict. After mode="json" the keys are
    # stringified ints; sort by their int value.
    raw_nld = payload.get("n_length_distribution")
    if isinstance(raw_nld, dict):
        sorted_items = sorted(raw_nld.items(), key=lambda kv: int(kv[0]))
        payload["n_length_distribution"] = dict(sorted_items)

    # Sort known_adapter_hits by name (stable alphabetical).
    raw_kah = payload.get("known_adapter_hits")
    if isinstance(raw_kah, dict):
        payload["known_adapter_hits"] = dict(sorted(raw_kah.items()))

    return payload


def write_library_report_json(report: LibraryReport, path: Path) -> str:
    """Write `report` to `path` as deterministic JSON; return sha256 hex.

    Output discipline (so the manifest's ``output_sha256`` is reproducible):

    - UTF-8 encoded, ``ensure_ascii=False`` (escapes only what JSON requires)
    - 2-space indent
    - trailing newline
    - top-level fields in pydantic declaration order
    - ``n_length_distribution`` in numeric key order
    - ``known_adapter_hits`` in alphabetical key order

    Note: this does NOT use ``json.dumps(sort_keys=True)`` because that
    would re-sort nested dicts lexically and undo the numeric sort above.
    """
    payload = _to_deterministic_payload(report)
    text = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    path.write_text(text, encoding="utf-8")
    return sha256(text.encode("utf-8")).hexdigest()


def read_library_report_json(path: Path) -> LibraryReport:
    """Load and validate a LibraryReport from disk.

    Pydantic v2 coerces the JSON's stringified int keys in
    ``n_length_distribution`` back to ``int``.
    """
    return LibraryReport.model_validate_json(path.read_text(encoding="utf-8"))
