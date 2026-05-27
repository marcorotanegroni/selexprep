"""Deterministic 5-level cascade for extracting SELEX round number from SRA/ENA metadata.

Priority order (first match wins):
  L1 — sample_attributes structured key/value     (confidence=HIGH)
  L2 — sample_title                                (confidence=HIGH or MEDIUM)
  L3 — library_name, experiment_title,
       design_description                          (confidence=MEDIUM)
  L4 — BioProject abstract (count only —
       NEVER assigns SRR→round)                    (confidence=LOW, informative)
  L5 — Unknown                                     (confidence=NONE)
       → ``round_number is None`` (caught by
         ``RoundRecord.is_unassigned``)
       → dump written to metadata/manual_review/

**Cardinal rule: never guess.** A catalog with 20% unknowns is better than one
with invented round assignments that contaminate enrichment signals downstream.
"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class RoundRecord:
    srr: str
    round_number: int | None
    confidence: str  # HIGH | MEDIUM | LOW | NONE
    source_field: str
    matched_pattern: str
    round_candidates: list[int] = field(default_factory=list)
    parser_notes: str = ""
    target_hint: str | None = None

    @property
    def is_unassigned(self) -> bool:
        """The run cannot safely contribute to ``rounds.tsv``.

        True iff either:

        - the cascade found no parseable round number
          (``round_number is None``), or
        - multiple distinct round numbers were parsed from the same text
          (``len(round_candidates) > 1``) — genuine ambiguity that
          needs a curator's eyeball before trusting any of them.

        Phase 6b.4 audit fix: replaces the older ``needs_manual_review``
        boolean field, which was set to ``True`` for every L3
        single-match parse (`base_confidence="MEDIUM"`) and caused
        ``FetchPlan.none_confidence_runs`` to refuse fetch on
        accessions where ``library_name="RAPT26-2R"`` parsed cleanly
        to round 2. The property makes the criterion explicit at
        every call site — there is no longer a per-record flag whose
        meaning could drift.
        """
        return self.round_number is None or len(self.round_candidates) > 1

    def to_dict(self) -> dict:
        d = asdict(self)
        d["round_candidates"] = json.dumps(d["round_candidates"])
        # Preserve the `needs_manual_review` column on serialized
        # ``rounds.csv`` rows so existing readers
        # (`selexprep.fetch.download.needs_manual_review` and curators
        # hand-editing `rounds_curated.csv`) keep working. Computed from
        # ``is_unassigned`` — the property is the single source of truth
        # on the dataclass; the CSV column is a derived view kept for
        # backward compat with the discovery → curation → download flow.
        d["needs_manual_review"] = self.is_unassigned
        return d


@dataclass
class BioprojectRoundSummary:
    bioproject_id: str
    n_srr: int
    n_high: int
    n_medium: int
    n_low: int
    n_none: int
    inconsistent_annotation: bool
    n_rounds_from_abstract: int | None


# ---------------------------------------------------------------------------
# Regex catalogue (ordered by specificity — first match wins at L2/L3)
# ---------------------------------------------------------------------------

_ROUND_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("round_word_digit", re.compile(r"(?:^|[\s_\-./])(?:[Rr]ound)[\s_\-]*(\d+)(?:[\s_\-./]|$)")),
    ("cycle_word_digit", re.compile(r"(?:^|[\s_\-./])[Cc]ycle[\s_\-]*(\d+)(?:[\s_\-./]|$)")),
    ("iteration_digit", re.compile(r"(?:^|[\s_\-./])[Ii]teration[\s_\-]*(\d+)(?:[\s_\-./]|$)")),
    ("pool_digit", re.compile(r"(?:^|[\s_\-./])[Pp]ool[\s_\-]*(\d+)(?:[\s_\-./]|$)")),
    ("R_digit_boundary", re.compile(r"(?:^|[\s_\-./])[Rr](\d+)(?:[\s_\-./]|$)")),
    ("C_digit_boundary", re.compile(r"(?:^|[\s_\-./])[Cc](\d+)(?:[\s_\-./]|$)")),
    ("digit_R_suffix", re.compile(r"(?:^|[\s_\-./])(\d+)[Rr](?:[\s_\-./]|$)")),
    # Phase 6b.5c — empirical patterns surfaced by the Phase 6b.4 audit
    # pilot's per-accession metadata inspection. Each was a SELEX deposit
    # the original cascade missed.
    #
    # ``RV01`` / ``RV02`` / ... (PRJEB51212 — Sall4 SELEX rounds). The R
    # and V are glued, so ``R_digit_boundary`` ('R immediately followed
    # by digit') doesn't match; this pattern handles the RV-glued case
    # explicitly.
    ("RV_digit", re.compile(r"(?:^|[\s_\-./])RV(\d+)(?:[\s_\-./]|$)")),
    # ``DNAFOXR00`` / ``DNAFOXR01`` / ... (PRJEB62756 — DNAFOX SELEX
    # series). A word prefix glued to ``R\d+`` with no separator.
    # Requires the prefix to be ≥3 letters (uppercase + ≥2 more letters)
    # to limit false positives from short protein abbreviations like
    # ``TFR1`` (transcription factor R1). Catches DNAFOX (6 letters) but
    # not TF (2).
    ("glued_word_R_digit", re.compile(r"[A-Z][A-Za-z]{2,}R(\d+)(?:[\s_\-./]|$)")),
    # ``Selex_7_cyc`` / ``7_cyc`` / ``7 cycles`` (PRJNA385825). Digit
    # BEFORE the cyc(le) token; complements ``cycle_word_digit`` which
    # only matches ``cycle 5`` (cycle BEFORE digit).
    (
        "digit_cyc_suffix",
        re.compile(
            # cyc / cycle / cycles — handles singular + plural
            r"(?:^|[\s_\-./])(\d+)[\s_\-]*cyc(?:les?)?(?:[\s_\-./]|$)",
            re.IGNORECASE,
        ),
    ),
]

# L4 — abstract round count (informative, never assigns SRR→round)
_ABSTRACT_COUNT_PATTERNS: list[re.Pattern[str]] = [
    re.compile(r"(\d+)\s+rounds?\s+of\s+(?:selection|SELEX)", re.IGNORECASE),
    re.compile(r"rounds?\s+(?:1|one)\s+(?:through|to|-)\s+(\d+)", re.IGNORECASE),
    re.compile(r"R0\s+(?:through|to|-)\s+R(\d+)", re.IGNORECASE),
    re.compile(r"(\d+)\s+(?:selection\s+)?cycles?", re.IGNORECASE),
]

# Keys identifying round number in sample_attributes (L1)
_ATTR_ROUND_KEYS: re.Pattern[str] = re.compile(
    r"^(?:selex_)?(?:round|cycle|iteration|r|selection_round|selection_cycle)(?:_num(?:ber)?)?$",
    re.IGNORECASE,
)

# Target hint — first capitalised word before a round indicator
_TARGET_HINT_PATTERN: re.Pattern[str] = re.compile(
    r"^([A-Z][A-Za-z0-9_\-]+)[\s_\-]+(?:[Rr]ound|[Rr]\d|[Cc]ycle)",
)


# ---------------------------------------------------------------------------
# Core parser
# ---------------------------------------------------------------------------


def parse_round(
    srr: str,
    sample_title: str = "",
    library_name: str = "",
    experiment_title: str = "",
    design_description: str = "",
    sample_attributes: dict[str, str] | None = None,
    abstract: str = "",
    manual_review_dir: Path | None = None,
    bioproject_id: str = "",
) -> RoundRecord:
    """Run the 5-level cascade and return a `RoundRecord`.

    Parameters
    ----------
    srr
        SRR accession string.
    sample_title, library_name, experiment_title, design_description
        SRA/ENA metadata text fields.
    sample_attributes
        Dict of custom BioSample attributes (key → value strings).
    abstract
        BioProject abstract — used for L4 informative count only.
    manual_review_dir
        If set, write a dump file for NONE-confidence records.
    bioproject_id
        Parent BioProject accession — groups manual-review dumps by BP
        (filename becomes ``{bp_id}_{srr}.txt``).
    """
    if sample_attributes is None:
        sample_attributes = {}

    all_metadata = {
        "srr": srr,
        "bioproject_id": bioproject_id,
        "sample_title": sample_title,
        "library_name": library_name,
        "experiment_title": experiment_title,
        "design_description": design_description,
        "sample_attributes": sample_attributes,
        "abstract_excerpt": abstract[:500],
    }

    # L1 — structured sample_attributes
    for key, val in sample_attributes.items():
        if _ATTR_ROUND_KEYS.match(key.strip()):
            stripped = val.strip()
            if re.fullmatch(r"\d+", stripped):
                round_num = int(stripped)
                return RoundRecord(
                    srr=srr,
                    round_number=round_num,
                    confidence="HIGH",
                    source_field="sample_attributes",
                    matched_pattern="structured_attr",
                    round_candidates=[round_num],
                    parser_notes=f"structured attribute key='{key}' value='{stripped}'",
                    target_hint=_extract_target_hint(sample_title),
                )

    # L2 — sample_title
    result = _match_text_field(srr, "sample_title", sample_title)
    if result is not None:
        result.target_hint = _extract_target_hint(sample_title)
        return result

    # L3 — library_name, experiment_title, design_description
    for field_name, text in [
        ("library_name", library_name),
        ("experiment_title", experiment_title),
        ("design_description", design_description),
    ]:
        result = _match_text_field(srr, field_name, text, base_confidence="MEDIUM")
        if result is not None:
            result.target_hint = _extract_target_hint(sample_title)
            return result

    # L4 — abstract count (informative only, NEVER assigns round)
    n_rounds_abstract = _extract_round_count_from_abstract(abstract)
    notes = ""
    if n_rounds_abstract is not None:
        notes = (
            f"abstract declares ~{n_rounds_abstract} rounds total (informative only, not assigned)"
        )

    # L5 — Unknown → manual review (round_number is None makes
    # ``RoundRecord.is_unassigned`` True; no per-record flag needed)
    record = RoundRecord(
        srr=srr,
        round_number=None,
        confidence="NONE",
        source_field="none",
        matched_pattern="none",
        round_candidates=[],
        parser_notes=notes or "no round indicator found in any metadata field",
        target_hint=_extract_target_hint(sample_title),
    )
    if manual_review_dir is not None:
        _write_manual_review_dump(record, all_metadata, manual_review_dir)
    return record


# ---------------------------------------------------------------------------
# L4 helper
# ---------------------------------------------------------------------------


def extract_round_count_from_abstract(abstract: str) -> int | None:
    """Public wrapper — used by `summarise_bioproject`."""
    return _extract_round_count_from_abstract(abstract)


def _extract_round_count_from_abstract(abstract: str) -> int | None:
    if not abstract:
        return None
    for pat in _ABSTRACT_COUNT_PATTERNS:
        m = pat.search(abstract)
        if m:
            return int(m.group(1))
    return None


# ---------------------------------------------------------------------------
# Text field matching
# ---------------------------------------------------------------------------


def _match_text_field(
    srr: str,
    field_name: str,
    text: str,
    base_confidence: str = "HIGH",
) -> RoundRecord | None:
    """Apply all patterns to `text`.

    Returns a `RoundRecord` if any pattern matches, `None` otherwise. Confidence
    downgrades to ``MEDIUM`` if multiple conflicting numbers are found.
    """
    if not text or not text.strip():
        return None

    found: list[tuple[str, int]] = []
    for pat_name, pat in _ROUND_PATTERNS:
        for m in pat.finditer(text):
            n = int(m.group(1))
            found.append((pat_name, n))

    if not found:
        return None

    unique_numbers = list(dict.fromkeys(n for _, n in found))
    unique_patterns = list(dict.fromkeys(p for p, _ in found))

    if len(unique_numbers) == 1:
        confidence = base_confidence
        notes = f"matched '{unique_patterns[0]}' in {field_name}: '{text[:120]}'"
    else:
        # Phase 6b.4 audit fix: the only legitimate "needs review" case
        # in this helper is genuine ambiguity (multiple distinct round
        # numbers from the same text). ``RoundRecord.is_unassigned``
        # picks this up via ``len(round_candidates) > 1`` — no per-record
        # boolean flag required. Pre-fix the field was set unconditionally
        # for MEDIUM and lumped single-match L3 parses into refusal
        # alongside genuine ambiguity, which is the bug this refactor
        # closes at source.
        confidence = "MEDIUM"
        notes = (
            f"conflicting candidates {unique_numbers} from patterns "
            f"{unique_patterns} in {field_name}: '{text[:120]}'"
        )

    return RoundRecord(
        srr=srr,
        round_number=unique_numbers[0],
        confidence=confidence,
        source_field=field_name,
        matched_pattern=unique_patterns[0],
        round_candidates=unique_numbers,
        parser_notes=notes,
    )


# ---------------------------------------------------------------------------
# Target hint
# ---------------------------------------------------------------------------


def _extract_target_hint(sample_title: str) -> str | None:
    if not sample_title:
        return None
    m = _TARGET_HINT_PATTERN.match(sample_title.strip())
    if m:
        return m.group(1)
    return None


# ---------------------------------------------------------------------------
# Manual review dump
# ---------------------------------------------------------------------------


def _write_manual_review_dump(
    record: RoundRecord,
    all_metadata: dict[str, Any],
    manual_review_dir: Path,
) -> None:
    manual_review_dir.mkdir(parents=True, exist_ok=True)
    bioproject = all_metadata.get("bioproject_id") or "unknown"
    out = manual_review_dir / f"{bioproject}_{record.srr}.txt"
    lines = [
        "=== MANUAL REVIEW REQUIRED ===",
        f"SRR: {record.srr}",
        f"BioProject: {bioproject}",
        f"Parser notes: {record.parser_notes}",
        "",
        "--- Metadata fields (search for round/cycle/R[N] indicators) ---",
    ]
    for k, v in all_metadata.items():
        if k == "sample_attributes":
            lines.append("sample_attributes:")
            for ak, av in v.items():
                lines.append(f"  {ak}: {av}")
        else:
            lines.append(f"{k}: {v}")
    lines += [
        "",
        "--- Action ---",
        "Fill in round_number below and change confidence to 'manual'.",
        "Then copy to rounds_curated.csv.",
        "",
        "round_number: ???",
        "confidence: manual",
        "source_field: manual",
        "parser_notes: <your notes>",
    ]
    out.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# BioProject-level summary
# ---------------------------------------------------------------------------


def summarise_bioproject(
    bioproject_id: str,
    records: list[RoundRecord],
    abstract: str = "",
) -> BioprojectRoundSummary:
    """Compute per-BioProject summary stats and flag inconsistent annotations."""
    n_high = sum(1 for r in records if r.confidence == "HIGH")
    n_medium = sum(1 for r in records if r.confidence == "MEDIUM")
    n_low = sum(1 for r in records if r.confidence == "LOW")
    n_none = sum(1 for r in records if r.confidence == "NONE")

    has_assigned = any(r.round_number is not None for r in records)
    has_missing = any(r.round_number is None for r in records)
    inconsistent = has_assigned and has_missing

    n_rounds_abstract = _extract_round_count_from_abstract(abstract)

    return BioprojectRoundSummary(
        bioproject_id=bioproject_id,
        n_srr=len(records),
        n_high=n_high,
        n_medium=n_medium,
        n_low=n_low,
        n_none=n_none,
        inconsistent_annotation=inconsistent,
        n_rounds_from_abstract=n_rounds_abstract,
    )


# ---------------------------------------------------------------------------
# Seed override support
# ---------------------------------------------------------------------------


def apply_seed_overrides(
    records: list[RoundRecord],
    manual_round_mapping: dict[str, int],
) -> list[RoundRecord]:
    """Apply manual `srr → round` overrides from a seed configuration.

    Overrides take precedence over all parser levels and produce ``HIGH``
    confidence with ``source_field='seed_override'``.
    """
    mapping = {str(k): int(v) for k, v in manual_round_mapping.items()}
    updated = []
    for r in records:
        if r.srr in mapping:
            r = RoundRecord(
                srr=r.srr,
                round_number=mapping[r.srr],
                confidence="HIGH",
                source_field="seed_override",
                matched_pattern="manual",
                round_candidates=[mapping[r.srr]],
                parser_notes=f"overridden by seed mapping: round={mapping[r.srr]}",
                target_hint=r.target_hint,
            )
        updated.append(r)
    return updated
