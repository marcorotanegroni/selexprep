"""Primer-equivalence rules for benchmarking inferred primers against paper-reported truth.

Implements the four locked-plan equivalence rules (line 364):

1. **Reverse-complement**  — many papers report the RC of what the
   sequencer sees; an observed primer that revcomp-matches truth is
   equivalent.
2. **U/T normalization** — RNA primers reported as U-containing sequences
   match T-containing observed (and vice versa).
3. **Barcode-prefix stripping** — sample-multiplex barcodes prepended to
   the 5' primer can mask the underlying primer; the curator may supply
   barcode strings to strip before comparison.
4. **IUPAC ambiguity rejection** — locked plan line 33 says "IUPAC
   ambiguous bases unsupported in v0.1 (counted separately)". A truth
   primer containing ``N``/``R``/``Y``/etc. is not silently fuzzy-matched;
   the result is ``IUPAC_UNSUPPORTED`` so the metric aggregator can
   count it in a dedicated bucket.

The public entry point :func:`primer_equivalent` returns a structured
:class:`EquivalenceResult` so callers (the metric aggregator + the Figure
A plot) can distinguish *how* a match was achieved — not just *whether*.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Literal

from selexprep.library.adapters import reverse_complement

EquivalenceKind = Literal[
    "EXACT",
    "REVCOMP",
    "U_T_NORMALIZED",
    "BARCODE_STRIPPED",
    "PARTIAL_5P",
    "PARTIAL_3P",
    "MISMATCH",
    "IUPAC_UNSUPPORTED",
]

# IUPAC ambiguity codes that v0.1 does not support.
_IUPAC_AMBIGUOUS: re.Pattern[str] = re.compile(r"[NRYSWKMBDHV]", re.IGNORECASE)


@dataclass(frozen=True)
class EquivalenceResult:
    """Outcome of a single observed-vs-truth primer comparison.

    ``matched`` is ``True`` for every kind except ``MISMATCH``,
    ``IUPAC_UNSUPPORTED``, and the (rare) all-empty input case.
    Callers should branch on ``equivalence_kind`` to render or
    aggregate the result, NOT on a string match against ``notes``.
    """

    matched: bool
    equivalence_kind: EquivalenceKind
    notes: str = ""


def _normalize_ut(seq: str) -> str:
    """Replace every U with T (case-preserving uppercase output)."""
    return seq.upper().replace("U", "T")


def _is_iupac_ambiguous(seq: str) -> bool:
    return bool(_IUPAC_AMBIGUOUS.search(seq))


def _strip_barcode_prefixes(observed: str, barcodes: tuple[str, ...]) -> tuple[str, str | None]:
    """Strip any matching barcode prefix from ``observed``.

    Returns ``(stripped, matched_barcode)``. ``matched_barcode`` is ``None``
    if no barcode prefix matched (in which case ``stripped == observed``).
    The longest matching barcode wins to avoid a 3-bp barcode prefix-
    matching when a 6-bp barcode is also present.
    """
    if not barcodes:
        return observed, None
    for bc in sorted(barcodes, key=len, reverse=True):
        if observed.upper().startswith(bc.upper()):
            return observed[len(bc) :], bc
    return observed, None


def _is_partial_5p(observed: str, truth: str) -> bool:
    """Observed is a (strict) prefix of truth, or truth is a strict prefix of observed.

    Either direction signals that one side captured only the 5' portion
    of the other — actionable for primer-recovery accounting.
    """
    o, t = observed.upper(), truth.upper()
    if o == t:
        return False
    return t.startswith(o) or o.startswith(t)


def _is_partial_3p(observed: str, truth: str) -> bool:
    """Observed and truth share a common 3' suffix without being equal."""
    o, t = observed.upper(), truth.upper()
    if o == t:
        return False
    return t.endswith(o) or o.endswith(t)


def primer_equivalent(
    observed: str | None,
    truth: str,
    *,
    allow_revcomp: bool = True,
    allow_ut: bool = True,
    strip_barcodes: tuple[str, ...] = (),
) -> EquivalenceResult:
    """Compare an observed primer against the paper-reported truth.

    Parameters
    ----------
    observed
        The primer ``selexprep`` inferred. ``None`` is treated as a
        ``MISMATCH`` (selexprep failed to recover the primer at all).
    truth
        The paper-reported primer sequence.
    allow_revcomp
        Accept reverse-complement equivalence. Default ``True`` (locked
        plan line 364).
    allow_ut
        Accept U↔T normalization. Default ``True``.
    strip_barcodes
        Tuple of barcode prefix sequences to strip from ``observed``
        before comparison. Empty by default.

    Order of checks (first match wins, so kinds are mutually exclusive):

    1. ``IUPAC_UNSUPPORTED`` — truth contains an IUPAC ambiguity code.
    2. ``MISMATCH`` — observed is ``None`` or truth is empty/blank.
    3. ``EXACT`` — case-folded exact match.
    4. ``BARCODE_STRIPPED`` — exact match after stripping a barcode prefix.
    5. ``U_T_NORMALIZED`` — exact match after U→T normalization.
    6. ``REVCOMP`` — observed (or its U-normalized form) matches the
       reverse complement of truth.
    7. ``PARTIAL_5P`` / ``PARTIAL_3P`` — overlapping prefix / suffix.
    8. ``MISMATCH`` — no rule fired.
    """
    if _is_iupac_ambiguous(truth):
        return EquivalenceResult(
            matched=False,
            equivalence_kind="IUPAC_UNSUPPORTED",
            notes=f"truth {truth!r} contains IUPAC ambiguity (v0.1 unsupported)",
        )

    truth_norm = truth.strip().upper()
    if not truth_norm:
        return EquivalenceResult(matched=False, equivalence_kind="MISMATCH", notes="empty truth")
    if observed is None:
        return EquivalenceResult(
            matched=False, equivalence_kind="MISMATCH", notes="observed is None"
        )

    observed_norm = observed.strip().upper()
    if not observed_norm:
        return EquivalenceResult(matched=False, equivalence_kind="MISMATCH", notes="empty observed")

    # 1. EXACT.
    if observed_norm == truth_norm:
        return EquivalenceResult(matched=True, equivalence_kind="EXACT")

    # 2. BARCODE_STRIPPED (exact match after removing a known barcode prefix).
    stripped, matched_bc = _strip_barcode_prefixes(observed_norm, strip_barcodes)
    if matched_bc is not None and stripped.upper() == truth_norm:
        return EquivalenceResult(
            matched=True,
            equivalence_kind="BARCODE_STRIPPED",
            notes=f"stripped barcode prefix {matched_bc!r}",
        )

    # 3. U_T_NORMALIZED.
    if allow_ut:
        obs_ut = _normalize_ut(observed_norm)
        truth_ut = _normalize_ut(truth_norm)
        if obs_ut == truth_ut and obs_ut != observed_norm:
            return EquivalenceResult(
                matched=True,
                equivalence_kind="U_T_NORMALIZED",
                notes="match after U→T normalization",
            )
        # Re-evaluate the EXACT check on the U-T normalized forms:
        # papers can report RNA primers with U, observed comes as T.
        if obs_ut == truth_ut:
            return EquivalenceResult(
                matched=True,
                equivalence_kind="U_T_NORMALIZED",
                notes="match after U→T normalization",
            )

    # 4. REVCOMP.
    if allow_revcomp:
        try:
            truth_rc = reverse_complement(truth_norm if not allow_ut else _normalize_ut(truth_norm))
        except ValueError:
            # Truth has chars revcomp can't handle (shouldn't happen — IUPAC
            # was caught upstream — but be defensive).
            truth_rc = ""
        obs_for_rc = _normalize_ut(observed_norm) if allow_ut else observed_norm
        if truth_rc and obs_for_rc == truth_rc:
            return EquivalenceResult(
                matched=True,
                equivalence_kind="REVCOMP",
                notes="match after reverse-complementing truth",
            )

    # 5. PARTIAL_5P / PARTIAL_3P. Compare on U-T normalized forms when allowed
    # so a partial-RNA-primer doesn't fall through to MISMATCH.
    obs_for_partial = _normalize_ut(observed_norm) if allow_ut else observed_norm
    truth_for_partial = _normalize_ut(truth_norm) if allow_ut else truth_norm
    if _is_partial_5p(obs_for_partial, truth_for_partial):
        return EquivalenceResult(
            matched=False,
            equivalence_kind="PARTIAL_5P",
            notes="observed and truth share a 5' prefix but are not equal",
        )
    if _is_partial_3p(obs_for_partial, truth_for_partial):
        return EquivalenceResult(
            matched=False,
            equivalence_kind="PARTIAL_3P",
            notes="observed and truth share a 3' suffix but are not equal",
        )

    return EquivalenceResult(matched=False, equivalence_kind="MISMATCH")
