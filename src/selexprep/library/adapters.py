"""Sequencing-adapter blacklist used by primer inference.

The locked plan (``~/.claude/plans/unified-seeking-treehouse.md`` line 291)
requires the LibraryReport pipeline to **record** how often known
sequencing adapters appear in the read pool, and to **exclude** those
adapters from primer candidates. It does NOT filter reads — under-trimmed
adapter sequences are diagnostic information, not error conditions.

**v0.1 set (conservative).** Locked plan line 291 says "TruSeq R1/R2,
Nextera, etc." — the v0.1 set is the two most common Illumina adapters
encountered in HT-SELEX deposits. Expanding to Illumina P5/P7, Small RNA,
and IonTorrent A is deferred to v0.2.

**Calibration status.** The exact composition is a `# CALIBRATION-TODO`:
the conservative set ships now; Codex review (rate-limited 2026-05-19 →
2026-05-26) gets the final list along with Phase 6 benchmark hit-rate
data.
"""

from __future__ import annotations

# CALIBRATION-TODO: locked plan line 291 ("TruSeq R1/R2, Nextera, etc.");
# Codex confirms or extends to the full Illumina set.
KNOWN_ADAPTERS: dict[str, str] = {
    # Illumina TruSeq Read 1 adapter prefix — the canonical contamination
    # probe in HT-SELEX deposits (see audit.py for the identical constant
    # used as the audit-module's standalone probe).
    "TRUSEQ_R1": "AGATCGGAAGAGC",
    # Nextera / Tn5 transposase adapter mosaic end — appears in Nextera
    # tagmentation libraries when extraction misses the constant end.
    "NEXTERA": "CTGTCTCTTATACACATCT",
}


# Watson-Crick complement table.
# U → A (RNA primers are reported as DNA per locked plan line 296; see
# `selexprep.library.detect._normalize_u_to_t`).
# IUPAC ambiguous bases (N, R, Y, ...) are explicitly unsupported in v0.1
# per locked plan line 33 ("IUPAC unsupported in v0.1").
_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "U": "A"}


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of an ACGTU sequence.

    Case-insensitive on input; returns uppercase. Raises ``ValueError`` on
    any non-ACGTU base (including ``N`` and other IUPAC codes — v0.1 does
    not support ambiguous bases).
    """
    upper = seq.upper()
    try:
        return "".join(_COMPLEMENT[b] for b in reversed(upper))
    except KeyError as e:
        bad = e.args[0]
        raise ValueError(
            f"reverse_complement: non-ACGTU base {bad!r} in {seq!r}; "
            "IUPAC ambiguity is unsupported in v0.1"
        ) from None


# Pre-compute the reverse complements once at module load so every call to
# count_adapter_hits is a cheap substring scan, not a per-call revcomp.
KNOWN_ADAPTERS_RC: dict[str, str] = {
    name: reverse_complement(adapter) for name, adapter in KNOWN_ADAPTERS.items()
}


# Default probe length: first ``k`` bp of an adapter is the canonical
# "matches the start of an adapter" signal used by both
# ``count_adapter_hits`` (substring scan over all reads) and
# ``matches_known_adapter_prefix`` (single-primer check). Colocated here
# so the two helpers stay in sync.
ADAPTER_PROBE_K = 13


def matches_known_adapter_prefix(
    primer: str | None,
    k: int = ADAPTER_PROBE_K,
) -> bool:
    """True if ``primer``'s first ``k`` bases match a known adapter or its RC.

    Used by primer-inference (``library/detect``) to drop adapter
    candidates after the auto-detection step, and by the extract runner
    (``extract/runner``) to warn when a manually-supplied
    ``--override-primer-*`` value matches a known sequencing adapter
    prefix. Returns False for ``None`` (no primer detected / no override
    given).
    """
    if not primer:
        return False
    probe = primer.upper()[:k]
    return any(adapter[:k] == probe for adapter in KNOWN_ADAPTERS.values()) or any(
        adapter_rc[:k] == probe for adapter_rc in KNOWN_ADAPTERS_RC.values()
    )


def count_adapter_hits(seqs: list[str], k: int = ADAPTER_PROBE_K) -> dict[str, int]:
    """Count sequences containing each known adapter (forward or RC) as a substring.

    For each adapter in :data:`KNOWN_ADAPTERS`, count how many sequences in
    ``seqs`` contain either the adapter's first ``k`` bases or its RC's
    first ``k`` bases as a substring. Sequences are scanned at full length
    (not just the flanks) because residual adapter readthrough can appear
    mid-read after low-quality trimming.

    Returns a dict mapping adapter name → number of sequences with at least
    one hit (forward OR RC; not summed). A sequence with both forward and
    RC of the same adapter counts as one hit for that adapter.

    Default ``k=13`` matches TruSeq R1's full length; shorter probes would
    cause false positives on random pools.
    """
    if k <= 0:
        raise ValueError(f"count_adapter_hits: k must be positive, got {k!r}")

    probes: dict[str, tuple[str, str]] = {}
    for name, adapter in KNOWN_ADAPTERS.items():
        fwd = adapter[:k]
        rev = KNOWN_ADAPTERS_RC[name][:k]
        probes[name] = (fwd, rev)

    hits: dict[str, int] = dict.fromkeys(KNOWN_ADAPTERS, 0)
    for seq in seqs:
        upper = seq.upper()
        for name, (fwd, rev) in probes.items():
            if fwd in upper or rev in upper:
                hits[name] += 1
    return hits
