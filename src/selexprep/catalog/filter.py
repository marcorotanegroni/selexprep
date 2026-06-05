"""Filter helpers for the discovery catalog DataFrame.

All filters are case-insensitive substring matches where applicable and
compose by stacking calls (each filter narrows the previous result). Callers
can also operate on the underlying DataFrame directly when these helpers
aren't expressive enough.
"""

from __future__ import annotations

import re

import pandas as pd

INSDC_PREFIX_RE = re.compile(r"^(PRJ[NDE][ABM]|SRP|ERP|DRP)")
"""Regex matching INSDC BioProject prefixes (PRJNA / PRJDB / PRJEB / SRP / ERP / DRP).

Public so the eligibility classifier and audit aggregator agree on what
counts as INSDC without each maintaining its own copy.
"""


def is_insdc_accession(accession: str) -> bool:
    """True iff ``accession`` matches the INSDC BioProject prefix scheme.

    Non-INSDC catalog rows (``figshare:*``, ``zenodo:*``, ``utexas:*``) are
    catalog-only references to published HT-SELEX data and cannot be
    fetched through ENA / SRA / DDBJ as raw FASTQ.
    """
    return bool(INSDC_PREFIX_RE.match(accession or ""))


_DISCOVERY_ONLY_PREFIXES = ("zenodo:", "figshare:", "utexas:")


def is_discovery_only(accession: str) -> bool:
    """True iff ``accession`` is a non-fetchable discovery-only catalog pointer.

    These rows (``zenodo:*``, ``figshare:*``, ``utexas:*``) reference published
    HT-SELEX data hosted outside INSDC, so ``fetch`` / ``run`` cannot retrieve
    them as raw FASTQ in v0.1. Unlike :func:`is_insdc_accession` (which only
    matches study-level prefixes), this is safe to call on run-level accessions
    (SRR/ERR/DRR) — they are not discovery-only.
    """
    return (accession or "").startswith(_DISCOVERY_ONLY_PREFIXES)


def filter_catalog(
    df: pd.DataFrame,
    target: str | None = None,
    organism: str | None = None,
    source_contains: str | None = None,
    min_rounds: int | None = None,
    insdc_only: bool = False,
) -> pd.DataFrame:
    """Apply common filters to a catalog DataFrame.

    Parameters
    ----------
    target
        Substring (case-insensitive) matched against ``protein_target``.
    organism
        Substring (case-insensitive) matched against ``target_organism``.
    source_contains
        Substring matched against ``source`` (e.g. ``"ena"`` keeps only
        ENA-discovered rows).
    min_rounds
        Drop rows where ``n_rounds_declared`` is missing or below this value.
    insdc_only
        If ``True``, keep only INSDC-format accessions
        (``PRJNA``/``PRJDB``/``PRJEB``/``SRP``/``ERP``/``DRP``) — drops
        processed-data deposits (``zenodo:*``, ``figshare:*``, ``utexas:*``)
        which can't be downloaded as raw FASTQ.
    """
    out = df

    if target:
        mask = out["protein_target"].fillna("").str.contains(target, case=False, na=False)
        out = out[mask]

    if organism:
        mask = out["target_organism"].fillna("").str.contains(organism, case=False, na=False)
        out = out[mask]

    if source_contains:
        mask = out["source"].fillna("").str.contains(source_contains, case=False, na=False)
        out = out[mask]

    if min_rounds is not None:
        nrounds = pd.to_numeric(out["n_rounds_declared"], errors="coerce").fillna(0)
        out = out[nrounds >= min_rounds]

    if insdc_only:
        mask = out["bioproject_id"].fillna("").str.match(INSDC_PREFIX_RE, na=False)
        out = out[mask]

    return out
