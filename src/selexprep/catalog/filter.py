"""Filter helpers for the discovery catalog DataFrame.

All filters are case-insensitive substring matches where applicable and
compose by stacking calls (each filter narrows the previous result). Callers
can also operate on the underlying DataFrame directly when these helpers
aren't expressive enough.
"""

from __future__ import annotations

import re

import pandas as pd

_INSDC_PREFIX_RE = re.compile(r"^(PRJ[NDE][ABM]|SRP|ERP|DRP)")


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
        mask = out["bioproject_id"].fillna("").str.match(_INSDC_PREFIX_RE, na=False)
        out = out[mask]

    return out
