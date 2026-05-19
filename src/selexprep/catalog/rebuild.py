"""Rebuild the bundled discovery catalog from broad ENA queries.

The catalog ships as package data — but it would go stale fast without a way
to refresh it. This module provides that path, in two layers:

- :func:`harvest_studies_from_ena` runs a deliberately broader set of queries
  than the original `selex_corpus.discover` (which was tuned for one
  researcher's thesis). The wider net catches studies the thesis-specific
  queries would have missed.
- :func:`rebuild_catalog` writes a fresh ``bioprojects.csv`` from the
  ENA result + any preserved non-INSDC entries (Zenodo / Figshare /
  processed-data deposits) from a previous catalog snapshot. Hand-enriched
  fields on known entries (``protein_target``, ``paper_doi``, ``paper_pmid``,
  ``n_rounds_declared``) are merged forward so refreshes never erase manual
  enrichment.

**Curation is intentionally NOT in scope.** The catalog reflects the public
archives, not any single researcher's "include/exclude" decisions. Per
the package's design philosophy, curation is the user's job downstream.
"""

from __future__ import annotations

import csv
import json
import logging
import time
import urllib.parse
import urllib.request
from pathlib import Path

logger = logging.getLogger(__name__)


ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/search"

#: Broader ENA queries than the thesis-specific selex_corpus.discover.py.
#: Maximises INSDC coverage of public SELEX/aptamer studies at the cost of
#: some noise (an "aptamer" mention in a study_title doesn't strictly mean
#: an HT-SELEX raw-reads deposit). Users filter further at the
#: `selexprep catalog list` stage.
ENA_QUERIES: tuple[str, ...] = (
    'study_title="HT-SELEX"',
    'study_title="SELEX-seq"',
    'study_title="SELEX"',
    'study_title="aptamer"',
    'study_title="Cell-SELEX"',
    'study_title="aptamer selection"',
    'study_title="systematic evolution"',
    'study_title="DNA aptamer"',
    'study_title="RNA aptamer"',
    'description="HT-SELEX"',
    'description="SELEX-seq"',
    'description="aptamer selection"',
    'description="SELEX rounds"',
)

_STUDY_FIELDS = (
    "study_accession,study_title,study_description,scientific_name,"
    "secondary_study_accession,first_public,last_updated"
)

PUBLIC_COLS = (
    "bioproject_id",
    "source",
    "study_title",
    "protein_target",
    "target_organism",
    "paper_doi",
    "paper_pmid",
    "n_rounds_declared",
    "abstract",
)


def _ena_get(url: str, timeout: int = 60) -> list[dict]:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            return json.loads(resp.read())
    except Exception as e:
        logger.warning("ENA query failed: %s", e)
        return []


def harvest_studies_from_ena(
    queries: tuple[str, ...] = ENA_QUERIES,
    request_pause_seconds: float = 0.5,
) -> dict[str, dict]:
    """Run each broad query against ENA Portal, return {study_accession: meta}.

    Deduplicates across queries; first-seen wins. Studies surfaced under
    multiple queries are listed once.
    """
    studies: dict[str, dict] = {}
    for q in queries:
        params = {
            "result": "study",
            "query": q,
            "fields": _STUDY_FIELDS,
            "limit": "10000",
            "format": "json",
        }
        url = f"{ENA_PORTAL}?{urllib.parse.urlencode(params)}"
        logger.info("[ENA] %s", q)
        data = _ena_get(url)
        new_here = 0
        for r in data:
            acc = (r.get("study_accession") or "").strip()
            if not acc or acc in studies:
                continue
            studies[acc] = r
            new_here += 1
        logger.info("[ENA] %s → %d studies (%d new)", q, len(data), new_here)
        time.sleep(request_pause_seconds)
    return studies


def _enrichment_index(catalog_path: Path) -> dict[str, dict]:
    """Build a {bioproject_id: row} index of hand-enriched fields from an old catalog.

    Only rows with at least one non-empty enrichable field are kept. Used by
    `rebuild_catalog` to merge enrichment forward across refreshes.
    """
    enriched: dict[str, dict] = {}
    if not catalog_path.exists():
        return enriched
    with open(catalog_path) as f:
        for row in csv.DictReader(f):
            bp = (row.get("bioproject_id") or "").strip()
            if not bp:
                continue
            if any(
                row.get(c)
                for c in ("protein_target", "paper_doi", "paper_pmid", "n_rounds_declared")
            ):
                enriched[bp] = row
    return enriched


def _passthrough_non_insdc(catalog_path: Path, exclude_ids: set[str]) -> list[dict]:
    """Carry over non-ENA-INSDC entries (Zenodo / Figshare / utexas processed
    deposits) from a previous catalog so they survive the refresh.
    """
    out: list[dict] = []
    if not catalog_path.exists():
        return out
    with open(catalog_path) as f:
        for row in csv.DictReader(f):
            bp = (row.get("bioproject_id") or "").strip()
            if not bp or bp in exclude_ids:
                continue
            if bp.startswith(("zenodo:", "figshare:", "utexas:")):
                out.append({k: row.get(k, "") for k in PUBLIC_COLS})
    return out


def rebuild_catalog(
    out_path: Path,
    preserve_from: Path | None = None,
    queries: tuple[str, ...] = ENA_QUERIES,
) -> int:
    """Refresh the catalog CSV in-place. Returns the count of bioprojects written.

    Behaviour:

    - Run every query in `queries` against ENA, union the resulting studies.
    - Map each ENA study row to the catalog schema.
    - When `preserve_from` is given AND the bioproject was hand-enriched in
      the old catalog, copy ``protein_target``, ``paper_doi``, ``paper_pmid``,
      ``n_rounds_declared`` forward into the new row.
    - Append non-INSDC entries (Zenodo / Figshare / utexas processed-data
      deposits) from the old catalog so refreshes don't lose them.

    No curation flags are added or preserved (``include`` /
    ``manual_curation_notes`` / etc.) — curation is the user's downstream
    job, not the package's.
    """
    studies = harvest_studies_from_ena(queries=queries)
    enriched = _enrichment_index(preserve_from) if preserve_from else {}

    rows: list[dict] = []

    # ENA studies → catalog rows
    for acc, meta in studies.items():
        enr = enriched.get(acc, {})
        rows.append(
            {
                "bioproject_id": acc,
                "source": "ena",
                "study_title": (meta.get("study_title") or "").strip(),
                "protein_target": (enr.get("protein_target") or "").strip(),
                "target_organism": (meta.get("scientific_name") or "").strip(),
                "paper_doi": (enr.get("paper_doi") or "").strip(),
                "paper_pmid": (enr.get("paper_pmid") or "").strip(),
                "n_rounds_declared": (enr.get("n_rounds_declared") or "").strip(),
                "abstract": (meta.get("study_description") or "").strip(),
            }
        )

    # Carry over non-INSDC deposits from the previous catalog
    if preserve_from:
        rows.extend(_passthrough_non_insdc(preserve_from, exclude_ids=set(studies)))

    # Stable ordering: INSDC studies first (alphabetical), then non-INSDC
    rows.sort(
        key=lambda r: (
            0 if r["source"] == "ena" else 1,
            r["bioproject_id"],
        )
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(PUBLIC_COLS))
        writer.writeheader()
        writer.writerows(rows)

    logger.info("Rebuilt catalog: %d rows → %s", len(rows), out_path)
    return len(rows)
