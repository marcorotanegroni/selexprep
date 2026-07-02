"""Discovery catalog of public SELEX bioprojects.

The catalog ships as package data so a fresh `pip install selexprep`
lets users immediately browse what is publicly available without
running live API queries against ENA / SRA / GEO / Zenodo / Figshare.

Two layers ship bundled:

- **Discovery layer** (:mod:`selexprep.catalog.reader`) — bioproject metadata
  only (title, target, organism, paper DOI, declared round count).
- **Curated metadata layer** (:mod:`selexprep.catalog.metadata`) — each
  deposit's experimental fields (study_type, target, target_class, chemistry,
  n_random, n_rounds, selection_format, counter_selection), curated by **dual
  independent LLM extraction** (Claude + Codex/GPT) and reconciled, with every
  value source-cited and genuine disagreements kept as both arms.

The LibraryReport enrichment (primer pair / extraction_mode *inferred* from the
reads by the preprocessing pipeline) remains future work.

Public API:

- :func:`load_catalog` — pandas DataFrame of every discovery-catalog row
- :func:`filter_catalog` — common filters (target, organism, INSDC-only, min rounds)
- :func:`catalog_version` — discovery snapshot identifier
- :func:`load_metadata` / :func:`load_metadata_records` — the curated metadata
  layer, flat (DataFrame) or provenance-rich (list of dicts)
- :func:`metadata_version` — curated-layer snapshot identifier
"""

from selexprep.catalog.filter import filter_catalog
from selexprep.catalog.metadata import (
    load_metadata,
    load_metadata_records,
    metadata_version,
)
from selexprep.catalog.reader import catalog_version, load_catalog

__all__ = [
    "catalog_version",
    "filter_catalog",
    "load_catalog",
    "load_metadata",
    "load_metadata_records",
    "metadata_version",
]
