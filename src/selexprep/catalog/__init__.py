"""Discovery catalog of public SELEX bioprojects.

The catalog ships as package data so a fresh `pip install selexprep`
lets users immediately browse what is publicly available without
running live API queries against ENA / SRA / GEO / Zenodo / Figshare.

This is the **discovery layer** of the catalog — bioproject metadata only
(title, target, organism, paper DOI, declared round count). The v0.2
**annotated layer** will enrich each row with the inferred LibraryReport
(primer pair, N-region length, extraction_mode, confidence) once the
full pipeline has been run over every entry.

Public API:

- :func:`load_catalog` — pandas DataFrame of every catalog row
- :func:`filter_catalog` — common filters (target, organism, INSDC-only,
  min rounds)
- :func:`catalog_version` — snapshot identifier
"""

from selexprep.catalog.filter import filter_catalog
from selexprep.catalog.reader import catalog_version, load_catalog

__all__ = ["catalog_version", "filter_catalog", "load_catalog"]
