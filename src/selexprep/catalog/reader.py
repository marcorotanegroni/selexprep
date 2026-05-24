"""Load the bundled public-SELEX discovery catalog."""

from __future__ import annotations

from importlib.resources import files
from pathlib import Path

import pandas as pd

#: Snapshot identifier for the bundled catalog. Bump this on every refresh.
#: Persona-1 users (ML researchers / atlas builders) rely on this to know
#: how stale the catalog is vs the upstream archives.
#:
#: Phase 6b.5a (2026-05-24): bumped from v0.1.5 to v0.1.6 because the
#: catalog now applies the per-run + per-BioProject library_strategy
#: filter (see ``selexprep.fetch.library_strategy``), which drops
#: 100%-blocklisted studies (RNA-Seq / ChIP-Seq / miRNA-Seq / etc.)
#: that the broad-recall text search would otherwise pick up. Mixed
#: studies are kept; their per-run mix is classified at audit time by
#: Phase 6b.5b.
CATALOG_VERSION = "v0.1.6-snapshot-2026-05-24"

_CATALOG_FILENAME = "bioprojects.csv"


def catalog_version() -> str:
    """Return the snapshot identifier of the bundled catalog."""
    return CATALOG_VERSION


def catalog_path() -> Path:
    """Filesystem path to the bundled catalog CSV.

    Uses ``importlib.resources`` so the path resolves correctly whether the
    package is installed from a wheel, an sdist, or a development checkout.
    """
    return Path(str(files("selexprep.catalog.data").joinpath(_CATALOG_FILENAME)))


def load_catalog() -> pd.DataFrame:
    """Load the bundled discovery catalog as a pandas DataFrame.

    Columns:
        - ``bioproject_id``: INSDC accession or processed-data identifier
          (PRJNA*, PRJDB*, PRJEB*, SRP*, ERP*, DRP*, zenodo:*, figshare:*)
        - ``source``: discovery adapter that surfaced this entry
        - ``study_title``, ``protein_target``, ``target_organism``
        - ``paper_doi``, ``paper_pmid``, ``n_rounds_declared``
        - ``abstract``: study abstract (may be long)
    """
    df = pd.read_csv(catalog_path(), dtype=str, keep_default_na=False)
    return df
