"""Load the bundled curated SELEX metadata — the **annotated layer**.

Where :mod:`selexprep.catalog.reader` is the *discovery* layer (which deposits
exist), this is the *annotated* layer: each deposit's experimental metadata
(``study_type``, ``target``, ``target_class``, ``chemistry``, ``n_random``,
``n_rounds``, ``selection_format``, ``counter_selection``), curated by **two
independent LLM extractions** (Claude + Codex/GPT) and reconciled.

Every value is source-cited (evidence quote + source + location). Where the two
independent extractions genuinely disagreed, **both are kept** with their
provenance rather than silently picking one.

Two shipped forms:

- ``curated_metadata.json`` — canonical, provenance-rich. Per cell: a
  ``status`` (``concordant`` / ``single_source`` / ``discordant`` /
  ``verified`` [a single-source cell adjudicated 3-way against the benchmark
  ground truth] / ``not_stated``), the value, and ``evidence_quote`` /
  ``source`` / ``location``.
  ``discordant`` cells carry both a ``claude`` and a ``codex`` sub-record.
  Loaded by :func:`load_metadata_records`.
- ``curated_metadata.csv`` — flat view: a ``<field>`` value column and a
  ``<field>_curation`` column per field; discordant shown as ``claude || codex``.
  Loaded by :func:`load_metadata`.

The extraction contract, the two raw arms, and the reconciliation method live in
``benchmarks/dual_extraction/`` in the source repository.
"""

from __future__ import annotations

import json
from importlib.resources import files
from pathlib import Path

import pandas as pd

#: Snapshot identifier for the bundled annotated layer. Bump on every rebuild.
METADATA_VERSION = "v0.2-dual-extraction-2026-07-01"

_JSON_FILENAME = "curated_metadata.json"
_CSV_FILENAME = "curated_metadata.csv"

#: The eight experimental fields curated by the dual-extraction pipeline.
METADATA_FIELDS = (
    "study_type",
    "target",
    "target_class",
    "chemistry",
    "n_random",
    "n_rounds",
    "selection_format",
    "counter_selection",
)


def metadata_version() -> str:
    """Return the snapshot identifier of the bundled annotated layer."""
    return METADATA_VERSION


def metadata_path(fmt: str = "json") -> Path:
    """Filesystem path to the bundled curated-metadata file.

    ``fmt`` is ``"json"`` (canonical, provenance-rich) or ``"csv"`` (flat view).
    Uses ``importlib.resources`` so it resolves in a wheel, sdist, or checkout.
    """
    name = _JSON_FILENAME if fmt == "json" else _CSV_FILENAME
    return Path(str(files("selexprep.catalog.data").joinpath(name)))


def load_metadata() -> pd.DataFrame:
    """Load the flat curated-metadata table as a pandas DataFrame.

    One row per deposit; per field a ``<field>`` value column and a
    ``<field>_curation`` column (``concordant`` / ``single_source:<arm>`` /
    ``discordant`` / ``not_stated``). Discordant values appear as
    ``claude || codex``; full provenance is in :func:`load_metadata_records`.
    """
    return pd.read_csv(metadata_path("csv"), dtype=str, keep_default_na=False)


def load_metadata_records() -> list[dict]:
    """Load the canonical, provenance-rich curated metadata as a list of dicts.

    Each record has ``accession``, ``source``, ``study_title`` and ``fields`` —
    a mapping of each field name to its reconciled cell (``status`` + ``value``
    + ``evidence_quote`` / ``source`` / ``location``; ``discordant`` cells carry
    both a ``claude`` and a ``codex`` sub-record).
    """
    return json.loads(metadata_path("json").read_text(encoding="utf-8"))
