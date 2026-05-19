"""SelexprepManifestV1 - reproducibility anchor for one extraction run.

Locked plan lines 162-175 spec the schema. One ``selexprep_manifest.json``
is emitted per dataset by ``selexprep extract``; downstream tooling (Phase 5
qc, future v0.2 AnnData export) consumes it as the single authoritative
record of what happened during extraction.

**Reproducibility discipline.** Output SHA256s are only guaranteed for
FASTA / TSV / JSON outputs (locked plan line 28). Parquet hashes are
pyarrow-version-dependent across releases; the manifest pins
``pyarrow_version`` instead, and Parquet hashes are intentionally
absent from ``output_sha256``.

Public API:

- :class:`SelexprepManifestV1` - the pydantic model (frozen).
- :func:`compute_sha256s` - hash FASTA/TSV/JSON outputs only.
- :func:`write_manifest_json` / :func:`read_manifest_json` -
  deterministic JSON I/O (matches the ``library/report.py`` discipline).
- :func:`build_manifest_from_extract_result` - convenience builder used
  by ``selexprep extract``.
"""

from __future__ import annotations

import json
import logging
import platform
import subprocess
from collections.abc import Iterable
from hashlib import sha256
from pathlib import Path
from typing import Literal

import dnaio
import pyarrow
from pydantic import BaseModel, ConfigDict

from selexprep import __version__ as _SELEXPREP_VERSION
from selexprep._io import sha256_file
from selexprep.library.report import (
    ExtractionMode,
    LibraryReport,
    ReadSource,
    RequiredAction,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Dependency-version capture
# ---------------------------------------------------------------------------


def _cutadapt_version() -> str:
    """Return the cutadapt version string by calling ``cutadapt --version``.

    Falls back to ``"unknown"`` if cutadapt is not on PATH (the runner
    refuses earlier in that case; this is defense in depth).
    """
    try:
        result = subprocess.run(
            ["cutadapt", "--version"],
            check=True,
            capture_output=True,
            text=True,
            timeout=10,
        )
        return result.stdout.strip() or "unknown"
    except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return "unknown"


# Reproducibility-tracked outputs only — Parquet hashes are advisory
# (locked plan line 28). The runner's outputs use these extensions.
_HASHABLE_SUFFIXES = {".fasta", ".fa", ".tsv", ".json"}


def _is_hashable_output(path: Path) -> bool:
    """True for FASTA/TSV/JSON outputs; False for Parquet/.cutadapt.json/etc."""
    # Strip .gz so foo.fasta.gz matches the same way foo.fasta would.
    name = path.name
    stem = name[:-3] if name.endswith(".gz") else name
    return any(stem.endswith(suf) for suf in _HASHABLE_SUFFIXES)


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------


class SelexprepManifestV1(BaseModel):
    """One reproducibility manifest per extraction run.

    Versioned schema: ``manifest_version`` carries the literal string
    so downstream tooling can detect v0.1 vs future revisions.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    manifest_version: Literal["selexprep_manifest_v1"] = "selexprep_manifest_v1"

    # Dependency versions
    selexprep_version: str
    python_version: str
    cutadapt_version: str
    dnaio_version: str
    pyarrow_version: str

    # Provenance
    accession: str | None
    bioproject_id: str | None
    runs: list[str]

    # Reproducibility hashes
    input_sha256: dict[str, str]
    output_sha256: dict[str, str]

    # Nested LibraryReport - the inference contract
    library_report: LibraryReport

    # Denormalized scan fields (locked plan line 170)
    extraction_mode: ExtractionMode
    read_source: ReadSource
    required_action: RequiredAction
    full_insert_recovered: bool

    # Run-level metadata
    parameters: dict[str, str]
    runtime_seconds_per_stage: dict[str, float]
    flags: list[str]
    sampling_seed: int


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def compute_sha256s(paths: Iterable[Path]) -> dict[str, str]:
    """SHA256 each FASTA/TSV/JSON path; skip Parquet & non-existent paths.

    Returns ``{basename: hex_digest}``. Locked plan line 28: Parquet
    hashes are advisory only; pyarrow_version pins their reproducibility.
    """
    out: dict[str, str] = {}
    for p in paths:
        if not p.exists():
            logger.debug("compute_sha256s: skipping missing path %s", p)
            continue
        if not _is_hashable_output(p):
            logger.debug("compute_sha256s: skipping non-hashable %s", p)
            continue
        out[p.name] = sha256_file(p)
    return out


def _to_deterministic_payload(manifest: SelexprepManifestV1) -> dict[str, object]:
    """Convert a manifest to a dict ready for deterministic JSON dump.

    - Top-level fields follow pydantic declaration order.
    - Nested ``library_report.n_length_distribution`` is sorted numerically
      (otherwise ``sort_keys`` would order "10", "100", "20" lexically).
    - ``input_sha256`` / ``output_sha256`` are sorted alphabetically by name.
    """
    payload: dict[str, object] = manifest.model_dump(mode="json")

    # Sort sha256 dicts by basename.
    for key in ("input_sha256", "output_sha256"):
        raw = payload.get(key)
        if isinstance(raw, dict):
            payload[key] = dict(sorted(raw.items()))

    # Sort the nested LibraryReport's int-keyed n_length_distribution numerically.
    lr_payload = payload.get("library_report")
    if isinstance(lr_payload, dict):
        nld = lr_payload.get("n_length_distribution")
        if isinstance(nld, dict):
            sorted_items = sorted(nld.items(), key=lambda kv: int(kv[0]))
            lr_payload["n_length_distribution"] = dict(sorted_items)
        kah = lr_payload.get("known_adapter_hits")
        if isinstance(kah, dict):
            lr_payload["known_adapter_hits"] = dict(sorted(kah.items()))

    return payload


def write_manifest_json(manifest: SelexprepManifestV1, path: Path) -> str:
    """Write `manifest` to `path` as deterministic JSON; return sha256 hex.

    Output discipline matches ``library.report.write_library_report_json``:

    - UTF-8 encoded, ``ensure_ascii=False``
    - 2-space indent
    - trailing newline
    - top-level fields in pydantic declaration order
    - nested ``n_length_distribution`` in numeric key order
    - ``input_sha256`` / ``output_sha256`` in alphabetical key order
    """
    payload = _to_deterministic_payload(manifest)
    text = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    path.write_text(text, encoding="utf-8")
    return sha256(text.encode("utf-8")).hexdigest()


def read_manifest_json(path: Path) -> SelexprepManifestV1:
    """Load and validate a manifest from disk."""
    return SelexprepManifestV1.model_validate_json(path.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------


def build_manifest_from_extract_result(
    *,
    library_report: LibraryReport,
    input_paths: Iterable[Path],
    output_paths: Iterable[Path],
    accession: str | None,
    bioproject_id: str | None,
    runs: list[str],
    parameters: dict[str, str],
    runtime_seconds_per_stage: dict[str, float] | None = None,
    flags: list[str] | None = None,
) -> SelexprepManifestV1:
    """Convenience builder used by `selexprep extract`.

    Computes input/output SHA256s (FASTA/TSV/JSON only), captures
    dependency versions, denormalizes the LibraryReport's classification
    fields, and assembles the manifest. ``flags`` defaults to ``[]`` —
    Phase 5's QC layer is where they get populated.
    """
    return SelexprepManifestV1(
        selexprep_version=_SELEXPREP_VERSION,
        python_version=platform.python_version(),
        cutadapt_version=_cutadapt_version(),
        dnaio_version=dnaio.__version__,
        pyarrow_version=pyarrow.__version__,
        accession=accession,
        bioproject_id=bioproject_id,
        runs=runs,
        input_sha256=compute_sha256s(input_paths),
        output_sha256=compute_sha256s(output_paths),
        library_report=library_report,
        extraction_mode=library_report.extraction_mode,
        read_source=library_report.read_source,
        required_action=library_report.required_action,
        full_insert_recovered=library_report.full_insert_recovered,
        parameters=parameters,
        runtime_seconds_per_stage=runtime_seconds_per_stage or {},
        flags=flags or [],
        sampling_seed=library_report.sampling_seed,
    )
