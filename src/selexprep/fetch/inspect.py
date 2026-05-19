"""Accession preview - ENA filereport REST without download.

Locked plan line 332: ``selexprep inspect <accession>`` fetches ENA/SRA
metadata only. Prints round count, ``library_strategy``,
``library_source``, per-run file sizes + MD5s. **The library_strategy
field is reported verbatim from SRA - it is NOT interpreted as a
DNA/RNA classification** (that's v0.2's `library-type-classifier`).

The ENA Portal filereport endpoint accepts study- (``PRJ*``/``SRP*``/
``ERP*``/``DRP*``) and run-level (``SRR*``/``ERR*``/``DRR*``) accessions
and returns one row per run.

Public API:

- :class:`InspectReport`, :class:`RunFileInfo` - typed result dataclasses.
- :func:`inspect_accession` - hit ENA, parse JSON, return ``InspectReport``.
- :func:`write_inspect_json` - emit a deterministic JSON sidecar when
  ``selexprep inspect --outdir`` is given.
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass
from pathlib import Path

import requests

logger = logging.getLogger(__name__)


ENA_FILEREPORT_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

# Fields requested in the ENA filereport response. Each is a column in the
# returned JSON; per-row, columns with multiple values (``fastq_md5``,
# ``fastq_bytes``, ``fastq_ftp``) are semicolon-delimited.
_ENA_FIELDS = (
    "run_accession,study_accession,study_title,library_strategy,"
    "library_source,read_count,base_count,fastq_md5,fastq_bytes,fastq_ftp"
)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RunFileInfo:
    """One sequencing run inside a study (one ENA filereport row)."""

    run_accession: str
    read_count: int
    base_count: int
    fastq_size_bytes: list[int]
    fastq_md5: list[str]


@dataclass(frozen=True)
class InspectReport:
    """Result of :func:`inspect_accession`.

    Study-level metadata + per-run file info. ``library_strategy`` /
    ``library_source`` are verbatim SRA strings; do NOT treat as a
    DNA/RNA hint (locked plan line 332).
    """

    accession: str
    bioproject_id: str | None
    study_title: str
    library_strategy: str
    library_source: str
    runs: list[RunFileInfo]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def _split_semicolon_list(value: object) -> list[str]:
    """Parse a semicolon-delimited string from ENA into a clean list.

    ENA returns columns like ``fastq_md5`` as ``"abc;def"`` for paired-end
    runs and ``""`` when absent. Some rows return numeric types; coerce
    defensively.
    """
    if not value:
        return []
    s = str(value)
    return [p for p in (chunk.strip() for chunk in s.split(";")) if p]


def inspect_accession(accession: str, *, timeout_s: int = 30) -> InspectReport:
    """Hit ENA filereport REST for `accession`; return parsed `InspectReport`.

    Accepts study- or run-level accessions. Raises:

    - ``requests.HTTPError`` on non-2xx (network / 404 unknown accession)
    - ``ValueError`` if ENA returns an empty record set
    """
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": _ENA_FIELDS,
        "format": "json",
    }
    logger.info("inspect_accession: GET %s params=%s", ENA_FILEREPORT_URL, params)
    response = requests.get(ENA_FILEREPORT_URL, params=params, timeout=timeout_s)
    response.raise_for_status()
    rows = response.json() or []
    if not rows:
        raise ValueError(f"ENA returned no records for accession {accession!r}")

    runs: list[RunFileInfo] = []
    for row in rows:
        size_strs = _split_semicolon_list(row.get("fastq_bytes"))
        sizes: list[int] = []
        for s in size_strs:
            try:
                sizes.append(int(s))
            except ValueError:
                logger.warning("inspect_accession: bad fastq_bytes entry %r in %s", s, accession)
        runs.append(
            RunFileInfo(
                run_accession=str(row.get("run_accession", "")),
                read_count=int(row.get("read_count") or 0),
                base_count=int(row.get("base_count") or 0),
                fastq_size_bytes=sizes,
                fastq_md5=_split_semicolon_list(row.get("fastq_md5")),
            )
        )

    first = rows[0]
    return InspectReport(
        accession=accession,
        bioproject_id=(str(first.get("study_accession")) if first.get("study_accession") else None),
        study_title=str(first.get("study_title") or ""),
        library_strategy=str(first.get("library_strategy") or ""),
        library_source=str(first.get("library_source") or ""),
        runs=runs,
    )


def write_inspect_json(report: InspectReport, path: Path) -> None:
    """Write an InspectReport as deterministic JSON (sort_keys + trailing newline).

    Useful for batch / scripted pipelines that consume the inspect
    output downstream of an interactive ``selexprep inspect`` call.
    """
    payload = {
        "accession": report.accession,
        "bioproject_id": report.bioproject_id,
        "study_title": report.study_title,
        "library_strategy": report.library_strategy,
        "library_source": report.library_source,
        "runs": [asdict(r) for r in report.runs],
    }
    text = json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    path.write_text(text, encoding="utf-8")
