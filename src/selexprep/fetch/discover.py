"""Multi-source discovery of public SELEX datasets.

Populates a three-table relational metadata schema (bioprojects, samples,
rounds) from up to nine sources, in priority order:

1. **Seed file** (caller-supplied YAML) — ground truth, always runs
2. **ENA REST API** — primary, full-text search
3. **NCBI SRA** (via ``pysradb``, optional dep) — complementary to ENA
4. **NCBI GEO** (via ``Bio.Entrez``, optional dep) — some SELEX studies here
5. **UTexas Aptamer Database** (Zenodo record 8387047) — meta-curated
6. **Zenodo API** — processed datasets
7. **Figshare API** — processed datasets
8. **Crossref** — DOI → paper metadata (enrichment)
9. **OpenAlex** — DOI → open-access URL

Each adapter inherits from ``SourceAdapter`` and returns
``(bioproject_rows, sample_rows)``. Optional dependencies (``pysradb``,
``Bio.Entrez``) soft-fail with a log warning when missing — discovery
still runs against the remaining sources.

**v0.1 packaging note:** the nine adapter classes live in this single
module. v0.2 will split them into ``selexprep.fetch.sources.*`` submodules
for testability (each adapter currently makes live HTTP calls in its
``search()`` method, which complicates unit testing).
"""

from __future__ import annotations

import contextlib
import csv
import io
import json
import logging
import re
import shutil
import time
from abc import ABC, abstractmethod
from collections import defaultdict
from pathlib import Path
from typing import Any, ClassVar

import requests
import yaml

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

BIOPROJECT_COLS = [
    "bioproject_id",
    "source",
    "study_title",
    "protein_target",
    "target_organism",
    "paper_doi",
    "paper_pmid",
    "n_rounds_declared",
    "library_type_verification",
    "library_type_evidence",
    "has_processed_counts",
    "abstract",
    "include",
    "manual_curation_notes",
]

SAMPLE_COLS = [
    "srr",
    "bioproject_id",
    "sample_title",
    "library_name",
    "experiment_title",
    "design_description",
    "sample_attributes",
    "total_bases",
    "n_reads",
    "target_hint",
    "raw_metadata",
]

ROUND_COLS = [
    "srr",
    "round_number",
    "confidence",
    "source_field",
    "matched_pattern",
    "round_candidates",
    "needs_manual_review",
    "parser_notes",
]


def empty_bioproject() -> dict:
    return {c: "" for c in BIOPROJECT_COLS}


def empty_sample() -> dict:
    return {c: "" for c in SAMPLE_COLS}


# ---------------------------------------------------------------------------
# HTTP helper
# ---------------------------------------------------------------------------


def _get(
    url: str,
    params: dict | None = None,
    retries: int = 3,
    timeout: int = 30,
) -> dict | list | None:
    """GET with exponential backoff. Returns parsed JSON or None on failure."""
    for attempt in range(retries):
        try:
            resp = requests.get(
                url,
                params=params,
                timeout=timeout,
                headers={"Accept": "application/json"},
            )
            if resp.status_code == 200:
                return resp.json()
            if resp.status_code == 429:
                wait = 2 ** (attempt + 2)
                logger.warning("Rate limited by %s, waiting %ds", url, wait)
                time.sleep(wait)
            else:
                logger.warning("HTTP %d from %s", resp.status_code, url)
                return None
        except requests.RequestException as e:
            logger.warning("Request error (attempt %d/%d): %s", attempt + 1, retries, e)
            time.sleep(2**attempt)
    return None


# ---------------------------------------------------------------------------
# Adapter base
# ---------------------------------------------------------------------------


class SourceAdapter(ABC):
    name: str = "base"

    @abstractmethod
    def search(self, query: str | None = None) -> tuple[list[dict], list[dict]]:
        """Return (bioproject_rows, sample_rows)."""
        ...

    def fetch_bioproject(self, bioproject_id: str) -> tuple[dict | None, list[dict]]:
        rows_bp, rows_s = self.search(query=bioproject_id)
        bp = next((r for r in rows_bp if r.get("bioproject_id") == bioproject_id), None)
        samples = [s for s in rows_s if s.get("bioproject_id") == bioproject_id]
        return bp, samples


# ---------------------------------------------------------------------------
# 1. Seed adapter (caller-supplied YAML)
# ---------------------------------------------------------------------------


class SeedAdapter(SourceAdapter):
    """Load a curated seed list from a user-supplied YAML file.

    Expected schema::

        entries:
          - bioproject_id: PRJNAxxxxxx
            protein_target: VEGF
            paper_doi: 10.xxxx/xxxxx
            paper_pmid: 12345678
            n_rounds_expected: 10
            library_type: RNA_confirmed   # or DNA_confirmed / ambiguous
            notes: |
              Optional curation notes.
            manual_round_mapping:         # optional
              SRRxxxxxxx: 0
              SRRxxxxxxy: 1
        blacklist:
          - bioproject_id: PRJxxxxxxx
        small_molecule_targets:
          - GTP
          - ATP
    """

    name = "seed"

    def __init__(self, seed_file: Path):
        self.seed_file = seed_file
        with open(seed_file, encoding="utf-8") as f:
            self.data = yaml.safe_load(f) or {}

    @property
    def blacklist(self) -> set[str]:
        return {e["bioproject_id"] for e in self.data.get("blacklist", [])}

    def search(self, query: str | None = None) -> tuple[list[dict], list[dict]]:
        bp_rows = []
        for entry in self.data.get("entries", []):
            if query and query not in (
                entry.get("bioproject_id", ""),
                entry.get("protein_target", ""),
            ):
                continue
            bp = empty_bioproject()
            bp["bioproject_id"] = entry["bioproject_id"]
            bp["source"] = "seed"
            bp["protein_target"] = entry.get("protein_target", "")
            bp["paper_doi"] = entry.get("paper_doi", "")
            bp["paper_pmid"] = str(entry.get("paper_pmid", ""))
            bp["n_rounds_declared"] = str(entry.get("n_rounds_expected", ""))
            bp["library_type_verification"] = entry.get("library_type", "ambiguous")
            bp["library_type_evidence"] = json.dumps({"source": "seed_file"})
            bp["include"] = "y" if entry.get("library_type") == "RNA_confirmed" else "n"
            bp["manual_curation_notes"] = entry.get("notes", "").strip()
            bp_rows.append(bp)
        return bp_rows, []


# ---------------------------------------------------------------------------
# 2. ENA adapter
# ---------------------------------------------------------------------------

ENA_PORTAL_URL = "https://www.ebi.ac.uk/ena/portal/api/search"

ENA_QUERIES = [
    'study_title="HT-SELEX"',
    'study_title="SELEX-seq"',
    'study_title="high-throughput SELEX"',
    'study_title="RNA aptamer"',
    'study_title="aptamer selection"',
    'study_title="systematic evolution of ligands"',
    'description="HT-SELEX"',
    'description="SELEX-seq"',
    'description="aptamer selection"',
    'description="RNA aptamer" AND description="round"',
]

ENA_RUN_FIELDS = (
    "run_accession,study_accession,study_title,sample_title,"
    "library_name,experiment_title,library_strategy,library_source,"
    "base_count,read_count,sample_alias,scientific_name"
)


class ENAAdapter(SourceAdapter):
    name = "ena"

    def search(self, query: str | None = None) -> tuple[list[dict], list[dict]]:
        from selexprep.fetch.library_strategy import (
            classify_study_by_library_strategies,
        )

        queries = [query] if query else ENA_QUERIES
        bp_rows: dict[str, dict] = {}
        sample_rows: list[dict] = []
        # Phase 6b.5a: collect per-BioProject library_strategy values so
        # we can apply the same per-run + per-study classifier as the
        # ``selexprep catalog refresh`` path (see
        # ``selexprep.catalog.rebuild``). Studies whose runs are 100%
        # blocklisted (RNA-Seq / ChIP-Seq / etc.) are filtered out here
        # too, keeping the discovery pipeline consistent with refresh.
        strategies_by_bp: dict[str, list[str]] = {}

        for q in queries:
            logger.info("[ENA] searching: %s", q)
            runs = _get(
                ENA_PORTAL_URL,
                params={
                    "result": "read_run",
                    "query": q,
                    "fields": ENA_RUN_FIELDS,
                    "limit": 5000,
                    "format": "json",
                },
            )
            if not runs:
                continue
            for run in runs:
                bp_id = run.get("study_accession", "")
                if not bp_id:
                    continue
                if bp_id not in bp_rows:
                    bp = empty_bioproject()
                    bp["bioproject_id"] = bp_id
                    bp["source"] = "ena"
                    bp["study_title"] = run.get("study_title", "")
                    bp_rows[bp_id] = bp
                    strategies_by_bp[bp_id] = []

                strategies_by_bp[bp_id].append((run.get("library_strategy") or "").strip())

                s = empty_sample()
                s["srr"] = run.get("run_accession", "")
                s["bioproject_id"] = bp_id
                s["sample_title"] = run.get("sample_title", "") or run.get("sample_alias", "")
                s["library_name"] = run.get("library_name", "")
                s["experiment_title"] = run.get("experiment_title", "")
                s["total_bases"] = str(run.get("base_count", ""))
                s["n_reads"] = str(run.get("read_count", ""))
                s["raw_metadata"] = json.dumps(run)
                sample_rows.append(s)
            time.sleep(0.4)

        # Apply per-study filter: drop BioProjects whose runs are 100%
        # blocklisted. Mixed studies (some compatible + some blocklisted)
        # are KEPT; the audit-eligibility layer (Phase 6b.5b) classifies
        # them downstream. See ``selexprep.fetch.library_strategy`` for
        # the decision rule.
        excluded_ids: set[str] = set()
        for bp_id, strategies in strategies_by_bp.items():
            classification = classify_study_by_library_strategies(bp_id, strategies)
            if classification.should_exclude:
                excluded_ids.add(bp_id)
                logger.info(
                    "[ENA] excluding %s: %s",
                    bp_id,
                    classification.exclusion_reason,
                )

        kept_bp_rows = [bp for bp_id, bp in bp_rows.items() if bp_id not in excluded_ids]
        kept_sample_rows = [s for s in sample_rows if s["bioproject_id"] not in excluded_ids]
        return kept_bp_rows, kept_sample_rows


# ---------------------------------------------------------------------------
# 3. SRA adapter (optional: pysradb)
# ---------------------------------------------------------------------------

SRA_QUERIES = [
    '"HT-SELEX"',
    '"high-throughput SELEX"',
    '"SELEX-seq"',
    '"aptamer" AND "selection cycles" AND "sequencing"',
    '"aptamer pool" AND "sequencing" AND "round"',
    '"aptamer" AND "NGS" AND "selection rounds"',
    '"RNA aptamer" AND "high-throughput sequencing"',
    '"in vitro selection" AND "aptamer" AND "RNA"',
    '"systematic evolution of ligands" AND "aptamer"',
]


class SRAAdapter(SourceAdapter):
    name = "sra"

    def __init__(self) -> None:
        try:
            from pysradb import SRAweb

            self._db = SRAweb()
        except ImportError:
            logger.warning(
                "pysradb not installed — SRA adapter disabled (install with `pip install selexprep[ncbi]`)"
            )
            self._db = None

    def search(self, query: str | None = None) -> tuple[list[dict], list[dict]]:
        if self._db is None:
            return [], []

        queries = [query] if query else SRA_QUERIES
        bp_rows: dict[str, dict] = {}
        sample_rows: list[dict] = []

        for q in queries:
            logger.info("[SRA] searching: %s", q)
            try:
                df = self._db.search_sra(q)
            except TypeError:
                try:
                    df = self._db.search_sra(q, return_max=500)
                except Exception as e:
                    logger.warning("[SRA] search failed: %s", e)
                    continue
            except Exception as e:
                logger.warning("[SRA] search failed: %s", e)
                continue

            if df is None or df.empty:
                continue

            cols = set(df.columns)
            logger.info("[SRA] %d rows", len(df))

            def _pick(row, *names, _cols=cols):
                for n in names:
                    if n in _cols:
                        v = row.get(n)
                        if v is not None and str(v).strip() not in ("", "nan"):
                            return str(v).strip()
                return ""

            for _, row in df.iterrows():
                bp_id = _pick(
                    row, "bioproject_accession", "study_accession", "bioproject", "study_alias"
                )
                if not bp_id:
                    continue
                if bp_id not in bp_rows:
                    bp = empty_bioproject()
                    bp["bioproject_id"] = bp_id
                    bp["source"] = "sra"
                    bp["study_title"] = _pick(row, "study_title")
                    bp["abstract"] = _pick(row, "study_abstract", "study_description")
                    bp_rows[bp_id] = bp

                srr = _pick(row, "run_accession", "run_1_accession")
                if not srr:
                    continue
                s = empty_sample()
                s["srr"] = srr
                s["bioproject_id"] = bp_id
                s["sample_title"] = _pick(row, "sample_title", "experiment_title")
                s["library_name"] = _pick(row, "library_name")
                s["experiment_title"] = _pick(row, "experiment_title")
                s["total_bases"] = _pick(row, "total_bases", "run_total_bases")
                s["n_reads"] = _pick(row, "total_spots", "run_total_spots")
                s["raw_metadata"] = json.dumps({k: str(v) for k, v in row.to_dict().items()})
                sample_rows.append(s)

        return list(bp_rows.values()), sample_rows


# ---------------------------------------------------------------------------
# 4. GEO adapter (optional: Bio.Entrez)
# ---------------------------------------------------------------------------


class GEOAdapter(SourceAdapter):
    name = "geo"

    GEO_QUERIES: ClassVar[list[str]] = [
        "HT-SELEX RNA aptamer",
        "SELEX-seq RNA aptamer protein",
        "aptamer selection RNA pool high-throughput sequencing",
    ]

    def __init__(self, email: str = "selexprep@local") -> None:
        try:
            from Bio import Entrez

            self.Entrez = Entrez
            Entrez.email = email
        except ImportError:
            logger.warning(
                "Biopython not installed — GEO adapter disabled (install with `pip install selexprep[ncbi]`)"
            )
            self.Entrez = None

    def search(self, query: str | None = None) -> tuple[list[dict], list[dict]]:
        if self.Entrez is None:
            return [], []

        queries = [query] if query else self.GEO_QUERIES
        bp_rows: dict[str, dict] = {}

        for q in queries:
            logger.info("[GEO] searching: %s", q)
            try:
                handle = self.Entrez.esearch(db="gds", term=q, retmax=200)
                record = self.Entrez.read(handle)
                handle.close()
                ids = record.get("IdList", [])
            except Exception as e:
                logger.warning("[GEO] search error: %s", e)
                continue

            for geo_id in ids:
                try:
                    handle = self.Entrez.esummary(db="gds", id=geo_id)
                    summary = self.Entrez.read(handle)
                    handle.close()
                except Exception:
                    continue

                for item in summary:
                    srp = item.get("Relations", {})
                    sra_accession = ""
                    for rel in srp if isinstance(srp, list) else []:
                        if rel.get("RelationshipType") == "SRA":
                            sra_accession = rel.get("TargetObject", "")
                            break

                    geo_acc = item.get("Accession", "")
                    title = item.get("title", "")
                    summary_text = item.get("summary", "")

                    key = sra_accession or geo_acc
                    if not key or key in bp_rows:
                        continue

                    bp = empty_bioproject()
                    bp["bioproject_id"] = key
                    bp["source"] = f"geo:{geo_acc}"
                    bp["study_title"] = title
                    bp["abstract"] = summary_text
                    bp_rows[key] = bp
                time.sleep(0.5)

        return list(bp_rows.values()), []


# ---------------------------------------------------------------------------
# 5. UTexas Aptamer Database (Zenodo record)
# ---------------------------------------------------------------------------

UTEXAS_ZENODO_RECORD = "8387047"
UTEXAS_ZENODO_URL = f"https://zenodo.org/api/records/{UTEXAS_ZENODO_RECORD}"


class UTexasDBAdapter(SourceAdapter):
    name = "utexas_db"

    def search(self, query: str | None = None) -> tuple[list[dict], list[dict]]:
        logger.info("[UTexasDB] fetching Zenodo record %s", UTEXAS_ZENODO_RECORD)
        meta = _get(UTEXAS_ZENODO_URL)
        if not meta:
            return [], []

        files = meta.get("files", [])
        data_url = next((f["links"]["self"] for f in files if f["key"].endswith(".csv")), None)
        is_xlsx = False
        if not data_url:
            data_url = next(
                (
                    f["links"]["self"]
                    for f in files
                    if f["key"].endswith(".xlsx") and "dataset" in f["key"].lower()
                ),
                None,
            )
            is_xlsx = True
        if not data_url:
            logger.warning("[UTexasDB] no CSV/XLSX file found in Zenodo record")
            return [], []

        logger.info("[UTexasDB] downloading %s from %s", "XLSX" if is_xlsx else "CSV", data_url)
        resp = requests.get(data_url, timeout=120)
        if resp.status_code != 200:
            logger.warning("[UTexasDB] failed to download: HTTP %d", resp.status_code)
            return [], []

        if is_xlsx:
            try:
                import pandas as pd

                df = pd.read_excel(io.BytesIO(resp.content), engine="openpyxl")
                rows_iter = df.to_dict(orient="records")
            except Exception as e:
                logger.warning("[UTexasDB] failed to parse XLSX: %s", e)
                return [], []
        else:
            rows_iter = list(csv.DictReader(io.StringIO(resp.text)))
        bp_rows: dict[str, dict] = {}

        for row in rows_iter:
            composition = str(row.get("Nucleic Acid Composition", "")).lower()
            target = str(row.get("Target Name", "")).strip()
            target_type = str(row.get("Target Type", "")).lower()

            if "rna" not in composition:
                continue
            if not target:
                continue
            if any(t in target_type for t in ("small molecule", "nucleic acid", "ion", "other")):
                continue

            doi = str(row.get("DOI", "")).strip()
            pmid = str(row.get("PMID", "")).strip()
            key = f"utexas:{doi or pmid or target}"

            if key not in bp_rows:
                bp = empty_bioproject()
                bp["bioproject_id"] = key
                bp["source"] = "utexas_db"
                bp["protein_target"] = target
                bp["paper_doi"] = doi
                bp["paper_pmid"] = pmid
                bp["library_type_verification"] = "RNA_confirmed"
                bp["library_type_evidence"] = json.dumps(
                    {"source": "UTexas Aptamer DB, RNA composition field"}
                )
                bp["include"] = "maybe"
                bp["manual_curation_notes"] = (
                    "From UTexas Aptamer Database — check if SRA raw SELEX data exists for this paper"
                )
                bp_rows[key] = bp

        logger.info("[UTexasDB] found %d RNA-vs-protein entries", len(bp_rows))
        return list(bp_rows.values()), []


# ---------------------------------------------------------------------------
# 6. Zenodo adapter
# ---------------------------------------------------------------------------

ZENODO_API = "https://zenodo.org/api/records"


class ZenodoAdapter(SourceAdapter):
    name = "zenodo"

    QUERIES: ClassVar[list[str]] = [
        "HT-SELEX aptamer",
        "SELEX RNA aptamer round",
        "aptamer selection sequencing",
    ]

    def search(self, query: str | None = None) -> tuple[list[dict], list[dict]]:
        queries = [query] if query else self.QUERIES
        bp_rows: dict[str, dict] = {}

        for q in queries:
            logger.info("[Zenodo] searching: %s", q)
            data = _get(ZENODO_API, params={"q": q, "size": 25})
            if not data:
                continue
            for hit in data.get("hits", {}).get("hits", []):
                zenodo_id = str(hit.get("id", ""))
                key = f"zenodo:{zenodo_id}"
                if key in bp_rows:
                    continue
                meta = hit.get("metadata", {})
                bp = empty_bioproject()
                bp["bioproject_id"] = key
                bp["source"] = f"zenodo:{zenodo_id}"
                bp["study_title"] = meta.get("title", "")
                bp["abstract"] = meta.get("description", "")
                bp["paper_doi"] = meta.get("doi", "")
                bp["has_processed_counts"] = "y"
                bp["include"] = "maybe"
                bp["manual_curation_notes"] = (
                    f"Zenodo record {zenodo_id} — verify if contains per-round SELEX counts"
                )
                bp_rows[key] = bp
            time.sleep(0.5)

        return list(bp_rows.values()), []


# ---------------------------------------------------------------------------
# 7. Figshare adapter
# ---------------------------------------------------------------------------

FIGSHARE_API = "https://api.figshare.com/v2/articles/search"


class FigshareAdapter(SourceAdapter):
    name = "figshare"

    QUERIES: ClassVar[list[str]] = ["HT-SELEX aptamer RNA", "aptamer SELEX rounds sequencing"]

    def search(self, query: str | None = None) -> tuple[list[dict], list[dict]]:
        queries = [query] if query else self.QUERIES
        bp_rows: dict[str, dict] = {}

        for q in queries:
            logger.info("[Figshare] searching: %s", q)
            payload = {"search_for": q, "page_size": 100}
            try:
                resp = requests.post(
                    FIGSHARE_API,
                    json=payload,
                    timeout=30,
                    headers={"Content-Type": "application/json"},
                )
                if resp.status_code != 200:
                    continue
                articles = resp.json()
            except Exception as e:
                logger.warning("[Figshare] error: %s", e)
                continue

            for art in articles:
                fig_id = str(art.get("id", ""))
                key = f"figshare:{fig_id}"
                if key in bp_rows:
                    continue
                bp = empty_bioproject()
                bp["bioproject_id"] = key
                bp["source"] = f"figshare:{fig_id}"
                bp["study_title"] = art.get("title", "")
                bp["paper_doi"] = art.get("doi", "")
                bp["has_processed_counts"] = "y"
                bp["include"] = "maybe"
                bp["manual_curation_notes"] = (
                    f"Figshare article {fig_id} — verify if contains per-round SELEX data"
                )
                bp_rows[key] = bp
            time.sleep(0.5)

        return list(bp_rows.values()), []


# ---------------------------------------------------------------------------
# 8. Crossref (DOI → paper metadata)
# ---------------------------------------------------------------------------

CROSSREF_API = "https://api.crossref.org/works"


class CrossrefAdapter(SourceAdapter):
    name = "crossref"

    def search(self, query: str | None = None) -> tuple[list[dict], list[dict]]:
        return [], []

    def enrich_doi(self, doi: str) -> dict:
        """Fetch structured metadata for a DOI."""
        url = f"{CROSSREF_API}/{doi}"
        data = _get(url)
        if not data:
            return {}
        msg = data.get("message", {})
        return {
            "title": " ".join(msg.get("title", [])),
            "authors": [
                f"{a.get('family', '')} {a.get('given', '')}".strip() for a in msg.get("author", [])
            ],
            "journal": msg.get("container-title", [""])[0] if msg.get("container-title") else "",
            "year": (msg.get("published-print") or msg.get("published-online") or {}).get(
                "date-parts", [[None]]
            )[0][0],
            "abstract": msg.get("abstract", ""),
        }


# ---------------------------------------------------------------------------
# 9. OpenAlex (DOI → OA URL)
# ---------------------------------------------------------------------------

OPENALEX_API = "https://api.openalex.org/works"


class OpenAlexAdapter(SourceAdapter):
    name = "openalex"

    def search(self, query: str | None = None) -> tuple[list[dict], list[dict]]:
        return [], []

    def get_oa_url(self, doi: str) -> str | None:
        """Return open-access PDF/HTML URL for a DOI, or None."""
        data = _get(f"{OPENALEX_API}/https://doi.org/{doi}")
        if not data:
            return None
        oa = data.get("open_access", {})
        return oa.get("oa_url")


# ---------------------------------------------------------------------------
# Blacklist filter
# ---------------------------------------------------------------------------


def load_blacklist(seed_file: Path) -> set[str]:
    """Load the `blacklist` block from a seed YAML, returning bioproject IDs."""
    with open(seed_file, encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    return {e["bioproject_id"] for e in data.get("blacklist", [])}


def load_small_molecule_targets(seed_file: Path) -> list[str]:
    """Load the `small_molecule_targets` list from a seed YAML (lowercased)."""
    with open(seed_file, encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    return [t.lower() for t in data.get("small_molecule_targets", [])]


_SELEX_CONTEXT_TOKENS = ("selex", "aptamer", "riboswitch")


def sm_mentioned_in_selex_context(abstract: str, sm_terms: list[str]) -> bool:
    """True iff a small-molecule term appears within ~200 chars of a SELEX keyword.

    Used to avoid false positives like "ATP-binding cassette" in unrelated
    protein-SELEX papers.
    """
    if not abstract:
        return False
    text = abstract.lower()
    ctx = "|".join(_SELEX_CONTEXT_TOKENS)
    for sm in sm_terms:
        sm_re = re.escape(sm.lower())
        pattern = re.compile(
            rf"\b{sm_re}\b.{{0,200}}(?:{ctx})|(?:{ctx}).{{0,200}}\b{sm_re}\b",
            re.IGNORECASE | re.DOTALL,
        )
        if pattern.search(text):
            return True
    return False


def is_blacklisted(bp: dict, blacklist: set[str], sm_targets: list[str]) -> bool:
    """Decide whether a discovered BioProject should be filtered out.

    Rejects:
    - explicit blacklist hits
    - small-molecule-target keywords in the protein_target field
    - small-molecule term co-occurring with a SELEX keyword in the abstract
    - whole-cell SELEX without a named protein target
    """
    if bp.get("bioproject_id") in blacklist:
        return True
    target = bp.get("protein_target", "").lower()
    abstract = bp.get("abstract", "") or ""
    if any(sm in target for sm in sm_targets):
        return True
    if sm_mentioned_in_selex_context(abstract, sm_targets[:5]):
        return True
    abstract_lc = abstract.lower()
    return ("whole cell" in abstract_lc or "differential cell selex" in abstract_lc) and not bp.get(
        "protein_target"
    )


# ---------------------------------------------------------------------------
# Library-type classification (optional — v0.2)
# ---------------------------------------------------------------------------


#: Sentinel written to `library_type_verification` when the v0.2 classifier
#: is not installed. Distinguishes "deliberately not assessed in v0.1" from a
#: bona fide empty string that would let downstream callers silently skip the
#: filter. Phase 2 will replace it with a real verdict when the classifier
#: lands.
NOT_ASSESSED_V0_1 = "not_assessed_v0.1"


def _classify_all(bp_rows: list[dict]) -> list[dict]:
    """Apply DNA/RNA library-type classification if the v0.2 classifier is
    available. If not, write the `not_assessed_v0.1` sentinel into every
    unclassified row's ``library_type_verification`` (plus an evidence JSON)
    so callers can detect the deferred-classifier path explicitly rather
    than treating empty cells as a successful "unclassified" outcome.
    """
    try:
        from selexprep.library.type_classifier import classify_bioproject  # type: ignore
    except ImportError:
        logger.info(
            "library_type_classifier not available (v0.2 feature); marking "
            "unclassified BioProjects with sentinel '%s'",
            NOT_ASSESSED_V0_1,
        )
        sentinel_evidence = json.dumps(
            {
                "source": "no_classifier",
                "reason": "library.type_classifier deferred to v0.2",
            }
        )
        for bp in bp_rows:
            if not bp.get("library_type_verification"):
                bp["library_type_verification"] = NOT_ASSESSED_V0_1
                bp["library_type_evidence"] = sentinel_evidence
        return bp_rows
    for bp in bp_rows:
        if bp.get("library_type_verification"):
            continue
        result = classify_bioproject(bp)
        bp["library_type_verification"] = result.verdict
        bp["library_type_evidence"] = result.to_json()
        if result.verdict == "DNA_confirmed":
            bp["include"] = "n"
        elif result.verdict == "ambiguous" and not bp.get("include"):
            bp["include"] = "maybe"
    return bp_rows


# ---------------------------------------------------------------------------
# Round parsing (via selexprep.fetch.metadata)
# ---------------------------------------------------------------------------


def load_seed_overrides(seed_file: Path) -> dict[str, dict[str, int]]:
    """Extract `manual_round_mapping` blocks from a seed YAML, keyed by BP id."""
    if not seed_file.exists():
        return {}
    with open(seed_file, encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    out: dict[str, dict[str, int]] = {}
    for entry in data.get("entries", []):
        mapping = entry.get("manual_round_mapping")
        if mapping:
            out[entry["bioproject_id"]] = {str(k): int(v) for k, v in mapping.items()}
    return out


def _parse_rounds(
    sample_rows: list[dict],
    manual_review_dir: Path,
    bp_abstracts: dict[str, str] | None = None,
    seed_overrides: dict[str, dict[str, int]] | None = None,
) -> list[dict]:
    from selexprep.fetch.metadata import apply_seed_overrides, parse_round

    bp_abstracts = bp_abstracts or {}
    seed_overrides = seed_overrides or {}

    records_by_bp: dict[str, list[Any]] = defaultdict(list)

    for s in sample_rows:
        attrs: dict[str, str] = {}
        raw_attrs = s.get("sample_attributes", "")
        if raw_attrs:
            with contextlib.suppress(json.JSONDecodeError):
                attrs = json.loads(raw_attrs)

        bp_id = s.get("bioproject_id", "")
        record = parse_round(
            srr=s["srr"],
            sample_title=s.get("sample_title", ""),
            library_name=s.get("library_name", ""),
            experiment_title=s.get("experiment_title", ""),
            design_description=s.get("design_description", ""),
            sample_attributes=attrs,
            abstract=bp_abstracts.get(bp_id, ""),
            manual_review_dir=manual_review_dir,
            bioproject_id=bp_id,
        )
        records_by_bp[bp_id].append(record)
        s["target_hint"] = record.target_hint or ""

    for bp_id, mapping in seed_overrides.items():
        if bp_id not in records_by_bp:
            continue
        records_by_bp[bp_id] = apply_seed_overrides(records_by_bp[bp_id], mapping)
        logger.info("[seed] applied manual_round_mapping for %s (%d SRRs)", bp_id, len(mapping))

    round_rows = []
    for _bp_id, recs in records_by_bp.items():
        for r in recs:
            round_rows.append(r.to_dict())
    return round_rows


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------


def deduplicate_bioprojects(rows: list[dict]) -> list[dict]:
    """Merge duplicate BioProject rows from multiple sources.

    First occurrence wins on bp_id; subsequent rows fill in any blank columns.
    Order in `rows` matters — pass seed-derived rows first to give them priority.
    """
    seen: dict[str, dict] = {}
    for bp in rows:
        key = bp["bioproject_id"]
        if key not in seen:
            seen[key] = bp
        else:
            existing = seen[key]
            for col in BIOPROJECT_COLS:
                if not existing.get(col) and bp.get(col):
                    existing[col] = bp[col]
    return list(seen.values())


def deduplicate_samples(rows: list[dict]) -> list[dict]:
    """First-wins dedup by `srr` accession."""
    seen: dict[str, dict] = {}
    for s in rows:
        key = s["srr"]
        if key and key not in seen:
            seen[key] = s
    return list(seen.values())


# ---------------------------------------------------------------------------
# CSV writer
# ---------------------------------------------------------------------------


def write_csv(rows: list[dict], cols: list[str], path: Path) -> None:
    """Write `rows` to a CSV at `path` with the given column order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Wrote %d rows → %s", len(rows), path)


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def run_discovery(
    output_dir: Path,
    seed_file: Path,
    sources: list[str] | None = None,
    query: str | None = None,
    seed_only: bool = False,
) -> tuple[list[dict], list[dict], list[dict]]:
    """Run discovery across all (or selected) sources.

    Returns ``(bioproject_rows, sample_rows, round_rows)`` and writes them to
    ``bioprojects.csv``, ``samples.csv``, ``rounds.csv`` under ``output_dir``.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    manual_review_dir = output_dir / "manual_review"

    blacklist = load_blacklist(seed_file)
    sm_targets = load_small_molecule_targets(seed_file)

    adapters: list[SourceAdapter] = [SeedAdapter(seed_file)]
    if not seed_only:
        all_adapters: list[SourceAdapter] = [
            ENAAdapter(),
            SRAAdapter(),
            GEOAdapter(),
            ZenodoAdapter(),
            FigshareAdapter(),
        ]
        if sources:
            all_adapters = [a for a in all_adapters if a.name in sources]
        adapters += all_adapters

    all_bp: list[dict] = []
    all_samples: list[dict] = []

    for adapter in adapters:
        logger.info("Running adapter: %s", adapter.name)
        try:
            bp_rows, sample_rows = adapter.search(query=query)
        except Exception as e:
            logger.error("[%s] adapter failed: %s", adapter.name, e)
            continue
        all_bp.extend(bp_rows)
        all_samples.extend(sample_rows)

    all_bp = deduplicate_bioprojects(all_bp)
    all_samples = deduplicate_samples(all_samples)

    # SRR lookup for BioProjects without samples (GEO/seed/Zenodo/Figshare)
    srr_bp_ids = {s["bioproject_id"] for s in all_samples if s.get("srr")}
    missing = [
        bp
        for bp in all_bp
        if bp["bioproject_id"] not in srr_bp_ids
        and not bp["bioproject_id"].startswith(("zenodo:", "figshare:", "utexas:"))
    ]
    if missing:
        logger.info("SRR lookup for %d BioProjects without samples ...", len(missing))
        for bp in missing:
            bp_id = bp["bioproject_id"]
            if bp_id.startswith(("SRP", "ERP", "DRP")):
                query_field = "secondary_study_accession"
            elif bp_id.startswith(("PRJNA", "PRJEB", "PRJDB")):
                query_field = "study_accession"
            elif bp_id.startswith("GSE"):
                query_field = "geo_accession"
            else:
                logger.info("  [skip] unknown accession format: %s", bp_id)
                continue

            runs = _get(
                ENA_PORTAL_URL,
                params={
                    "result": "read_run",
                    "query": f'{query_field}="{bp_id}"',
                    "fields": ENA_RUN_FIELDS,
                    "limit": 5000,
                    "format": "json",
                },
            )
            if not runs:
                logger.info("  [%s] no runs found via %s", bp_id, query_field)
                continue
            logger.info("  [%s] found %d runs", bp_id, len(runs))
            for run in runs:
                s = empty_sample()
                s["srr"] = run.get("run_accession", "")
                s["bioproject_id"] = bp_id
                s["sample_title"] = run.get("sample_title", "") or run.get("sample_alias", "")
                s["library_name"] = run.get("library_name", "")
                s["experiment_title"] = run.get("experiment_title", "")
                s["total_bases"] = str(run.get("base_count", ""))
                s["n_reads"] = str(run.get("read_count", ""))
                s["raw_metadata"] = json.dumps(run)
                all_samples.append(s)
            time.sleep(0.3)
        all_samples = deduplicate_samples(all_samples)
        logger.info("After SRR lookup: %d total samples", len(all_samples))

    n_before = len(all_bp)
    all_bp = [bp for bp in all_bp if not is_blacklisted(bp, blacklist, sm_targets)]
    logger.info("Blacklist/SM filter: %d → %d BioProjects", n_before, len(all_bp))

    all_bp = _classify_all(all_bp)

    bp_abstracts = {bp["bioproject_id"]: bp.get("abstract", "") for bp in all_bp}
    seed_overrides = load_seed_overrides(seed_file)
    all_rounds = _parse_rounds(
        all_samples,
        manual_review_dir,
        bp_abstracts=bp_abstracts,
        seed_overrides=seed_overrides,
    )

    write_csv(all_bp, BIOPROJECT_COLS, output_dir / "bioprojects.csv")
    write_csv(all_samples, SAMPLE_COLS, output_dir / "samples.csv")
    write_csv(all_rounds, ROUND_COLS, output_dir / "rounds.csv")

    rounds_curated = output_dir / "rounds_curated.csv"
    if not rounds_curated.exists():
        shutil.copyfile(output_dir / "rounds.csv", rounds_curated)
        logger.info(
            "Bootstrapped %s from rounds.csv; edit for manual round curation",
            rounds_curated,
        )

    n_rna = sum(1 for bp in all_bp if bp.get("library_type_verification") == "RNA_confirmed")
    n_ambig = sum(1 for bp in all_bp if bp.get("library_type_verification") == "ambiguous")
    n_manual = sum(1 for r in all_rounds if r.get("needs_manual_review"))
    logger.info(
        "Discovery complete: %d BioProjects (%d RNA_confirmed, %d ambiguous), "
        "%d samples, %d SRRs need manual round review",
        len(all_bp),
        n_rna,
        n_ambig,
        len(all_samples),
        n_manual,
    )

    return all_bp, all_samples, all_rounds
