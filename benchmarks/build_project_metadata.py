"""Build the static per-project SELEX metadata table shipped with the benchmark.

For every SELEX deposit in the bundled catalog this emits one structured row of
*experiment characteristics* — the kind of thing a reader would otherwise have
to dig out of the paper or the raw ENA record by hand. It is a **build
artifact**, joined from sources that each own exactly one slice of the truth so
they never drift:

- ``ground_truth.tsv`` + ``project_annotations.tsv`` — the **verified** benchmark
  layer (chemistry, target class, random-region length, target identity,
  selection format, counter-selection). Hand-checked against each paper's M&M.
- ``catalog_annotations.tsv`` *(optional)* — target / selection format / rounds
  **extracted from the fetched study or paper abstract**, each value traceable
  to its source text. Curated by hand because a small local LLM is not reliable
  enough; this file is also the eval-set for that LLM later.
- ``bioprojects.csv`` (the catalog) — deterministic record fields: study title,
  target organism, declared rounds, paper DOI/PMID.
- OpenAlex (DOI -> open access) — the paper-obtainability flag.

**Every row carries two honesty signals:**

1. ``curation_level`` — ``verified`` (benchmark, paper-checked) /
   ``extracted`` (abstract-derived, review-grade) / ``none`` (only the
   deterministic record layer). A reader always knows how a value was obtained.
2. ``metadata_tier`` — RECORD_ONLY / ABSTRACT / FULL_TEXT. The paper is the axis
   along which most of this metadata is even obtainable, so the gap is a
   labelled, filterable field — not a silent hole.

*Cardinal rule, same as the round-assignment cascade:* a blank cell means "not
stated in the source", never a guessed value.

The per-round enrichment trajectory (n_reads / n_unique / singleton_frac per
round) is populated into each deposit's ``rounds`` slot from ``--results-dir``
(a ``selexprep run`` output tree); deposits with no count run keep
``rounds: null``. The trajectory is nested, so it lives in the JSON only — the
CSV stays flat.

**Output:**

- ``benchmarks/project_metadata.csv``  — one flat row per deposit.
- ``benchmarks/project_metadata.json`` — same rows, plus a per-round ``rounds``
  trajectory (null until ``--results-dir`` is supplied).

Usage::

    python -m benchmarks.build_project_metadata                  # bundled defaults, hits OpenAlex
    python -m benchmarks.build_project_metadata --no-network     # offline; tier capped at ABSTRACT
    python -m benchmarks.build_project_metadata --results-dir out/  # + per-round trajectory
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

_HERE = Path(__file__).resolve().parent
_CATALOG = _HERE.parent / "src" / "selexprep" / "catalog" / "data" / "bioprojects.csv"
DEFAULT_GROUND_TRUTH = _HERE / "ground_truth.tsv"
DEFAULT_ANNOTATIONS = _HERE / "project_annotations.tsv"
DEFAULT_CATALOG = _CATALOG
DEFAULT_CATALOG_ANNOTATIONS = _HERE / "catalog_annotations.tsv"
DEFAULT_OUT_DIR = _HERE
DEFAULT_CACHE = _HERE / ".oa_cache.json"
OPENALEX_API = "https://api.openalex.org/works"

# Flat CSV column order (scalar fields only; trajectory lives in the JSON).
COLUMNS = [
    "accession",
    "study_type",
    "target",
    "target_class",
    "chemistry",
    "n_random",
    "selection_format",
    "counter_selection",
    "target_organism",
    "n_rounds",
    "study_title",
    "paper_doi",
    "paper_pmid",
    "paper_linked",
    "paper_oa",
    "paper_oa_url",
    "metadata_tier",
    "curation_level",
]


# ---------------------------------------------------------------------------
# Source readers
# ---------------------------------------------------------------------------
def _read_delimited(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def load_descriptors(path: Path) -> dict[str, dict[str, str]]:
    """ground_truth.tsv -> verified library descriptors, keyed by accession."""
    out: dict[str, dict[str, str]] = {}
    for row in _read_delimited(path, "\t"):
        acc = row["accession"].strip()
        out[acc] = {
            "chemistry": row.get("library_kind", "").strip(),
            "target_class": row.get("target_kind", "").strip(),
            "n_random": row.get("n_length_truth", "").strip(),
            "paper_doi": row.get("paper_doi", "").strip(),
            "paper_pmid": row.get("paper_pmid", "").strip(),
        }
    return out


def load_annotations(path: Path) -> dict[str, dict[str, str]]:
    """project_annotations.tsv / catalog_annotations.tsv -> prose-derived fields."""
    if not path.exists():
        return {}
    out: dict[str, dict[str, str]] = {}
    for row in _read_delimited(path, "\t"):
        acc = (row.get("accession") or row.get("bioproject_id") or "").strip()
        if not acc:
            continue
        out[acc] = {
            "study_type": row.get("study_type", "").strip(),
            "target": row.get("target", "").strip(),
            "selection_format": row.get("selection_format", "").strip(),
            "counter_selection": row.get("counter_selection", "").strip(),
            "n_rounds": row.get("n_rounds", "").strip(),
            # discovery-only deposits (figshare/zenodo) carry the paper DOI resolved
            # from their host API — the catalog never populated it for them.
            "paper_doi": row.get("paper_doi", "").strip(),
        }
    return out


def load_catalog(path: Path) -> dict[str, dict[str, str]]:
    """bioprojects.csv -> deterministic record fields, keyed by bioproject id."""
    out: dict[str, dict[str, str]] = {}
    for row in _read_delimited(path, ","):
        acc = row["bioproject_id"].strip()
        out[acc] = {
            "study_title": row.get("study_title", "").strip(),
            "protein_target": row.get("protein_target", "").strip(),
            "target_organism": row.get("target_organism", "").strip(),
            "n_rounds_declared": row.get("n_rounds_declared", "").strip(),
            "paper_doi": row.get("paper_doi", "").strip(),
            "paper_pmid": row.get("paper_pmid", "").strip(),
        }
    return out


# ---------------------------------------------------------------------------
# Paper obtainability (DOI -> open access), with an on-disk cache
# ---------------------------------------------------------------------------
def _fetch_oa(doi: str, *, mailto: str | None, timeout: float) -> tuple[bool | None, str | None]:
    params = {"mailto": mailto} if mailto else {}
    url = f"{OPENALEX_API}/doi:{urllib.parse.quote(doi)}"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "selexprep-benchmark"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.load(resp)
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        logger.warning("OpenAlex lookup failed for %s: %s", doi, exc)
        return None, None
    oa = data.get("open_access") or {}
    return bool(oa.get("is_oa")), oa.get("oa_url")


def metadata_tier(paper_doi: str, is_oa: bool | None) -> str:
    """RECORD_ONLY (no paper) < ABSTRACT (paper, not OA/unknown) < FULL_TEXT (OA)."""
    if not paper_doi:
        return "RECORD_ONLY"
    if is_oa:
        return "FULL_TEXT"
    return "ABSTRACT"


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------
def _first(*values: str) -> str:
    """First non-empty value (precedence resolution)."""
    for v in values:
        if v:
            return v
    return ""


def build_rows(
    *,
    catalog: dict[str, dict[str, str]],
    descriptors: dict[str, dict[str, str]],
    project_ann: dict[str, dict[str, str]],
    catalog_ann: dict[str, dict[str, str]],
    oa: dict[str, tuple[bool | None, str | None]],
) -> list[dict[str, Any]]:
    accessions = sorted(set(catalog) | set(descriptors))
    rows: list[dict[str, Any]] = []
    for acc in accessions:
        cat = catalog.get(acc, {})
        gt = descriptors.get(acc, {})
        verified = acc in descriptors
        ann = project_ann.get(acc, {}) if verified else catalog_ann.get(acc, {})
        if verified:
            curation = "verified"
        elif acc in catalog_ann:
            curation = "extracted"
        else:
            curation = "none"

        doi = _first(ann.get("paper_doi", ""), gt.get("paper_doi", ""), cat.get("paper_doi", ""))
        is_oa, oa_url = oa.get(doi, (None, None))
        rows.append(
            {
                "accession": acc,
                "study_type": ann.get("study_type", ""),
                "target": _first(ann.get("target", ""), cat.get("protein_target", "")),
                "target_class": gt.get("target_class", ""),
                "chemistry": gt.get("chemistry", ""),
                "n_random": gt.get("n_random", ""),
                "selection_format": ann.get("selection_format", ""),
                "counter_selection": ann.get("counter_selection", ""),
                "target_organism": cat.get("target_organism", ""),
                "n_rounds": _first(ann.get("n_rounds", ""), cat.get("n_rounds_declared", "")),
                "study_title": cat.get("study_title", ""),
                "paper_doi": doi,
                "paper_pmid": _first(gt.get("paper_pmid", ""), cat.get("paper_pmid", "")),
                "paper_linked": bool(doi),
                "paper_oa": is_oa,
                "paper_oa_url": oa_url or "",
                "metadata_tier": metadata_tier(doi, is_oa),
                "curation_level": curation,
            }
        )
    return rows


def resolve_all_oa(
    dois: set[str],
    *,
    use_network: bool,
    mailto: str | None,
    timeout: float,
    cache_path: Path,
) -> dict[str, tuple[bool | None, str | None]]:
    """Resolve OA status for a set of DOIs, reusing/refreshing an on-disk cache."""
    cache: dict[str, list[Any]] = {}
    if cache_path.exists():
        cache = json.loads(cache_path.read_text(encoding="utf-8"))
    resolved: dict[str, tuple[bool | None, str | None]] = {
        doi: (val[0], val[1]) for doi, val in cache.items()
    }
    if use_network:
        missing = sorted(d for d in dois if d and d not in resolved)
        for i, doi in enumerate(missing, 1):
            resolved[doi] = _fetch_oa(doi, mailto=mailto, timeout=timeout)
            if i % 25 == 0:
                logger.info("  OA lookup %d/%d", i, len(missing))
        if missing:
            cache_path.write_text(
                json.dumps({d: list(v) for d, v in resolved.items()}, indent=0),
                encoding="utf-8",
            )
            logger.info(
                "OA: resolved %d new, %d cached", len(missing), len(resolved) - len(missing)
            )
    return resolved


def _csv_cell(value: Any) -> str:
    """Render bool as true/false, None as empty — stable, greppable CSV."""
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(COLUMNS)
        for row in rows:
            writer.writerow([_csv_cell(row[c]) for c in COLUMNS])


def load_trajectory(results_dir: Path, accession: str) -> list[dict[str, Any]] | None:
    """Per-round enrichment trajectory for one deposit from ``selexprep run`` outputs.

    Reads ``<results_dir>/<accession>/round_*/counts.parquet`` (the layout the
    ``run`` driver writes) and recomputes per-round depth/diversity from the
    ``reads`` column — the same derivation the counter uses for cached rounds.
    Returns ``None`` when no per-round counts exist for the accession, so the
    JSON keeps ``rounds: null`` for deposits that were never counted.
    """
    acc_dir = results_dir / accession
    if not acc_dir.is_dir():
        return None
    import pandas as pd  # lazy: only needed when --results-dir is supplied

    rounds: list[dict[str, Any]] = []
    for round_dir in sorted(acc_dir.glob("round_*")):
        parquet = round_dir / "counts.parquet"
        if not parquet.is_file():
            continue
        reads = pd.read_parquet(parquet, columns=["reads"])["reads"].to_numpy()
        n_unique = int(reads.size)
        n_singletons = int((reads == 1).sum())
        label = round_dir.name.split("_", 1)[1]  # "00".."NN" or "unknown"
        rounds.append(
            {
                "round": int(label) if label.isdigit() else label,
                "n_reads": int(reads.sum()),
                "n_unique": n_unique,
                "singleton_frac": (n_singletons / n_unique) if n_unique else 0.0,
            }
        )
    if not rounds:
        return None
    # numeric rounds first (ascending), then any string labels (e.g. "unknown")
    rounds.sort(key=lambda r: (isinstance(r["round"], str), r["round"]))
    return rounds


def write_json(rows: list[dict[str, Any]], path: Path) -> None:
    # ``rounds`` carries the per-round trajectory, or null for un-counted deposits.
    payload = [{**row, "rounds": row.get("rounds")} for row in rows]
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    parser.add_argument("--ground-truth", type=Path, default=DEFAULT_GROUND_TRUTH)
    parser.add_argument("--annotations", type=Path, default=DEFAULT_ANNOTATIONS)
    parser.add_argument("--catalog-annotations", type=Path, default=DEFAULT_CATALOG_ANNOTATIONS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument(
        "--no-network",
        action="store_true",
        help="skip the OpenAlex OA lookup; metadata_tier is then capped at ABSTRACT for linked papers",
    )
    parser.add_argument("--mailto", default=None, help="contact email for the OpenAlex polite pool")
    parser.add_argument("--timeout", type=float, default=15.0)
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=None,
        help="a `selexprep run` output dir; populates each deposit's per-round trajectory "
        "from <results-dir>/<accession>/round_*/counts.parquet (else rounds stays null)",
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    catalog = load_catalog(args.catalog)
    descriptors = load_descriptors(args.ground_truth)
    project_ann = load_annotations(args.annotations)
    catalog_ann = load_annotations(args.catalog_annotations)

    dois = {
        d
        for src in (catalog.values(), descriptors.values())
        for r in src
        if (d := r.get("paper_doi", ""))
    }
    oa = resolve_all_oa(
        dois,
        use_network=not args.no_network,
        mailto=args.mailto,
        timeout=args.timeout,
        cache_path=args.cache,
    )

    rows = build_rows(
        catalog=catalog,
        descriptors=descriptors,
        project_ann=project_ann,
        catalog_ann=catalog_ann,
        oa=oa,
    )

    if args.results_dir is not None:
        n_with = 0
        for row in rows:
            row["rounds"] = load_trajectory(args.results_dir, row["accession"])
            n_with += row["rounds"] is not None
        logger.info("trajectory: per-round counts attached for %d/%d deposits", n_with, len(rows))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.out_dir / "project_metadata.csv"
    json_path = args.out_dir / "project_metadata.json"
    write_csv(rows, csv_path)
    write_json(rows, json_path)

    def _counts(field: str, keys: tuple[str, ...]) -> str:
        c = {k: sum(r[field] == k for r in rows) for k in keys}
        return ", ".join(f"{k}={v}" for k, v in c.items())

    logger.info("wrote %s and %s (%d deposits)", csv_path.name, json_path.name, len(rows))
    logger.info("curation_level: %s", _counts("curation_level", ("verified", "extracted", "none")))
    logger.info(
        "metadata_tier:  %s", _counts("metadata_tier", ("FULL_TEXT", "ABSTRACT", "RECORD_ONLY"))
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
