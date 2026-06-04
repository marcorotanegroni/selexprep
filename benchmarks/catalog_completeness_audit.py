"""one-off catalog discovery completeness audit.

Measures what fraction of ENA studies explicitly typed as
``library_strategy="SELEX"`` are accounted for by selexprep's bundled
catalog (either present in ``bioprojects.csv`` or recorded in
``bioprojects_excluded.csv`` with a documented reason). This is the
audit referenced from Figure B's narrative — the discovery-side
counterpart to 's fetcher / inference audit.

**Not a recurring pipeline.** Run this when:

- You've just refreshed the catalog and want to verify the discovery
  layer caught everything ENA types as ``library_strategy="SELEX"``.
- A paper reviewer asks for a coverage number that can be re-derived
  from public APIs.
- A future ENA refresh adds new SELEX-tagged deposits and you want to
  know which ones slipped through.

**Output:**

- ``benchmarks/catalog_completeness_audit.json`` — structured audit
  report with the ENA count, present / excluded / unaccounted-for
  breakdowns, missing accession list, and coverage percentages.
- ``benchmarks/catalog_completeness_audit.tsv`` — per-accession
  diff (one row per ENA-SELEX accession with its catalog status).

**Methodology**:

1. ENA Portal API at the data-type level:
   ``result=read_run&query=library_strategy="SELEX"``. This is the
   structured-metadata query, orthogonal to the text-pattern queries
   selexprep's discovery layer normally uses.
2. Deduplicate to unique ``study_accession`` values.
3. Diff against ``bioprojects.csv`` (kept) and
   ``bioprojects_excluded.csv`` (filtered).
4. Anything in the ENA set but neither file is an unaccounted-for
   discovery miss — surface in the missing list.

**Note on ``AMPLICON``.** A second arm of the original audit also
checked ``library_strategy="AMPLICON"`` intersected with SELEX/aptamer
text queries; that check produced 26/26 already-covered and is
therefore not a discovery gap. This script doesn't repeat it (the
catalog refresh's text-pattern queries already catch that subset).

Usage:

    python -m benchmarks.catalog_completeness_audit \\
        --catalog src/selexprep/catalog/data/bioprojects.csv \\
        --excluded src/selexprep/catalog/data/bioprojects_excluded.csv \\
        --out-dir benchmarks/

Or with all defaults (uses the bundled catalog snapshot):

    python -m benchmarks.catalog_completeness_audit
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

from selexprep.catalog.filter import is_insdc_accession
from selexprep.catalog.reader import catalog_path, catalog_version

ENA_URL = (
    "https://www.ebi.ac.uk/ena/portal/api/search"
    "?result=read_run"
    "&query=" + urllib.parse.quote('library_strategy="SELEX"') + "&fields=study_accession"
    "&format=tsv"
    "&limit=100000"
)


def fetch_ena_selex_studies(timeout_s: int = 60) -> set[str]:
    """Hit ENA at the structured-metadata level and return unique study_accessions.

    Uses ``library_strategy="SELEX"`` rather than text patterns so the
    result is orthogonal to selexprep's text-based discovery layer —
    that orthogonality is the whole point of this audit.
    """
    req = urllib.request.Request(ENA_URL, headers={"User-Agent": "selexprep-audit/0.1"})
    with urllib.request.urlopen(req, timeout=timeout_s) as fh:
        text = fh.read().decode("utf-8")
    studies: set[str] = set()
    for line in text.strip().splitlines()[1:]:  # skip header
        parts = line.split("\t")
        if len(parts) >= 2 and parts[1]:
            studies.add(parts[1])
    return studies


def audit_coverage(
    catalog_csv: Path,
    excluded_csv: Path,
    *,
    timeout_s: int = 60,
) -> dict:
    """Diff ENA's library_strategy=SELEX set against the bundled catalog.

    Returns a structured report with the ENA count, present / excluded /
    unaccounted breakdowns, the missing accession list, and percentages.
    """
    ena_studies = fetch_ena_selex_studies(timeout_s=timeout_s)

    cat = pd.read_csv(catalog_csv)
    catalog_ids = set(cat["bioproject_id"].dropna().astype(str).tolist())

    if excluded_csv.exists():
        exc = pd.read_csv(excluded_csv)
        excluded_ids = set(exc["bioproject_id"].dropna().astype(str).tolist())
    else:
        excluded_ids = set()

    present = ena_studies & catalog_ids
    in_excluded = ena_studies & excluded_ids
    unaccounted = ena_studies - catalog_ids - excluded_ids

    # All ENA-SELEX studies should be INSDC by construction; sanity-check.
    non_insdc = [a for a in ena_studies if not is_insdc_accession(a)]

    return {
        "catalog_version": catalog_version(),
        "ena_query_url": ENA_URL,
        "n_ena_selex_studies": len(ena_studies),
        "n_present_in_catalog": len(present),
        "n_present_in_excluded_sidecar": len(in_excluded),
        "n_unaccounted_for": len(unaccounted),
        "discovery_coverage_pct": round(
            100.0 * (len(present) + len(in_excluded)) / max(1, len(ena_studies)), 2
        ),
        "auditable_coverage_pct": round(100.0 * len(present) / max(1, len(ena_studies)), 2),
        "missing_accessions": sorted(unaccounted),
        "ena_non_insdc_count": len(non_insdc),
    }


def write_per_accession_tsv(
    catalog_csv: Path,
    excluded_csv: Path,
    out_path: Path,
    *,
    timeout_s: int = 60,
) -> None:
    """Per-accession diff: one row per ENA-SELEX study with catalog status.

    Columns: ``accession``, ``in_catalog``, ``in_excluded``, ``in_neither``,
    ``catalog_study_title``, ``exclusion_reason``. Sorted by accession for
    deterministic diffs.
    """
    ena_studies = fetch_ena_selex_studies(timeout_s=timeout_s)
    cat = pd.read_csv(catalog_csv).set_index("bioproject_id")
    exc = (
        pd.read_csv(excluded_csv).set_index("bioproject_id")
        if excluded_csv.exists()
        else pd.DataFrame()
    )

    rows = []
    for acc in sorted(ena_studies):
        in_catalog = acc in cat.index
        in_excluded = acc in exc.index
        rows.append(
            {
                "accession": acc,
                "in_catalog": in_catalog,
                "in_excluded": in_excluded,
                "in_neither": not (in_catalog or in_excluded),
                "catalog_study_title": (str(cat.at[acc, "study_title"]) if in_catalog else ""),
                "exclusion_reason": (str(exc.at[acc, "exclusion_reason"]) if in_excluded else ""),
            }
        )
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    p.add_argument(
        "--catalog",
        type=Path,
        default=None,
        help="Path to bioprojects.csv. Defaults to the bundled snapshot.",
    )
    p.add_argument(
        "--excluded",
        type=Path,
        default=None,
        help=("Path to bioprojects_excluded.csv. Defaults to the sidecar next to the catalog."),
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path("benchmarks"),
        help="Where to write the JSON + per-accession TSV.",
    )
    p.add_argument(
        "--timeout-s",
        type=int,
        default=60,
        help="HTTP timeout for the ENA query.",
    )
    args = p.parse_args(argv)

    catalog_csv = args.catalog if args.catalog is not None else catalog_path()
    excluded_csv = (
        args.excluded
        if args.excluded is not None
        else catalog_csv.parent / "bioprojects_excluded.csv"
    )

    print(f"Catalog:  {catalog_csv}", file=sys.stderr)
    print(f"Excluded: {excluded_csv}", file=sys.stderr)
    print(f"ENA URL:  {ENA_URL}", file=sys.stderr)

    report = audit_coverage(catalog_csv, excluded_csv, timeout_s=args.timeout_s)
    print(
        f"\nENA library_strategy=SELEX studies: {report['n_ena_selex_studies']}\n"
        f"  in bioprojects.csv (kept):   {report['n_present_in_catalog']}\n"
        f"  in bioprojects_excluded.csv: {report['n_present_in_excluded_sidecar']}\n"
        f"  unaccounted (still missing): {report['n_unaccounted_for']}\n"
        f"\nDiscovery coverage: {report['discovery_coverage_pct']}%\n"
        f"Auditable coverage: {report['auditable_coverage_pct']}%",
        file=sys.stderr,
    )
    if report["missing_accessions"]:
        print(
            f"\nMissing accessions: {' '.join(report['missing_accessions'])}",
            file=sys.stderr,
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    json_out = args.out_dir / "catalog_completeness_audit.json"
    tsv_out = args.out_dir / "catalog_completeness_audit.tsv"

    json_out.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"\nwrote {json_out}", file=sys.stderr)

    write_per_accession_tsv(catalog_csv, excluded_csv, tsv_out, timeout_s=args.timeout_s)
    print(f"wrote {tsv_out}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
