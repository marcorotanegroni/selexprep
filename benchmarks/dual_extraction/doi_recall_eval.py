"""Internal evaluation: did source-grounded dual extraction recover the API DOIs?

Compares the paper DOIs that the two independent extraction arms
(``claude_extractions.json`` / ``codex_extractions.json``) identified per deposit
against the deterministic accession->DOI links the discovery catalog got from the
APIs (``bioprojects.csv:paper_doi``). This is an INTERNAL evaluation of whether
reading the sources paid off — NOT a dataset field, and nothing here is
back-filled into the curated metadata.

Metric, over the deposits that were dual-extracted:
  - recovered_same_doi : an API-DOI deposit where >=1 arm's paper DOI equals the
                         API DOI (normalized) -> recall of the API links.
  - beyond_api         : a deposit with NO API DOI where >=1 arm found a paper.

Run:  python benchmarks/dual_extraction/doi_recall_eval.py
Reads only committed files; writes doi_recall_eval.json next to this script.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
CATALOG = HERE.parent.parent / "src" / "selexprep" / "catalog" / "data" / "bioprojects.csv"


# figshare/institutional-repository DOI prefixes: these identify a *dataset/asset*
# deposit, not a journal article. When the catalog's paper_doi is one of these, the
# API never had a real paper link — so an extraction that finds an article DOI here
# is ADDING a paper, not disagreeing.
_REPO_PREFIXES = (
    "10.6084/m9.figshare",  # figshare
    "10.26686/wgtn",  # Victoria U. Wellington repository
    "10.48546",
    "10.5281/zenodo",  # Zenodo
    "10.5061/dryad",  # Dryad data repository
    "10.17811",  # RUO / institutional dataset repository
)


def is_repository_doi(raw: str | None) -> bool:
    n = norm_doi(raw, strip_asset=False)
    return n.startswith(_REPO_PREFIXES)


def norm_doi(raw: str | None, strip_asset: bool = True) -> str:
    """Normalize a DOI for comparison: lowercase bare ``10.xxxx/...`` form.

    ``strip_asset`` removes trailing figshare/publisher sub-asset suffixes so a
    figure/table/supplement DOI (``...0097574.g002``, ``...t001``, ``.s001``,
    ``.v1``) collapses to its parent *article* DOI — the two then compare equal.
    """
    if not raw:
        return ""
    s = raw.strip().lower()
    s = re.sub(r"^(https?://)?(dx\.)?doi\.org/", "", s)
    s = s.strip().rstrip(".").strip()
    if strip_asset:
        # iteratively peel .v1 (version) and .g2/.t1/.s1/.d1/.e1/.f1 (figure/table/
        # supplement/data/extended/file) asset suffixes appended to an article DOI.
        while True:
            new = re.sub(r"\.(v|g|t|s|d|e|f)\d+$", "", s)
            if new == s:
                break
            s = new
    return s


def paper_doi_of(rec: dict) -> str:
    for sd in rec.get("source_documents", []):
        if sd.get("type") == "paper":
            return norm_doi(sd.get("id"))
    return ""


def load_arm(path: Path) -> dict[str, str]:
    return {r["accession"]: paper_doi_of(r) for r in json.loads(path.read_text())}


def main() -> int:
    import csv

    claude = load_arm(HERE / "claude_extractions.json")
    codex = load_arm(HERE / "codex_extractions.json")
    extracted = set(claude) | set(codex)

    api_raw: dict[str, str] = {}
    with CATALOG.open() as fh:
        for row in csv.DictReader(fh):
            acc = row["bioproject_id"]
            if acc in extracted:
                api_raw[acc] = (row.get("paper_doi") or "").strip()

    # split the API paper_doi values into real article DOIs vs repository/dataset DOIs
    api_paper = {a: norm_doi(d) for a, d in api_raw.items() if d and not is_repository_doi(d)}
    api_repo = {a for a, d in api_raw.items() if d and is_repository_doi(d)}
    api_none = extracted - set(api_paper) - api_repo

    def arm_doi(a: str) -> set[str]:
        return {d for d in (claude.get(a, ""), codex.get(a, "")) if d}

    # recall of the API's *article* DOIs (after asset-suffix normalization)
    recovered_same = {a for a in api_paper if api_paper[a] in arm_doi(a)}
    # extraction upgraded a repository/dataset DOI to a real article
    repo_upgraded = {a for a in api_repo if arm_doi(a)}
    # deposits the API had no DOI for at all, where extraction found a paper
    beyond = {a for a in api_none if arm_doi(a)}

    report = {
        "n_dual_extracted": len(extracted),
        "n_api_article_doi": len(api_paper),
        "recovered_same_article_doi": len(recovered_same),
        "recovered_same_article_doi_pct": round(
            100.0 * len(recovered_same) / max(1, len(api_paper)), 1
        ),
        "n_api_repository_doi": len(api_repo),
        "repo_doi_upgraded_to_article": len(repo_upgraded),
        "n_api_no_doi": len(api_none),
        "beyond_api_found_paper": len(beyond),
        "papers_added_by_extraction": len(repo_upgraded) + len(beyond),
        "not_recovered_same_article_doi": sorted(set(api_paper) - recovered_same),
    }
    (HERE / "doi_recall_eval.json").write_text(json.dumps(report, indent=2) + "\n")
    for k, v in report.items():
        if k != "not_recovered_same_article_doi":
            print(f"{k:32s}: {v}")
    print(
        f"\n(genuinely-different / not-matched article DOIs: "
        f"{len(report['not_recovered_same_article_doi'])} accessions)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
