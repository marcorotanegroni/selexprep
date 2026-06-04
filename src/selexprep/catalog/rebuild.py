"""Rebuild the bundled discovery catalog from broad ENA queries.

The catalog ships as package data — but it would go stale fast without a way
to refresh it. This module provides that path, in two layers:

- :func:`harvest_runs_from_ena` runs a deliberately broader set of queries
  than the original `selex_corpus.discover` (which was tuned for one
  researcher's thesis). The wider net catches studies the thesis-specific
  queries would have missed.
- :func:`rebuild_catalog` writes a fresh ``bioprojects.csv`` from the
  ENA result + any preserved non-INSDC entries (Zenodo / Figshare /
  processed-data deposits) from a previous catalog snapshot. Hand-enriched
  fields on known entries (``protein_target``, ``paper_doi``, ``paper_pmid``,
  ``n_rounds_declared``) are merged forward so refreshes never erase manual
  enrichment.

**library_strategy hygiene.** Queries now run at the **run
level** (``result=read_run``) instead of the study level so the per-run
``library_strategy`` field is available. Runs are grouped by study and
each study is classified via
:func:`selexprep.fetch.library_strategy.classify_study_by_library_strategies`.
Studies whose runs are 100% in
:data:`~selexprep.fetch.library_strategy.LIBRARY_STRATEGY_BLOCKLIST` are
dropped from the catalog and their exclusion reason is recorded in a
sidecar ``bioprojects_excluded.csv`` next to the main catalog file.
Mixed studies (some compatible runs + some blocklisted) are kept; the
audit-eligibility layer classifies those at audit time.

**Curation is intentionally NOT in scope.** The catalog reflects the public
archives, not any single researcher's "include/exclude" decisions. Per
the package's design philosophy, curation is the user's job downstream.
"""

from __future__ import annotations

import csv
import json
import logging
import time
import urllib.parse
import urllib.request
from pathlib import Path

from selexprep.fetch.library_strategy import (
    StudyStrategyClassification,
    classify_study_by_library_strategies,
)

logger = logging.getLogger(__name__)


ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/search"

#: Broader ENA queries than the thesis-specific selex_corpus.discover.py.
#: Maximises INSDC coverage of public SELEX/aptamer studies at the cost of
#: some noise (an "aptamer" mention in a study_title doesn't strictly mean
#: an HT-SELEX raw-reads deposit). The library_strategy
#: filter (per-run + per-study aggregation) drops obvious non-SELEX
#: false positives at refresh time; the audit-eligibility
#: layer classifies the rest. Users filter further at the
#: `selexprep catalog list` stage.
#:
#: ENA's quoted-phrase syntax is **exact-token**:
#: ``study_title="aptamer"`` does NOT match a title containing
#: ``aptamers`` (the plural is a distinct token). Both singular and
#: plural variants are listed below so studies like PRJNA935703 ("DNA
#: aptamers for detection of pyoverdine pf5") aren't silently dropped
#: from the catalog. Verified empirically against ENA Portal on
#: 2026-05-24.
ENA_QUERIES: tuple[str, ...] = (
    'study_title="HT-SELEX"',
    'study_title="SELEX-seq"',
    'study_title="SELEX"',
    'study_title="aptamer"',
    'study_title="aptamers"',
    'study_title="Cell-SELEX"',
    'study_title="aptamer selection"',
    'study_title="systematic evolution"',
    'study_title="DNA aptamer"',
    'study_title="DNA aptamers"',
    'study_title="RNA aptamer"',
    'study_title="RNA aptamers"',
    'description="HT-SELEX"',
    'description="SELEX-seq"',
    'description="aptamer selection"',
    'description="SELEX rounds"',
    # data-type positive query. The catalog
    # completeness audit found that 49.5% (51/103) of ENA studies
    # explicitly typed as ``library_strategy="SELEX"`` were missing
    # from the text-pattern-only catalog. Adding the strategy as a
    # positive query recovers them; obvious mis-labels (e.g. Onion
    # GBS, ATAC-seq deposits with library_strategy=SELEX) are filtered
    # by ``MANUAL_EXCLUSIONS`` below. The empirical AMPLICON+text-query
    # intersection found 26/26 already in the
    # catalog, so AMPLICON does NOT need to be promoted to a positive
    # query — text matching catches that subset.
    'library_strategy="SELEX"',
)


#: accessions ENA labels as ``library_strategy="SELEX"`` but
#: that are NOT aptamer HT-SELEX deposits per per-study metadata triage.
#: These pass the library_strategy positive query above (so they'd land
#: in the catalog) but per-study evidence (library_source, library_name,
#: study_title, run-level metadata) shows they're submission-metadata
#: mis-labels of genuine non-SELEX experiments. Forced into
#: ``bioprojects_excluded.csv`` at refresh time with the reason string
#: surfaced verbatim. Adding to this dict is the v0.2 path for any new
#: mis-labels surfaced by future catalog completeness audits — the
#: the design's "record reasons, not silent deletion" discipline.
#:
#: Each entry was verified against ENA's read_run + per-study metadata
#: on 2026-05-28 via the catalog completeness audit (see
#: ``benchmarks/catalog_completeness_audit.py``).
MANUAL_EXCLUSIONS: dict[str, str] = {
    # ---- Originally identified mis-labels ----
    "PRJDB19386": (
        "Onion GBS genomic sequencing (KAP NGS support project); "
        "library_source=GENOMIC across all runs (Onion_GBS_*_GENOMIC "
        "library_names); library_strategy=SELEX is a submission-"
        "metadata error"
    ),
    "PRJDB4462": (
        "Rice (Oryza sativa) root developmental transcriptomics; "
        "library_source=TRANSCRIPTOMIC across all runs; "
        "library_strategy=SELEX is a submission-metadata error"
    ),
    "PRJEB108136": (
        "Human milk protein profiling using pre-enriched RNA-sequence "
        "libraries (custom 'Kklib' libraries); library_source=SYNTHETIC "
        "but library_selection=RT-PCR and study scope is protein "
        "profiling, not aptamer HT-SELEX"
    ),
    "PRJNA1104196": (
        "AEGIS (Artificially Expanded Genetic Information System) DNA "
        "sequencing methodology (ESEGA-seq); 50 runs profiling synthetic "
        "expanded-alphabet DNA, not aptamer HT-SELEX selection rounds"
    ),
    "PRJNA1338232": (
        "Arabidopsis thaliana NLP7 transcription-factor regulatory "
        "network study (mixed ChIP-Seq + RNA-Seq + WGS + SELEX-tagged "
        "runs across 38 samples); the SELEX-tagged subset is TF-binding "
        "selection, not aptamer HT-SELEX"
    ),
    "PRJNA577206": (
        "Small RNA-seq under carbon/nitrogen/sulfur limitation in "
        "Arabidopsis thaliana; library_source=TRANSCRIPTOMIC, "
        "library_selection='size fractionation' across all runs; "
        "library_strategy=SELEX is a submission-metadata error"
    ),
    "PRJNA608749": (
        "In vitro CRISPR enzyme binding-affinity and endonuclease-"
        "activity profiling (Cas variant 'SELEX-like' selections "
        "against substrate libraries); library_source=SYNTHETIC + "
        "library_selection=RANDOM but the selection is CRISPR off-"
        "target profiling, not aptamer HT-SELEX"
    ),
    "PRJNA751745": (
        "Glycine max (soybean) mixed ATAC-seq + RNA-seq deposit; "
        "SELEX-tagged runs have library_name='Flower ATAC-seq rep*' "
        "and library_source=GENOMIC — mis-labeled ATAC-seq libraries, "
        "not aptamer HT-SELEX"
    ),
    # ---- gSELEX / genomic-fragment SELEX ----
    # All have library_source=GENOMIC (or mixed GENOMIC+SYNTHETIC) on
    # their SELEX-tagged runs, indicating the variable region is genomic-
    # fragment-derived rather than primer-flanked synthetic random
    # oligos. selexprep's primer-inference + random-region-extraction
    # logic is built for the synthetic-random-library shape; gSELEX-
    # style libraries have a different shape (genomic fragments with
    # flanking adapters, not a fixed N-region) and would be silently
    # mis-processed.
    "PRJDB16474": (
        "Aspergillus oryzae KojR target genes via gSELEX-Seq (literal "
        "'gSELEX-Seq' in title); library_source=GENOMIC; uses genomic "
        "DNA fragments as the library rather than synthetic random "
        "oligos — out of scope for selexprep's primer-flanked-N-region "
        "extraction"
    ),
    "PRJDB4820": (
        "Genome-wide fungal transcriptional regulation (Aspergillus "
        "nidulans) via gSELEX; library_source=GENOMIC across all 4 "
        "runs; genomic-fragment SELEX library, not synthetic-random "
        "aptamer SELEX"
    ),
    "PRJDB6696": (
        "Genome-wide fungal transcriptional regulation (Aspergillus "
        "oryzae RIB40) via gSELEX; library_source=GENOMIC across all "
        "8 runs; genomic-fragment SELEX library, not synthetic-random "
        "aptamer SELEX"
    ),
    "PRJEB50170": (
        "Helicase-SELEX for E. coli Rho natural substrates; "
        "library_source=GENOMIC; selection on natural transcripts (not "
        "a synthetic random library) — out of scope for selexprep's "
        "extraction"
    ),
    "PRJEB15639": (
        "Active transcription factor identification (ATI) assay (Mus "
        "musculus); 59 SELEX-tagged runs with mixed GENOMIC + "
        "SYNTHETIC library_source — gSELEX-style genome-wide TF "
        "binding survey rather than synthetic-random aptamer SELEX"
    ),
    "PRJEB7934": (
        "Heterodimeric transcription factor complex specificities "
        "(human); 199 SELEX-tagged runs with mixed GENOMIC + SYNTHETIC "
        "library_source + 1 ChIP-Exo control — gSELEX/TF-binding "
        "survey on genomic substrates, not synthetic-random aptamer "
        "SELEX"
    ),
    "PRJEB9797": (
        "CpG methylation effects on TF binding (mouse); 194 SELEX-"
        "tagged runs with mixed GENOMIC + SYNTHETIC library_source + "
        "ChIP-Seq runs — gSELEX/methyl-SELEX on genomic substrates"
    ),
    "PRJNA272858": (
        "DNA shape recognition vs sequence (Drosophila melanogaster); "
        "32 SELEX-tagged runs, library_source=GENOMIC across all — "
        "gSELEX on Drosophila genomic fragments, not synthetic-random "
        "aptamer SELEX"
    ),
    "PRJNA378233": (
        "Arabidopsis MADS-box homeotic protein complex DNA-binding "
        "specificity; 28 SELEX runs, library_source=GENOMIC across all "
        "— gSELEX-style TF binding survey, not synthetic-random "
        "aptamer SELEX"
    ),
    "PRJNA486548": (
        "Systematic TF binding to noncoding variants (human); 197 "
        "SELEX-tagged runs ALL library_source=GENOMIC — gSELEX-style "
        "screening of genomic noncoding regions, not synthetic-random "
        "aptamer SELEX"
    ),
    "PRJNA611637": (
        "Compendium of DNA-Binding Specificities of TFs in Pseudomonas "
        "savastanoi; 184 SELEX-tagged runs ALL library_source=GENOMIC "
        "— gSELEX on bacterial genomic substrates, not synthetic-"
        "random aptamer SELEX"
    ),
    "PRJNA820961": (
        "FRUITFULL MADS-domain TF target gene selection (Arabidopsis); "
        "28 SELEX-tagged runs, library_source=GENOMIC — gSELEX-style "
        "TF binding survey, not synthetic-random aptamer SELEX"
    ),
    "PRJNA1113781": (
        "High-throughput discovery of functional RNA domains; 5 SELEX-"
        "tagged runs with mixed GENOMIC + SYNTHETIC library_source "
        "(library_names like 'FGD3_Circ', 'Were1_Pool') — selection "
        "on natural-RNA-derived substrates, not synthetic-random "
        "aptamer SELEX"
    ),
}

#: Run-level fields needed for catalog row + per-study library_strategy
#: classification. Querying at ``result=read_run`` means one row per run
#: (some studies have many); we dedupe by run_accession across queries
#: and aggregate by study_accession.
#:
#: NOTE: ENA's `result=read_run` does NOT expose
#: ``study_description`` / ``first_public`` / ``last_updated`` — those
#: are study-level fields only available with ``result=study``. As a
#: trade-off, refresh runs without abstracts for *new* entries; existing
#: entries keep their abstracts via the enrichment-preserve path
#: (:func:`_enrichment_index`). Adding a second study-level query pass
#: just for abstracts would double the ENA call count for marginal
#: benefit and is deferred.
_RUN_FIELDS = (
    "run_accession,study_accession,study_title,scientific_name,library_strategy,library_source"
)

PUBLIC_COLS = (
    "bioproject_id",
    "source",
    "study_title",
    "protein_target",
    "target_organism",
    "paper_doi",
    "paper_pmid",
    "n_rounds_declared",
    "abstract",
)

#: Schema of the exclusion sidecar
#: (``bioprojects_excluded.csv``). Emitted alongside ``bioprojects.csv``
#: so users can audit "what got filtered and why" — the design
#: forbids silent deletion.
EXCLUDED_COLS = (
    "bioproject_id",
    "source",
    "study_title",
    "n_runs_total",
    "n_runs_blocklisted",
    "blocklisted_strategies",
    "exclusion_reason",
)


def _ena_get(url: str, timeout: int = 60) -> list[dict]:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            return json.loads(resp.read())
    except Exception as e:
        logger.warning("ENA query failed: %s", e)
        return []


def harvest_runs_from_ena(
    queries: tuple[str, ...] = ENA_QUERIES,
    request_pause_seconds: float = 0.5,
) -> dict[str, dict]:
    """Run each broad query at the run level, group by study.

    Returns ``{study_accession: meta}`` where ``meta`` carries the
    study-level fields plus a ``"library_strategies"`` list with one
    entry per observed run. Run-accession-level deduplication is done
    across queries (first occurrence wins for the study-level fields).
    """
    studies: dict[str, dict] = {}
    seen_runs: set[str] = set()

    for q in queries:
        # bumped from 10000 to 100000. ENA returns up to
        # ~15000 runs for ``study_title="SELEX"`` alone (verified
        # 2026-05-24), so a 10k limit truncated and silently dropped
        # later-paginated studies (e.g. PRJEB70964 — alpha-synuclein SELEX,
        # one of the Tier 1 ground-truth rows). ``limit=0`` also means
        # "unlimited" for ENA Portal, but 100k is explicit and survives
        # any future ENA semantic change.
        params = {
            "result": "read_run",
            "query": q,
            "fields": _RUN_FIELDS,
            "limit": "100000",
            "format": "json",
        }
        url = f"{ENA_PORTAL}?{urllib.parse.urlencode(params)}"
        logger.info("[ENA] %s", q)
        data = _ena_get(url)

        new_runs_here = 0
        for r in data:
            run_acc = (r.get("run_accession") or "").strip()
            if not run_acc or run_acc in seen_runs:
                continue
            seen_runs.add(run_acc)
            study_acc = (r.get("study_accession") or "").strip()
            if not study_acc:
                continue

            if study_acc not in studies:
                studies[study_acc] = {
                    "study_accession": study_acc,
                    "study_title": (r.get("study_title") or "").strip(),
                    # ``study_description`` not available at result=read_run;
                    # preserved via _enrichment_index for previously-known
                    # entries. New entries land with empty abstract.
                    "study_description": "",
                    "scientific_name": (r.get("scientific_name") or "").strip(),
                    "library_strategies": [],
                }
            studies[study_acc]["library_strategies"].append(
                (r.get("library_strategy") or "").strip()
            )
            new_runs_here += 1

        logger.info(
            "[ENA] %s → %d runs (%d new across %d studies)",
            q,
            len(data),
            new_runs_here,
            len(studies),
        )
        time.sleep(request_pause_seconds)
    return studies


def _enrichment_index(catalog_path: Path) -> dict[str, dict]:
    """Build a {bioproject_id: row} index of preservable fields from an old catalog.

    Used by `rebuild_catalog` to merge fields forward across refreshes:

    - **Hand-enrichment** (``protein_target`` / ``paper_doi`` /
      ``paper_pmid`` / ``n_rounds_declared``) — manually curated, must
      survive a refresh.
    - **Abstract** — ENA's ``result=read_run`` query
      doesn't expose ``study_description``, so abstracts can only be
      carried forward from a previous catalog snapshot. Including
      ``abstract`` in the trigger keeps entries indexed even when
      they have no hand-enrichment but do have a useful abstract from
      a prior refresh.
    """
    enriched: dict[str, dict] = {}
    if not catalog_path.exists():
        return enriched
    with open(catalog_path) as f:
        for row in csv.DictReader(f):
            bp = (row.get("bioproject_id") or "").strip()
            if not bp:
                continue
            if any(
                row.get(c)
                for c in (
                    "protein_target",
                    "paper_doi",
                    "paper_pmid",
                    "n_rounds_declared",
                    "abstract",
                )
            ):
                enriched[bp] = row
    return enriched


def _passthrough_non_insdc(catalog_path: Path, exclude_ids: set[str]) -> list[dict]:
    """Carry over non-ENA-INSDC entries (Zenodo / Figshare / utexas processed
    deposits) from a previous catalog so they survive the refresh.
    """
    out: list[dict] = []
    if not catalog_path.exists():
        return out
    with open(catalog_path) as f:
        for row in csv.DictReader(f):
            bp = (row.get("bioproject_id") or "").strip()
            if not bp or bp in exclude_ids:
                continue
            if bp.startswith(("zenodo:", "figshare:", "utexas:")):
                out.append({k: row.get(k, "") for k in PUBLIC_COLS})
    return out


def _classify_all_studies(
    studies: dict[str, dict],
) -> dict[str, StudyStrategyClassification]:
    """Apply the per-run + per-study library_strategy classifier to every study.

    after the library_strategy classifier runs, accessions
    listed in :data:`MANUAL_EXCLUSIONS` are forced into the excluded
    bucket regardless of the classifier's verdict — these are studies
    that ENA labels as ``library_strategy="SELEX"`` but per-study
    evidence (library_source, library_name, study_title, run-level
    metadata) shows are mis-labels of non-SELEX experiments. The
    classifier alone can't catch them because ``SELEX`` is a compatible
    strategy by definition.
    """
    classifications: dict[str, StudyStrategyClassification] = {}
    for acc, meta in studies.items():
        if acc in MANUAL_EXCLUSIONS:
            n_runs = len(meta.get("library_strategies", []))
            classifications[acc] = StudyStrategyClassification(
                bioproject_id=acc,
                n_runs_total=n_runs,
                # Manual exclusion is at the study-context level, not
                # the library_strategy level — leave the strategy
                # bookkeeping at 0 so a reader can tell at a glance
                # this isn't a strategy-blocklist exclusion.
                n_runs_compatible=n_runs,
                n_runs_blocklisted=0,
                blocklisted_strategies={},
                should_exclude=True,
                is_mixed_strategy=False,
                exclusion_reason=f"manual exclusion: {MANUAL_EXCLUSIONS[acc]}",
            )
            continue
        classifications[acc] = classify_study_by_library_strategies(
            bioproject_id=acc,
            library_strategies=meta.get("library_strategies", []),
        )
    return classifications


def _excluded_sidecar_path(catalog_path: Path) -> Path:
    """Where ``bioprojects_excluded.csv`` lands relative to the main catalog."""
    return catalog_path.parent / "bioprojects_excluded.csv"


def _write_excluded_sidecar(
    path: Path,
    studies: dict[str, dict],
    classifications: dict[str, StudyStrategyClassification],
) -> int:
    """Emit ``bioprojects_excluded.csv``. Returns the number of rows written.

    Sorted by ``bioproject_id`` for deterministic diffs. Writes even when
    there are zero exclusions (header-only file) so the schema is
    discoverable.
    """
    rows: list[dict] = []
    for acc, classification in classifications.items():
        if not classification.should_exclude:
            continue
        meta = studies.get(acc, {})
        rows.append(
            {
                "bioproject_id": acc,
                "source": "ena",
                "study_title": meta.get("study_title", ""),
                "n_runs_total": classification.n_runs_total,
                "n_runs_blocklisted": classification.n_runs_blocklisted,
                "blocklisted_strategies": json.dumps(
                    classification.blocklisted_strategies,
                    sort_keys=True,
                ),
                "exclusion_reason": classification.exclusion_reason,
            }
        )
    rows.sort(key=lambda r: str(r["bioproject_id"]))

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(EXCLUDED_COLS))
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def rebuild_catalog(
    out_path: Path,
    preserve_from: Path | None = None,
    queries: tuple[str, ...] = ENA_QUERIES,
) -> int:
    """Refresh the catalog CSV in-place. Returns the count of bioprojects written.

    Behaviour:

    - Run every query in `queries` against ENA at the run level
      (``result=read_run``) so per-run ``library_strategy`` is available.
    - Group runs by study and apply the per-run + per-BioProject
      classifier
      (:func:`selexprep.fetch.library_strategy.classify_study_by_library_strategies`):
      drop studies whose runs are 100% blocklisted; keep all others.
    - Emit ``bioprojects_excluded.csv`` (sidecar next to ``out_path``)
      with exclusion reasons.
    - Map each kept ENA study row to the catalog schema.
    - When `preserve_from` is given AND the bioproject was hand-enriched in
      the old catalog, copy ``protein_target``, ``paper_doi``, ``paper_pmid``,
      ``n_rounds_declared`` forward into the new row.
    - Append non-INSDC entries (Zenodo / Figshare / utexas processed-data
      deposits) from the old catalog so refreshes don't lose them.

    No curation flags are added or preserved (``include`` /
    ``manual_curation_notes`` / etc.) — curation is the user's downstream
    job, not the package's.
    """
    studies = harvest_runs_from_ena(queries=queries)
    classifications = _classify_all_studies(studies)
    enriched = _enrichment_index(preserve_from) if preserve_from else {}

    rows: list[dict] = []
    n_kept_mixed = 0

    # Apply per-study library_strategy decision: keep all but the
    # 100%-blocklisted studies.
    kept_ids: set[str] = set()
    for acc, meta in studies.items():
        classification = classifications[acc]
        if classification.should_exclude:
            continue
        kept_ids.add(acc)
        if classification.is_mixed_strategy:
            n_kept_mixed += 1

        enr = enriched.get(acc, {})
        # abstract fallback: ENA's read_run query doesn't
        # expose study_description, so newly-discovered entries get a
        # blank abstract from ``meta``. Previously-known entries fall
        # back to the abstract preserved in the old catalog via
        # _enrichment_index so refresh doesn't strip them.
        abstract = (meta.get("study_description") or enr.get("abstract") or "").strip()
        rows.append(
            {
                "bioproject_id": acc,
                "source": "ena",
                "study_title": (meta.get("study_title") or "").strip(),
                "protein_target": (enr.get("protein_target") or "").strip(),
                "target_organism": (meta.get("scientific_name") or "").strip(),
                "paper_doi": (enr.get("paper_doi") or "").strip(),
                "paper_pmid": (enr.get("paper_pmid") or "").strip(),
                "n_rounds_declared": (enr.get("n_rounds_declared") or "").strip(),
                "abstract": abstract,
            }
        )

    # Carry over non-INSDC deposits from the previous catalog
    if preserve_from:
        rows.extend(_passthrough_non_insdc(preserve_from, exclude_ids=kept_ids))

    # Stable ordering: INSDC studies first (alphabetical), then non-INSDC
    rows.sort(
        key=lambda r: (
            0 if r["source"] == "ena" else 1,
            r["bioproject_id"],
        )
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(PUBLIC_COLS))
        writer.writeheader()
        writer.writerows(rows)

    # sidecar: ``bioprojects_excluded.csv`` lives next to
    # the main catalog file. Emitted even when empty so the schema is
    # discoverable and downstream consumers can rely on its presence.
    excluded_path = _excluded_sidecar_path(out_path)
    n_excluded = _write_excluded_sidecar(excluded_path, studies, classifications)

    logger.info(
        "Rebuilt catalog: %d rows → %s (kept %d, %d mixed-strategy; excluded %d → %s)",
        len(rows),
        out_path,
        len(kept_ids),
        n_kept_mixed,
        n_excluded,
        excluded_path,
    )
    return len(rows)


__all__ = [
    "ENA_PORTAL",
    "ENA_QUERIES",
    "EXCLUDED_COLS",
    "MANUAL_EXCLUSIONS",
    "PUBLIC_COLS",
    "harvest_runs_from_ena",
    "rebuild_catalog",
]
