"""Build a fetch plan for an accession: ENA metadata + round inference.

Single source of truth for "what does this accession look like and how do
its runs map to rounds?" — consumed by both ``selexprep fetch`` (single-
accession CLI) and ``selexprep run`` (batch driver).

Calls ENA filereport with **extended fields** (adds ``sample_title``,
``library_name``, ``experiment_title``, ``sample_accession``) so the
5-level cascade in :mod:`selexprep.fetch.metadata` can run on real
metadata, not just the bare-bones ``inspect`` field set.

Public API:

- :class:`FetchRun` — one run inside the plan (SRR + round assignment +
  FASTQ URLs/MD5s/bytes + paired-end flag).
- :class:`FetchPlan` — accession + study-level metadata + list of
  :class:`FetchRun`.
- :func:`build_fetch_plan` — fetch + parse.
- :func:`write_fetch_metadata_json` — deterministic JSON sidecar.

The plan is the audit trail. Resume oracles (downstream) verify FASTQ
files on disk against the plan; do NOT treat the plan's presence as a
"download complete" sentinel.
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path

from selexprep.fetch.inspect import query_ena_filereport, split_semicolon_list
from selexprep.fetch.metadata import RoundRecord, parse_round

logger = logging.getLogger(__name__)


_ENA_FETCH_FIELDS = (
    "run_accession,study_accession,study_title,library_strategy,"
    "library_source,library_name,experiment_title,sample_title,"
    "sample_accession,read_count,base_count,fastq_md5,fastq_bytes,fastq_ftp"
)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FetchRun:
    """One sequencing run inside a fetch plan.

    Adds ``library_strategy`` per-run. Previously only the
    study-level value was preserved on :class:`FetchPlan`; per-run
    granularity is required for the audit eligibility layer
    (:mod:`selexprep.benchmark.eligibility`) to detect mixed BioProjects
    where some runs are SELEX-compatible and others are blocklisted
    (RNA-Seq controls, ChIP-Seq adjacent assays, etc.). Defaults to
    empty string for backward compatibility with old fetch_metadata.json
    files written before the field existed.
    """

    srr: str
    sample_accession: str
    sample_title: str
    library_name: str
    experiment_title: str
    read_count: int
    base_count: int
    fastq_urls: list[str]
    fastq_md5s: list[str]
    fastq_bytes: list[int]
    paired_end: bool
    round_record: RoundRecord
    library_strategy: str = ""


@dataclass(frozen=True)
class FetchPlan:
    """Complete fetch plan for one accession.

    ``runs`` is sorted by ``srr`` for determinism.
    """

    accession: str
    bioproject_id: str | None
    study_title: str
    library_strategy: str
    library_source: str
    runs: list[FetchRun] = field(default_factory=list)

    @property
    def has_any_assigned_rounds(self) -> bool:
        """True iff at least one run has a HIGH/MEDIUM round assignment."""
        return any(
            run.round_record.round_number is not None
            and run.round_record.confidence in ("HIGH", "MEDIUM")
            for run in self.runs
        )

    @property
    def none_confidence_runs(self) -> list[FetchRun]:
        """Runs the cascade could not safely assign a round to.

        Delegates to ``RoundRecord.is_unassigned`` so the criterion is named once: a run is in this list
        iff its parsed `round_number` is None OR multiple distinct round
        numbers were parsed from the same metadata text (genuine
        ambiguity). MEDIUM records with a single unambiguous parse
        — e.g. ``library_name=RAPT26-2R`` → round 2 — are NOT covered
        here; they're trusted as L3 assignments per the cascade
        docstring.

        Both cases legitimately block fetch (the trusted-assignments
        contract that ``detect`` / ``extract`` consume requires a
        single unambiguous round per run).
        """
        return [run for run in self.runs if run.round_record.is_unassigned]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def build_fetch_plan(accession: str, *, timeout_s: int = 30) -> FetchPlan:
    """Hit ENA filereport with extended fields + parse rounds per run.

    Single ENA round-trip; no per-run sample-attribute calls (deferred to
    a follow-up if benchmarking shows the filereport fields are
    insufficient).

    Raises:
        requests.HTTPError on non-2xx.
        ValueError if ENA returns no records.
    """
    rows = query_ena_filereport(accession, fields=_ENA_FETCH_FIELDS, timeout_s=timeout_s)

    runs: list[FetchRun] = []
    for row in rows:
        srr = str(row.get("run_accession", "")).strip()
        if not srr:
            logger.warning("build_fetch_plan: skipping row with empty run_accession")
            continue

        urls = split_semicolon_list(row.get("fastq_ftp"))
        md5s = split_semicolon_list(row.get("fastq_md5"))
        size_strs = split_semicolon_list(row.get("fastq_bytes"))
        sizes: list[int] = []
        for s in size_strs:
            try:
                sizes.append(int(s))
            except ValueError:
                logger.warning("build_fetch_plan: bad fastq_bytes %r in %s", s, srr)

        sample_title = str(row.get("sample_title") or "")
        library_name = str(row.get("library_name") or "")
        experiment_title = str(row.get("experiment_title") or "")

        round_record = parse_round(
            srr=srr,
            sample_title=sample_title,
            library_name=library_name,
            experiment_title=experiment_title,
        )

        runs.append(
            FetchRun(
                srr=srr,
                sample_accession=str(row.get("sample_accession") or ""),
                sample_title=sample_title,
                library_name=library_name,
                experiment_title=experiment_title,
                read_count=int(row.get("read_count") or 0),
                base_count=int(row.get("base_count") or 0),
                fastq_urls=urls,
                fastq_md5s=md5s,
                fastq_bytes=sizes,
                paired_end=len(urls) == 2,
                round_record=round_record,
                library_strategy=str(row.get("library_strategy") or "").strip(),
            )
        )

    runs.sort(key=lambda r: r.srr)

    first = rows[0]
    return FetchPlan(
        accession=accession,
        bioproject_id=(str(first.get("study_accession")) if first.get("study_accession") else None),
        study_title=str(first.get("study_title") or ""),
        library_strategy=str(first.get("library_strategy") or ""),
        library_source=str(first.get("library_source") or ""),
        runs=runs,
    )


def fastq_filenames_for_run(run: FetchRun) -> list[str]:
    """Canonical FASTQ basenames the download will produce for this run.

    derive directly from ``run.fastq_urls`` (the same
    URLs ``download_srr_ena_direct`` uses at ``fetch/download.py:351``
    where ``dest = output_dir / Path(url_path).name``). Synthesizing
    ``{SRR}.fastq.gz`` is right for ENA's canonical naming but brittle
    if URLs ever use a different scheme — deriving from the URLs is the
    source of truth.

    Fallback: when ``fastq_urls`` is empty (defensive only — a real ENA
    response always has them when the run is downloadable), synthesize
    from ``paired_end`` so the resume oracle still returns a defined
    answer instead of an empty list.
    """
    if run.fastq_urls:
        return [Path(url).name for url in run.fastq_urls]
    if run.paired_end:
        return [f"{run.srr}_1.fastq.gz", f"{run.srr}_2.fastq.gz"]
    return [f"{run.srr}.fastq.gz"]


def write_fetch_metadata_json(plan: FetchPlan, path: Path) -> None:
    """Emit a deterministic JSON serialization of a FetchPlan.

    Sorted top-level keys; runs already sorted by SRR (build_fetch_plan
    invariant). Used as the audit trail next to ``rounds.tsv``.
    """
    payload = {
        "accession": plan.accession,
        "bioproject_id": plan.bioproject_id,
        "study_title": plan.study_title,
        "library_strategy": plan.library_strategy,
        "library_source": plan.library_source,
        "runs": [_run_to_dict(r) for r in plan.runs],
    }
    text = json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    path.write_text(text, encoding="utf-8")


def _run_to_dict(run: FetchRun) -> dict:
    d = asdict(run)
    rr = d["round_record"]
    rr["round_candidates"] = list(run.round_record.round_candidates)
    return d
