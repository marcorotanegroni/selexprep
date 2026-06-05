"""Fetch orchestrator: accession → per-round FASTQs + round map + audit trail.

Wraps :func:`selexprep.fetch.plan.build_fetch_plan` (metadata) with
:func:`selexprep.fetch.download.download_srr` (transfer) and emits the
artefacts ``selexprep detect`` / ``selexprep extract`` consume:

::

    outdir/
    ├── round_NN/SRRxxxxx[_1|_2].fastq.gz   (HIGH/MEDIUM-confidence runs)
    ├── round_unknown/SRRzzzzz.fastq.gz     (only with --allow-manual-review)
    ├── rounds.tsv                          (HIGH/MEDIUM only; trusted contract)
    ├── manual_review.tsv                   (only if --allow-manual-review used)
    └── fetch_metadata.json                 (audit trail; not a resume oracle)

**Cardinal rule:** never guess a round
assignment. NONE-confidence runs do NOT appear in ``rounds.tsv`` —
downstream tools see only trusted assignments. The user opts in to
``allow_manual_review`` to download them into ``round_unknown/`` AND
surface them in ``manual_review.tsv``.

**contract relaxation.** When some (but not all) runs are
NONE-confidence and ``allow_manual_review`` is False, the runner
*skips* those runs (logging the SRRs at WARNING level) and proceeds
with the HIGH/MEDIUM ones. Previously the runner hard-refused the
whole accession, which threw away clean SELEX time series when even a
single run lacked a parseable round (audit cohort: PRJNA1244796 with
260/268 parseable, PRJNA809588 with 10/12). The trusted-assignments
contract is preserved because ``rounds.tsv`` is still HIGH/MEDIUM-only;
the relaxation is purely about *which* runs we attempt to download.

If every run is NONE-confidence (no HIGH/MEDIUM rounds), the runner
fails fast *unless* ``allow_manual_review`` is passed. Without the
flag, refusing avoids wasting bandwidth on a download whose rounds
can't be trusted. With it, the user has explicitly opted
in to manual curation, so all runs download into ``round_unknown/``,
``rounds.tsv`` is emitted empty (no trusted assignments), and a
hand-supplied round-map drives ``detect`` downstream. This is what the
Tier-1 benchmark's ``curated`` rows rely on for deposits whose round
structure is absent from ENA metadata.

Public API:

- :class:`FetchResult` — summary of one :func:`run_fetch` invocation.
- :func:`run_fetch` — single-accession orchestrator.
- :func:`check_fetch_inventory` — resume oracle used by both ``fetch``
  (skip already-present files) and ``selexprep run`` (re-run fetch if
  any expected FASTQ is missing).
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path

from selexprep.fetch.download import DownloadBackend, download_srr, validate_fastq_gz
from selexprep.fetch.plan import (
    FetchPlan,
    FetchRun,
    build_fetch_plan,
    fastq_filenames_for_run,
    write_fetch_metadata_json,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FetchResult:
    """Outcome of one :func:`run_fetch` invocation."""

    plan: FetchPlan
    outdir: Path
    rounds_tsv: Path | None
    manual_review_tsv: Path | None
    fetch_metadata_json: Path
    downloaded_srrs: list[str] = field(default_factory=list)
    skipped_srrs: list[str] = field(default_factory=list)
    failed_srrs: list[str] = field(default_factory=list)
    manual_review_srrs: list[str] = field(default_factory=list)
    refused_reason: str | None = None


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def run_fetch(
    accession: str,
    outdir: Path,
    *,
    backend: DownloadBackend = "ena",
    allow_manual_review: bool = False,
    dry_run: bool = False,
    timeout_s: int = 30,
) -> FetchResult:
    """Build a plan + download per-run FASTQs + emit round map.

     Default ``backend="ena"`` matches the CLI's paper-grade default
    . Library-
     level ``download_srr`` keeps its ``auto`` default for Python API
     callers who want the convenience chain.
    """
    from selexprep.catalog.filter import is_discovery_only

    if is_discovery_only(accession):
        raise ValueError(
            f"{accession!r} is a discovery-only catalog pointer (non-INSDC: "
            "figshare/zenodo/utexas) and cannot be fetched as raw FASTQ in v0.1. "
            "Only INSDC accessions (PRJNA/PRJDB/PRJEB/SRP/ERP/DRP) are retrievable."
        )

    outdir.mkdir(parents=True, exist_ok=True)

    plan = build_fetch_plan(accession, timeout_s=timeout_s)

    # Refusal: no run produced a HIGH/MEDIUM round assignment at all.
    # gated by ``not allow_manual_review`` (consistent with the
    # some-unassigned block below). Without the flag we refuse — a download
    # whose rounds can't be trusted wastes bandwidth. With it, the user has
    # explicitly opted in to manual curation, so we fall through: every run
    # downloads into round_unknown/, rounds.tsv comes out empty, and a
    # curated round-map drives detect downstream (the Tier-1 benchmark's
    # ``curated`` path for metadata-less deposits).
    if not plan.has_any_assigned_rounds and not allow_manual_review:
        reason = (
            f"All {len(plan.runs)} run(s) in {accession} are unassigned "
            "(no HIGH/MEDIUM round could be parsed from sample_attributes, "
            "sample_title, library_name, or experiment_title). detect's "
            "cross-round persistence is unusable without HIGH/MEDIUM "
            "rounds; refusing to download. Re-run with --allow-manual-review "
            "to download all runs into round_unknown/ and drive detect with "
            "a hand-curated round-map, or open the accession's metadata to "
            "confirm round indicators are missing."
        )
        logger.error("run_fetch[%s]: %s", accession, reason)
        metadata_path = outdir / "fetch_metadata.json"
        write_fetch_metadata_json(plan, metadata_path)
        return FetchResult(
            plan=plan,
            outdir=outdir,
            rounds_tsv=None,
            manual_review_tsv=None,
            fetch_metadata_json=metadata_path,
            refused_reason=reason,
        )

    # all-unassigned BUT --allow-manual-review is set → proceed.
    # The download loop routes every (unassigned) run to round_unknown/ and
    # emits manual_review.tsv; rounds.tsv will be empty (no trusted rounds).
    if not plan.has_any_assigned_rounds:
        logger.warning(
            "run_fetch[%s]: all %d run(s) unassigned; downloading to "
            "round_unknown/ under --allow-manual-review. rounds.tsv will be "
            "empty — supply a curated round-map to detect.",
            accession,
            len(plan.runs),
        )

    # partial-parseability is no longer a hard refusal.
    # When some runs are unassigned and --allow-manual-review is NOT set,
    # log them at WARNING level and skip them in the download loop; the
    # rest of the accession proceeds normally. The hard refusal at line
    # ~102 (all-unassigned) still fires when there's nothing to fetch.
    none_runs = plan.none_confidence_runs
    if none_runs and not allow_manual_review:
        srrs = ", ".join(r.srr for r in none_runs)
        logger.warning(
            "run_fetch[%s]: skipping %d unassigned run(s) (no parseable "
            "round, or conflicting round numbers): %s. Pass "
            "--allow-manual-review to download them into round_unknown/ "
            "and surface them in manual_review.tsv. They will never "
            "enter rounds.tsv (trusted-assignments contract).",
            accession,
            len(none_runs),
            srrs,
        )

    # Download per run.
    downloaded: list[str] = []
    skipped: list[str] = []
    failed: list[str] = []
    manual_review_srrs: list[str] = []

    for run in plan.runs:
        if run.round_record.is_unassigned and not allow_manual_review:
            # skip unassigned runs unless the user opted in
            # to manual-review downloading. The WARNING above already
            # listed the SRRs and the reason; recording in ``skipped``
            # keeps the result shape stable.
            skipped.append(run.srr)
            continue
        target_dir = _target_dir_for_run(outdir, run)
        target_dir.mkdir(parents=True, exist_ok=True)
        if run.round_record.is_unassigned:
            manual_review_srrs.append(run.srr)

        # skip-already-present must validate the gzip stream,
        # not just check Path.exists. A SIGKILL'd download leaves a corrupt
        # .fastq.gz on disk; bare exists() would skip it and downstream
        # detect/extract would hit gzip errors deep in the pipeline.
        # validate_fastq_gz mirrors the same check download_srr_ena_direct
        # performs at fetch/download.py:288.
        expected_files = [target_dir / name for name in fastq_filenames_for_run(run)]
        if expected_files and all(validate_fastq_gz(p) for p in expected_files):
            logger.info("run_fetch[%s/%s]: FASTQs already present — skipping", accession, run.srr)
            skipped.append(run.srr)
            continue

        if dry_run:
            logger.info(
                "run_fetch[%s/%s]: dry-run; would download to %s", accession, run.srr, target_dir
            )
            skipped.append(run.srr)
            continue

        ok = download_srr(run.srr, output_dir=target_dir, backend=backend)
        if ok:
            downloaded.append(run.srr)
        else:
            failed.append(run.srr)

    # Emit rounds.tsv (HIGH/MEDIUM only).
    rounds_tsv = outdir / "rounds.tsv"
    _write_rounds_tsv(rounds_tsv, plan, outdir)

    # Emit manual_review.tsv only if any manual-review runs were processed.
    manual_review_tsv: Path | None = None
    if manual_review_srrs:
        manual_review_tsv = outdir / "manual_review.tsv"
        _write_manual_review_tsv(manual_review_tsv, plan, outdir)

    # Emit fetch_metadata.json (audit trail).
    metadata_path = outdir / "fetch_metadata.json"
    write_fetch_metadata_json(plan, metadata_path)

    return FetchResult(
        plan=plan,
        outdir=outdir,
        rounds_tsv=rounds_tsv,
        manual_review_tsv=manual_review_tsv,
        fetch_metadata_json=metadata_path,
        downloaded_srrs=downloaded,
        skipped_srrs=skipped,
        failed_srrs=failed,
        manual_review_srrs=manual_review_srrs,
    )


def check_fetch_inventory(plan: FetchPlan, outdir: Path) -> list[str]:
    """Return basenames of expected FASTQs from `plan` that are missing or corrupt.

    Used by the run-batch driver as the file-inventory resume oracle for
    the fetch stage. Returns an
    empty list iff every expected FASTQ is on disk in its expected
    ``round_NN/`` / ``round_unknown/`` directory AND passes
    :func:`validate_fastq_gz` (gzip-stream integrity, min-size, ``gzip -t``).

    previous version used bare ``Path.exists()`` which
    would skip a SIGKILL'd corrupt partial download.
    """
    missing: list[str] = []
    for run in plan.runs:
        target_dir = _target_dir_for_run(outdir, run)
        for name in fastq_filenames_for_run(run):
            path = target_dir / name
            if not validate_fastq_gz(path):
                missing.append(name)
    return missing


def read_fetch_metadata_json(path: Path) -> FetchPlan:
    """Reconstruct a FetchPlan from the on-disk audit trail.

    Used by the run-batch driver during resume to compare expected
    FASTQs (from the audit trail) against what's on disk.
    """
    from selexprep.fetch.metadata import RoundRecord  # local import: rare path

    payload = json.loads(path.read_text(encoding="utf-8"))
    runs: list[FetchRun] = []
    for r in payload.get("runs", []):
        rr_data = r["round_record"]
        # audit refactor: ``needs_manual_review`` was removed
        # from ``RoundRecord`` (replaced by the ``is_unassigned`` property
        # computed from ``round_number`` + ``round_candidates``). Old
        # ``fetch_metadata.json`` files on disk still have the field;
        # silently drop it. ``round_candidates`` carries the information
        # needed to reconstruct ``is_unassigned`` correctly.
        round_record = RoundRecord(
            srr=rr_data["srr"],
            round_number=rr_data["round_number"],
            confidence=rr_data["confidence"],
            source_field=rr_data["source_field"],
            matched_pattern=rr_data["matched_pattern"],
            round_candidates=rr_data.get("round_candidates", []) or [],
            parser_notes=rr_data.get("parser_notes", ""),
            target_hint=rr_data.get("target_hint"),
        )
        runs.append(
            FetchRun(
                srr=r["srr"],
                sample_accession=r.get("sample_accession", ""),
                sample_title=r.get("sample_title", ""),
                library_name=r.get("library_name", ""),
                experiment_title=r.get("experiment_title", ""),
                read_count=int(r.get("read_count", 0)),
                base_count=int(r.get("base_count", 0)),
                fastq_urls=list(r.get("fastq_urls", [])),
                fastq_md5s=list(r.get("fastq_md5s", [])),
                fastq_bytes=list(r.get("fastq_bytes", [])),
                paired_end=bool(r.get("paired_end", False)),
                round_record=round_record,
                # per-run library_strategy. Old
                # fetch_metadata.json files written before the field
                # existed default to empty string.
                library_strategy=str(r.get("library_strategy", "")),
            )
        )
    return FetchPlan(
        accession=payload["accession"],
        bioproject_id=payload.get("bioproject_id"),
        study_title=payload.get("study_title", ""),
        library_strategy=payload.get("library_strategy", ""),
        library_source=payload.get("library_source", ""),
        runs=runs,
    )


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _target_dir_for_run(outdir: Path, run: FetchRun) -> Path:
    """Where this run's FASTQs land — `round_NN/` or `round_unknown/`."""
    rr = run.round_record
    if rr.is_unassigned:
        return outdir / "round_unknown"
    # is_unassigned False implies round_number is not None — narrow for mypy.
    assert rr.round_number is not None
    return outdir / f"round_{rr.round_number:02d}"


def _write_rounds_tsv(path: Path, plan: FetchPlan, outdir: Path) -> None:
    """Emit `file<TAB>round_number` for HIGH/MEDIUM-confidence runs only.

    Matches ``cli._load_round_map`` format verbatim. Sorted
    by filename for determinism.
    """
    rows: list[tuple[str, int]] = []
    for run in plan.runs:
        rr = run.round_record
        if rr.is_unassigned:
            continue
        if rr.confidence not in ("HIGH", "MEDIUM"):
            continue
        assert rr.round_number is not None  # narrow for mypy (is_unassigned False)
        for name in fastq_filenames_for_run(run):
            rows.append((name, rr.round_number))
    rows.sort(key=lambda x: x[0])

    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("file\tround_number\n")
        for name, rn in rows:
            fh.write(f"{name}\t{rn}\n")


def _write_manual_review_tsv(path: Path, plan: FetchPlan, outdir: Path) -> None:
    """Emit per-FASTQ manual-review records.

    Header: ``file<TAB>srr<TAB>parser_notes<TAB>sample_title<TAB>library_name``.
    Sorted by SRR for determinism.
    """
    rows: list[tuple[str, str, str, str, str]] = []
    for run in plan.runs:
        rr = run.round_record
        if not rr.is_unassigned:
            continue
        for name in fastq_filenames_for_run(run):
            rows.append(
                (
                    name,
                    run.srr,
                    rr.parser_notes,
                    run.sample_title,
                    run.library_name,
                )
            )
    rows.sort(key=lambda x: (x[1], x[0]))

    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("file\tsrr\tparser_notes\tsample_title\tlibrary_name\n")
        for name, srr, notes, sample_title, library_name in rows:
            fh.write(f"{name}\t{srr}\t{notes}\t{sample_title}\t{library_name}\n")
