"""Batch driver: fetch → detect → extract → count → qc per accession.

Closes the v0.1 single-dataset CLI surface so ``selexprep run
accessions.tsv`` can process multiple accessions in one invocation.
Each accession runs through the full pipeline as a self-contained
sub-directory under ``outdir/<accession>/`` so a SIGKILL mid-batch
leaves a partially-completed state ``--resume`` can pick up from.

**Stage-by-stage resume oracles** (locked plan + user peer-review
point 5; file-inventory based, not sentinel-flag based):

| Stage   | Oracle                                                              |
|---------|---------------------------------------------------------------------|
| fetch   | ``fetch_metadata.json`` exists AND every expected FASTQ is on disk  |
| detect  | ``library_report.json`` exists                                      |
| extract | ``selexprep_manifest.json`` exists                                  |
| count   | per-round ``round_NN/counts.parquet`` (granular: only missing rounds re-run) |
| qc      | ``qc/flags.yaml`` exists                                            |

**Paired-end handling** (user peer-review point 1): FASTQs are grouped
by ENA naming convention — ``<srr>_1.fastq.gz`` = R1,
``<srr>_2.fastq.gz`` = R2. R1-only sequences feed
``compute_library_report``'s primary stream; R2 sequences are passed
as ``paired_mate_streams``; both are threaded through ``run_extract``
as ``paired_r2_inputs``. No commingled detect inputs.

**Split-primer skip** (user peer-review point 2): if the inferred
``LibraryReport.required_action == "READ_MERGING_RECOMMENDED"`` (i.e.,
``PAIRED_END_SPLIT_PRIMERS``, ``full_insert_recovered=False``), count
and qc are skipped and status records
``SKIPPED_READ_MERGING_RECOMMENDED``. Joining R1+R2 by read ID alone
would produce misleading half-insert counts (locked plan line 325).

**Manual-review separation** (user peer-review point 3): inherited
from ``selexprep.fetch.runner.run_fetch`` — NONE-confidence runs land
in ``round_unknown/`` and ``manual_review.tsv``; ``rounds.tsv``
contains HIGH/MEDIUM only.
"""

from __future__ import annotations

import gzip
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import pandas as pd
import requests
import yaml

from selexprep.count.counter import count_fasta
from selexprep.extract import run_extract
from selexprep.fetch import (
    FetchPlan,
    check_fetch_inventory,
    read_fetch_metadata_json,
    run_fetch,
)
from selexprep.fetch.download import DownloadBackend
from selexprep.library import (
    compute_library_report,
    read_library_report_json,
    write_library_report_json,
)
from selexprep.library.report import LibraryReport, ReadSource
from selexprep.qc.runner import run_qc

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------


RunStatus = Literal[
    "OK",
    "SKIPPED_READ_MERGING_RECOMMENDED",
    "FETCH_REFUSED",
    "FETCH_FAILED",
    "DETECT_FAILED",
    "EXTRACT_REFUSED",
    "EXTRACT_FAILED",
    "COUNT_FAILED",
    "QC_FAILED",
    # Codex pass 1 fix: a truly-unexpected exception that escaped the
    # per-stage try/except blocks should not be misreported as
    # EXTRACT_FAILED. Used only by the outer except in run_batch.
    "UNEXPECTED_FAILURE",
]


@dataclass
class RunRowReport:
    """Status of one accession in a batch."""

    accession: str
    status: RunStatus
    last_stage_completed: str
    extraction_mode: str | None = None
    required_action: str | None = None
    confidence: float | None = None
    library_report_status: str | None = None
    flags_raised: int | None = None
    notes: str = ""


@dataclass
class RunReport:
    """Outcome of one :func:`run_batch` invocation."""

    accessions_tsv: Path
    outdir: Path
    rows: list[RunRowReport] = field(default_factory=list)
    summary_tsv: Path | None = None


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def run_batch(
    accessions_tsv: Path,
    outdir: Path,
    *,
    resume: bool = False,
    stop_on_error: bool = False,
    backend: DownloadBackend = "ena",
    allow_manual_review: bool = False,
    timeout_s: int = 30,
) -> RunReport:
    """Run the full pipeline for every accession in ``accessions_tsv``.

    Args:
        accessions_tsv: TSV with at least an ``accession`` column.
        outdir: parent directory; one sub-directory created per accession.
        resume: skip stages whose resume oracle passes (per-stage table
            above). Default False = re-run everything.
        stop_on_error: fail-fast on the first per-accession failure
            instead of recording status and continuing.
        backend: download backend passed to ``run_fetch`` (default
            ``"ena"`` for paper-grade reproducibility).
        allow_manual_review: thread through to ``run_fetch`` —
            NONE-confidence runs go to ``round_unknown/`` and
            ``manual_review.tsv``, not ``rounds.tsv``.
        timeout_s: HTTP timeout for the ENA metadata call.

    Returns:
        :class:`RunReport` with per-accession :class:`RunRowReport`
        rows and the path to the emitted ``run_summary.tsv``.
    """
    accessions = _load_accessions(accessions_tsv)
    outdir.mkdir(parents=True, exist_ok=True)
    rows: list[RunRowReport] = []

    for accession in accessions:
        logger.info("run_batch[%s]: starting", accession)
        acc_dir = outdir / accession
        acc_dir.mkdir(parents=True, exist_ok=True)
        try:
            row = _process_one_accession(
                accession=accession,
                acc_dir=acc_dir,
                resume=resume,
                backend=backend,
                allow_manual_review=allow_manual_review,
                timeout_s=timeout_s,
            )
        except Exception as e:
            # Codex pass 1 fix: this is the outer safety net for exceptions
            # that escaped EVERY per-stage try/except (i.e., a truly-
            # unexpected programming error or OS condition). It must NOT
            # claim EXTRACT_FAILED — that misreports an unrelated failure
            # (e.g., an ENA HTTPError during fetch) as an extract bug.
            # Stage-specific exceptions get caught + classified inside
            # _process_one_accession.
            logger.exception("run_batch[%s]: unexpected error", accession)
            row = RunRowReport(
                accession=accession,
                status="UNEXPECTED_FAILURE",
                last_stage_completed="unknown",
                notes=f"unexpected: {type(e).__name__}: {e}",
            )

        rows.append(row)
        if stop_on_error and row.status not in ("OK", "SKIPPED_READ_MERGING_RECOMMENDED"):
            logger.error("run_batch: --stop-on-error tripped on %s (%s)", accession, row.status)
            break

    summary_tsv = outdir / "run_summary.tsv"
    _write_summary_tsv(summary_tsv, rows)
    return RunReport(
        accessions_tsv=accessions_tsv,
        outdir=outdir,
        rows=rows,
        summary_tsv=summary_tsv,
    )


# ---------------------------------------------------------------------------
# TSV loader
# ---------------------------------------------------------------------------


def _load_accessions(tsv_path: Path) -> list[str]:
    df = pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False)
    if "accession" not in df.columns:
        raise ValueError(
            f"accessions TSV {tsv_path} must have an 'accession' column; found {list(df.columns)}"
        )
    accs = [a.strip() for a in df["accession"].tolist() if a.strip()]
    dups = sorted({a for a in accs if accs.count(a) > 1})
    if dups:
        raise ValueError(f"accessions TSV has duplicate accessions: {dups}")
    return accs


# ---------------------------------------------------------------------------
# Per-accession pipeline
# ---------------------------------------------------------------------------


def _process_one_accession(
    *,
    accession: str,
    acc_dir: Path,
    resume: bool,
    backend: DownloadBackend,
    allow_manual_review: bool,
    timeout_s: int,
) -> RunRowReport:
    """Run fetch → detect → extract → count → qc for one accession."""

    # ----------------------------------------------------------- fetch stage
    plan: FetchPlan | None = None
    fetch_metadata_path = acc_dir / "fetch_metadata.json"
    rounds_tsv_path = acc_dir / "rounds.tsv"

    if resume and _fetch_inventory_complete(acc_dir):
        logger.info("run_batch[%s]: fetch resume oracle passed", accession)
        try:
            plan = read_fetch_metadata_json(fetch_metadata_path)
        except (ValueError, KeyError, OSError) as e:
            logger.exception("run_batch[%s]: failed to read resumed fetch_metadata.json", accession)
            return RunRowReport(
                accession=accession,
                status="FETCH_FAILED",
                last_stage_completed="fetch",
                notes=f"resume read failed: {type(e).__name__}: {e}",
            )
    else:
        # Codex pass 1 fix: HTTPError / RequestException / ValueError from
        # build_fetch_plan or download_srr must be classified as FETCH_FAILED,
        # not bubble up to the outer except (which would mislabel them).
        try:
            fetch_result = run_fetch(
                accession,
                acc_dir,
                backend=backend,
                allow_manual_review=allow_manual_review,
                timeout_s=timeout_s,
            )
        except (requests.HTTPError, requests.RequestException, ValueError) as e:
            logger.exception("run_batch[%s]: fetch stage raised", accession)
            return RunRowReport(
                accession=accession,
                status="FETCH_FAILED",
                last_stage_completed="fetch",
                notes=f"{type(e).__name__}: {e}",
            )
        if fetch_result.refused_reason is not None:
            return RunRowReport(
                accession=accession,
                status="FETCH_REFUSED",
                last_stage_completed="fetch",
                notes=fetch_result.refused_reason,
            )
        if fetch_result.failed_srrs:
            return RunRowReport(
                accession=accession,
                status="FETCH_FAILED",
                last_stage_completed="fetch",
                notes=f"failed SRRs: {', '.join(fetch_result.failed_srrs)}",
            )
        plan = fetch_result.plan

    # ----------------------------------------------------------- pair R1/R2
    r1_by_round, r2_by_round = _collect_round_inputs(acc_dir)
    is_paired = r2_by_round is not None and any(r2_by_round.values())
    if not r1_by_round:
        return RunRowReport(
            accession=accession,
            status="FETCH_FAILED",
            last_stage_completed="fetch",
            notes="No R1 FASTQs discovered after fetch (round_*/)",
        )

    # ---------------------------------------------------------- detect stage
    library_report_path = acc_dir / "library_report.json"
    if resume and library_report_path.exists():
        logger.info("run_batch[%s]: detect resume oracle passed", accession)
        lr = read_library_report_json(library_report_path)
    else:
        try:
            sequences_by_round = _read_sequences_by_round(r1_by_round)
            paired_mate_streams: dict[int, list[str]] | None = None
            read_source: ReadSource = "R1"
            if is_paired and r2_by_round is not None:
                paired_mate_streams = _read_sequences_by_round(r2_by_round)
                read_source = "R1_AND_R2"
            lr = compute_library_report(
                sequences_by_round,
                read_source=read_source,
                paired_mate_streams=paired_mate_streams,
            )
            write_library_report_json(lr, library_report_path)
        except Exception as e:
            logger.exception("run_batch[%s]: detect failed", accession)
            return RunRowReport(
                accession=accession,
                status="DETECT_FAILED",
                last_stage_completed="fetch",
                notes=f"{type(e).__name__}: {e}",
            )

    # ---------------------------------------------------------- extract stage
    manifest_path = acc_dir / "selexprep_manifest.json"
    if resume and manifest_path.exists():
        logger.info("run_batch[%s]: extract resume oracle passed", accession)
    else:
        round_map = (
            _build_round_map(rounds_tsv_path)
            if rounds_tsv_path.exists()
            else _build_round_map_from_inputs(r1_by_round, r2_by_round)
        )
        all_r1: list[Path] = sorted(p for paths in r1_by_round.values() for p in paths)
        try:
            extract_result = run_extract(
                lr,
                fastq_inputs=all_r1,
                outdir=acc_dir,
                round_map=round_map,
                paired_r2_inputs=r2_by_round if is_paired else None,
                accession=accession,
                bioproject_id=plan.bioproject_id,
                runs=[run.srr for run in plan.runs],
                parameters={"backend": backend, "resume": str(resume).lower()},
            )
        except Exception as e:
            logger.exception("run_batch[%s]: extract failed", accession)
            return RunRowReport(
                accession=accession,
                status="EXTRACT_FAILED",
                last_stage_completed="detect",
                notes=f"{type(e).__name__}: {e}",
                extraction_mode=lr.extraction_mode,
                required_action=lr.required_action,
                confidence=lr.confidence,
                library_report_status=lr.status,
            )
        if extract_result.skipped_reason is not None:
            return RunRowReport(
                accession=accession,
                status="EXTRACT_REFUSED",
                last_stage_completed="detect",
                notes=extract_result.skipped_reason,
                extraction_mode=lr.extraction_mode,
                required_action=lr.required_action,
                confidence=lr.confidence,
                library_report_status=lr.status,
            )

    # ----------------------------------- Split-primer guard (point 2)
    if lr.required_action == "READ_MERGING_RECOMMENDED":
        return RunRowReport(
            accession=accession,
            status="SKIPPED_READ_MERGING_RECOMMENDED",
            last_stage_completed="extract",
            notes=(
                "PAIRED_END_SPLIT_PRIMERS: count + qc skipped. R1 and R2 "
                "cover different parts of the insert; joining by read ID "
                "alone is biologically wrong without merging (v0.2)."
            ),
            extraction_mode=lr.extraction_mode,
            required_action=lr.required_action,
            confidence=lr.confidence,
            library_report_status=lr.status,
        )

    # ----------------------------------------------------------- count stage
    try:
        _run_count_stage(acc_dir, lr, resume=resume)
    except Exception as e:
        logger.exception("run_batch[%s]: count failed", accession)
        return RunRowReport(
            accession=accession,
            status="COUNT_FAILED",
            last_stage_completed="extract",
            notes=f"{type(e).__name__}: {e}",
            extraction_mode=lr.extraction_mode,
            required_action=lr.required_action,
            confidence=lr.confidence,
            library_report_status=lr.status,
        )

    # -------------------------------------------------------------- qc stage
    qc_flags_path = acc_dir / "qc" / "flags.yaml"
    if resume and qc_flags_path.exists():
        logger.info("run_batch[%s]: qc resume oracle passed", accession)
        flags_raised = _count_yaml_flags(qc_flags_path)
    else:
        try:
            qc_result = run_qc(manifest_path)
            flags_raised = qc_result.n_flags_raised
        except Exception as e:
            logger.exception("run_batch[%s]: qc failed", accession)
            return RunRowReport(
                accession=accession,
                status="QC_FAILED",
                last_stage_completed="count",
                notes=f"{type(e).__name__}: {e}",
                extraction_mode=lr.extraction_mode,
                required_action=lr.required_action,
                confidence=lr.confidence,
                library_report_status=lr.status,
            )

    return RunRowReport(
        accession=accession,
        status="OK",
        last_stage_completed="qc",
        extraction_mode=lr.extraction_mode,
        required_action=lr.required_action,
        confidence=lr.confidence,
        library_report_status=lr.status,
        flags_raised=flags_raised,
    )


# ---------------------------------------------------------------------------
# Stage helpers
# ---------------------------------------------------------------------------


def _fetch_inventory_complete(acc_dir: Path) -> bool:
    """Resume oracle for fetch (user peer-review point 5)."""
    metadata_path = acc_dir / "fetch_metadata.json"
    rounds_tsv = acc_dir / "rounds.tsv"
    if not metadata_path.exists():
        return False
    try:
        plan = read_fetch_metadata_json(metadata_path)
    except (ValueError, KeyError) as e:
        logger.warning("fetch resume oracle: cannot read %s — %s", metadata_path, e)
        return False
    missing = check_fetch_inventory(plan, acc_dir)
    if missing:
        logger.info("fetch resume oracle: missing FASTQs %s — will re-run fetch", missing)
        return False
    if plan.has_any_assigned_rounds and not rounds_tsv.exists():
        logger.info("fetch resume oracle: rounds.tsv missing — will re-run fetch")
        return False
    return True


_R1_R2_PATTERN = re.compile(r"^(?P<srr>[^.]+?)(?P<suffix>_[12])?\.fastq\.gz$")


def _collect_round_inputs(
    acc_dir: Path,
) -> tuple[dict[int, list[Path]], dict[int, list[Path]] | None]:
    """Group ``round_NN/*.fastq.gz`` into R1/R2 buckets by ``_1``/``_2`` suffix."""
    r1_by_round: dict[int, list[Path]] = {}
    r2_by_round: dict[int, list[Path]] = {}
    any_r2 = False

    for round_dir in sorted(acc_dir.glob("round_*")):
        if round_dir.name == "round_unknown":
            continue
        try:
            r = int(round_dir.name.split("_")[1])
        except (IndexError, ValueError):
            continue
        for fq in sorted(round_dir.glob("*.fastq.gz")):
            match = _R1_R2_PATTERN.match(fq.name)
            if match and match.group("suffix") == "_2":
                r2_by_round.setdefault(r, []).append(fq)
                any_r2 = True
            else:
                r1_by_round.setdefault(r, []).append(fq)

    return r1_by_round, (r2_by_round if any_r2 else None)


def _read_sequences_by_round(by_round: dict[int, list[Path]]) -> dict[int, list[str]]:
    """Read sequence lines from per-round FASTQ inputs."""
    out: dict[int, list[str]] = {}
    for r, paths in by_round.items():
        seqs: list[str] = []
        for fq in paths:
            seqs.extend(_read_fastq_sequences(fq))
        out[r] = seqs
    return out


def _read_fastq_sequences(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    seqs: list[str] = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seqs.append(line.strip())
    return seqs


def _build_round_map(rounds_tsv: Path) -> dict[str, int]:
    df = pd.read_csv(rounds_tsv, sep="\t")
    return {Path(str(row["file"])).name: int(row["round_number"]) for _, row in df.iterrows()}


def _build_round_map_from_inputs(
    r1_by_round: dict[int, list[Path]],
    r2_by_round: dict[int, list[Path]] | None,
) -> dict[str, int]:
    """Fallback round map for resume paths where rounds.tsv is missing."""
    rm: dict[str, int] = {}
    for r, paths in r1_by_round.items():
        for p in paths:
            rm[p.name] = r
    if r2_by_round:
        for r, paths in r2_by_round.items():
            for p in paths:
                rm[p.name] = r
    return rm


def _run_count_stage(acc_dir: Path, lr: LibraryReport, *, resume: bool) -> None:
    """Per-round count: granular resume — skip rounds whose parquet exists."""
    extracted_name = _extracted_filename_for_mode(lr.extraction_mode)
    if extracted_name is None:
        logger.warning(
            "count stage: no countable extracted FASTA for mode=%s — skipping",
            lr.extraction_mode,
        )
        return

    for round_dir in sorted(acc_dir.glob("round_*")):
        if round_dir.name == "round_unknown":
            continue
        extracted = round_dir / extracted_name
        if not extracted.exists():
            logger.info("count stage: no %s in %s — skipping round", extracted_name, round_dir)
            continue
        out_parquet = round_dir / "counts.parquet"
        if resume and out_parquet.exists():
            logger.info("count stage: %s already present — skipping (resume)", out_parquet)
            continue
        count_fasta(extracted, out_parquet)


def _extracted_filename_for_mode(mode: str) -> str | None:
    """Map ``extraction_mode`` → the single canonical countable FASTA filename.

    Returns None for modes that should NOT be counted (the split-primer
    guard catches PAIRED_END_SPLIT_PRIMERS before this is reached, but
    the lookup is kept defensive).
    """
    return {
        "BOTH_PRIMERS_SINGLE_READ": "extracted.fasta.gz",
        "FIVE_PRIME_ONLY": "partial_5p_extracted.fasta.gz",
        "THREE_PRIME_ONLY": "partial_3p_extracted.fasta.gz",
    }.get(mode)


def _count_yaml_flags(yaml_path: Path) -> int:
    """Count flag entries in a resumed flags.yaml.

    Codex pass 1 fix: ``write_flags_yaml`` uses ``safe_dump(sort_keys=True)``
    which sorts dict keys alphabetically per entry, so the first line of an
    entry can be ``- evidence:`` rather than ``- name:`` whenever an
    ``evidence`` key is non-empty. The previous ``startswith("- name:")``
    grep undercounted resumed runs that had any non-trivial evidence. Parse
    the YAML as the source of truth.
    """
    if not yaml_path.exists():
        return 0
    try:
        payload = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
    except yaml.YAMLError:
        logger.warning("flags.yaml at %s is malformed; reporting 0 flags", yaml_path)
        return 0
    if not isinstance(payload, list):
        return 0
    return sum(1 for entry in payload if isinstance(entry, dict) and "name" in entry)


# ---------------------------------------------------------------------------
# Summary writer
# ---------------------------------------------------------------------------


_SUMMARY_COLUMNS = (
    "accession",
    "status",
    "last_stage_completed",
    "library_report_status",
    "extraction_mode",
    "required_action",
    "confidence",
    "flags_raised",
    "notes",
)


def _write_summary_tsv(path: Path, rows: list[RunRowReport]) -> None:
    """Emit a deterministic per-accession summary TSV (sorted by accession)."""
    sorted_rows = sorted(rows, key=lambda r: r.accession)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(_SUMMARY_COLUMNS) + "\n")
        for row in sorted_rows:
            values = [
                row.accession,
                row.status,
                row.last_stage_completed,
                row.library_report_status or "",
                row.extraction_mode or "",
                row.required_action or "",
                f"{row.confidence:.4f}" if row.confidence is not None else "",
                str(row.flags_raised) if row.flags_raised is not None else "",
                row.notes.replace("\t", " ").replace("\n", " "),
            ]
            fh.write("\t".join(values) + "\n")


# Re-export for type-checker callers
__all__ = ["RunReport", "RunRowReport", "RunStatus", "run_batch"]
