"""Extract orchestrator: ``LibraryReport`` -> per-round FASTAs.

Dispatches per ``LibraryReport.extraction_mode``, wiring the optional
sample-sheet demux pre-step, the strand-orientation pre-step, and the
per-mode cutadapt invocation in :mod:`selexprep.extract.trim`.

**Refusal rule**: if the LibraryReport's status is
``UNABLE_TO_INFER`` or its ``extraction_mode`` is ``UNABLE_TO_EXTRACT``,
the runner refuses without an explicit override. No silent miscalls.

**Output layout**:

::

    outdir/
    +- round_00/
    |  +- extracted.fasta.gz           (BOTH_PRIMERS_SINGLE_READ)
    |     OR partial_5p_extracted.fasta.gz             (FIVE_PRIME_ONLY)
    |     OR partial_3p_extracted.fasta.gz             (THREE_PRIME_ONLY)
    |     OR partial_5p_extracted_R1.fasta.gz +        (PAIRED_END_SPLIT_PRIMERS)
    |        partial_3p_extracted_R2.fasta.gz
    +- round_NN/...
    +- strand_report.tsv               (only if orientation in {MIXED, REVERSE})
    +- trim_reports.json               (per-round cutadapt argv + counts;
                                        precursor to the manifest)

Joined-counts (``joined_counts.tsv``) is NOT emitted in v0.1 - joining R1
+R2 by read ID alone without read merging is biologically wrong.
"""

from __future__ import annotations

import gzip
import json
import logging
from dataclasses import dataclass, field, replace
from pathlib import Path

from selexprep._io import open_gzip_text_deterministic
from selexprep.extract import strand as strand_module
from selexprep.extract import trim as trim_module
from selexprep.extract.demux import demux_sample_sheet
from selexprep.library.adapters import matches_known_adapter_prefix, reverse_complement
from selexprep.library.report import LibraryReport
from selexprep.manifest import (
    build_manifest_from_extract_result,
    read_manifest_json,
    write_manifest_json,
)

logger = logging.getLogger(__name__)


# How many reads per round to sample when computing the strand_report
# distribution. Cheap QC summary; not used for any classification call
# .
_STRAND_REPORT_SAMPLE_PER_ROUND = 10_000


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ExtractResult:
    """Summary of one ``run_extract`` invocation."""

    outputs: list[Path] = field(default_factory=list)
    strand_report_path: Path | None = None
    trim_reports_path: Path | None = None
    manifest_path: Path | None = None
    extract_diff_path: Path | None = None
    skipped_reason: str | None = None
    trim_reports: list[trim_module.TrimReport] = field(default_factory=list)

    @property
    def skipped(self) -> bool:
        return self.skipped_reason is not None


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _refusal_reason(library_report: LibraryReport) -> str | None:
    """Return a refusal reason if LR signals UNABLE; otherwise None.

    the design: refuse without an explicit override; no silent
    miscalls.
    """
    if library_report.status == "UNABLE_TO_INFER":
        return (
            f"LibraryReport status is UNABLE_TO_INFER "
            f"(reason: {library_report.failure_reason or 'unspecified'}). "
            "Pass --override-primer-5p / --override-primer-3p (and the "
            "runner will promote extraction_mode + status), or hand-edit "
            "the LibraryReport JSON to proceed."
        )
    if library_report.extraction_mode == "UNABLE_TO_EXTRACT":
        return (
            f"LibraryReport extraction_mode is UNABLE_TO_EXTRACT "
            f"(reason: {library_report.failure_reason or 'unspecified'}). "
            "Pass --override-primer-5p / --override-primer-3p (and the "
            "runner will promote extraction_mode + status), or hand-edit "
            "the LibraryReport JSON to proceed."
        )
    return None


def _group_by_round(fastq_inputs: list[Path], round_map: dict[str, int]) -> dict[int, list[Path]]:
    """Group FASTQ paths by round number using basename lookup in round_map."""
    out: dict[int, list[Path]] = {}
    for fq in fastq_inputs:
        r = round_map.get(fq.name)
        if r is None:
            raise ValueError(f"FASTQ {fq.name!r} not in round map")
        out.setdefault(r, []).append(fq)
    return out


def _check_no_clobber(targets: list[Path], rebuild: bool) -> None:
    """Refuse to overwrite existing outputs unless ``rebuild`` is True."""
    existing = [t for t in targets if t.exists()]
    if existing and not rebuild:
        names = ", ".join(str(p) for p in existing)
        raise FileExistsError(f"Output already exists: {names}. Pass --rebuild to overwrite.")


def _sample_sequences(fq: Path, max_reads: int) -> list[str]:
    """Read up to ``max_reads`` sequences from a FASTQ.gz."""
    seqs: list[str] = []
    with gzip.open(fq, "rt", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seqs.append(line.strip())
                if len(seqs) >= max_reads:
                    break
    return seqs


def _maybe_reorient(
    round_inputs: dict[int, list[Path]],
    reorient_dir: Path,
    *,
    do_revcomp: bool,
) -> dict[int, list[Path]]:
    """Optionally revcomp every read in every input; return new path mapping.

    When ``do_revcomp`` is False this is a no-op (returns the input map
    unchanged). When True, every FASTQ.gz is re-emitted to
    ``reorient_dir/round_NN/<basename>`` with each read revcomp'd.
    """
    if not do_revcomp:
        return round_inputs
    out: dict[int, list[Path]] = {}
    for r, paths in round_inputs.items():
        round_dir = reorient_dir / f"round_{r:02d}"
        round_dir.mkdir(parents=True, exist_ok=True)
        out[r] = []
        for fq in paths:
            target = round_dir / fq.name
            n_records = strand_module.reorient_fastq_gz(fq, target)
            logger.info("reoriented %d records: %s -> %s", n_records, fq, target)
            out[r].append(target)
    return out


def _compute_strand_distributions(
    round_inputs: dict[int, list[Path]],
    primer_5p: str | None,
    primer_3p: str | None,
) -> dict[int, dict[str, int]]:
    """Sample reads per round and compute forward/reverse/ambiguous counts."""
    dists: dict[int, dict[str, int]] = {}
    for r, paths in round_inputs.items():
        seqs: list[str] = []
        for fq in paths:
            seqs.extend(_sample_sequences(fq, _STRAND_REPORT_SAMPLE_PER_ROUND))
        dists[r] = strand_module.detect_strand_distribution(seqs, primer_5p, primer_3p)
    return dists


def _write_trim_reports(reports: list[trim_module.TrimReport], path: Path) -> None:
    """Write per-round cutadapt argv + counts as JSON (manifest precursor)."""
    payload = [
        {
            "cutadapt_cmd": r.cutadapt_cmd,
            "n_in": r.n_in,
            "n_out": r.n_out,
            "return_code": r.return_code,
            "output_paths": [str(p) for p in r.output_paths],
        }
        for r in reports
    ]
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _concat_fasta_gz_deterministic(parts: list[Path], target: Path) -> None:
    """Concatenate `parts` (gzipped FASTAs) into `target` deterministically.

    Reads each part in order, decompresses, and re-emits through
    ``open_gzip_text_deterministic`` (mtime=0 header). Caller is responsible
    for the input order: pass the temp files in the same order as the
    per-input trim calls. Temp inputs are NOT cleaned up here — the caller
    does that.
    """
    target.parent.mkdir(parents=True, exist_ok=True)
    with open_gzip_text_deterministic(target) as out:
        for part in parts:
            with gzip.open(part, "rt", encoding="utf-8") as fh:
                for line in fh:
                    out.write(line)


def _trim_round_single_end(
    round_inputs: list[Path],
    final_target: Path,
    trim_fn,  # type: ignore[no-untyped-def]
    **trim_kwargs,
) -> list[trim_module.TrimReport]:
    """Trim N single-end FASTQs of one round to `final_target`.

    For one input, calls ``trim_fn`` directly (fast path). For multiple
    inputs, trims each to a per-input temp target, then concatenates
    deterministically into ``final_target`` (previously
    each iteration overwrote the same target, dropping all but the last).

    Returns one TrimReport per input — their ``output_paths`` are rewritten
    to ``[final_target]`` so the downstream per-round aggregation in
    ``_trim_reports_by_round_from_objects`` groups them correctly.
    """
    if len(round_inputs) == 1:
        return [trim_fn(round_inputs[0], final_target, **trim_kwargs)]

    final_target.parent.mkdir(parents=True, exist_ok=True)
    temps: list[Path] = []
    reports: list[trim_module.TrimReport] = []
    try:
        for i, fq in enumerate(round_inputs):
            tmp = final_target.parent / f".{final_target.name}.part_{i:03d}.fasta.gz"
            report = trim_fn(fq, tmp, **trim_kwargs)
            temps.append(tmp)
            reports.append(replace(report, output_paths=[final_target]))
        _concat_fasta_gz_deterministic(temps, final_target)
    finally:
        for tmp in temps:
            tmp.unlink(missing_ok=True)
    return reports


def _trim_round_paired_split(
    r1_inputs: list[Path],
    r2_inputs: list[Path],
    final_r1_target: Path,
    final_r2_target: Path,
    *,
    primer_5p_r1: str,
    primer_5p_r2: str,
) -> list[trim_module.TrimReport]:
    """Same as ``_trim_round_single_end`` but for paired-end split mode.

    R1 and R2 are processed in lock-step (one cutadapt call per (R1, R2)
    pair). Multiple pairs in the same round are concatenated into the
    two final targets deterministically.
    """
    if len(r1_inputs) != len(r2_inputs):
        raise ValueError(f"R1/R2 pair-count mismatch: {len(r1_inputs)} R1 vs {len(r2_inputs)} R2")
    if len(r1_inputs) == 1:
        report = trim_module.trim_paired_split(
            r1_inputs[0],
            r2_inputs[0],
            final_r1_target,
            final_r2_target,
            primer_5p_r1=primer_5p_r1,
            primer_5p_r2=primer_5p_r2,
        )
        return [report]

    final_r1_target.parent.mkdir(parents=True, exist_ok=True)
    final_r2_target.parent.mkdir(parents=True, exist_ok=True)
    temps_r1: list[Path] = []
    temps_r2: list[Path] = []
    reports: list[trim_module.TrimReport] = []
    try:
        for i, (r1, r2) in enumerate(zip(r1_inputs, r2_inputs, strict=True)):
            tmp_r1 = final_r1_target.parent / f".{final_r1_target.name}.part_{i:03d}.fasta.gz"
            tmp_r2 = final_r2_target.parent / f".{final_r2_target.name}.part_{i:03d}.fasta.gz"
            report = trim_module.trim_paired_split(
                r1,
                r2,
                tmp_r1,
                tmp_r2,
                primer_5p_r1=primer_5p_r1,
                primer_5p_r2=primer_5p_r2,
            )
            temps_r1.append(tmp_r1)
            temps_r2.append(tmp_r2)
            reports.append(replace(report, output_paths=[final_r1_target, final_r2_target]))
        _concat_fasta_gz_deterministic(temps_r1, final_r1_target)
        _concat_fasta_gz_deterministic(temps_r2, final_r2_target)
    finally:
        for tmp in temps_r1 + temps_r2:
            tmp.unlink(missing_ok=True)
    return reports


def _trim_reports_by_round_from_objects(
    reports: list[trim_module.TrimReport],
) -> dict[int, dict[str, int]]:
    """Aggregate live TrimReport objects by round (read from output path's parent)."""
    by_round: dict[int, dict[str, int]] = {}
    for r in reports:
        if not r.output_paths:
            continue
        round_dir = r.output_paths[0].parent.name
        if not round_dir.startswith("round_"):
            continue
        round_num = int(round_dir.split("_")[1])
        agg = by_round.setdefault(round_num, {"n_in": 0, "n_out": 0})
        agg["n_in"] += r.n_in
        agg["n_out"] += r.n_out
    return by_round


def _trim_reports_by_round_from_json(payload: list[dict]) -> dict[int, dict[str, int]]:
    """Aggregate trim_reports.json entries by round."""
    by_round: dict[int, dict[str, int]] = {}
    for entry in payload:
        paths = entry.get("output_paths") or []
        if not paths:
            continue
        round_dir = Path(paths[0]).parent.name
        if not round_dir.startswith("round_"):
            continue
        round_num = int(round_dir.split("_")[1])
        agg = by_round.setdefault(round_num, {"n_in": 0, "n_out": 0})
        agg["n_in"] += int(entry.get("n_in", 0))
        agg["n_out"] += int(entry.get("n_out", 0))
    return by_round


def _read_baseline_for_diff(
    outdir: Path,
) -> tuple[LibraryReport, dict[int, dict[str, int]]] | None:
    """Read the previous run's manifest + trim_reports.json from `outdir`.

    Returns ``(baseline_library_report, baseline_by_round)`` or ``None``
    if either artifact is missing or unparseable (caller skips diff
    emission in that case).
    """
    manifest_path = outdir / "selexprep_manifest.json"
    trim_path = outdir / "trim_reports.json"
    if not (manifest_path.exists() and trim_path.exists()):
        return None
    try:
        baseline_manifest = read_manifest_json(manifest_path)
        baseline_lr = baseline_manifest.library_report
        baseline_payload = json.loads(trim_path.read_text(encoding="utf-8"))
        baseline_by_round = _trim_reports_by_round_from_json(baseline_payload)
    except Exception as e:
        # Defensive: baseline artifacts may be corrupt in many ways; the
        # diff is informational, so we degrade gracefully rather than fail
        # the whole rebuild on a malformed baseline.
        logger.warning("Could not read baseline for diff: %s", e)
        return None
    return baseline_lr, baseline_by_round


def _write_extract_diff(
    *,
    baseline_lr: LibraryReport,
    baseline_by_round: dict[int, dict[str, int]],
    new_lr: LibraryReport,
    new_by_round: dict[int, dict[str, int]],
    path: Path,
) -> None:
    """Write per-round comparison TSV.

    Columns: ``round\\tprimer_5p_baseline\\tprimer_5p_new\\tprimer_3p_baseline\\t
    primer_3p_new\\tn_in\\tn_out_baseline\\tn_out_new\\tdelta_n_out``.
    Sorted by round.
    """
    rounds_all = sorted(set(baseline_by_round) | set(new_by_round))
    header = (
        "round\tprimer_5p_baseline\tprimer_5p_new\tprimer_3p_baseline\t"
        "primer_3p_new\tn_in\tn_out_baseline\tn_out_new\tdelta_n_out\n"
    )
    lines = [header]
    for r in rounds_all:
        b = baseline_by_round.get(r, {"n_in": 0, "n_out": 0})
        n = new_by_round.get(r, {"n_in": 0, "n_out": 0})
        delta = n["n_out"] - b["n_out"]
        lines.append(
            f"{r}\t"
            f"{baseline_lr.primer_5p or ''}\t{new_lr.primer_5p or ''}\t"
            f"{baseline_lr.primer_3p or ''}\t{new_lr.primer_3p or ''}\t"
            f"{n['n_in']}\t{b['n_out']}\t{n['n_out']}\t{delta}\n"
        )
    path.write_text("".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def run_extract(
    library_report: LibraryReport,
    fastq_inputs: list[Path],
    outdir: Path,
    *,
    round_map: dict[str, int],
    sample_sheet: Path | None = None,
    rebuild: bool = False,
    paired_r2_inputs: dict[int, list[Path]] | None = None,
    override_primer_5p: str | None = None,
    override_primer_3p: str | None = None,
    accession: str | None = None,
    bioproject_id: str | None = None,
    runs: list[str] | None = None,
    parameters: dict[str, str] | None = None,
) -> ExtractResult:
    """Run the extraction pipeline for one dataset.

    owns the cutadapt-driven trimming. Adds:

    - ``override_primer_{5p,3p}``: clone the LibraryReport with swapped
      primer fields. Without ``rebuild`` the override outputs go to a
      separate ``<outdir>/overridden/`` subtree (avoids clobbering
      baseline). With ``rebuild`` the baseline is overwritten and an
      ``extract_diff.tsv`` is emitted comparing baseline vs new
      per-round read counts.
    - ``selexprep_manifest.json`` is now emitted automatically.

    Args:
        library_report: typed contract from 's ``compute_library_report``.
        fastq_inputs: R1 FASTQ(.gz) files (or single-end inputs).
        outdir: where outputs are written. Override-only path routes
            outputs to ``<outdir>/overridden/``.
        round_map: basename -> round number.
        sample_sheet: optional TSV for demux pre-step.
        rebuild: overwrite existing outputs.
        paired_r2_inputs: optional ``{round: [r2_fastq, ...]}``.
        override_primer_5p / override_primer_3p: optional primer
            overrides applied via ``LibraryReport.model_copy``.
        accession / bioproject_id / runs / parameters: provenance fields
            captured into the emitted ``selexprep_manifest.json``.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # apply override BEFORE the refusal check.
    # The design says ``--override-primer-{5p,3p}`` is an
    # explicit user bypass for UNABLE_TO_INFER / UNABLE_TO_EXTRACT. To honor
    # that, override has to (1) replace the primer values AND (2) promote the
    # classification fields the refusal check inspects (status,
    # extraction_mode), otherwise the override is a no-op for the cases it
    # was designed to handle.
    has_override = override_primer_5p is not None or override_primer_3p is not None
    if has_override:
        updates: dict[str, object] = {}
        if override_primer_5p is not None:
            updates["primer_5p"] = override_primer_5p
        if override_primer_3p is not None:
            updates["primer_3p"] = override_primer_3p

        # Warn (but don't refuse) if an override matches a known sequencing
        # adapter prefix — the user typed it explicitly, so respect the
        # escape-hatch semantics, but surface the foot-gun (
        # non-blocking).
        for label, primer in (("5p", override_primer_5p), ("3p", override_primer_3p)):
            if primer and matches_known_adapter_prefix(primer):
                logger.warning(
                    "Override primer %s matches a known sequencing adapter prefix; "
                    "this is allowed as an explicit escape hatch but is usually a "
                    "mistake (override value: %r)",
                    label,
                    primer,
                )

        # Promote extraction_mode when the baseline was UNABLE_TO_EXTRACT and
        # the user is declaring biology via override. The resulting mode is
        # inferred from which primer sides are now defined.
        has_5p_after = override_primer_5p is not None or library_report.primer_5p is not None
        has_3p_after = override_primer_3p is not None or library_report.primer_3p is not None
        if library_report.extraction_mode == "UNABLE_TO_EXTRACT":
            if has_5p_after and has_3p_after:
                updates["extraction_mode"] = "BOTH_PRIMERS_SINGLE_READ"
                updates["full_insert_recovered"] = True
                updates["required_action"] = "NONE"
            elif has_5p_after:
                updates["extraction_mode"] = "FIVE_PRIME_ONLY"
                updates["full_insert_recovered"] = False
                updates["required_action"] = "NONE"
            elif has_3p_after:
                updates["extraction_mode"] = "THREE_PRIME_ONLY"
                updates["full_insert_recovered"] = False
                updates["required_action"] = "NONE"
            # else: no primers at all even after override — UNABLE stays,
            # refusal will fire below with the hand-edit-LR hint.

        # Manual override = MEDIUM confidence by convention; the user must
        # hand-edit the LibraryReport JSON to claim HIGH.
        if library_report.status == "UNABLE_TO_INFER":
            updates["status"] = "MEDIUM"
            updates["failure_reason"] = None

        library_report = library_report.model_copy(update=updates)
        logger.info(
            "Applied primer overrides: 5p=%r 3p=%r -> mode=%s status=%s",
            override_primer_5p,
            override_primer_3p,
            library_report.extraction_mode,
            library_report.status,
        )

    # Refusal check happens AFTER override so the user's explicit escape
    # hatch can clear the UNABLE state.
    refusal = _refusal_reason(library_report)
    if refusal is not None:
        logger.warning("Refusing to extract: %s", refusal)
        return ExtractResult(skipped_reason=refusal)

    # When override is set but rebuild is NOT, redirect outputs to a
    # subtree to avoid clobbering the baseline (the design: rebuild
    # required to overwrite).
    if has_override and not rebuild:
        outdir = outdir / "overridden"
        outdir.mkdir(parents=True, exist_ok=True)

    # When rebuild + override: read baseline manifest + trim_reports BEFORE
    # the per-round outputs get overwritten, so we can emit extract_diff.tsv.
    baseline_for_diff: tuple[LibraryReport, dict[int, dict[str, int]]] | None = None
    if rebuild and has_override:
        baseline_for_diff = _read_baseline_for_diff(outdir)

    # Optional demux pre-step (multiplexed input).
    round_inputs: dict[int, list[Path]] | None = None
    input_root: Path | None = None  # set in sample-sheet mode so manifest
    # input_sha256 keys are relative to demux dir (demuxed basenames
    # collide across rounds — follow-up fix).
    if sample_sheet is not None:
        demux_dir = outdir / "demux"
        demux_dir.mkdir(parents=True, exist_ok=True)
        demux_sample_sheet(sample_sheet, out_root=demux_dir)
        # After demux, per-round files live under demux/round_NN/<srr>.fastq.gz
        # (single-end) or demux/round_NN/<srr>_1.fastq.gz (paired-end).
        # We collect BOTH R1 and R2 paths here (previously
        # paired_r2_inputs was never rebuilt from demux outputs, so
        # PAIRED_END_SPLIT_PRIMERS + --sample-sheet failed with "requires
        # --paired-r2 inputs"). And we build round_inputs path-aware
        # directly from the demuxed parent dirs — going via the basename-
        # keyed round_map would collapse across rounds because the demuxed
        # files share basenames (e.g. ``srr_1.fastq.gz`` in every
        # ``round_NN/`` folder).
        all_demuxed = sorted(demux_dir.glob("round_*/*.fastq.gz"))
        r2_demuxed = [p for p in all_demuxed if p.name.endswith("_2.fastq.gz")]
        r1_demuxed = [p for p in all_demuxed if not p.name.endswith("_2.fastq.gz")]
        if not r1_demuxed:
            return ExtractResult(
                skipped_reason=f"Sample-sheet demux produced no FASTQs under {demux_dir}",
            )
        fastq_inputs = r1_demuxed
        input_root = demux_dir
        round_inputs = {}
        for p in r1_demuxed:
            r = int(p.parent.name.split("_")[1])
            round_inputs.setdefault(r, []).append(p)
        if r2_demuxed:
            paired_r2_inputs = {}
            for p in r2_demuxed:
                r = int(p.parent.name.split("_")[1])
                paired_r2_inputs.setdefault(r, []).append(p)

    # Path-aware sample-sheet path already populated round_inputs; otherwise
    # group via the basename-keyed round_map (the CLI guarantees basenames
    # are unique on this code path).
    if round_inputs is None:
        round_inputs = _group_by_round(fastq_inputs, round_map)

    # Strand orientation pre-step.
    primer_5p = library_report.primer_5p
    primer_3p = library_report.primer_3p
    orientation = library_report.orientation
    strand_distributions: dict[int, dict[str, int]] = {}
    if orientation in ("REVERSE", "MIXED"):
        strand_distributions = _compute_strand_distributions(round_inputs, primer_5p, primer_3p)

    do_revcomp = orientation == "REVERSE"
    reorient_dir = outdir / "reoriented"
    round_inputs = _maybe_reorient(round_inputs, reorient_dir, do_revcomp=do_revcomp)

    # Per-mode trim dispatch.
    extraction_mode = library_report.extraction_mode
    trim_reports: list[trim_module.TrimReport] = []
    all_outputs: list[Path] = []

    rounds_sorted = sorted(round_inputs.keys())
    if extraction_mode == "BOTH_PRIMERS_SINGLE_READ":
        if not primer_5p or not primer_3p:
            return ExtractResult(
                skipped_reason="BOTH_PRIMERS_SINGLE_READ requires both primers in LibraryReport",
            )
        targets = [outdir / f"round_{r:02d}" / "extracted.fasta.gz" for r in rounds_sorted]
        _check_no_clobber(targets, rebuild)
        for r, target in zip(rounds_sorted, targets, strict=True):
            reports = _trim_round_single_end(
                round_inputs[r],
                target,
                trim_module.trim_single_end_linked,
                primer_5p=primer_5p,
                primer_3p=primer_3p,
            )
            trim_reports.extend(reports)
            all_outputs.append(target)

    elif extraction_mode == "FIVE_PRIME_ONLY":
        if not primer_5p:
            return ExtractResult(skipped_reason="FIVE_PRIME_ONLY requires primer_5p")
        targets = [
            outdir / f"round_{r:02d}" / "partial_5p_extracted.fasta.gz" for r in rounds_sorted
        ]
        _check_no_clobber(targets, rebuild)
        for r, target in zip(rounds_sorted, targets, strict=True):
            reports = _trim_round_single_end(
                round_inputs[r],
                target,
                trim_module.trim_single_end_5p,
                primer_5p=primer_5p,
            )
            trim_reports.extend(reports)
            all_outputs.append(target)

    elif extraction_mode == "THREE_PRIME_ONLY":
        if not primer_3p:
            return ExtractResult(skipped_reason="THREE_PRIME_ONLY requires primer_3p")
        targets = [
            outdir / f"round_{r:02d}" / "partial_3p_extracted.fasta.gz" for r in rounds_sorted
        ]
        _check_no_clobber(targets, rebuild)
        for r, target in zip(rounds_sorted, targets, strict=True):
            reports = _trim_round_single_end(
                round_inputs[r],
                target,
                trim_module.trim_single_end_3p,
                primer_3p=primer_3p,
            )
            trim_reports.extend(reports)
            all_outputs.append(target)

    elif extraction_mode == "PAIRED_END_SPLIT_PRIMERS":
        if not primer_5p or not primer_3p:
            return ExtractResult(
                skipped_reason="PAIRED_END_SPLIT_PRIMERS requires both primers in LibraryReport",
            )
        if paired_r2_inputs is None:
            return ExtractResult(
                skipped_reason="PAIRED_END_SPLIT_PRIMERS requires --paired-r2 inputs",
            )
        primer_5p_r2 = reverse_complement(primer_3p)
        r1_targets = [
            outdir / f"round_{r:02d}" / "partial_5p_extracted_R1.fasta.gz" for r in rounds_sorted
        ]
        r2_targets = [
            outdir / f"round_{r:02d}" / "partial_3p_extracted_R2.fasta.gz" for r in rounds_sorted
        ]
        _check_no_clobber(r1_targets + r2_targets, rebuild)
        for r, r1_target, r2_target in zip(rounds_sorted, r1_targets, r2_targets, strict=True):
            r1_paths = round_inputs[r]
            r2_paths = paired_r2_inputs.get(r, [])
            if len(r1_paths) != len(r2_paths):
                return ExtractResult(
                    skipped_reason=(
                        f"round {r}: R1/R2 pair-count mismatch "
                        f"({len(r1_paths)} R1 vs {len(r2_paths)} R2)"
                    ),
                )
            reports = _trim_round_paired_split(
                r1_paths,
                r2_paths,
                r1_target,
                r2_target,
                primer_5p_r1=primer_5p,
                primer_5p_r2=primer_5p_r2,
            )
            trim_reports.extend(reports)
            all_outputs.append(r1_target)
            all_outputs.append(r2_target)
    else:
        return ExtractResult(
            skipped_reason=f"Unknown extraction_mode: {extraction_mode!r}",
        )

    # Emit strand report (if applicable).
    strand_report_path: Path | None = None
    if strand_distributions:
        strand_report_path = outdir / "strand_report.tsv"
        strand_module.write_strand_report(strand_distributions, strand_report_path)
        all_outputs.append(strand_report_path)

    # Emit trim reports JSON.
    trim_reports_path = outdir / "trim_reports.json"
    _write_trim_reports(trim_reports, trim_reports_path)
    all_outputs.append(trim_reports_path)

    # emit extract_diff.tsv when we have a baseline to compare against.
    extract_diff_path: Path | None = None
    if baseline_for_diff is not None:
        baseline_lr, baseline_by_round = baseline_for_diff
        new_by_round = _trim_reports_by_round_from_objects(trim_reports)
        extract_diff_path = outdir / "extract_diff.tsv"
        _write_extract_diff(
            baseline_lr=baseline_lr,
            baseline_by_round=baseline_by_round,
            new_lr=library_report,
            new_by_round=new_by_round,
            path=extract_diff_path,
        )
        all_outputs.append(extract_diff_path)

    # emit the reproducibility manifest.
    # Collate input paths (R1 + R2 if paired).
    input_paths: list[Path] = list(fastq_inputs)
    if paired_r2_inputs is not None:
        for paths in paired_r2_inputs.values():
            input_paths.extend(paths)

    manifest = build_manifest_from_extract_result(
        library_report=library_report,
        input_paths=input_paths,
        output_paths=all_outputs,
        output_root=outdir,
        input_root=input_root,
        accession=accession,
        bioproject_id=bioproject_id,
        runs=runs or [],
        parameters=parameters or {},
    )
    manifest_path = outdir / "selexprep_manifest.json"
    write_manifest_json(manifest, manifest_path)
    all_outputs.append(manifest_path)

    return ExtractResult(
        outputs=all_outputs,
        strand_report_path=strand_report_path,
        trim_reports_path=trim_reports_path,
        manifest_path=manifest_path,
        extract_diff_path=extract_diff_path,
        skipped_reason=None,
        trim_reports=trim_reports,
    )
