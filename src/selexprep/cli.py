"""Typer CLI dispatcher for selexprep.

Phase 0 scaffold + progressive wiring: subcommands not yet implemented exit
with code 2. ``catalog`` (Phase 1.5) and ``detect`` (Phase 2) are live;
``inspect``/``fetch``/``extract``/``count``/``qc``/``run`` arrive in later
phases.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import pandas as pd
import typer

from selexprep import __version__
from selexprep.catalog.cli import app as catalog_app
from selexprep.count.counter import count_fasta, count_fastq_pretrimmed
from selexprep.extract import run_extract
from selexprep.library import (
    compute_library_report,
    read_library_report_json,
    write_library_report_json,
)
from selexprep.qc.runner import run_qc

logger = logging.getLogger(__name__)

app = typer.Typer(
    name="selexprep",
    help="Accession-first preprocessing for public HT-SELEX with primer auto-inference.",
    no_args_is_help=True,
    add_completion=False,
)
app.add_typer(catalog_app, name="catalog")


def _not_implemented(name: str) -> None:
    typer.secho(
        f"selexprep {name}: not yet implemented (Phase 0 scaffold).",
        fg=typer.colors.YELLOW,
        err=True,
    )
    raise typer.Exit(code=2)


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"selexprep {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        None,
        "--version",
        "-V",
        callback=_version_callback,
        is_eager=True,
        help="Show version and exit.",
    ),
) -> None:
    """Accession-first preprocessing for public HT-SELEX with primer auto-inference."""


@app.command()
def inspect(
    accession: str = typer.Argument(..., help="ENA / SRA / DDBJ accession."),
    outdir: Path | None = typer.Option(
        None, "--outdir", help="If given, also write inspect.json to this directory."
    ),
    timeout_s: int = typer.Option(30, "--timeout-s", help="HTTP timeout (seconds)."),
) -> None:
    """Preview accession metadata (round count, library_strategy, files) without downloading."""
    from selexprep.fetch.inspect import inspect_accession, write_inspect_json

    try:
        report = inspect_accession(accession, timeout_s=timeout_s)
    except ValueError as e:
        typer.secho(f"inspect: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=2) from e

    typer.echo(f"Accession:        {report.accession}")
    typer.echo(f"BioProject:       {report.bioproject_id or '-'}")
    typer.echo(f"Study title:      {report.study_title or '-'}")
    typer.echo(
        f"library_strategy: {report.library_strategy or '-'}  (SRA verbatim; not classified)"
    )
    typer.echo(f"library_source:   {report.library_source or '-'}")
    typer.echo(f"Runs ({len(report.runs)}):")
    for run in report.runs:
        size_str = ", ".join(f"{n} B" for n in run.fastq_size_bytes) or "-"
        typer.echo(
            f"  {run.run_accession}  reads={run.read_count}  bases={run.base_count}  "
            f"files=[{size_str}]"
        )

    if outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        out_path = outdir / "inspect.json"
        write_inspect_json(report, out_path)
        typer.echo(f"\nwrote {out_path}")


@app.command()
def fetch(
    accession: str = typer.Argument(..., help="ENA / SRA / DDBJ accession."),
    outdir: Path = typer.Option(..., "--outdir", help="Output directory."),
    backend: str = typer.Option(
        "ena",
        "--backend",
        help=(
            "Download backend. Default 'ena' (paper-grade reproducibility; "
            "fail-fast if ENA can't serve). Use 'auto' for the convenience "
            "fallback chain (ENA → kingfisher → sra-toolkit; GPL-3.0 "
            "tools opt-in)."
        ),
    ),
    allow_manual_review: bool = typer.Option(
        False,
        "--allow-manual-review",
        help=(
            "NONE-confidence runs are downloaded to round_unknown/ and "
            "surfaced in manual_review.tsv; they NEVER enter rounds.tsv."
        ),
    ),
    dry_run: bool = typer.Option(False, "--dry-run", help="Print plan; no FASTQ download."),
    timeout_s: int = typer.Option(30, "--timeout-s", help="HTTP timeout (seconds)."),
) -> None:
    """Download FASTQ + metadata for an accession; round map auto-populated.

    Emits ``outdir/rounds.tsv`` (the trusted-assignments contract for
    detect/extract), per-round FASTQs under ``round_NN/``, plus
    ``fetch_metadata.json`` as the audit trail. Refuses upfront if any
    run has a NONE-confidence round assignment unless
    ``--allow-manual-review`` is passed.
    """
    from selexprep.fetch import run_fetch as run_fetch_fn

    allowed_backends = {"auto", "ena", "kingfisher", "sra"}
    if backend not in allowed_backends:
        typer.secho(
            f"fetch: invalid --backend {backend!r}; allowed: {sorted(allowed_backends)}",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=2)

    try:
        result = run_fetch_fn(
            accession,
            outdir,
            backend=backend,  # type: ignore[arg-type]
            allow_manual_review=allow_manual_review,
            dry_run=dry_run,
            timeout_s=timeout_s,
        )
    except ValueError as e:
        typer.secho(f"fetch: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=2) from e

    if result.refused_reason is not None:
        typer.secho(f"fetch: refused — {result.refused_reason}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=2)

    typer.echo(f"Accession:        {result.plan.accession}")
    typer.echo(f"BioProject:       {result.plan.bioproject_id or '-'}")
    typer.echo(f"Study title:      {result.plan.study_title or '-'}")
    typer.echo(f"Runs:             {len(result.plan.runs)}")
    typer.echo(f"  downloaded:     {len(result.downloaded_srrs)}")
    typer.echo(f"  already-present:{len(result.skipped_srrs)}")
    typer.echo(f"  failed:         {len(result.failed_srrs)}")
    typer.echo(f"  manual-review:  {len(result.manual_review_srrs)}")
    if result.rounds_tsv is not None:
        typer.echo(f"rounds.tsv -> {result.rounds_tsv}")
    if result.manual_review_tsv is not None:
        typer.echo(f"manual_review.tsv -> {result.manual_review_tsv}")
    typer.echo(f"fetch_metadata.json -> {result.fetch_metadata_json}")

    if result.failed_srrs:
        typer.secho(
            f"fetch: {len(result.failed_srrs)} run(s) failed: {', '.join(result.failed_srrs)}",
            fg=typer.colors.YELLOW,
            err=True,
        )
        raise typer.Exit(code=1)


def _load_round_map(path: Path) -> dict[str, int]:
    """Parse a 2-column TSV (``file<TAB>round_number``) → basename → round.

    Matching is by basename so relative/absolute path differences between the
    round map and the FASTQ arguments are tolerated. The TSV must have a
    header row with columns ``file`` and ``round_number``.
    """
    df = pd.read_csv(path, sep="\t")
    if "file" not in df.columns or "round_number" not in df.columns:
        raise typer.BadParameter(
            f"round map {path} must have columns 'file' and 'round_number'; "
            f"found {list(df.columns)}"
        )
    out: dict[str, int] = {}
    for _, row in df.iterrows():
        out[Path(str(row["file"])).name] = int(row["round_number"])
    return out


def _read_fastq_sequences(path: Path) -> list[str]:
    """Read sequence lines (every 4th line, offset 1) from a FASTQ(.gz)."""
    opener = gzip.open if path.suffix == ".gz" else open
    seqs: list[str] = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seqs.append(line.strip())
    return seqs


@app.command()
def detect(
    fastq: list[Path] = typer.Argument(..., help="Input FASTQ files."),
    round_map: Path | None = typer.Option(
        None,
        "--round-map",
        help="TSV mapping FASTQ file → round number (required for local FASTQs).",
    ),
    paired_r2: list[Path] | None = typer.Option(
        None,
        "--paired-r2",
        help="R2 FASTQ files (one or more); enables paired split-primer detection.",
    ),
    outdir: Path = typer.Option(..., "--outdir", help="Output directory."),
    sampling_seed: int = typer.Option(
        42, "--sampling-seed", help="Seed for primer-inference subsampling RNG."
    ),
    max_reads_per_round: int | None = typer.Option(
        None,
        "--max-reads-per-round",
        help="Subsample each round to at most N reads (default: use all).",
    ),
) -> None:
    """Auto-infer primers + library structure from FASTQ; emit LibraryReport JSON."""
    if round_map is None:
        typer.secho(
            "detect: --round-map is required for local FASTQs (cross-round "
            "persistence is a core inference signal).",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=2)

    all_input_paths: list[Path] = list(fastq) + list(paired_r2 or [])
    all_basenames = [p.name for p in all_input_paths]
    if len(set(all_basenames)) != len(all_basenames):
        from collections import Counter as _Counter

        dups = sorted(name for name, count in _Counter(all_basenames).items() if count > 1)
        typer.secho(
            "detect: duplicate FASTQ basenames in inputs: "
            f"{dups}. Round-map matching is by basename.",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=2)

    round_by_basename = _load_round_map(round_map)
    sequences_by_round: dict[int, list[str]] = {}
    for fq in fastq:
        r = round_by_basename.get(fq.name)
        if r is None:
            typer.secho(
                f"detect: FASTQ {fq.name!r} not in round map {round_map}",
                fg=typer.colors.RED,
                err=True,
            )
            raise typer.Exit(code=2)
        sequences_by_round.setdefault(r, []).extend(_read_fastq_sequences(fq))

    paired_mate_streams: dict[int, list[str]] | None = None
    read_source = "R1"
    if paired_r2:
        paired_mate_streams = {}
        read_source = "R1_AND_R2"
        for fq in paired_r2:
            r = round_by_basename.get(fq.name)
            if r is None:
                typer.secho(
                    f"detect: R2 FASTQ {fq.name!r} not in round map {round_map}",
                    fg=typer.colors.RED,
                    err=True,
                )
                raise typer.Exit(code=2)
            paired_mate_streams.setdefault(r, []).extend(_read_fastq_sequences(fq))

    report = compute_library_report(
        sequences_by_round,
        read_source=read_source,  # type: ignore[arg-type]
        paired_mate_streams=paired_mate_streams,
        sampling_seed=sampling_seed,
        max_reads_per_round=max_reads_per_round,
    )

    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / "library_report.json"
    sha = write_library_report_json(report, out_path)

    typer.echo(f"library_report.json -> {out_path}")
    typer.echo(f"  sha256:          {sha[:12]}...")
    typer.echo(f"  extraction_mode: {report.extraction_mode}")
    typer.echo(f"  required_action: {report.required_action}")
    typer.echo(f"  status:          {report.status}")
    if report.failure_reason:
        typer.echo(f"  failure_reason:  {report.failure_reason}")


@app.command()
def extract(
    fastq: list[Path] = typer.Argument(..., help="Input FASTQ files (R1 or single-end)."),
    library_report: Path = typer.Option(
        ..., "--library-report", help="Path to LibraryReport JSON from `detect`."
    ),
    round_map: Path | None = typer.Option(
        None,
        "--round-map",
        help="TSV mapping FASTQ file -> round number (required unless --sample-sheet).",
    ),
    sample_sheet: Path | None = typer.Option(
        None, "--sample-sheet", help="Optional sample sheet for demultiplexing."
    ),
    paired_r2: list[Path] | None = typer.Option(
        None,
        "--paired-r2",
        help="R2 FASTQ files (one per R1 in --fastq); required for PAIRED_END_SPLIT_PRIMERS.",
    ),
    override_primer_5p: str | None = typer.Option(None, "--override-primer-5p"),
    override_primer_3p: str | None = typer.Option(None, "--override-primer-3p"),
    rebuild: bool = typer.Option(
        False, "--rebuild", help="Force re-extraction; overwrite existing outputs."
    ),
    outdir: Path = typer.Option(..., "--outdir", help="Output directory."),
) -> None:
    """Trim primers + extract random region per round.

    Phase 4: ``--override-primer-{5p,3p}`` now applies primer overrides
    via ``LibraryReport.model_copy``. Pair with ``--rebuild`` to overwrite
    the baseline outputs in place AND emit ``extract_diff.tsv`` comparing
    baseline vs override per-round read counts.
    """
    # Basename collision check (Codex pass 1 fix): the round-map lookup is
    # basename-keyed, so two input FASTQs with the same name in different
    # directories would silently overwrite each other in the lookup dict.
    # Refuse early with a clear error rather than producing a wrong result.
    all_input_paths: list[Path] = list(fastq) + list(paired_r2 or [])
    all_basenames = [p.name for p in all_input_paths]
    if len(set(all_basenames)) != len(all_basenames):
        from collections import Counter as _Counter

        dups = sorted(name for name, count in _Counter(all_basenames).items() if count > 1)
        typer.secho(
            "extract: duplicate FASTQ basenames in --fastq + --paired-r2: "
            f"{dups}. Round-map matching is by basename — rename or "
            "stage inputs in non-colliding paths.",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=2)

    if sample_sheet is None and round_map is None:
        typer.secho(
            "extract: either --round-map or --sample-sheet is required "
            "(per-round routing of input FASTQs).",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=2)

    lr = read_library_report_json(library_report)

    # Round-map assignment when not using sample-sheet demux.
    rm: dict[str, int] = {}
    if round_map is not None:
        rm = _load_round_map(round_map)
        # Sanity-check inputs are covered.
        for fq in fastq:
            if fq.name not in rm:
                typer.secho(
                    f"extract: FASTQ {fq.name!r} not in round map {round_map}",
                    fg=typer.colors.RED,
                    err=True,
                )
                raise typer.Exit(code=2)

    # Paired R2 grouping.
    paired_r2_inputs: dict[int, list[Path]] | None = None
    if paired_r2:
        if not round_map:
            typer.secho(
                "extract: --paired-r2 requires --round-map for round assignment.",
                fg=typer.colors.RED,
                err=True,
            )
            raise typer.Exit(code=2)
        paired_r2_inputs = {}
        for r2_fq in paired_r2:
            r = rm.get(r2_fq.name)
            if r is None:
                typer.secho(
                    f"extract: R2 FASTQ {r2_fq.name!r} not in round map.",
                    fg=typer.colors.RED,
                    err=True,
                )
                raise typer.Exit(code=2)
            paired_r2_inputs.setdefault(r, []).append(r2_fq)

    # Capture argv into the manifest's `parameters` field.
    parameters: dict[str, str] = {
        "fastq": ",".join(str(p) for p in fastq),
        "library_report": str(library_report),
        "round_map": str(round_map) if round_map else "",
        "sample_sheet": str(sample_sheet) if sample_sheet else "",
        "paired_r2": ",".join(str(p) for p in (paired_r2 or [])),
        "override_primer_5p": override_primer_5p or "",
        "override_primer_3p": override_primer_3p or "",
        "rebuild": str(rebuild),
        "outdir": str(outdir),
    }

    result = run_extract(
        lr,
        fastq,
        outdir,
        round_map=rm,
        sample_sheet=sample_sheet,
        rebuild=rebuild,
        paired_r2_inputs=paired_r2_inputs,
        override_primer_5p=override_primer_5p,
        override_primer_3p=override_primer_3p,
        parameters=parameters,
    )

    if result.skipped:
        typer.secho(f"extract: skipped - {result.skipped_reason}", fg=typer.colors.YELLOW, err=True)
        raise typer.Exit(code=2)

    typer.echo(f"extract: wrote {len(result.outputs)} files under {outdir}")
    for p in result.outputs:
        typer.echo(f"  {p}")


def _parse_round_label(label: str) -> int:
    """Parse ``--round R0`` / ``--round 0`` / ``--round 00`` -> integer 0.

    Accepts an optional leading ``R``/``r``/``round_`` prefix. Raises
    ``typer.BadParameter`` on anything else.
    """
    s = label.strip()
    for prefix in ("round_", "Round_", "R", "r"):
        if s.startswith(prefix):
            s = s[len(prefix) :]
            break
    try:
        r = int(s)
    except ValueError as e:
        raise typer.BadParameter(
            f"--round expects a non-negative integer (optionally R-prefixed); got {label!r}"
        ) from e
    if r < 0:
        raise typer.BadParameter(f"--round must be >= 0; got {label!r} (parsed as {r})")
    return r


@app.command()
def count(
    extracted: Path = typer.Argument(..., help="Extracted FASTA(.gz) file (output of `extract`)."),
    round_label: str = typer.Option(..., "--round", help="Round label, e.g. R0 or 0."),
    outdir: Path = typer.Option(..., "--outdir", help="Output directory (round_NN/ subdir)."),
    from_pretrimmed_fastq: bool = typer.Option(
        False,
        "--from-pretrimmed-fastq",
        help=(
            "Opt-in: count a pre-trimmed FASTQ (primers already stripped, "
            "e.g. by AptaPLEX / EasyDIVER+ / external cutadapt). selexprep "
            "cannot verify the trimming state and will warn at invocation. "
            "Default behavior (no flag) hard-rejects FASTQ inputs to keep "
            "the v0.1 pipeline contract `extract -> count` unambiguous."
        ),
    ),
) -> None:
    """Count unique sequences from an extracted FASTA -> counts.parquet."""
    # Phase 5 Codex pass 1: hard-reject FASTQ by default (selexprep count
    # accepts only EXTRACTED FASTA from `selexprep extract`, primer-stripped
    # and random-region-only). FASTQ inputs would silently get parsed as if
    # every other line is a sequence header. Power users with externally
    # pre-trimmed FASTQ can opt in via --from-pretrimmed-fastq, which routes
    # to count_fastq_pretrimmed and surfaces a loud "cannot verify trimming"
    # warning. FASTQ counting otherwise stays a library API
    # (selexprep.count.counter.count_round, used in the thesis pipeline).
    name_lower = extracted.name.lower()
    is_fastq = name_lower.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))

    if is_fastq and not from_pretrimmed_fastq:
        typer.secho(
            "count: expected extracted FASTA(.gz) from `selexprep extract`, "
            "got FASTQ.\n"
            "  - If you want selexprep's primer-aware extraction, run "
            "`selexprep extract` first, then count the emitted "
            "extracted.fasta.gz / partial_*_extracted.fasta.gz file.\n"
            "  - If your FASTQ is ALREADY primer-stripped by another tool "
            "(AptaPLEX, EasyDIVER+, external cutadapt), pass "
            "--from-pretrimmed-fastq to opt in. selexprep will not verify "
            "the trimming state.",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=2)

    if from_pretrimmed_fastq and not is_fastq:
        typer.secho(
            f"count: --from-pretrimmed-fastq passed but input is not FASTQ. got: {extracted.name}",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=2)

    r = _parse_round_label(round_label)
    round_dir = outdir / f"round_{r:02d}"
    round_dir.mkdir(parents=True, exist_ok=True)
    out_parquet = round_dir / "counts.parquet"

    if from_pretrimmed_fastq:
        # Print AND log: the typer.secho is the user-facing safety footnote
        # (Typer apps don't configure stdlib logging by default, so a bare
        # logger.warning would be invisible to most CLI users). The
        # logger.warning still goes out for downstream tooling that captures
        # the selexprep.cli logger explicitly (e.g., qc/runner audit trails).
        typer.secho(
            "WARNING: --from-pretrimmed-fastq accepts the input as-is. "
            "selexprep cannot verify that primers/adapters have been "
            "stripped; if they remain in the input, unique-sequence "
            "counts will be inflated. For primer-aware extraction, run "
            "`selexprep extract` first and count the emitted FASTA.",
            fg=typer.colors.YELLOW,
            err=True,
        )
        logger.warning(
            "Counting pre-trimmed FASTQ %s via --from-pretrimmed-fastq. "
            "selexprep cannot verify that primers/adapters have been stripped; "
            "if they remain in the input, unique-sequence counts will be "
            "inflated. For selexprep's primer-aware extraction, use "
            "`selexprep extract` then count the emitted FASTA.",
            extracted,
        )
        stats = count_fastq_pretrimmed(extracted, out_parquet)
    else:
        stats = count_fasta(extracted, out_parquet)

    typer.echo(f"counts.parquet -> {out_parquet}")
    typer.echo(f"  unique sequences: {stats['n_unique']:,}")
    typer.echo(f"  total reads:      {stats['n_reads']:,}")
    typer.echo(
        f"  top sequence:     {stats['top_seq_reads']:,} reads ({stats['top_seq_rpm']:.0f} RPM)"
    )
    typer.echo(f"  Shannon entropy:  {stats['shannon_entropy_bits']:.2f} bits")


@app.command()
def qc(
    manifest: Path = typer.Argument(..., help="Path to selexprep_manifest.json."),
    counts_dir: Path | None = typer.Option(
        None,
        "--counts-dir",
        help="Directory containing round_*/counts.parquet (default: manifest's parent).",
    ),
    outdir: Path | None = typer.Option(
        None,
        "--outdir",
        help="Where to write flags.yaml + PNG plots (default: <manifest_parent>/qc/).",
    ),
) -> None:
    """Produce QC plots + depth-aware suspicion flags from a manifest."""
    result = run_qc(manifest, counts_dir=counts_dir, outdir=outdir)
    typer.echo(
        f"qc: {result.n_flags_raised} flag(s) raised; {len(result.plot_paths)} plot(s) written"
    )
    if result.flags_yaml_path is not None:
        typer.echo(f"  flags.yaml: {result.flags_yaml_path}")
    for p in result.plot_paths:
        typer.echo(f"  plot: {p}")
    for f in result.flags:
        typer.echo(f"  [{f.severity.upper()}] {f.name}")


@app.command()
def run(
    accessions: Path = typer.Argument(..., help="TSV of accessions to batch-process."),
    outdir: Path = typer.Option(..., "--outdir", help="Output directory."),
    resume: bool = typer.Option(False, "--resume", help="Resume an interrupted batch."),
    stop_on_error: bool = typer.Option(
        False,
        "--stop-on-error",
        help="Halt the batch on first per-accession failure instead of logging+continuing.",
    ),
    backend: str = typer.Option(
        "ena",
        "--backend",
        help=(
            "Download backend (default 'ena' for paper-grade reproducibility; "
            "use 'auto' for the convenience fallback chain)."
        ),
    ),
    allow_manual_review: bool = typer.Option(
        False,
        "--allow-manual-review",
        help=(
            "NONE-confidence runs are downloaded to round_unknown/ and "
            "surfaced in manual_review.tsv; they NEVER enter rounds.tsv."
        ),
    ),
    timeout_s: int = typer.Option(30, "--timeout-s", help="HTTP timeout (seconds)."),
) -> None:
    """Batch-process a list of accessions; emits per-dataset + corpus-level outputs."""
    from selexprep.run import run_batch

    allowed_backends = {"auto", "ena", "kingfisher", "sra"}
    if backend not in allowed_backends:
        typer.secho(
            f"run: invalid --backend {backend!r}; allowed: {sorted(allowed_backends)}",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=2)

    try:
        report = run_batch(
            accessions,
            outdir,
            resume=resume,
            stop_on_error=stop_on_error,
            backend=backend,  # type: ignore[arg-type]
            allow_manual_review=allow_manual_review,
            timeout_s=timeout_s,
        )
    except ValueError as e:
        typer.secho(f"run: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=2) from e

    n_ok = sum(1 for r in report.rows if r.status == "OK")
    n_skipped = sum(1 for r in report.rows if r.status == "SKIPPED_READ_MERGING_RECOMMENDED")
    n_failed = sum(
        1 for r in report.rows if r.status not in ("OK", "SKIPPED_READ_MERGING_RECOMMENDED")
    )
    typer.echo(f"Accessions processed: {len(report.rows)}")
    typer.echo(f"  OK:                 {n_ok}")
    typer.echo(f"  skipped (merging):  {n_skipped}")
    typer.echo(f"  failed:             {n_failed}")
    if report.summary_tsv is not None:
        typer.echo(f"run_summary.tsv -> {report.summary_tsv}")

    for row in report.rows:
        if row.status not in ("OK", "SKIPPED_READ_MERGING_RECOMMENDED"):
            typer.secho(
                f"  [{row.status}] {row.accession}: {row.notes}",
                fg=typer.colors.YELLOW,
                err=True,
            )

    # NOTE (Phase 6b.4 HPC audit fix): ``selexprep run`` is a batch driver —
    # per-accession failures are first-class data captured in
    # ``run_summary.tsv``, which is the report. A non-zero exit here would
    # conflate "the runner did its job and recorded failures" (a normal
    # operational outcome for a noisy public corpus) with "the runner
    # itself crashed" (which is already handled separately by the outer
    # ``except ValueError`` → exit 2). The audit Snakefile + any
    # downstream automation (CI, monitoring, the Tier 2 ``rule run_corpus``
    # under ``set -e``) depends on a clean exit when the summary was
    # written. Users who want fail-fast can pass ``--stop-on-error``,
    # which also writes the summary before halting.
    #
    # An empty input TSV (zero accessions) DOES still need to fail loudly
    # — that's an operator error, not "the audited corpus is messy".
    if not report.rows:
        typer.secho(
            "run: accessions TSV produced zero rows after parsing — refusing "
            "to emit an empty summary.",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=2)


if __name__ == "__main__":
    app()
