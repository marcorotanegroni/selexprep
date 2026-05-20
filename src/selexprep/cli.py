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
from selexprep.count.counter import count_fasta
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
) -> None:
    """Download FASTQ + metadata for an accession; round map auto-populated."""
    _not_implemented("fetch")


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

    # Phase 2 CLI is single-end only; paired-end (R2 stream) arrives in Phase 3.
    report = compute_library_report(
        sequences_by_round,
        read_source="R1",
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
        return int(s)
    except ValueError as e:
        raise typer.BadParameter(
            f"--round expects an integer (optionally R-prefixed); got {label!r}"
        ) from e


@app.command()
def count(
    extracted: Path = typer.Argument(..., help="Extracted FASTA(.gz) file (output of `extract`)."),
    round_label: str = typer.Option(..., "--round", help="Round label, e.g. R0 or 0."),
    outdir: Path = typer.Option(..., "--outdir", help="Output directory (round_NN/ subdir)."),
) -> None:
    """Count unique sequences from an extracted FASTA -> counts.parquet."""
    r = _parse_round_label(round_label)
    round_dir = outdir / f"round_{r:02d}"
    round_dir.mkdir(parents=True, exist_ok=True)
    out_parquet = round_dir / "counts.parquet"

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
) -> None:
    """Batch-process a list of accessions; emits per-dataset + corpus-level outputs."""
    _not_implemented("run")


if __name__ == "__main__":
    app()
