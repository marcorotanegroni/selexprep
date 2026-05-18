"""Typer CLI dispatcher for selexprep.

Phase 0 scaffold: all subcommands are stubs that exit with code 2 (not yet
implemented). They exist so `--help` lists the full command surface and so the
package can be installed and tested end-to-end on TestPyPI before Phase 1
ports the real logic from `selex_corpus/`.
"""

from __future__ import annotations

from pathlib import Path

import typer

from selexprep import __version__

app = typer.Typer(
    name="selexprep",
    help="Accession-first preprocessing for public HT-SELEX with primer auto-inference.",
    no_args_is_help=True,
    add_completion=False,
)


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
) -> None:
    """Preview accession metadata (round count, library_strategy, files) without downloading."""
    _not_implemented("inspect")


@app.command()
def fetch(
    accession: str = typer.Argument(..., help="ENA / SRA / DDBJ accession."),
    outdir: Path = typer.Option(..., "--outdir", help="Output directory."),
) -> None:
    """Download FASTQ + metadata for an accession; round map auto-populated."""
    _not_implemented("fetch")


@app.command()
def detect(
    fastq: list[Path] = typer.Argument(..., help="Input FASTQ files."),
    round_map: Path | None = typer.Option(
        None, "--round-map", help="TSV mapping FASTQ file → round number (required for local FASTQs)."
    ),
    outdir: Path = typer.Option(..., "--outdir", help="Output directory."),
) -> None:
    """Auto-infer primers + library structure from FASTQ; emit LibraryReport JSON."""
    _not_implemented("detect")


@app.command()
def extract(
    fastq: list[Path] = typer.Argument(..., help="Input FASTQ files."),
    library_report: Path = typer.Option(
        ..., "--library-report", help="Path to LibraryReport JSON from `detect`."
    ),
    sample_sheet: Path | None = typer.Option(
        None, "--sample-sheet", help="Optional sample sheet for demultiplexing."
    ),
    override_primer_5p: str | None = typer.Option(None, "--override-primer-5p"),
    override_primer_3p: str | None = typer.Option(None, "--override-primer-3p"),
    rebuild: bool = typer.Option(False, "--rebuild", help="Force re-extraction; emit diff report."),
    outdir: Path = typer.Option(..., "--outdir", help="Output directory."),
) -> None:
    """Trim primers + extract random region per round."""
    _not_implemented("extract")


@app.command()
def count(
    extracted: Path = typer.Argument(..., help="Extracted FASTA file."),
    round: str = typer.Option(..., "--round", help="Round label (e.g. R1)."),
    outdir: Path = typer.Option(..., "--outdir", help="Output directory."),
) -> None:
    """Count unique sequences per round (raw + RPM/frequency)."""
    _not_implemented("count")


@app.command()
def qc(
    manifest: Path = typer.Argument(..., help="Path to selexprep manifest JSON."),
) -> None:
    """Produce QC plots + depth-aware suspicion flags from a manifest."""
    _not_implemented("qc")


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
