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
from selexprep.library import compute_library_report, write_library_report_json

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
