"""Typer subapp for the ``selexprep catalog`` CLI verb."""

from __future__ import annotations

from pathlib import Path

import typer

from selexprep.catalog.filter import filter_catalog
from selexprep.catalog.reader import catalog_version, load_catalog

app = typer.Typer(
    name="catalog",
    help="Browse the bundled public-SELEX discovery catalog.",
    no_args_is_help=True,
)


_DEFAULT_LIST_COLS = (
    "bioproject_id",
    "protein_target",
    "target_organism",
    "n_rounds_declared",
    "source",
)


@app.command("list")
def list_catalog(
    target: str = typer.Option(
        None, "--target", help="Filter by protein target (case-insensitive substring)."
    ),
    organism: str = typer.Option(None, "--organism", help="Filter by target organism (substring)."),
    source: str = typer.Option(
        None, "--source", help="Filter by discovery source (substring, e.g. 'ena')."
    ),
    min_rounds: int = typer.Option(
        None, "--min-rounds", help="Drop bioprojects with fewer than N declared rounds."
    ),
    insdc_only: bool = typer.Option(
        False,
        "--insdc-only",
        help=(
            "Keep only PRJ*/SRP/ERP/DRP accessions (drops Zenodo / Figshare "
            "processed-data deposits)."
        ),
    ),
    limit: int = typer.Option(50, "--limit", help="Maximum rows to print."),
) -> None:
    """List public-SELEX bioprojects in the bundled catalog."""
    df = load_catalog()
    filtered = filter_catalog(
        df,
        target=target,
        organism=organism,
        source_contains=source,
        min_rounds=min_rounds,
        insdc_only=insdc_only,
    )
    n_total = len(filtered)
    head = filtered.head(limit)

    cols = list(_DEFAULT_LIST_COLS)
    # Truncate long titles in the listing view
    if not head.empty:
        typer.echo(head[cols].to_string(index=False))
    else:
        typer.echo("(no matches)")
    typer.echo("")
    suffix = "" if n_total <= limit else f"  (showing first {limit}; use --limit)"
    typer.echo(f"Matched: {n_total}{suffix}    Catalog: {catalog_version()}")


@app.command("show")
def show_bp(
    accession: str = typer.Argument(..., help="Bioproject accession (e.g. PRJNA315881)."),
) -> None:
    """Show full detail for one bioproject in the catalog."""
    df = load_catalog()
    rows = df[df["bioproject_id"] == accession]
    if rows.empty:
        typer.secho(
            f"Accession {accession!r} not found in catalog ({catalog_version()})",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=1)
    row = rows.iloc[0]
    for col in df.columns:
        val = str(row[col])
        if not val:
            continue
        typer.echo(f"{col}: {val}")


@app.command("version")
def version() -> None:
    """Print the catalog snapshot identifier."""
    typer.echo(catalog_version())


@app.command("refresh")
def refresh(
    out: Path = typer.Option(
        None,
        "--out",
        help="Output CSV path. Defaults to the bundled package data file (in-place refresh).",
    ),
    preserve_enrichment: bool = typer.Option(
        True,
        "--preserve-enrichment/--no-preserve-enrichment",
        help=(
            "Carry hand-enriched fields (protein_target, paper_doi, "
            "paper_pmid, n_rounds_declared) forward from the existing "
            "catalog when an accession is still found upstream."
        ),
    ),
) -> None:
    """Refresh the catalog by re-running broad ENA queries.

    Hits ENA's public Portal API with ~13 broad SELEX/aptamer keyword
    queries (see `selexprep.catalog.rebuild.ENA_QUERIES`), unions the
    resulting INSDC studies, merges hand-enriched fields forward from the
    current catalog (when `--preserve-enrichment` is on), and writes the
    result to `--out` (or the bundled package data file by default).

    Takes ~30-60 seconds depending on ENA latency. No NCBI / Zenodo /
    Figshare calls here — ENA mirrors INSDC, so SRA + DDBJ deposits are
    already covered.
    """
    from selexprep.catalog.reader import catalog_path
    from selexprep.catalog.rebuild import rebuild_catalog

    current = catalog_path()
    out_path = out if out is not None else current
    preserve = current if preserve_enrichment else None

    typer.echo(f"Refreshing catalog → {out_path}")
    if preserve:
        typer.echo(f"  preserving enrichment from {preserve}")

    n = rebuild_catalog(out_path=out_path, preserve_from=preserve)
    typer.echo(f"Wrote {n} bioprojects.")
