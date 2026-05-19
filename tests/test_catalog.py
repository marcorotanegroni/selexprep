"""Unit tests for selexprep.catalog (discovery catalog + filters + CLI)."""

from __future__ import annotations

import pandas as pd
from typer.testing import CliRunner

from selexprep.catalog import catalog_version, filter_catalog, load_catalog
from selexprep.catalog.cli import app as catalog_app
from selexprep.cli import app as root_app

runner = CliRunner()


# ----- reader / shape -----


def test_load_catalog_returns_dataframe_with_expected_columns() -> None:
    df = load_catalog()
    assert isinstance(df, pd.DataFrame)
    expected = {
        "bioproject_id",
        "source",
        "study_title",
        "protein_target",
        "target_organism",
        "paper_doi",
        "paper_pmid",
        "n_rounds_declared",
        "abstract",
    }
    assert set(df.columns) == expected


def test_load_catalog_drops_thesis_specific_columns() -> None:
    """Thesis-only fields must NOT be in the public v0.1 catalog."""
    df = load_catalog()
    forbidden = {
        "include",
        "manual_curation_notes",
        "library_type_verification",
        "library_type_evidence",
        "has_processed_counts",
    }
    assert not (set(df.columns) & forbidden), "thesis-specific columns leaked into public catalog"


def test_load_catalog_has_seed_anchor_entries() -> None:
    """At minimum the 4 seed bioprojects (HT-SELEX benchmarks) must be present
    so users can sanity-check the catalog ships correctly."""
    df = load_catalog()
    ids = set(df["bioproject_id"])
    assert "PRJNA315881" in ids  # Hoinka et al. IL-10RA — AptaSUITE benchmark
    assert "PRJNA321551" in ids  # Dao et al. CCR7 — AptaTRACE benchmark


def test_catalog_version_is_a_nonempty_string() -> None:
    assert isinstance(catalog_version(), str)
    assert catalog_version()  # non-empty


# ----- filter -----


def test_filter_catalog_target_substring_match() -> None:
    df = load_catalog()
    out = filter_catalog(df, target="IL-10")
    assert len(out) >= 1
    assert all("il-10" in t.lower() for t in out["protein_target"] if t)


def test_filter_catalog_target_case_insensitive() -> None:
    df = load_catalog()
    upper = filter_catalog(df, target="CCR7")
    lower = filter_catalog(df, target="ccr7")
    assert len(upper) == len(lower)


def test_filter_catalog_insdc_only_drops_processed_data() -> None:
    df = load_catalog()
    out = filter_catalog(df, insdc_only=True)
    # No zenodo/figshare/utexas accessions remain
    assert not out["bioproject_id"].str.startswith("zenodo:").any()
    assert not out["bioproject_id"].str.startswith("figshare:").any()
    assert not out["bioproject_id"].str.startswith("utexas:").any()
    # And every retained ID is INSDC-formatted
    assert all(bp.startswith(("PRJ", "SRP", "ERP", "DRP")) for bp in out["bioproject_id"])


def test_filter_catalog_min_rounds() -> None:
    df = load_catalog()
    out = filter_catalog(df, min_rounds=5)
    # Every retained row has n_rounds_declared >= 5 (or, since we dropped any
    # row whose n_rounds wasn't a valid int, none of the survivors fall below)
    nrounds = pd.to_numeric(out["n_rounds_declared"], errors="coerce")
    assert (nrounds.fillna(0) >= 5).all()


def test_filter_catalog_source() -> None:
    df = load_catalog()
    out = filter_catalog(df, source_contains="ena")
    assert all("ena" in src.lower() for src in out["source"] if src)


def test_filter_catalog_filters_compose() -> None:
    df = load_catalog()
    chained = filter_catalog(df, insdc_only=True, min_rounds=3)
    indep = filter_catalog(filter_catalog(df, insdc_only=True), min_rounds=3)
    assert len(chained) == len(indep)


def test_filter_catalog_no_filters_returns_full_catalog() -> None:
    df = load_catalog()
    out = filter_catalog(df)
    assert len(out) == len(df)


# ----- CLI smoke tests (subapp + root app integration) -----


def test_cli_catalog_list_runs_clean() -> None:
    result = runner.invoke(catalog_app, ["list", "--limit", "5"])
    assert result.exit_code == 0
    assert "bioproject_id" in result.stdout
    assert "Catalog:" in result.stdout


def test_cli_catalog_list_filters_by_target() -> None:
    result = runner.invoke(catalog_app, ["list", "--target", "IL-10", "--limit", "20"])
    assert result.exit_code == 0
    # IL-10RA should appear; row count should be small
    assert "IL-10" in result.stdout or "il-10" in result.stdout.lower()


def test_cli_catalog_show_known_accession() -> None:
    result = runner.invoke(catalog_app, ["show", "PRJNA315881"])
    assert result.exit_code == 0
    assert "PRJNA315881" in result.stdout
    assert "IL-10" in result.stdout  # in protein_target / title / abstract


def test_cli_catalog_show_unknown_accession_exits_nonzero() -> None:
    result = runner.invoke(catalog_app, ["show", "PRJNOTREAL"])
    assert result.exit_code != 0


def test_cli_catalog_version() -> None:
    result = runner.invoke(catalog_app, ["version"])
    assert result.exit_code == 0
    assert result.stdout.strip()  # non-empty


def test_cli_root_app_exposes_catalog_subapp() -> None:
    """`selexprep catalog --help` must work from the root dispatcher."""
    result = runner.invoke(root_app, ["catalog", "--help"])
    assert result.exit_code == 0
    assert "list" in result.stdout
    assert "show" in result.stdout
