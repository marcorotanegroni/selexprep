"""Unit tests for selexprep.catalog (discovery catalog + filters + CLI + rebuild)."""

from __future__ import annotations

import csv
from pathlib import Path
from unittest import mock

import pandas as pd
import pytest
from typer.testing import CliRunner

from selexprep.catalog import catalog_version, filter_catalog, load_catalog
from selexprep.catalog.cli import app as catalog_app
from selexprep.catalog.rebuild import (
    PUBLIC_COLS,
    _enrichment_index,
    _passthrough_non_insdc,
    rebuild_catalog,
)
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


def test_cli_catalog_list_marks_fetchability() -> None:
    """`list` surfaces a fetchable column; non-INSDC rows read 'discovery-only'."""
    result = runner.invoke(catalog_app, ["list", "--source", "zenodo", "--limit", "5"])
    assert result.exit_code == 0
    assert "fetchable" in result.stdout
    assert "discovery-only" in result.stdout


def test_cli_catalog_show_marks_discovery_only() -> None:
    """`show` on a non-INSDC pointer flags it as discovery-only."""
    df = load_catalog()
    disc = df[df["bioproject_id"].str.startswith(("zenodo:", "figshare:"))]
    assert not disc.empty, "expected at least one discovery-only catalog row"
    result = runner.invoke(catalog_app, ["show", disc.iloc[0]["bioproject_id"]])
    assert result.exit_code == 0
    assert "discovery-only" in result.stdout


def test_inspect_accession_refuses_discovery_only() -> None:
    """Non-INSDC discovery pointers are refused before any network call."""
    from selexprep.fetch.inspect import inspect_accession

    with pytest.raises(ValueError, match="discovery-only"):
        inspect_accession("zenodo:99999")


def test_run_fetch_refuses_discovery_only(tmp_path: Path) -> None:
    from selexprep.fetch import run_fetch

    with pytest.raises(ValueError, match="discovery-only"):
        run_fetch("zenodo:99999", tmp_path)


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
    assert "refresh" in result.stdout


# ----- Rebuild (uses fake ENA results — no network) -----


def _write_catalog_csv(path: Path, rows: list[dict]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(PUBLIC_COLS))
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in PUBLIC_COLS})


def test_enrichment_index_picks_rows_with_any_hand_field(tmp_path: Path) -> None:
    p = tmp_path / "cat.csv"
    _write_catalog_csv(
        p,
        [
            {"bioproject_id": "PRJ1", "protein_target": "VEGF"},  # enriched
            {"bioproject_id": "PRJ2", "paper_doi": "10.1/x"},  # enriched
            {"bioproject_id": "PRJ3", "n_rounds_declared": "8"},  # enriched
            {"bioproject_id": "PRJ4"},  # bare
        ],
    )
    idx = _enrichment_index(p)
    assert set(idx) == {"PRJ1", "PRJ2", "PRJ3"}


def test_passthrough_non_insdc_carries_only_zenodo_figshare_utexas(tmp_path: Path) -> None:
    p = tmp_path / "cat.csv"
    _write_catalog_csv(
        p,
        [
            {"bioproject_id": "PRJ1", "source": "ena"},
            {"bioproject_id": "zenodo:123", "source": "zenodo:123", "study_title": "Z"},
            {"bioproject_id": "figshare:456", "source": "figshare:456", "study_title": "F"},
            {"bioproject_id": "utexas:doi:x", "source": "utexas_db", "study_title": "U"},
        ],
    )
    carried = _passthrough_non_insdc(p, exclude_ids=set())
    bps = {r["bioproject_id"] for r in carried}
    assert bps == {"zenodo:123", "figshare:456", "utexas:doi:x"}


def test_passthrough_excludes_ids_already_seen(tmp_path: Path) -> None:
    p = tmp_path / "cat.csv"
    _write_catalog_csv(
        p,
        [
            {"bioproject_id": "zenodo:123", "source": "zenodo:123", "study_title": "Z"},
            {"bioproject_id": "zenodo:456", "source": "zenodo:456", "study_title": "ZZ"},
        ],
    )
    carried = _passthrough_non_insdc(p, exclude_ids={"zenodo:123"})
    bps = {r["bioproject_id"] for r in carried}
    assert bps == {"zenodo:456"}


def test_rebuild_catalog_merges_enrichment_and_preserves_non_insdc(tmp_path: Path) -> None:
    """End-to-end: rebuild_catalog with mocked ENA harvest merges hand-enriched
    fields forward and carries non-INSDC deposits across."""
    old_path = tmp_path / "old.csv"
    _write_catalog_csv(
        old_path,
        [
            # Seed entry with full enrichment — must be preserved when ENA
            # re-surfaces this accession.
            {
                "bioproject_id": "PRJTEST1",
                "source": "seed",
                "study_title": "stale title from old catalog",
                "protein_target": "IL-10RA",
                "target_organism": "Homo sapiens",
                "paper_doi": "10.1/x",
                "paper_pmid": "12345",
                "n_rounds_declared": "5",
                "abstract": "stale abstract",
            },
            # Non-INSDC deposit — must be carried forward
            {
                "bioproject_id": "zenodo:99",
                "source": "zenodo:99",
                "study_title": "Z processed data",
            },
        ],
    )

    fake_ena = {
        "PRJTEST1": {
            "study_accession": "PRJTEST1",
            "study_title": "FRESH ENA TITLE",
            "study_description": "Fresh abstract from ENA",
            "scientific_name": "Homo sapiens",
            "library_strategies": ["OTHER"],
        },
        "PRJTEST2": {
            "study_accession": "PRJTEST2",
            "study_title": "New aptamer study",
            "study_description": "Discovered via broader query",
            "scientific_name": "Mus musculus",
            "library_strategies": ["OTHER", "OTHER"],
        },
    }

    new_path = tmp_path / "new.csv"
    with mock.patch(
        "selexprep.catalog.rebuild.harvest_runs_from_ena",
        return_value=fake_ena,
    ):
        n = rebuild_catalog(out_path=new_path, preserve_from=old_path)

    assert n == 3  # 2 ENA + 1 zenodo

    with open(new_path) as f:
        rebuilt = {r["bioproject_id"]: r for r in csv.DictReader(f)}

    # PRJTEST1: enrichment preserved (target/doi/pmid/rounds carried)
    # but the title + abstract come from fresh ENA
    assert rebuilt["PRJTEST1"]["protein_target"] == "IL-10RA"
    assert rebuilt["PRJTEST1"]["paper_doi"] == "10.1/x"
    assert rebuilt["PRJTEST1"]["paper_pmid"] == "12345"
    assert rebuilt["PRJTEST1"]["n_rounds_declared"] == "5"
    assert rebuilt["PRJTEST1"]["study_title"] == "FRESH ENA TITLE"
    assert rebuilt["PRJTEST1"]["abstract"] == "Fresh abstract from ENA"

    # PRJTEST2: brand new from ENA, no enrichment available
    assert rebuilt["PRJTEST2"]["protein_target"] == ""
    assert rebuilt["PRJTEST2"]["study_title"] == "New aptamer study"

    # zenodo:99: passed through from old catalog
    assert rebuilt["zenodo:99"]["study_title"] == "Z processed data"


def test_rebuild_catalog_no_preserve_drops_enrichment_and_non_insdc(tmp_path: Path) -> None:
    """Without `preserve_from`, the rebuild is a clean slate: ENA-only."""
    fake_ena = {
        "PRJTEST1": {
            "study_accession": "PRJTEST1",
            "study_title": "Fresh title",
            "study_description": "Fresh abstract",
            "library_strategies": ["OTHER"],
        },
    }
    new_path = tmp_path / "new.csv"
    with mock.patch(
        "selexprep.catalog.rebuild.harvest_runs_from_ena",
        return_value=fake_ena,
    ):
        n = rebuild_catalog(out_path=new_path, preserve_from=None)

    assert n == 1
    with open(new_path) as f:
        rebuilt = list(csv.DictReader(f))
    assert rebuilt[0]["bioproject_id"] == "PRJTEST1"
    assert rebuilt[0]["protein_target"] == ""


def test_rebuild_catalog_excludes_all_blocklisted_studies(tmp_path: Path) -> None:
    """a study whose runs are 100% RNA-Seq must be dropped
    from the catalog AND recorded in ``bioprojects_excluded.csv`` with
    a human-readable reason. Mixed studies (some OTHER + some RNA-Seq)
    must be KEPT — they contain real SELEX data alongside controls.
    """
    fake_ena = {
        "PRJ_REAL_SELEX": {
            "study_accession": "PRJ_REAL_SELEX",
            "study_title": "Real HT-SELEX study",
            "study_description": "...",
            "scientific_name": "synthetic construct",
            "library_strategies": ["OTHER", "OTHER", "OTHER"],
        },
        "PRJ_RNASEQ_FALSE_POSITIVE": {
            "study_accession": "PRJ_RNASEQ_FALSE_POSITIVE",
            "study_title": "Bladder cancer transcriptome (mentions aptamer)",
            "study_description": "...",
            "scientific_name": "Homo sapiens",
            "library_strategies": ["RNA-Seq"] * 5,
        },
        "PRJ_MIXED_SELEX_AND_CONTROLS": {
            "study_accession": "PRJ_MIXED_SELEX_AND_CONTROLS",
            "study_title": "SELEX with RNA-seq controls",
            "study_description": "...",
            "scientific_name": "Mus musculus",
            "library_strategies": ["OTHER", "OTHER", "RNA-Seq"],
        },
    }

    new_path = tmp_path / "new.csv"
    with mock.patch(
        "selexprep.catalog.rebuild.harvest_runs_from_ena",
        return_value=fake_ena,
    ):
        n = rebuild_catalog(out_path=new_path, preserve_from=None)

    # 2 kept (real SELEX + mixed), 1 excluded (RNA-Seq-only)
    assert n == 2
    with open(new_path) as f:
        rebuilt = {r["bioproject_id"] for r in csv.DictReader(f)}
    assert rebuilt == {"PRJ_REAL_SELEX", "PRJ_MIXED_SELEX_AND_CONTROLS"}
    assert "PRJ_RNASEQ_FALSE_POSITIVE" not in rebuilt

    # Sidecar present + populated
    excluded_path = tmp_path / "bioprojects_excluded.csv"
    assert excluded_path.exists()
    with open(excluded_path) as f:
        excluded_rows = list(csv.DictReader(f))
    assert len(excluded_rows) == 1
    row = excluded_rows[0]
    assert row["bioproject_id"] == "PRJ_RNASEQ_FALSE_POSITIVE"
    assert row["n_runs_total"] == "5"
    assert row["n_runs_blocklisted"] == "5"
    assert "RNA-Seq" in row["exclusion_reason"]
    assert "5 runs" in row["exclusion_reason"]
    # blocklisted_strategies is a JSON-encoded dict for downstream parsing
    assert "RNA-Seq" in row["blocklisted_strategies"]


def test_rebuild_catalog_emits_empty_exclusion_sidecar(tmp_path: Path) -> None:
    """When every study is compatible, the sidecar still exists (header-only).
    Downstream consumers can rely on its presence."""
    fake_ena = {
        "PRJ_OK": {
            "study_accession": "PRJ_OK",
            "study_title": "...",
            "study_description": "...",
            "scientific_name": "",
            "library_strategies": ["OTHER"],
        },
    }
    new_path = tmp_path / "new.csv"
    with mock.patch(
        "selexprep.catalog.rebuild.harvest_runs_from_ena",
        return_value=fake_ena,
    ):
        rebuild_catalog(out_path=new_path, preserve_from=None)

    excluded_path = tmp_path / "bioprojects_excluded.csv"
    assert excluded_path.exists()
    with open(excluded_path) as f:
        rows = list(csv.DictReader(f))
    assert rows == []  # header-only
