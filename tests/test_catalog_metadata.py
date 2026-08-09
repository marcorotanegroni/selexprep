"""Unit tests for the curated (annotated) metadata layer."""

from __future__ import annotations

import pandas as pd

from selexprep.catalog import load_metadata, load_metadata_records, metadata_version
from selexprep.catalog.metadata import METADATA_FIELDS

_N_DEPOSITS = 240


def test_metadata_version_is_nonempty_str() -> None:
    v = metadata_version()
    assert isinstance(v, str) and v


def test_load_metadata_flat_shape() -> None:
    df = load_metadata()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == _N_DEPOSITS
    assert "bioproject_id" in df.columns
    for f in METADATA_FIELDS:
        assert f in df.columns, f
        assert f"{f}_curation" in df.columns, f


def test_load_metadata_records_shape_and_provenance() -> None:
    recs = load_metadata_records()
    assert isinstance(recs, list) and len(recs) == _N_DEPOSITS
    r = recs[0]
    assert "accession" in r and "fields" in r
    for f in METADATA_FIELDS:
        assert f in r["fields"]
        assert "status" in r["fields"][f]
    # at least some concordant cells carry a citable source
    assert any(
        fd.get("status") == "concordant" and fd.get("source")
        for rec in recs
        for fd in rec["fields"].values()
    )


def test_cells_where_arms_disagreed_keep_both_arms() -> None:
    """Adjudication resolves a disagreement; it never erases it.

    Cells the two extractions disagreed on are now ``adjudicated`` rather than
    left ``discordant``, but both arms stay on record so any reader can audit
    the call that was made.
    """
    recs = load_metadata_records()
    disagreed = [
        fd
        for rec in recs
        for fd in rec["fields"].values()
        if fd.get("status") in ("discordant", "adjudicated")
    ]
    assert disagreed, "expected cells where the two arms disagreed"
    for fd in disagreed:
        assert "claude" in fd and "codex" in fd
        assert "value" in fd["claude"] and "value" in fd["codex"]


def test_adjudicated_cells_carry_value_and_rule() -> None:
    """An adjudicated cell must state what was chosen and on what grounds."""
    recs = load_metadata_records()
    adjudicated = [
        fd for rec in recs for fd in rec["fields"].values() if fd.get("status") == "adjudicated"
    ]
    assert adjudicated, "expected adjudicated cells"
    for fd in adjudicated:
        assert fd.get("value"), "adjudicated cell must carry the resolved value"
        assert fd.get("adjudication_rule"), "adjudicated cell must name the rule applied"
        assert fd.get("adjudication_note"), "adjudicated cell must record the reasoning"


def test_flat_and_records_agree_on_accessions() -> None:
    df = load_metadata()
    recs = load_metadata_records()
    assert set(df["bioproject_id"]) == {r["accession"] for r in recs}
