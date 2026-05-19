"""Unit tests for selexprep.fetch.discover (offline only — network adapters
are tested separately under @pytest.mark.network)."""

from __future__ import annotations

import csv
from pathlib import Path

import yaml

from selexprep.fetch.discover import (
    BIOPROJECT_COLS,
    NOT_ASSESSED_V0_1,
    ROUND_COLS,
    SAMPLE_COLS,
    SeedAdapter,
    _classify_all,
    deduplicate_bioprojects,
    deduplicate_samples,
    empty_bioproject,
    empty_sample,
    is_blacklisted,
    load_blacklist,
    load_seed_overrides,
    load_small_molecule_targets,
    sm_mentioned_in_selex_context,
    write_csv,
)

# ----- Schema helpers -----


def test_empty_bioproject_has_all_columns() -> None:
    bp = empty_bioproject()
    assert set(bp.keys()) == set(BIOPROJECT_COLS)
    assert all(v == "" for v in bp.values())


def test_empty_sample_has_all_columns() -> None:
    s = empty_sample()
    assert set(s.keys()) == set(SAMPLE_COLS)
    assert all(v == "" for v in s.values())


# ----- SeedAdapter -----


def _write_seed(
    path: Path,
    entries: list[dict],
    blacklist: list[dict] | None = None,
    sm_targets: list[str] | None = None,
) -> None:
    payload = {"entries": entries}
    if blacklist is not None:
        payload["blacklist"] = blacklist
    if sm_targets is not None:
        payload["small_molecule_targets"] = sm_targets
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")


def test_seed_adapter_loads_entries(tmp_path: Path) -> None:
    seed = tmp_path / "seed.yaml"
    _write_seed(
        seed,
        [
            {
                "bioproject_id": "PRJ1",
                "protein_target": "VEGF",
                "library_type": "RNA_confirmed",
                "paper_doi": "10.1/x",
                "n_rounds_expected": 8,
                "notes": "test",
            },
            {"bioproject_id": "PRJ2", "protein_target": "TNF", "library_type": "DNA_confirmed"},
        ],
    )
    bp_rows, sample_rows = SeedAdapter(seed).search()
    assert len(bp_rows) == 2
    assert sample_rows == []
    rna = next(r for r in bp_rows if r["bioproject_id"] == "PRJ1")
    assert rna["include"] == "y"
    assert rna["protein_target"] == "VEGF"
    assert rna["n_rounds_declared"] == "8"
    dna = next(r for r in bp_rows if r["bioproject_id"] == "PRJ2")
    assert dna["include"] == "n"


def test_seed_adapter_filters_by_query(tmp_path: Path) -> None:
    seed = tmp_path / "seed.yaml"
    _write_seed(
        seed,
        [
            {"bioproject_id": "PRJ1", "library_type": "RNA_confirmed"},
            {"bioproject_id": "PRJ2", "library_type": "RNA_confirmed"},
        ],
    )
    bp_rows, _ = SeedAdapter(seed).search(query="PRJ1")
    assert len(bp_rows) == 1
    assert bp_rows[0]["bioproject_id"] == "PRJ1"


def test_seed_adapter_blacklist_accessor(tmp_path: Path) -> None:
    seed = tmp_path / "seed.yaml"
    _write_seed(
        seed,
        entries=[],
        blacklist=[{"bioproject_id": "PRJBAD1"}, {"bioproject_id": "PRJBAD2"}],
    )
    assert SeedAdapter(seed).blacklist == {"PRJBAD1", "PRJBAD2"}


# ----- Blacklist / SM helpers -----


def test_load_blacklist(tmp_path: Path) -> None:
    seed = tmp_path / "seed.yaml"
    _write_seed(seed, [], blacklist=[{"bioproject_id": "X"}])
    assert load_blacklist(seed) == {"X"}


def test_load_small_molecule_targets_lowercased(tmp_path: Path) -> None:
    seed = tmp_path / "seed.yaml"
    _write_seed(seed, [], sm_targets=["ATP", "GTP", "Theophylline"])
    assert load_small_molecule_targets(seed) == ["atp", "gtp", "theophylline"]


def test_is_blacklisted_explicit() -> None:
    bp = empty_bioproject()
    bp["bioproject_id"] = "PRJBAD"
    assert is_blacklisted(bp, {"PRJBAD"}, [])


def test_is_blacklisted_small_molecule_target() -> None:
    bp = empty_bioproject()
    bp["bioproject_id"] = "PRJ1"
    bp["protein_target"] = "ATP riboswitch"
    assert is_blacklisted(bp, set(), ["atp"])


def test_is_blacklisted_whole_cell_no_target() -> None:
    bp = empty_bioproject()
    bp["bioproject_id"] = "PRJ1"
    bp["abstract"] = "We used whole cell SELEX against pancreatic cells..."
    bp["protein_target"] = ""
    assert is_blacklisted(bp, set(), [])


def test_is_blacklisted_passes_clean_bp() -> None:
    bp = empty_bioproject()
    bp["bioproject_id"] = "PRJ1"
    bp["protein_target"] = "VEGF"
    bp["abstract"] = "Selection of RNA aptamers against VEGF protein."
    assert not is_blacklisted(bp, set(), ["atp", "gtp"])


# ----- SM context check -----


def test_sm_context_matches_with_selex_keyword() -> None:
    abstract = "We performed SELEX against theophylline binding RNAs."
    assert sm_mentioned_in_selex_context(abstract, ["theophylline"])


def test_sm_context_misses_when_no_selex_keyword() -> None:
    abstract = "ATP-binding cassette transporters are membrane proteins."
    assert not sm_mentioned_in_selex_context(abstract, ["atp"])


def test_sm_context_empty_abstract() -> None:
    assert not sm_mentioned_in_selex_context("", ["atp"])


# ----- Seed overrides -----


def test_load_seed_overrides(tmp_path: Path) -> None:
    seed = tmp_path / "seed.yaml"
    _write_seed(
        seed,
        [
            {"bioproject_id": "PRJ1", "manual_round_mapping": {"SRR1": 0, "SRR2": 1}},
            {"bioproject_id": "PRJ2"},  # no mapping
        ],
    )
    overrides = load_seed_overrides(seed)
    assert overrides == {"PRJ1": {"SRR1": 0, "SRR2": 1}}


def test_load_seed_overrides_missing_file_returns_empty(tmp_path: Path) -> None:
    assert load_seed_overrides(tmp_path / "nope.yaml") == {}


# ----- Deduplication -----


def test_deduplicate_bioprojects_first_wins() -> None:
    rows = [
        {**empty_bioproject(), "bioproject_id": "PRJ1", "source": "seed", "protein_target": "VEGF"},
        {
            **empty_bioproject(),
            "bioproject_id": "PRJ1",
            "source": "ena",
            "study_title": "VEGF aptamers",
        },
    ]
    out = deduplicate_bioprojects(rows)
    assert len(out) == 1
    # First row wins on bp_id; subsequent fills in blank columns
    assert out[0]["source"] == "seed"
    assert out[0]["protein_target"] == "VEGF"
    assert out[0]["study_title"] == "VEGF aptamers"  # filled from second row


def test_deduplicate_samples_by_srr() -> None:
    rows = [
        {**empty_sample(), "srr": "SRR1", "bioproject_id": "PRJ1"},
        {**empty_sample(), "srr": "SRR1", "bioproject_id": "PRJ_dup"},
        {**empty_sample(), "srr": "SRR2", "bioproject_id": "PRJ2"},
    ]
    out = deduplicate_samples(rows)
    assert len(out) == 2
    assert next(r for r in out if r["srr"] == "SRR1")["bioproject_id"] == "PRJ1"


def test_deduplicate_samples_skips_empty_srr() -> None:
    rows = [{**empty_sample(), "srr": "", "bioproject_id": "PRJ1"}]
    out = deduplicate_samples(rows)
    assert out == []


# ----- CSV writer -----


def test_write_csv_round_trip(tmp_path: Path) -> None:
    rows = [
        {**empty_bioproject(), "bioproject_id": "PRJ1", "source": "seed"},
        {**empty_bioproject(), "bioproject_id": "PRJ2", "source": "ena"},
    ]
    path = tmp_path / "bioprojects.csv"
    write_csv(rows, BIOPROJECT_COLS, path)

    with open(path, encoding="utf-8") as f:
        read = list(csv.DictReader(f))
    assert len(read) == 2
    assert read[0]["bioproject_id"] == "PRJ1"
    assert read[0]["source"] == "seed"
    # Confirm column order is preserved
    with open(path, encoding="utf-8") as f:
        header = f.readline().strip().split(",")
    assert header == BIOPROJECT_COLS


def test_write_csv_extras_action_ignore(tmp_path: Path) -> None:
    """Rows with extra keys should not raise; the extras are silently dropped."""
    rows = [{**empty_sample(), "srr": "SRR1", "extra_column": "should be dropped"}]
    path = tmp_path / "samples.csv"
    write_csv(rows, SAMPLE_COLS, path)
    with open(path, encoding="utf-8") as f:
        first = next(csv.DictReader(f))
    assert "extra_column" not in first
    assert first["srr"] == "SRR1"


# ----- Sentinel for deferred library-type classifier -----


def test_classify_all_writes_sentinel_when_classifier_absent() -> None:
    """When library_type_classifier is not importable (v0.2 deferred), every
    unclassified BioProject must get the `not_assessed_v0.1` sentinel + an
    evidence JSON — not an empty string that downstream callers would treat
    as a successful 'unclassified' verdict."""
    rows = [
        {**empty_bioproject(), "bioproject_id": "PRJ1"},
        {
            **empty_bioproject(),
            "bioproject_id": "PRJ2",
            "library_type_verification": "RNA_confirmed",
        },
    ]
    out = _classify_all(rows)
    assert out[0]["library_type_verification"] == NOT_ASSESSED_V0_1
    assert "no_classifier" in out[0]["library_type_evidence"]
    # Pre-existing verdict (e.g. from seed) is preserved
    assert out[1]["library_type_verification"] == "RNA_confirmed"


# ----- Round columns sanity -----


def test_round_cols_matches_round_record_to_dict() -> None:
    """RoundRecord.to_dict() keys must match ROUND_COLS so CSV writes don't drop fields."""
    from selexprep.fetch.metadata import RoundRecord

    rr = RoundRecord(
        srr="SRR1",
        round_number=1,
        confidence="HIGH",
        source_field="sample_title",
        matched_pattern="round_word_digit",
    )
    d = rr.to_dict()
    for col in ROUND_COLS:
        assert col in d, f"ROUND_COLS includes {col!r} but RoundRecord.to_dict() does not"
