"""Phase 6b.9 + 6b.10 integrity tests for the Tier-1 benchmark curation.

These don't execute the Snakefile (not a CI workload) — they pin the
*curation contract* so it can't regress silently:

- the curated round-maps exist, are valid TSVs, and are referenced
  correctly from ground_truth.tsv;
- PRJNA315881 stays ``auto`` (it's a multiplexed deposit, deliberately
  not curated to a single round);
- the Snakefile manifest enforces the R1-only paired-end policy, and the
  fetch rule has the Phase 6b.10 partial-fetch fix;
- Phase 6b.10 two-arm reframe: every row carries a valid ``read_state``
  arm label + boolean flags; out-of-scope rows live in
  ``excluded_datasets.tsv`` (not ground_truth); read_state labels are
  mirrored in ``read_state_evidence.tsv``; the pre-detect screening
  decision for the whole original-11 pool is recorded in
  ``screening_log.tsv`` (non-circularity audit trail).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

_BENCH = Path(__file__).resolve().parents[1] / "benchmarks"
_GROUND_TRUTH = _BENCH / "ground_truth.tsv"
_SNAKEFILE = _BENCH / "Snakefile"
_READ_STATE_EVIDENCE = _BENCH / "read_state_evidence.tsv"
_EXCLUDED = _BENCH / "excluded_datasets.tsv"
_SCREENING_LOG = _BENCH / "screening_log.tsv"

# Curated rows remaining in the benchmark after the 6b.10 cleanup
# (PRJNA935703 + PRJNA975735 were removed to excluded_datasets.tsv;
# PRJNA315881 stays auto — it's multiplexed).
_CURATED = {
    "PRJEB28411",
    "PRJNA883192",
}

# Phase 6b.10 benchmark set (post-cleanup), labeled by read-state arm.
_RECOVERY = {"PRJDB9110", "PRJDB9111", "PRJNA883192", "PRJNA315881", "PRJEB70964"}
_SPECIFICITY = {"PRJEB28411", "PRJEB22637", "PRJNA990511"}
# Out-of-scope rows removed from ground_truth.tsv → excluded_datasets.tsv.
_EXCLUDED_ACCESSIONS = {"PRJNA728693", "PRJNA935703", "PRJNA975735"}
# Evidence-based, pre-detect reason vocabulary (NEVER "primers_unrecoverable").
_REASON_VOCAB = {
    "primer_truth_unverified",
    "inaccessible",
    "unsupported_architecture",
    "not_raw_reads",
    "insufficient_depth",
    "multiplexed_without_demux",
}


def _gt() -> pd.DataFrame:
    return pd.read_csv(_GROUND_TRUTH, sep="\t", dtype=str).fillna("")


def test_exactly_these_rows_are_curated() -> None:
    gt = _gt()
    curated = set(gt.loc[gt["round_map_source"] == "curated", "accession"])
    assert curated == _CURATED


def test_prjna315881_stays_auto_multiplex() -> None:
    """PRJNA315881's SRR3279660 is a rounds-1-4 multiplexed FASTQ; mapping it
    to one round would merge 4 rounds, so it is deliberately left auto."""
    gt = _gt()
    row = gt.loc[gt["accession"] == "PRJNA315881"].iloc[0]
    assert row["round_map_source"] == "auto"


def test_curated_round_maps_exist_and_are_valid() -> None:
    gt = _gt().set_index("accession")
    for acc in _CURATED:
        rel = gt.at[acc, "round_map_path"]
        assert rel, f"{acc}: curated row missing round_map_path"
        path = _BENCH / rel
        assert path.exists(), f"{acc}: round-map {path} missing"
        rm = pd.read_csv(path, sep="\t")
        assert list(rm.columns) == ["file", "round_number"], f"{acc}: bad columns"
        assert len(rm) >= 1
        # round_number column is integer-valued
        rm["round_number"].astype(int)
        # every listed file is a .fastq.gz basename, and R2 mates are NOT
        # listed (R1-only policy — paired rows use _1 only)
        for fn in rm["file"]:
            assert fn.endswith(".fastq.gz"), f"{acc}: {fn} not a fastq.gz"
            assert not fn.endswith("_2.fastq.gz"), f"{acc}: {fn} is an R2 mate"


def test_curated_notes_carry_inference_caveat() -> None:
    gt = _gt().set_index("accession")
    for acc in _CURATED:
        notes = gt.at[acc, "notes"]
        assert "inferred from sample" in notes, f"{acc}: missing round-inference caveat"
        assert "not per-round biological claims" in notes, f"{acc}: missing claim caveat"


def test_snakefile_enforces_r1_only_manifest() -> None:
    """The manifest filter excludes R2 mates (paired-end policy). Guards
    against a silent regression that would re-introduce R1/R2 mixing."""
    text = _SNAKEFILE.read_text(encoding="utf-8")
    assert "fastqs.manifest" in text
    assert "! -name '*_2.fastq.gz'" in text


def test_snakefile_detect_guards_empty_manifest() -> None:
    """detect uses an if/fi guard (not a brace group, which collides with
    Snakemake's {} placeholders) so a failed fetch fails fast + readably."""
    text = _SNAKEFILE.read_text(encoding="utf-8")
    assert "if [ ! -s {input.manifest} ]; then" in text


# ===========================================================================
# Phase 6b.10 — two-arm reframe integrity
# ===========================================================================


def test_snakefile_has_partial_fetch_fix() -> None:
    """The fetch rule scopes errexit off around fetch (set +e / set -e) so a
    partial download (e.g. PRJEB70964 17/27) proceeds while a total crash
    still fails. Guards against re-introducing a bare ``&&`` chain that would
    abort the rule on any non-zero fetch exit under set -euo pipefail."""
    text = _SNAKEFILE.read_text(encoding="utf-8")
    assert "set +e" in text
    assert "set -e" in text
    assert "rc=$?" in text
    # requires BOTH a non-empty manifest AND a non-empty fetch_metadata.json
    assert "{output.marker}" in text
    assert "partial download; proceeding" in text


def test_excluded_rows_removed_from_ground_truth() -> None:
    """Out-of-scope deposits are removed from ground_truth.tsv entirely
    (they are NOT figure arms — they live in excluded_datasets.tsv)."""
    accessions = set(_gt()["accession"])
    assert accessions.isdisjoint(_EXCLUDED_ACCESSIONS)
    # and the benchmark set is exactly recovery + specificity
    assert accessions == _RECOVERY | _SPECIFICITY


def test_every_row_has_valid_read_state_and_flags() -> None:
    gt = _gt()
    valid_states = {"raw_standard", "pre_trimmed"}
    flag_cols = ["mono_round", "partial_fetch", "paired_end_r1_only", "demultiplexed"]
    for _, row in gt.iterrows():
        acc = row["accession"]
        assert row["read_state"] in valid_states, f"{acc}: bad read_state {row['read_state']!r}"
        for col in flag_cols:
            assert row[col] in {"true", "false"}, f"{acc}: {col}={row[col]!r} not true/false"


def test_read_state_arms_match_expected_membership() -> None:
    gt = _gt().set_index("accession")
    for acc in _RECOVERY:
        assert gt.at[acc, "read_state"] == "raw_standard", f"{acc} should be recovery arm"
    for acc in _SPECIFICITY:
        assert gt.at[acc, "read_state"] == "pre_trimmed", f"{acc} should be specificity arm"


def test_flag_assignments_match_plan() -> None:
    """The plan-named orthogonal flags are set on exactly the named rows."""
    gt = _gt().set_index("accession")
    assert gt.at["PRJNA315881", "mono_round"] == "true"
    assert gt.at["PRJEB70964", "partial_fetch"] == "true"
    assert gt.at["PRJNA883192", "paired_end_r1_only"] == "true"
    # no row in the current set is included via reproducible demux
    assert (gt["demultiplexed"] == "true").sum() == 0


# --- read_state_evidence.tsv -----------------------------------------------


def test_read_state_evidence_schema_and_coverage() -> None:
    assert _READ_STATE_EVIDENCE.exists()
    ev = pd.read_csv(_READ_STATE_EVIDENCE, sep="\t", dtype=str).fillna("")
    assert list(ev.columns) == [
        "accession",
        "expected_n",
        "expected_full_length",
        "observed_read_len_median",
        "observed_read_len_range",
        "n_reads_checked",
        "architecture_evidence",
        "read_state",
        "evidence_note",
    ]
    # covers exactly the benchmark set, and read_state mirrors ground_truth
    gt = _gt().set_index("accession")
    assert set(ev["accession"]) == _RECOVERY | _SPECIFICITY
    for _, row in ev.iterrows():
        acc = row["accession"]
        assert row["read_state"] in {"raw_standard", "pre_trimmed"}
        assert row["read_state"] == gt.at[acc, "read_state"], (
            f"{acc}: evidence/gt read_state mismatch"
        )


# --- excluded_datasets.tsv -------------------------------------------------


def test_excluded_datasets_schema_and_reasons() -> None:
    assert _EXCLUDED.exists()
    ex = pd.read_csv(_EXCLUDED, sep="\t", dtype=str).fillna("")
    assert "accession" in ex.columns
    assert "reason" in ex.columns
    assert set(ex["accession"]) == _EXCLUDED_ACCESSIONS
    for reason in ex["reason"]:
        assert reason in _REASON_VOCAB, f"excluded reason {reason!r} not in vocab"
        assert reason != "primers_unrecoverable"


# --- screening_log.tsv (non-circularity audit trail) -----------------------


def test_screening_log_covers_full_original_pool() -> None:
    """The pre-detect screening decision is recorded for ALL original-11
    candidates — proof selection wasn't conditioned on detect success."""
    assert _SCREENING_LOG.exists()
    sl = pd.read_csv(_SCREENING_LOG, sep="\t", dtype=str).fillna("")
    assert list(sl.columns) == ["accession", "read_state", "included", "arm", "reason", "evidence"]
    assert set(sl["accession"]) == _RECOVERY | _SPECIFICITY | _EXCLUDED_ACCESSIONS


def test_screening_log_included_flags_match_arms() -> None:
    sl = pd.read_csv(_SCREENING_LOG, sep="\t", dtype=str).fillna("").set_index("accession")
    for acc in _RECOVERY | _SPECIFICITY:
        assert sl.at[acc, "included"] == "true", f"{acc} should be included"
    for acc in _EXCLUDED_ACCESSIONS:
        assert sl.at[acc, "included"] == "false", f"{acc} should be excluded"


def test_screening_log_excluded_reasons_are_evidence_based() -> None:
    """Excluded rows carry an evidence-based, pre-detect reason — NEVER the
    post-hoc / detect-dependent ``primers_unrecoverable``."""
    sl = pd.read_csv(_SCREENING_LOG, sep="\t", dtype=str).fillna("")
    excluded = sl[sl["included"] == "false"]
    for _, row in excluded.iterrows():
        assert row["reason"] in _REASON_VOCAB, f"{row['accession']}: reason {row['reason']!r}"
        assert row["reason"] != "primers_unrecoverable"
