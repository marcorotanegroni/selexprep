"""+ integrity tests for the Tier-1 benchmark curation.

These don't execute the Snakefile (not a CI workload) — they pin the
*curation contract* so it can't regress silently:

- the curated round-maps exist, are valid TSVs, and are referenced
  correctly from ground_truth.tsv;
- PRJNA315881 is excluded until a reproducible demultiplexing step exists;
- the Snakefile writes separate R1/R2 manifests and passes R2 through
  ``selexprep detect --paired-r2`` instead of silently forcing paired
  datasets into R1-only mode;
- two-arm reframe: every row carries a valid ``read_state``
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

# Curated rows remaining in the benchmark after the cleanup
# (PRJNA935703 + PRJNA975735 + PRJNA315881 were moved to excluded_datasets.tsv).
_CURATED = {
    "PRJEB28411",
    "PRJNA883192",
}

# benchmark set (post-cleanup), labeled by read-state role.
# PRJNA615076 = first recovery addition (Kolm 2020, E. faecalis DNA whole-cell SELEX)
_RECOVERY = {
    "PRJDB9110",
    "PRJDB9111",
    # Ground truth from RIBOMIC patent WO2020204151 (SEQ ID NO 39-41) — external
    # to both the deposited reads and any selexprep output.
    "PRJDB19098",
    # Retained as a row but scored out (verified=false): its 3' truth is
    # read-derived, so the pair has no fully external ground truth, and the
    # deposit is paired-end split-primer, so its extraction cannot be verified
    # without read merging. Reported under skipped_unverified_accessions.
    "PRJNA883192",
    "PRJNA615076",
    "PRJNA809588",
    "PRJNA1395820",
    "PRJEB62495",
}
_SPECIFICITY = {"PRJEB28411", "PRJEB22637", "PRJNA990511"}
# Adapter-collision negative control (5' constant = revcomp of a known adapter);
# excluded from the recovery denominator.
_ADAPTER_CONTROL = {"PRJEB70964"}
# All in-benchmark rows (any of the three roles).
_BENCHMARK = _RECOVERY | _SPECIFICITY | _ADAPTER_CONTROL
# Out-of-scope rows removed from ground_truth.tsv → excluded_datasets.tsv.
# PRJNA315881: round-5 pool multiplexed-by-condition (undocumented barcodes) →
# not a fair primer-inference test (multiplexing/offset/truncation);
# re-entrant after reproducible demux that does not use primer info.
_EXCLUDED_ACCESSIONS = {
    "PRJNA728693",
    "PRJNA935703",
    "PRJNA975735",
    "PRJNA315881",
    "PRJNA1068659",
}
# Evidence-based, pre-detect reason vocabulary (NEVER "primers_unrecoverable").
_REASON_VOCAB = {
    "primer_truth_unverified",
    "primer_truth_conflict",
    "inaccessible",
    "unsupported_architecture",
    "not_raw_reads",
    "insufficient_depth",
    "multiplexed_without_demux",
    "boundary_ambiguous",
}


def _gt() -> pd.DataFrame:
    return pd.read_csv(_GROUND_TRUTH, sep="\t", dtype=str).fillna("")


def test_exactly_these_rows_are_curated() -> None:
    gt = _gt()
    curated = set(gt.loc[gt["round_map_source"] == "curated", "accession"])
    assert curated == _CURATED


def test_prjna315881_excluded_multiplexed_without_demux() -> None:
    """PRJNA315881's only fetchable pool is multiplexed-by-condition (undocumented
    barcodes), so it's not a fair primer-inference test — removed
    from ground_truth and recorded in excluded_datasets with reason
    multiplexed_without_demux; re-entrant after reproducible demux."""
    assert "PRJNA315881" not in set(_gt()["accession"])
    ex = pd.read_csv(_EXCLUDED, sep="\t", dtype=str).fillna("").set_index("accession")
    assert "PRJNA315881" in ex.index
    assert ex.at["PRJNA315881", "reason"] == "multiplexed_without_demux"


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
        # every listed file is a .fastq.gz basename
        for fn in rm["file"]:
            assert fn.endswith(".fastq.gz"), f"{acc}: {fn} not a fastq.gz"
        if acc == "PRJNA883192":
            # Paired-end row: curated map must route both mates so the
            # benchmark can pass R2 through detect instead of doing R1-only.
            assert set(rm["file"]) == {
                "SRR21674163_1.fastq.gz",
                "SRR21674163_2.fastq.gz",
                "SRR21674162_1.fastq.gz",
                "SRR21674162_2.fastq.gz",
            }


def test_curated_notes_carry_inference_caveat() -> None:
    gt = _gt().set_index("accession")
    for acc in _CURATED:
        notes = gt.at[acc, "notes"]
        assert "inferred from sample" in notes, f"{acc}: missing round-inference caveat"
        assert "not per-round biological claims" in notes, f"{acc}: missing claim caveat"


def test_snakefile_writes_separate_r1_r2_manifests() -> None:
    """The R1 manifest excludes R2 mates, but R2 is captured separately and
    passed via ``--paired-r2``. Guards against both old failure modes:
    R1/R2 commingling and R2 silently discarded from the Tier 1 benchmark."""
    text = _SNAKEFILE.read_text(encoding="utf-8")
    assert "fastqs.manifest" in text
    assert "fastqs.r2.manifest" in text
    assert "! -name '*_2.fastq.gz'" in text
    assert "-name '*_2.fastq.gz'" in text
    assert "--paired-r2" in text


def test_snakefile_detect_guards_empty_manifest() -> None:
    """detect uses an if/fi guard (not a brace group, which collides with
    Snakemake's {} placeholders) so a failed fetch fails fast + readably."""
    text = _SNAKEFILE.read_text(encoding="utf-8")
    assert "if [ ! -s {input.manifest} ]; then" in text


# ===========================================================================
# two-arm reframe integrity
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
    # and the benchmark set is exactly the three roles
    assert accessions == _BENCHMARK


def test_every_row_has_valid_read_state_and_flags() -> None:
    gt = _gt()
    valid_states = {"raw_standard", "pre_trimmed", "adapter_control"}
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
    for acc in _ADAPTER_CONTROL:
        assert gt.at[acc, "read_state"] == "adapter_control", f"{acc} should be adapter-control"


def test_flag_assignments_match_plan() -> None:
    """The plan-named orthogonal flags are set on exactly the named rows."""
    gt = _gt().set_index("accession")
    assert gt.at["PRJNA990511", "mono_round"] == "true"
    assert gt.at["PRJEB70964", "partial_fetch"] == "true"
    assert gt.at["PRJEB28411", "paired_end_r1_only"] == "false"
    assert gt.at["PRJNA883192", "paired_end_r1_only"] == "false"
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
    assert set(ev["accession"]) == _BENCHMARK
    for _, row in ev.iterrows():
        acc = row["accession"]
        assert row["read_state"] in {"raw_standard", "pre_trimmed", "adapter_control"}
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
    assert set(sl["accession"]) == _BENCHMARK | _EXCLUDED_ACCESSIONS


def test_screening_log_included_flags_match_arms() -> None:
    sl = pd.read_csv(_SCREENING_LOG, sep="\t", dtype=str).fillna("").set_index("accession")
    for acc in _BENCHMARK:
        assert sl.at[acc, "included"] == "true", f"{acc} should be included"
    for acc in _EXCLUDED_ACCESSIONS:
        assert sl.at[acc, "included"] == "false", f"{acc} should be excluded"


# detect-inference OUTPUTS that must not leak into pre-detect evidence files
# (inclusion-protocol hard rule 4 — non-circularity). Bare "detect" is allowed
# (evidence may say "...NOT detect" to assert independence); these are the
# inference outputs whose presence would be leakage.
_FORBIDDEN_DETECT_STRINGS = (
    "detect match",
    "match_rate",
    "extraction_mode",
    "library_report",
    "unable_to_infer",
    "unable_to_extract",
    "primers_unrecoverable",
    "n_length_mode",
)


def test_evidence_files_have_no_detect_derived_strings() -> None:
    """Non-circularity guard: pre-detect evidence files must not cite
    detect's inference outputs. Pins inclusion-protocol hard rule 4 against
    regression."""
    for name in ("screening_log.tsv", "excluded_datasets.tsv", "read_state_evidence.tsv"):
        text = (_BENCH / name).read_text(encoding="utf-8").lower()
        for token in _FORBIDDEN_DETECT_STRINGS:
            assert token not in text, f"{name}: detect-derived string {token!r} present (leakage)"


def test_screening_log_excluded_reasons_are_evidence_based() -> None:
    """Excluded rows carry an evidence-based, pre-detect reason — NEVER the
    post-hoc / detect-dependent ``primers_unrecoverable``."""
    sl = pd.read_csv(_SCREENING_LOG, sep="\t", dtype=str).fillna("")
    excluded = sl[sl["included"] == "false"]
    for _, row in excluded.iterrows():
        assert row["reason"] in _REASON_VOCAB, f"{row['accession']}: reason {row['reason']!r}"
        assert row["reason"] != "primers_unrecoverable"
