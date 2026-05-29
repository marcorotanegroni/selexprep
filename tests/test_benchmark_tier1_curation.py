"""Phase 6b.9 integrity tests for the Tier-1 benchmark curation.

These don't execute the Snakefile (not a CI workload) — they pin the
*curation contract* so it can't regress silently:

- the 4 curated round-maps exist, are valid TSVs, and are referenced
  correctly from ground_truth.tsv;
- PRJNA315881 stays ``auto`` (it's a multiplexed deposit, deliberately
  not curated to a single round);
- the Snakefile manifest enforces the R1-only paired-end policy.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

_BENCH = Path(__file__).resolve().parents[1] / "benchmarks"
_GROUND_TRUTH = _BENCH / "ground_truth.tsv"
_SNAKEFILE = _BENCH / "Snakefile"

# The 4 zero-round blockers curated in Phase 6b.9 (PRJNA315881 is NOT here).
_CURATED = {
    "PRJEB28411",
    "PRJNA935703",
    "PRJNA975735",
    "PRJNA883192",
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
