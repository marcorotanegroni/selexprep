"""End-to-end test for ``selexprep.qc.runner.run_qc``."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml

from selexprep.library.report import LibraryReport
from selexprep.manifest import build_manifest_from_extract_result, write_manifest_json
from selexprep.qc.runner import run_qc


def _make_library_report(**overrides: object) -> LibraryReport:
    base = {
        "primer_5p": "GGTAATACGACTCACTATAGGG",
        "primer_3p": "CCATGCATGCATGCATGCAT",
        "variants_5p": [],
        "variants_3p": [],
        "known_adapter_hits": {},
        "extraction_mode": "BOTH_PRIMERS_SINGLE_READ",
        "full_insert_recovered": True,
        "read_source": "R1",
        "required_action": "NONE",
        "orientation": "FORWARD",
        "n_length_mode": 30,
        "n_length_distribution": {30: 1000},
        "n_length_confidence": 1.0,
        "match_rate_5p": 0.95,
        "match_rate_3p": 0.92,
        "position_consistency_5p": 0.95,
        "position_consistency_3p": 0.92,
        "read_fraction_used_for_inference": 1.0,
        "sampling_seed": 42,
        "confidence": 0.85,
        "status": "HIGH",
        "failure_reason": None,
    }
    base.update(overrides)
    return LibraryReport(**base)  # type: ignore[arg-type]


def _build_outdir(tmp_path: Path) -> Path:
    """Build a synthetic outdir matching the post-extract layout.

    Includes:
    - selexprep_manifest.json
    - trim_reports.json
    - round_NN/counts.parquet
    - (optional) strand_report.tsv
    """
    outdir = tmp_path / "ds"
    outdir.mkdir()

    # Per-round counts.parquet
    for r in range(3):
        round_dir = outdir / f"round_{r:02d}"
        round_dir.mkdir()
        df = pd.DataFrame(
            {
                "sequence": [f"r{r}_seq_{i}" for i in range(50)],
                "reads": [(50 - i) for i in range(50)],
                "rank": list(range(1, 51)),
                "rpm": [(50 - i) * 1000.0 for i in range(50)],
            }
        )
        df.to_parquet(round_dir / "counts.parquet", index=False, compression="zstd")

    # trim_reports.json
    trim_reports = [
        {
            "cutadapt_cmd": ["cutadapt", "..."],
            "n_in": 1000 * (3 - r),
            "n_out": 950 * (3 - r),
            "return_code": 0,
            "output_paths": [str(outdir / f"round_{r:02d}" / "extracted.fasta.gz")],
        }
        for r in range(3)
    ]
    (outdir / "trim_reports.json").write_text(json.dumps(trim_reports), encoding="utf-8")

    # Manifest
    manifest = build_manifest_from_extract_result(
        library_report=_make_library_report(),
        input_paths=[],
        output_paths=[],
        accession="SRR_TEST",
        bioproject_id="PRJ_TEST",
        runs=["SRR_TEST"],
        parameters={},
    )
    write_manifest_json(manifest, outdir / "selexprep_manifest.json")
    return outdir


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------


def test_run_qc_emits_flags_yaml_and_plots(tmp_path: Path) -> None:
    outdir = _build_outdir(tmp_path)
    result = run_qc(outdir / "selexprep_manifest.json")

    assert result.flags_yaml_path is not None
    assert result.flags_yaml_path.exists()
    assert len(result.plot_paths) == 4  # all 4 plots
    for p in result.plot_paths:
        assert p.exists()


def test_run_qc_writes_into_qc_subdir_by_default(tmp_path: Path) -> None:
    outdir = _build_outdir(tmp_path)
    result = run_qc(outdir / "selexprep_manifest.json")
    assert result.flags_yaml_path is not None
    assert result.flags_yaml_path.parent == outdir / "qc"


def test_run_qc_writes_into_custom_outdir(tmp_path: Path) -> None:
    outdir = _build_outdir(tmp_path)
    custom = tmp_path / "custom_qc"
    result = run_qc(outdir / "selexprep_manifest.json", outdir=custom)
    assert result.flags_yaml_path is not None
    assert result.flags_yaml_path.parent == custom


def test_run_qc_flags_yaml_is_valid_yaml(tmp_path: Path) -> None:
    outdir = _build_outdir(tmp_path)
    result = run_qc(outdir / "selexprep_manifest.json")
    assert result.flags_yaml_path is not None
    parsed = yaml.safe_load(result.flags_yaml_path.read_text())
    # Either empty list or list of flag dicts.
    assert isinstance(parsed, list)
    for entry in parsed:
        assert "name" in entry
        assert "severity" in entry
        assert "evidence" in entry


# ---------------------------------------------------------------------------
# Flag firing
# ---------------------------------------------------------------------------


def test_run_qc_fires_low_total_reads_for_tiny_pool(tmp_path: Path) -> None:
    outdir = tmp_path / "ds"
    outdir.mkdir()
    # One round, only 50 reads total -> below the 10k threshold.
    round_dir = outdir / "round_00"
    round_dir.mkdir()
    df = pd.DataFrame(
        {
            "sequence": [f"s{i}" for i in range(5)],
            "reads": [10, 10, 10, 10, 10],
            "rank": [1, 2, 3, 4, 5],
            "rpm": [200_000.0, 200_000.0, 200_000.0, 200_000.0, 200_000.0],
        }
    )
    df.to_parquet(round_dir / "counts.parquet", index=False)

    manifest = build_manifest_from_extract_result(
        library_report=_make_library_report(),
        input_paths=[],
        output_paths=[],
        accession=None,
        bioproject_id=None,
        runs=[],
        parameters={},
    )
    write_manifest_json(manifest, outdir / "selexprep_manifest.json")

    result = run_qc(outdir / "selexprep_manifest.json")
    flag_names = {f.name for f in result.flags}
    assert "low_total_reads" in flag_names


def test_run_qc_fires_read_merging_info_for_paired_end_split(tmp_path: Path) -> None:
    outdir = tmp_path / "ds"
    outdir.mkdir()
    # No counts.parquet needed for this manifest-only flag.
    manifest = build_manifest_from_extract_result(
        library_report=_make_library_report(
            extraction_mode="PAIRED_END_SPLIT_PRIMERS",
            required_action="READ_MERGING_RECOMMENDED",
            full_insert_recovered=False,
        ),
        input_paths=[],
        output_paths=[],
        accession=None,
        bioproject_id=None,
        runs=[],
        parameters={},
    )
    write_manifest_json(manifest, outdir / "selexprep_manifest.json")

    result = run_qc(outdir / "selexprep_manifest.json")
    flag_names = {f.name for f in result.flags}
    assert "requires_read_merging_for_full_insert" in flag_names
