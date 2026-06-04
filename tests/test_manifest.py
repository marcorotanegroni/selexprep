"""Unit tests for ``selexprep.manifest``."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from selexprep._io import sha256_file
from selexprep.library.report import LibraryReport
from selexprep.manifest import (
    SelexprepManifestV1,
    build_manifest_from_extract_result,
    compute_sha256s,
    read_manifest_json,
    write_manifest_json,
)


def _make_library_report() -> LibraryReport:
    """Minimal LibraryReport for manifest tests."""
    return LibraryReport(
        primer_5p="ACGTACGTACGTACGT",
        primer_3p="TGCATGCATGCATGCA",
        variants_5p=[],
        variants_3p=[],
        known_adapter_hits={"TRUSEQ_R1": 0},
        extraction_mode="BOTH_PRIMERS_SINGLE_READ",
        full_insert_recovered=True,
        read_source="R1",
        required_action="NONE",
        orientation="FORWARD",
        n_length_mode=30,
        n_length_distribution={30: 1000},
        n_length_confidence=1.0,
        match_rate_5p=0.95,
        match_rate_3p=0.95,
        position_consistency_5p=0.95,
        position_consistency_3p=0.95,
        read_fraction_used_for_inference=1.0,
        sampling_seed=42,
        confidence=0.85,
        status="HIGH",
        failure_reason=None,
    )


# ---------------------------------------------------------------------------
# Schema + roundtrip
# ---------------------------------------------------------------------------


def test_manifest_version_pinned_to_v1() -> None:
    """manifest_version is a Literal — defaults to the v1 string."""
    m = build_manifest_from_extract_result(
        library_report=_make_library_report(),
        input_paths=[],
        output_paths=[],
        accession="SRR000000",
        bioproject_id="PRJEB00000",
        runs=["SRR000000"],
        parameters={},
    )
    assert m.manifest_version == "selexprep_manifest_v1"


def test_manifest_denormalizes_classification_fields() -> None:
    """extraction_mode/read_source/required_action/full_insert_recovered
    are copied from the nested LibraryReport for quick scan."""
    lr = _make_library_report()
    m = build_manifest_from_extract_result(
        library_report=lr,
        input_paths=[],
        output_paths=[],
        accession=None,
        bioproject_id=None,
        runs=[],
        parameters={},
    )
    assert m.extraction_mode == lr.extraction_mode
    assert m.read_source == lr.read_source
    assert m.required_action == lr.required_action
    assert m.full_insert_recovered == lr.full_insert_recovered
    assert m.sampling_seed == lr.sampling_seed


def test_manifest_captures_dep_versions() -> None:
    m = build_manifest_from_extract_result(
        library_report=_make_library_report(),
        input_paths=[],
        output_paths=[],
        accession=None,
        bioproject_id=None,
        runs=[],
        parameters={},
    )
    # Non-empty strings; exact values depend on the install environment.
    assert m.selexprep_version
    assert m.python_version
    assert m.dnaio_version
    assert m.pyarrow_version
    # cutadapt may report "unknown" if not on PATH; verify it's a string.
    assert isinstance(m.cutadapt_version, str)


def test_manifest_is_frozen() -> None:
    m = build_manifest_from_extract_result(
        library_report=_make_library_report(),
        input_paths=[],
        output_paths=[],
        accession=None,
        bioproject_id=None,
        runs=[],
        parameters={},
    )
    with pytest.raises(ValidationError):
        m.sampling_seed = 0  # type: ignore[misc]


def test_manifest_rejects_extra_fields() -> None:
    with pytest.raises(ValidationError):
        SelexprepManifestV1(
            selexprep_version="x",
            python_version="x",
            cutadapt_version="x",
            dnaio_version="x",
            pyarrow_version="x",
            accession=None,
            bioproject_id=None,
            runs=[],
            input_sha256={},
            output_sha256={},
            library_report=_make_library_report(),
            extraction_mode="BOTH_PRIMERS_SINGLE_READ",
            read_source="R1",
            required_action="NONE",
            full_insert_recovered=True,
            parameters={},
            runtime_seconds_per_stage={},
            flags=[],
            sampling_seed=42,
            mystery_field="oops",  # type: ignore[call-arg]
        )


def test_manifest_json_roundtrip_preserves_data(tmp_path: Path) -> None:
    original = build_manifest_from_extract_result(
        library_report=_make_library_report(),
        input_paths=[],
        output_paths=[],
        accession="SRR123",
        bioproject_id="PRJ123",
        runs=["SRR123"],
        parameters={"foo": "bar"},
    )
    path = tmp_path / "manifest.json"
    write_manifest_json(original, path)
    loaded = read_manifest_json(path)
    assert loaded == original


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------


def test_manifest_json_is_deterministic_across_writes(tmp_path: Path) -> None:
    m = build_manifest_from_extract_result(
        library_report=_make_library_report(),
        input_paths=[],
        output_paths=[],
        accession="SRR123",
        bioproject_id="PRJ123",
        runs=["SRR123"],
        parameters={"foo": "bar"},
    )
    path_a = tmp_path / "a.json"
    path_b = tmp_path / "b.json"
    write_manifest_json(m, path_a)
    write_manifest_json(m, path_b)
    assert path_a.read_bytes() == path_b.read_bytes()
    assert sha256_file(path_a) == sha256_file(path_b)


def test_manifest_int_keys_serialized_numerically(tmp_path: Path) -> None:
    """The nested LibraryReport's int-keyed n_length_distribution should
    be sorted numerically, not lexically ('10', '100', '20')."""
    lr = _make_library_report().model_copy(
        update={"n_length_distribution": {10: 1, 20: 2, 100: 3, 30: 4}}
    )
    m = build_manifest_from_extract_result(
        library_report=lr,
        input_paths=[],
        output_paths=[],
        accession=None,
        bioproject_id=None,
        runs=[],
        parameters={},
    )
    path = tmp_path / "manifest.json"
    write_manifest_json(m, path)
    payload = json.loads(path.read_text())
    keys = list(payload["library_report"]["n_length_distribution"].keys())
    assert keys == sorted(keys, key=int)


# ---------------------------------------------------------------------------
# compute_sha256s — FASTA/TSV/JSON only
# ---------------------------------------------------------------------------


def test_compute_sha256s_hashes_fasta_tsv_json(tmp_path: Path) -> None:
    fa = tmp_path / "x.fasta.gz"
    fa.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x00fakefasta")
    tsv = tmp_path / "y.tsv"
    tsv.write_text("col\nrow1\n", encoding="utf-8")
    js = tmp_path / "z.json"
    js.write_text("{}\n", encoding="utf-8")

    hashes = compute_sha256s([fa, tsv, js])
    assert set(hashes) == {"x.fasta.gz", "y.tsv", "z.json"}
    assert all(len(h) == 64 for h in hashes.values())


def test_compute_sha256s_skips_parquet(tmp_path: Path) -> None:
    """Parquet hashes are advisory per the design; not in
    output_sha256."""
    pq = tmp_path / "counts.parquet"
    pq.write_bytes(b"PAR1\x00")  # not real parquet, but enough for the suffix test
    hashes = compute_sha256s([pq])
    assert hashes == {}


def test_compute_sha256s_skips_missing_paths(tmp_path: Path) -> None:
    nonexistent = tmp_path / "missing.json"
    hashes = compute_sha256s([nonexistent])
    assert hashes == {}


def test_compute_sha256s_distinct_keys_per_round_with_root(tmp_path: Path) -> None:
    """regression: per-round outputs that share the
    same basename (round_00/extracted.fasta.gz vs round_01/extracted.fasta.gz)
    must NOT collide in output_sha256. Pass ``root=outdir`` to key by
    relative path; basename-only keying would last-write-wins one of them."""
    outdir = tmp_path / "out"
    (outdir / "round_00").mkdir(parents=True)
    (outdir / "round_01").mkdir(parents=True)
    p1 = outdir / "round_00" / "extracted.fasta.gz"
    p2 = outdir / "round_01" / "extracted.fasta.gz"
    # Distinct content so the hashes themselves differ — the test verifies
    # both end up in the dict, not just that the hash function works.
    p1.write_bytes(b"round_00 distinct content " + b"A" * 100)
    p2.write_bytes(b"round_01 distinct content " + b"B" * 100)

    hashes = compute_sha256s([p1, p2], root=outdir)

    assert "round_00/extracted.fasta.gz" in hashes
    assert "round_01/extracted.fasta.gz" in hashes
    assert hashes["round_00/extracted.fasta.gz"] != hashes["round_01/extracted.fasta.gz"]


def test_compute_sha256s_without_root_falls_back_to_basename(tmp_path: Path) -> None:
    """Backward-compat: without ``root``, the key is the basename (suitable
    for inputs where filenames are unique by construction)."""
    fa = tmp_path / "file.fasta.gz"
    fa.write_bytes(b"x")
    hashes = compute_sha256s([fa])
    assert "file.fasta.gz" in hashes


def test_compute_sha256s_handles_mixed_input(tmp_path: Path) -> None:
    fa = tmp_path / "good.fasta.gz"
    fa.write_bytes(b"data")
    pq = tmp_path / "skipme.parquet"
    pq.write_bytes(b"data")
    cut = tmp_path / "x.cutadapt.json"  # intermediate JSON we did NOT name as an output
    cut.write_text("{}", encoding="utf-8")

    hashes = compute_sha256s([fa, pq, cut])
    assert "good.fasta.gz" in hashes
    assert "skipme.parquet" not in hashes
    # x.cutadapt.json ends with .json so the basename-extension rule
    # currently includes it; this is intentional - the runner cleans up
    # cutadapt's intermediate JSONs (see extract/trim.py) so it never
    # reaches compute_sha256s in real flows.
    assert "x.cutadapt.json" in hashes
