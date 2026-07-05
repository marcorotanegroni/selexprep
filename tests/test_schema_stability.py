"""Schema-stability regression guards for the public data contracts.

selexprep is Beta: the ``library_report.json`` and ``selexprep_manifest.json``
schemas, their enumerations, and the ``counts.parquet`` columns are a stable
public surface (see ``STABILITY.md``). These tests pin the exact field/column
sets so that removing or renaming any of them fails CI, forcing a deliberate
schema-version bump + CHANGELOG entry rather than a silent breaking change.

Additive changes (a new optional field, a new enum value) are a minor-version
concern; if you add one on purpose, update the expected set here in the same
commit — that edit is the audit trail.
"""

from __future__ import annotations

from typing import get_args

from selexprep.library.report import (
    ExtractionMode,
    LibraryReport,
    Orientation,
    ReadSource,
    RequiredAction,
    Status,
)
from selexprep.manifest import SelexprepManifestV1
from selexprep.run.runner import _SUMMARY_COLUMNS

# --- library_report.json (LibraryReport) ------------------------------------

_LIBRARY_REPORT_FIELDS = {
    "primer_5p",
    "primer_3p",
    "variants_5p",
    "variants_3p",
    "known_adapter_hits",
    "extraction_mode",
    "full_insert_recovered",
    "read_source",
    "required_action",
    "orientation",
    "n_length_mode",
    "n_length_distribution",
    "n_length_confidence",
    "match_rate_5p",
    "match_rate_3p",
    "position_consistency_5p",
    "position_consistency_3p",
    "read_fraction_used_for_inference",
    "sampling_seed",
    "confidence",
    "status",
    "failure_reason",
}

# --- selexprep_manifest.json (SelexprepManifestV1) --------------------------

_MANIFEST_FIELDS = {
    "manifest_version",
    "selexprep_version",
    "python_version",
    "cutadapt_version",
    "dnaio_version",
    "pyarrow_version",
    "accession",
    "bioproject_id",
    "runs",
    "input_sha256",
    "output_sha256",
    "library_report",
    "extraction_mode",
    "read_source",
    "required_action",
    "full_insert_recovered",
    "parameters",
    "runtime_seconds_per_stage",
    "flags",
    "sampling_seed",
}

# --- stable enumerations ----------------------------------------------------

_ENUM_VALUES = {
    "Status": (Status, {"HIGH", "MEDIUM", "LOW", "UNABLE_TO_INFER"}),
    "ExtractionMode": (
        ExtractionMode,
        {
            "BOTH_PRIMERS_SINGLE_READ",
            "FIVE_PRIME_ONLY",
            "THREE_PRIME_ONLY",
            "PAIRED_END_SPLIT_PRIMERS",
            "UNABLE_TO_EXTRACT",
        },
    ),
    "ReadSource": (ReadSource, {"R1", "R2", "R1_AND_R2", "INTERLEAVED", "UNKNOWN"}),
    "RequiredAction": (
        RequiredAction,
        {"NONE", "MANUAL_PRIMERS_REQUIRED", "READ_MERGING_RECOMMENDED"},
    ),
    "Orientation": (Orientation, {"FORWARD", "REVERSE", "MIXED"}),
}

# ``counts.parquet`` columns (written by count.counter._counter_to_parquet).
_COUNTS_PARQUET_COLUMNS = ["sequence", "reads", "rank", "rpm"]

# ``run_summary.tsv`` columns — the corpus-level output of ``selexprep run``.
_RUN_SUMMARY_COLUMNS = (
    "accession",
    "status",
    "last_stage_completed",
    "library_report_status",
    "extraction_mode",
    "required_action",
    "confidence",
    "flags_raised",
    "notes",
)


def test_library_report_field_set_is_stable() -> None:
    assert set(LibraryReport.model_fields) == _LIBRARY_REPORT_FIELDS


def test_manifest_field_set_is_stable() -> None:
    assert set(SelexprepManifestV1.model_fields) == _MANIFEST_FIELDS


def test_manifest_version_tag_is_stable() -> None:
    # the literal that downstream tooling keys on to detect the schema revision
    assert get_args(SelexprepManifestV1.model_fields["manifest_version"].annotation) == (
        "selexprep_manifest_v1",
    )


def test_enum_value_sets_are_stable() -> None:
    for name, (alias, expected) in _ENUM_VALUES.items():
        assert set(get_args(alias)) == expected, name


def test_public_schemas_reject_unknown_fields() -> None:
    # the mechanism that makes the contract enforceable: frozen + extra="forbid"
    for model in (LibraryReport, SelexprepManifestV1):
        assert model.model_config.get("frozen") is True
        assert model.model_config.get("extra") == "forbid"


def test_counts_parquet_columns_are_documented() -> None:
    # guards the documented column contract against a silent rename in the docstring
    from selexprep.count import counter

    assert all(col in counter.__doc__ for col in _COUNTS_PARQUET_COLUMNS)


def test_run_summary_columns_are_stable() -> None:
    assert _SUMMARY_COLUMNS == _RUN_SUMMARY_COLUMNS
