# Stability policy

`selexprep` is **Beta**. The core accession/local-FASTQ preprocessing pipeline is
feature-complete and covered by continuous integration, and the public surfaces
listed below are treated as **stable within the 0.x series**: they will not be
removed or changed incompatibly without a **minor-version bump** and a
`CHANGELOG` entry. Additive changes (a new optional field, a new CLI flag) may
land in a patch release; a **new enumeration value is a minor-version change**,
because consumers may be matching the value set exhaustively (see *Stable*
below).

The stable data contracts are enforced in code — `LibraryReport` and
`SelexprepManifestV1` are frozen Pydantic models with `extra="forbid"`, and
`tests/test_schema_stability.py` pins their exact field sets, the enumerations,
and the `counts.parquet` columns, so a breaking change fails CI rather than
shipping silently.

## Stable

- **CLI commands** and their core behaviour: `inspect`, `fetch`, `detect`,
  `extract`, `count`, `qc`, `run`, `catalog`.
- **`library_report.json`** — the `LibraryReport` schema: `primer_5p`/`primer_3p`,
  `extraction_mode`, `read_source`, `required_action`, `orientation`, `status`,
  `confidence`, the `n_length_*` and `match_rate_*` fields, `sampling_seed`, and
  `failure_reason` (full field set in `library/report.py`).
- **`selexprep_manifest.json`** — the `SelexprepManifestV1` schema, carrying the
  literal `manifest_version = "selexprep_manifest_v1"` so downstream tooling can
  detect the revision; tool/dependency versions, provenance, `input_sha256` /
  `output_sha256`, the nested `library_report`, run parameters, and `sampling_seed`.
- **`counts.parquet`** — columns `sequence`, `reads`, `rank`, `rpm`.
- **`run_summary.tsv`** — the corpus-level output of `run`; columns `accession`,
  `status`, `last_stage_completed`, `library_report_status`, `extraction_mode`,
  `required_action`, `confidence`, `flags_raised`, `notes`.
- **`rounds.tsv`** — columns `file`, `round_number`.
- **Curated metadata layer** — the `load_metadata` / `load_metadata_records` API
  and its eight fields (`study_type`, `target`, `target_class`, `chemistry`,
  `n_random`, `n_rounds`, `selection_format`, `counter_selection`) plus the
  per-cell provenance/status structure (pinned by `tests/test_catalog_metadata.py`).
- **Enumerations** — the value sets of `status` (HIGH/MEDIUM/LOW/UNABLE_TO_INFER),
  `extraction_mode`, `read_source`, `required_action`, and `orientation`. New
  values are a minor-version change.
- **Determinism** — for a fixed `--sampling-seed`, FASTA/TSV/JSON outputs are
  byte-identical across runs; the manifest's `output_sha256` is the audit anchor.

## Experimental / not yet implemented

These are documented boundaries of the current version — the operating envelope,
not hidden gaps. They are **not** a commitment that every item will be
implemented; some are deliberate scope boundaries (handled by downstream tools),
others possible future directions:

- **Paired-end read merging** — split-primer rounds are trimmed on each side and
  reported, but not stitched into a full-length insert, so `count`/`qc` are
  skipped for those rounds.
- **Multiplex auto-detection** — multiplexed deposits currently need a
  user-supplied `--sample-sheet`.
- **figshare/Zenodo fetch backends** — those deposits are discovery-only pointers
  in the catalog; only INSDC accessions are fetchable.
- **AnnData export**, **BibTeX auto-citation**, **library-type classification**.
- **Inference thresholds** may be retuned as the benchmark expands. Tests assert
  on *behaviour* (e.g. the returned `status`/`extraction_mode`), never on the
  numeric threshold values, so retuning stays within the stable contract.

## Versioning

`selexprep` follows semantic versioning within the 0.x series. A change to any
stable surface above requires a minor-version bump; the schema literals
(`manifest_version`) and this document are updated in the same change.
