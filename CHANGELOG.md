# Changelog

All notable changes to `selexprep` are documented here.

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versioning: [SemVer](https://semver.org/).

## [Unreleased]

### Added

**Phase 1 — library modules ported (166 tests, all green):**

- `selexprep._common` — shared utilities: `iter_srr_files` (exact-name match, no SRR1234↔SRR12345 collisions), `load_csv`, `parse_round_number`, `setup_logging`.
- `selexprep.count.counter` — FASTQ.gz → parquet sequence counting with anchored linked-adapter cutadapt trimming; paired-end pair-sync; pyfastx fast path with gzip fallback; Shannon-entropy + singleton-fraction pool stats; multi-target BioProject layout support.
- `selexprep.fetch.metadata` — deterministic 5-level cascade for round inference (sample_attributes → sample_title → library_name/experiment_title/design_description → abstract count → manual review). Never guesses; 20% unknowns preferred over silent miscalls.
- `selexprep.fetch.discover` — multi-source SELEX discovery across nine sources (Seed YAML, ENA, NCBI SRA, GEO, UTexas Aptamer DB, Zenodo, Figshare, Crossref, OpenAlex). Optional dependencies (pysradb, Bio.Entrez) soft-fail. Library-type classification is conditional (deferred to v0.2).
- `selexprep.fetch.download` — MIT-compatible-first download dispatcher: **ENA-direct (default)** → kingfisher (optional, GPL-3.0) → sra-toolkit (optional). Pure-`requests` ENA filereport API with Range-resumable streaming + MD5 verification; no external tools required for the default install.
- `selexprep.library.detect` — empirical primer/constant-region inference. **Default: scans every unique sequence in the parquet (no top-N subsampling)** so the long tail of rare unique sequences also confirms the consensus.
- `selexprep.library.audit` — pre/post-extraction structure audits: raw-FASTQ sampling with 3'-aligned positional base frequencies (variable-length-tolerant), trimmed-parquet length distribution + TruSeq R1 contamination probe. Pure dataclasses; no CLI.
- `selexprep.extract.demux` — sample-sheet-driven barcode demultiplexer for pooled multi-round SELEX runs. Validates barcodes for Hamming distance, keeps paired-end R1/R2 in lockstep, never trims R2.
- `selexprep.qc.coverage` — per-BioProject round-coverage classification (all_rounds_public / partial / multiplexed_unrecoverable / unknown). Filter is a configurable callable (no thesis-specific hard-coding).
- `selexprep.qc.consistency` — k-mer Jaccard distance + monotonicity check across rounds. Strictly diagnostic; never reassigns rounds from enrichment signal.
- `selexprep.qc.readiness` — eight-section sequence-level readiness review (`pre / alphabet / lengths / trim_seq / composition / diversity / selection / consistency`). Tag-aware composition/diversity thresholds; tag is a per-call parameter (not a hard-coded BP map).

**CLI:** seven stub subcommands (`inspect`, `fetch`, `detect`, `extract`, `count`, `qc`, `run`) ship in the Typer dispatcher; wiring to the ported library modules lands in Phase 2.

### Notes
- v0.1 packaging note: `discover.py` keeps its nine adapter classes in one file. Splitting into a `selexprep.fetch.sources.*` subpackage is a v0.2 cleanup.
- `kingfisher` (GPL-3.0) was dropped from `pyproject.toml` — it remains a runtime-detected optional subprocess backend so a default `pip install selexprep` stays MIT-only.

## [0.1.0] — TBD

First public release. See implementation plan for the locked v0.1 feature set.
