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

### Known v0.1 follow-ups (documented Codex peer-review findings)

These were flagged during the Phase 0/1 peer-review and are **not** blocking Phase 2; capturing them so they don't get lost.

- **`count.counter` still trims raw FASTQs inline.** The Phase 1 port preserved the original `selex_corpus` behavior — `count_round()` can run cutadapt on raw FASTQ inputs. The plan's final shape is *`extract` produces primer-stripped FASTAs and `count` only counts those*. The separation lands when the `extract` step is wired in Phase 2/3; until then, `count.counter` is dual-purpose.
- **`qc.readiness` requires clusters / enrich parquets.** The module is a faithful port and still expects `round_*.clusters.parquet`, `enrich_*.parquet`, `summary.json`, and `cluster_stats.json` — artifacts that v0.1 does *not* produce (clustering / enrichment are out of v0.1 scope). It remains exposed as a library API so the thesis pipeline can use it, but it is **not** wired into the `selexprep qc` CLI verb. The v0.1 `qc` verb will get a thinner, manifest-driven implementation when `extract`/`count` are fully separated.
- **Mocked-HTTP coverage gap.** The nine network adapters in `fetch.discover` and `download_srr_*` paths beyond ENA-direct don't yet have offline mocked tests (only their parsing helpers + dispatcher + SeedAdapter are covered). To be addressed before PyPI release.

### Phase 1.5 — discovery catalog (new)

The biggest "where do I even start?" UX gap in v0.1 was: a user installs the
package, knows nothing about which public SELEX accessions exist, and is
expected to run a multi-API discovery scan before anything useful happens.
Phase 1.5 fills this gap by shipping a **bundled discovery catalog** as
package data.

- **`selexprep.catalog`** — new subpackage with `load_catalog()`,
  `filter_catalog()` (target / organism / source / min-rounds / INSDC-only),
  and a `catalog_version()` snapshot identifier.
- **`selexprep catalog list/show/version`** — new Typer subapp wired into the
  root CLI. `list` supports the same filters as the Python API; `show
  <accession>` prints the full row including the study abstract.
- **Catalog content (v0.1 snapshot):** 219 public SELEX bioprojects with
  bioproject_id / source / study_title / protein_target / target_organism /
  paper_doi / paper_pmid / n_rounds_declared / abstract. Thesis-specific
  columns (`include`, `manual_curation_notes`,
  `library_type_verification`, `library_type_evidence`, `has_processed_counts`)
  are intentionally stripped so the catalog reflects the public archives,
  not any single researcher's curation.
- **v0.2 plan:** enrich each catalog row with the inferred `LibraryReport`
  (primer pair, N-region length, extraction_mode, confidence) once the full
  Phase 2 pipeline runs end-to-end. The enriched catalog will be deposited
  to Zenodo with a DOI and unlocks the NAR Database Issue paper venue
  alongside the planned Bioinformatics Advances Application Note.

### Improvements landed alongside the Phase 0/1 wrap

- `selexprep._io` — new module with `open_gzip_text_deterministic()` and `sha256_file()`. All `.gz` writes now produce bit-identical bytes across reruns (gzip header `mtime=0`), making manifest SHA256 hashes reproducible.
- `extract.demux` switched to the deterministic gzip writer; new regression test (`test_demux_output_is_deterministic_across_reruns`) compares byte-for-byte.
- `download_srr(backend="auto|ena|kingfisher|sra")` — explicit backend selector. `auto` (default) preserves the ENA-first dispatch; `ena` is the pure-MIT path with no fallback; `kingfisher` and `sra` force a specific backend. Any kingfisher invocation now logs an explicit GPL-3.0 notice.
- `fetch.discover._classify_all` writes a `NOT_ASSESSED_V0_1` sentinel + evidence JSON when the v0.2 classifier is absent. Distinguishes "deferred to v0.2" from an empty-string verdict downstream callers might silently treat as success.
- `qc.readiness.review_bioproject(tag=None)` now WARNs at call time when callers omit the tag (silently defaulting to `"untagged"` is the most likely call-site bug).
- `numpy` added as an explicit core dependency (it was being imported directly while only present transitively via `pandas`).

## [0.1.0] — TBD

First public release. See implementation plan for the locked v0.1 feature set.
