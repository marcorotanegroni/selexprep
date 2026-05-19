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

### Phase 2 — LibraryReport schema + cross-round inference (2026-05-19)

Adds the `LibraryReport` pydantic schema and the cross-round primer
inference pipeline that turns Phase 1's single-pool flank detector into
the typed contract every downstream stage (`extract`, `count`, `qc`,
`manifest`) consumes.

- **NEW: `selexprep.library.report`** — `LibraryReport(BaseModel)`
  with the locked schema (plan lines 233-285), `Literal` aliases for the
  five categorical fields (`ExtractionMode`, `ReadSource`,
  `RequiredAction`, `Orientation`, `Status`), the `_classify` pure
  function implementing the locked decision table (plan lines 300-309)
  with the no-round-map status cap (line 289), and deterministic JSON
  I/O (`write_library_report_json` / `read_library_report_json` —
  bit-identical output across reruns, numeric ordering for int-keyed
  dicts). Strict-mypy clean (pydantic plugin enabled).
- **NEW: `selexprep.library.adapters`** — conservative v0.1 blacklist
  (TruSeq R1 + Nextera) with auto-computed reverse complements,
  `reverse_complement()` helper (rejects IUPAC ambiguity in v0.1), and
  `count_adapter_hits()` substring scanner (records hits; does NOT
  filter reads).
- **EXTENDED: `selexprep.library.detect`** — `compute_library_report()`
  orchestrator. Cross-round persistence as `1 - clip(stdev/mean, 0, 1)`,
  position consistency with ±2 nt tolerance, U→T normalization of RNA
  primers, paired-end split detection (R1 5' + R2 5' = revcomp(3'
  primer)), MIXED/FORWARD/REVERSE orientation diagnostic, composite
  confidence via weighted sum (two regimes: with vs without round map).
  Phase 1 functions (`detect_flank`, `detect_primers`,
  `detect_from_parquet`, `earliest_round_parquet`) are unchanged and
  remain the algorithmic primitives the new orchestrator consumes.
- **WIRED: `selexprep detect`** — CLI command parses `--round-map` TSV
  (columns `file<TAB>round_number`), groups FASTQs by round, runs
  `compute_library_report`, writes `library_report.json` to `--outdir`.
  Refuses to run without `--round-map` (cross-round persistence is a
  core inference signal). Single-end only in Phase 2; paired-end via
  `compute_library_report`'s `paired_mate_streams` kwarg awaits CLI
  surface in Phase 3.
- **INFRASTRUCTURE:** `pyproject.toml` declares
  `plugins = ["pydantic.mypy"]` under `[tool.mypy]` (required for
  strict-mypy to resolve `BaseModel` field types — without it every
  field decays to `Any`).
- **Tests added (32 new, 234 + 1 xfailed total):** `tests/test_adapters.py`
  (14 tests covering revcomp + blacklist + substring scan);
  `tests/test_report.py` (19 tests covering every row of the locked
  classification table plus edge cases — status cap, adapter
  demotion, orientation, U→T, deterministic serialization, schema
  immutability, numeric int-key ordering, empty input, sub-floor
  input); `tests/test_cli.py` (3 new tests for `detect` CLI:
  missing-round-map, end-to-end round-trip, FASTQ-not-in-map).

#### Calibration status — placeholder pending Codex peer review

Codex usage was rate-limited 2026-05-19 → 2026-05-26 when Phase 2
shipped, so calibration numbers ship as locked-plan literals (match
rates `> 0.7`, n_length confidence `> 0.8`, UNABLE floor `< 0.4`) plus
placeholders for everything the locked plan does not pin down
(composite weights, status cutoffs, position-consistency tolerance,
adapter list exact composition, persistence formula). Every numeric
placeholder carries a `# CALIBRATION-TODO` comment naming the locked
plan line (where applicable) or "not in locked plan - placeholder
pending Codex" otherwise.

Test discipline: **all tests assert on behavior, never on threshold
values** — e.g. `assert report.extraction_mode == "BOTH_PRIMERS_SINGLE_READ"`,
never `assert HIGH_CONFIDENCE_CUTOFF == 0.80`. So when Codex tuning
lands (or Phase 6 benchmark numbers update the constants), the test
suite stays green by construction.

Recovery list: `grep -rn "CALIBRATION-TODO" src/` returns the full
inventory.

### Phase 1.5.1 — catalog refresh against broad ENA queries

A Codex / sanity-check pass after Phase 1.5 revealed that the initial
bundled catalog (219 rows, sourced from the thesis-specific
`selex_corpus.discover` run) under-counted INSDC studies by ~50%: 94
INSDC accessions vs ~120 unique studies surfaced by broader keyword
queries against ENA. The thesis queries combined keywords with AND
clauses to maximize precision; a generic-tool catalog wants the
broader OR-style net.

- **NEW: `selexprep.catalog.rebuild`** — reproducible refresh script
  that runs 13 broad ENA queries (HT-SELEX, SELEX-seq, SELEX, aptamer,
  Cell-SELEX, RNA aptamer, DNA aptamer, systematic evolution, SELEX
  rounds, …), unions the studies, merges hand-enriched fields
  (`protein_target` / `paper_doi` / `paper_pmid` /
  `n_rounds_declared`) forward from the previous catalog when an
  accession is still upstream, and carries non-INSDC deposits
  (zenodo/figshare/utexas processed-data entries) across refreshes
  unchanged.
- **NEW: `selexprep catalog refresh [--out PATH --no-preserve-enrichment]`**
  — CLI verb that wraps `rebuild_catalog`. Lets users (or CI) refresh
  the catalog on demand without touching the package source.
- **Catalog refreshed in-place** for v0.1.5: 273 bioprojects (148
  ENA-discovered INSDC studies + 125 carried-forward
  Zenodo/Figshare deposits). Snapshot bumped to
  `v0.1.5-snapshot-2026-05-19`. The 4 seed entries (Hoinka IL-10RA,
  Dao CCR7, …) keep their hand-curated enrichment.
- **No curation flags.** Confirmed with the PI: the package never
  ships `include` / `manual_curation_notes` columns. Curation is the
  user's downstream job; the catalog reflects the public archives
  only.

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
