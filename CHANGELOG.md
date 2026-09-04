# Changelog

All notable changes to `selexprep` are documented here.

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versioning: [SemVer](https://semver.org/).

## [Unreleased]

_No unreleased changes yet._

## [0.4.1] - 2026-09-04

### Fixed

- **The curated metadata now reads in English.** The 47 adjudication notes -- the
  record of how each disagreement between the two independent extractions was
  settled -- were written in Italian, as was one `adjudication_rule` value and,
  more consequentially, one `n_random` value exposed by both the JSON and the
  flat CSV. They ship inside the wheel and the Zenodo archive, so a reader could
  not check an entry against its source, which is the whole point of shipping
  the provenance. Translated literally: accessions, field names, quoted evidence
  and method names are carried over verbatim, and the data itself is unchanged --
  240 records, 1920 cells, 1206 concordant / 236 single-source / 47 adjudicated /
  2 verified / 429 not stated, with both extraction arms retained on every
  adjudicated cell. `METADATA_VERSION` is bumped to
  `v0.3.2-dual-extraction-adjudicated-en-2026-09-04`.

## [0.4.0] - 2026-08-10

### Fixed

- **`detect` no longer reports a library where there is none.** When one
  sequence dominates a pool — a monoclonal small-RNA library, an adapter-dimer
  pool, an already-collapsed deposit — every position is conserved, so the
  positional-consensus walk finds no constant/random boundary, consumes the
  whole read from both ends, and infers a random region of zero nucleotides.
  `detect` used to emit that as a healthy two-sided library
  (`full_insert_recovered=True`, `required_action="NONE"`), which carried `run`
  through to `count` and wrote a table of zero-length sequences without raising
  anything. It now returns `UNABLE_TO_INFER` / `MANUAL_PRIMERS_REQUIRED` with a
  `failure_reason` that explains the cause and points at
  `--override-primer-5p/-3p`. Paired split-primer deposits are unaffected: they
  report `n_length_mode=None`, not `0`, and keep asking for read merging.
  Found by the benchmark's adapter-control arm on PRJDB7022.

  The refusal deliberately claims only that no variable region is present, not
  that the input is not SELEX: a late-round pool that has converged onto a
  single winner has a randomised region as conserved as its primers and
  produces the same signal, and `detect` cannot separate the two from reads
  alone. `MANUAL_PRIMERS_REQUIRED` is the right instruction in both cases —
  the user who knows the construct supplies the primers and extraction
  proceeds.
- **Outer-edge core rescue in flank detection.** Constant technical sequence
  sitting between the read edge and the library constant (a truncated
  sequencing adapter, an index remnant) was absorbed into the called flank,
  collapsing the whole-primer match rate below the primer-found threshold and
  downgrading a two-sided library to one-sided extraction. `detect` now retries
  with the well-supported core of the flank when — and only when — the
  full-length match rate has already failed. The trim span is unchanged, so the
  random-region boundary cannot move. On PRJEB62495 this turns
  `FIVE_PRIME_ONLY` (65 nt output) into `BOTH_PRIMERS_SINGLE_READ` (40 nt, the
  published length).

### Changed

- **Tier-1 benchmark rebalanced to three arms of seven** (21 source-verified
  deposits, up from 11): four pre-trimmed deposits added to the specificity arm
  and six non-SELEX small-RNA deposits to the adapter-control arm, both selected
  from archive metadata before any inference was run. PRJDB19098 (ground truth
  from patent WO2020204151) joins the recovery arm; PRJNA883192 leaves the
  scored set because its 3′ constant is resolvable only from the deposited
  reads, which makes it circular to score.
- **Curated metadata: the 47 cells where the two independent extractions
  disagreed are now adjudicated** (`v0.3.1-dual-extraction-adjudicated`), each
  carrying the resolved value, the rule applied, and the reasoning, with both
  arms preserved on record.

## [0.3.0] - 2026-07-03

### Changed

- **Development status is now Beta** (`Development Status :: 4 - Beta`). The core
  accession/local-FASTQ preprocessing pipeline is feature-complete and CI-tested,
  and the CLI commands and primary output schemas are treated as stable within
  the 0.x series.
- README status reads "beta" and points to the new stability policy; refreshed
  stale catalog/version wording.

### Added

- **`STABILITY.md`** — a stability policy declaring the stable public surface
  (the `inspect`/`fetch`/`detect`/`extract`/`count`/`qc`/`run`/`catalog` commands;
  the `library_report.json`, `selexprep_manifest.json`, `counts.parquet`, and
  `rounds.tsv` schemas; the enumerations; and the determinism guarantee) versus
  the experimental / not-yet-implemented features.
- **Schema-stability regression tests** (`tests/test_schema_stability.py`) that
  pin the `LibraryReport` and `SelexprepManifestV1` field sets, the enumerations,
  and the `counts.parquet` columns, so a breaking change to a public data
  contract fails CI and forces a deliberate schema-version bump.

## [0.2.1] - 2026-07-03

### Changed

- Refreshed the discovery catalog and curated metadata layer from **238 to 240
  deposits**. Two newly-deposited public ENA SELEX studies were curated in (each
  by the same two-independent-extraction method): `PRJEB114397` (an aptamer
  selection against perfluorooctanoic acid) and `PRJNA1481083` (automated SELEX
  against 96 protein targets). Two further new deposits were classified out of
  scope and recorded in the exclusion sidecar: `PRJEB88669` (genomic
  Helicase-SELEX) and `PRJNA860038` (a transcription-factor binding-motif
  SELEX-seq).
- Bumped `METADATA_VERSION` and `CATALOG_VERSION` to the `2026-07-03` snapshot.

### Fixed

- `catalog_version()` no longer returns a stale `v0.1.7-snapshot-2026-05-28`
  identifier — it had not been bumped through the catalog's 250 → 238 → 240
  changes.

## [0.2.0] - 2026-07-02

### Added

- **Curated metadata layer** (`selexprep.catalog.metadata`) — the annotated
  layer anticipated in v0.1. Each of the 238 bundled SELEX deposits now ships
  with curated experimental metadata: `study_type`, `target`, `target_class`,
  `chemistry`, `n_random`, `n_rounds`, `selection_format`, `counter_selection`.
  Built by **two independent LLM extractions** (Claude + Codex/GPT), reconciled,
  with per-value provenance (evidence quote + source + location). Where the two
  extractions genuinely disagreed, **both are kept** rather than silently
  picking one. This raises experimental-field coverage from ~4 to ~1479 filled,
  source-cited cells versus the discovery catalog alone.
- **New public API**: `load_metadata()` (flat `DataFrame`),
  `load_metadata_records()` (provenance-rich list of dicts), `metadata_version()`.
  Bundled as `curated_metadata.json` (canonical) + `curated_metadata.csv`
  (flat view). The extraction contract, both raw arms, and the reconciliation
  method live under `benchmarks/dual_extraction/` (not shipped on PyPI).

## [0.1.1] - 2026-06-16

### Fixed

- **cutadapt discovery**: `extract`, `count`, and the manifest's version
  capture now locate cutadapt next to the running Python interpreter when it is
  not on `$PATH`. This fixes "cutadapt not found on PATH" under `pipx install`
  (which exposes only selexprep's own entry point), absolute-path invocation,
  or a workflow runner with a sanitized PATH — cases where cutadapt is installed
  alongside selexprep but the environment isn't "activated".

## [0.1.0] - 2026-06-13

First public release: accession-first preprocessing for high-throughput SELEX
(HT-SELEX) sequencing deposits, with automatic primer / constant-region
inference. Give it an INSDC accession (ENA / SRA / DDBJ) and it fetches the
runs, infers the library's flanking constants from the reads, and extracts the
random regions — no manual primer entry required.

### Added

- **Command-line interface** (`selexprep <verb>`):
  - `inspect` — summarize an accession's runs and metadata.
  - `fetch` — download FASTQs for an accession (ENA-direct by default).
  - `detect` — infer the 5′/3′ constant regions (primers) from the reads.
  - `extract` — strip the inferred constants and emit the random-region reads.
  - `count` — collapse extracted reads to unique-sequence counts per round.
  - `qc` — quality-control flags and plots.
  - `run` — end-to-end fetch → detect → extract → count, with `--resume`.
  - `catalog` — browse the bundled discovery catalog of SELEX deposits.
- **Primer / constant-region inference** (`detect`): position-anchored
  consensus over the read pool, with a typed `LibraryReport` describing the
  inferred 5′/3′ constants, random-region length, match rates, read state
  (raw vs. pre-trimmed), and a confidence-graded status. Cross-round inference
  reconciles evidence across selection rounds.
- **Adapter awareness**: known Illumina sequencing adapters (e.g. TruSeq) are
  recorded where present and excluded from primer candidates, so adapter
  read-through is reported as diagnostic information rather than mistaken for a
  library constant.
- **`extract`**: paired-end handling, strand-orientation detection, and
  per-mode adapter handling. cutadapt is invoked as a subprocess (its CLI is
  the stable contract). Supports `--override-primer-{5p,3p}` to bypass
  inference, and a rebuild path for manually corrected primers.
- **Discovery catalog** (`catalog`): a bundled, refreshable index of public
  SELEX deposits built from INSDC `library_strategy="SELEX"` queries, with
  per-run and per-BioProject strategy filtering and manual exclusions for
  mislabelled deposits.
- **Quality control** (`qc`): diversity / rarefaction helpers and depth-aware
  flags (e.g. unexpected rarefied-diversity increase, modal-length spread,
  orientation skew, low read depth, adapter contamination).
- **Deterministic outputs**: all gzip writes are byte-identical across reruns
  (gzip header `mtime=0`) and JSON is written with sorted keys, so a
  `SelexprepManifestV1` run manifest carries reproducible `sha256` hashes.
- **Benchmark suite** (under `benchmarks/`, not shipped on PyPI): a Tier 1
  primer-recovery benchmark against paper-grounded ground truth and a Tier 2
  corpus-audit pipeline over the discovery catalog.

### Known limitations (v0.2 carry-forward)

- **Multiplexed (inline-barcoded) deposits** need a user-supplied sample sheet;
  automatic demultiplex detection is deferred.
- **Read merging** of overlapping mates is not implemented.
- **`qc.readiness`** (clustering / enrichment review) is a faithful library
  API but expects clustering artifacts that v0.1 does not produce; it is not
  wired into the `qc` verb.
- **`count.counter`** can still trim raw FASTQs inline; the clean split
  (`extract` strips, `count` only counts) is partial.
- **`--from-pretrimmed-fastq`** validates record completeness but not per-line
  FASTQ conformance — adequate for the power-user opt-in.
- **Network coverage**: the non-ENA fetch backends still lack offline mocked
  tests (carried into v0.2).

### Packaging

- MIT-licensed. A default `pip install selexprep` pulls only MIT-compatible
  dependencies (pydantic v2, Typer, pandas, numpy); cutadapt is invoked as a
  subprocess.
- `kingfisher` (GPL-3.0) is an optional, runtime-detected subprocess backend —
  not a declared dependency — so the default install stays MIT-only.

[Unreleased]: https://github.com/marcorotanegroni/selexprep/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/marcorotanegroni/selexprep/compare/v0.2.1...v0.3.0
[0.2.1]: https://github.com/marcorotanegroni/selexprep/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/marcorotanegroni/selexprep/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/marcorotanegroni/selexprep/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/marcorotanegroni/selexprep/releases/tag/v0.1.0
