# Changelog

All notable changes to `selexprep` are documented here.

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versioning: [SemVer](https://semver.org/).

## [Unreleased]

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
  tests (to be closed before the PyPI release).

### Packaging

- MIT-licensed. A default `pip install selexprep` pulls only MIT-compatible
  dependencies (pydantic v2, Typer, pandas, numpy); cutadapt is invoked as a
  subprocess.
- `kingfisher` (GPL-3.0) is an optional, runtime-detected subprocess backend —
  not a declared dependency — so the default install stays MIT-only.
