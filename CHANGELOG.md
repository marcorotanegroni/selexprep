# Changelog

All notable changes to `selexprep` are documented here.

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versioning: [SemVer](https://semver.org/).

## [Unreleased]

### Fixed

- **Catalogue: PRJEB114397 `target_class` is `small molecule`.** The layer held
  the source's descriptive phrase ("toxic, bioaccumulative chemical") instead of
  a class. Both extractions agree in substance (Claude quoted the phrase, Codex
  read "chemical"); PFOA is a small molecule, the class the layer uses for
  chemical targets. The cell stays concordant and carries a note; no count in
  the layer changes.
- **Catalogue provenance: both raw arms for the two deposits added in 0.2.1.**
  PRJEB114397 and PRJNA1481083 were extracted by both arms on 2026-07-03, but the
  raw outputs were not in the repository and their 12 concordant cells lacked
  `codex_concurs`. The outputs are now in `benchmarks/dual_extraction/`
  (`*_extractions_new4.json`, unchanged), and the 12 cells carry the value the
  second arm read; every concordant cell now has it.

## [0.4.5] - 2026-09-30

### Added

- **Benchmark: round-assignment coverage of the release**
  (`benchmarks/round_validation/coverage.{py,sh}`). For every run of the 127
  INSDC deposits in the catalogue it records the round `fetch` assigns, with its
  confidence, source field and notes, and reports how many runs receive a round
  and how many are left unassigned, and the agreement with the curated round
  maps on the runs both cover. Descriptive and in-sample, as fixed in
  `benchmarks/round_validation/PROTOCOL.md` (Amendment 1): an assignment is never
  counted as correct.
- **Benchmark: the deposit excluded after inference is scored as a sensitivity
  analysis.** PRJNA1244358 was excluded from the recovery arm after `detect`
  had refused it (v0.4.0), which `benchmarks/README.md` now discloses. It is back
  in `ground_truth.tsv` with `read_state` `excluded_after_inference`: out of
  every arm and every distribution, and scored only in a new
  `excluded_after_inference` block of `metrics.json`, which recomputes the
  recovery arm with it counted in, whatever its result with this version.
  `table_1.md` shows it in a separate sensitivity section, never as a control.

- **`detect` and `run` say what `status` means.** Whenever a status is shown
  next to inferred primers, a note says that it rates how well the reads
  support them, not a verified boundary, and points to the new "Known limits"
  section of the LibraryReport reference. The semi-synthetic benchmark found
  the case: one sequence making up 80% of every round, the earliest included,
  is taken for constant with status HIGH, and the call reaches 1 nt to almost
  the whole random region into it (19 of 30 test runs). The inference is
  unchanged; the README no longer promises "no silent miscalls".
- **Semi-synthetic benchmark** (`benchmarks/synthetic/`): a pre-registered
  design, a donor screen, a generator that builds libraries around real
  random regions with the truth in every read name, an evaluator independent
  of `selexprep`, and development and test results.

- **`fetch` now says when one round label covers runs of different
  BioSamples.** Deposits with parallel selections often number each one's
  rounds from 1 (`SELEX S1 Round 04`, `SELEX S2 Round 04`), so both runs get
  round 4 and `count` merges them into one pool. That is right for replicates
  of one selection and wrong for separate selections, and the metadata cannot
  tell which. The assignment is unchanged; `fetch` logs a warning naming the
  rounds and runs, and `fetch_metadata.json` records them under
  `rounds_with_several_samples`. Lanes of one BioSample are not flagged.
- **Benchmark: all-read verification of every deposit's read state**
  (`benchmarks/verify_read_state.{py,sh}`, results in
  `benchmarks/read_state_full.tsv`). It streams every read of every fetched
  FASTQ and reports, per file, the length distribution, where the published
  constants sit, adapter probes and a positional conservation profile.
  `read_state_evidence.tsv` is now built from it instead of spot-checks of the
  leading reads.
- **Benchmark: per-run extraction yield for the recovery arm.** `detect`
  infers the constants from the earliest round; the Tier-1 workflow now also
  runs `extract` on every round and writes `results/extraction_yield.tsv`
  (reads in / kept per input FASTQ, modal extracted length per round), so a
  run the called constants do not fit is visible instead of hidden behind a
  correct random-region length.

### Fixed

- **`fetch` reads the round from the sample's attributes, and no longer reads
  ArrayExpress sample numbers as rounds.** The round parser has always ranked
  structured sample attributes first, but `fetch` never passed any: the ENA
  filereport carries none. It now reads them from ENA's sample records, where
  ArrayExpress writes the SDRF's `selex cycle` and GEO its `round`. An
  attribute key names a round when it has one round word (round, cycle,
  iteration, r) and its other words come from a fixed vocabulary (selex,
  selection, enrichment, of, number, pcr); this covers the spellings found in
  the catalogue (`selex cycle`, `Characteristics[selex cycle]`,
  `Selection cycle` with values such as `Cycle3`, `selex round`,
  `round of selex enrichment`, `round_selection`) and keeps out keys such as
  `cell cycle` or `cycle threshold`. `pcr cycle` is read at MEDIUM confidence:
  PRJEB38961 uses it for the selection cycle, but it can also mean
  amplification cycles. The pattern that read `RV01`…`RV39` as rounds 1–39 is
  removed: in PRJEB51212, PRJEB51473 and PRJEB38961 those are sample numbers,
  and the attributes give cycles 0, 1, 3 and 6. An attribute outranks a title
  that reads as another round, and the conflict is named in the parser notes.
  Every round attribute is weighed before deciding, so the result does not
  depend on the order the archive lists them: unambiguous attributes that
  disagree leave the run unassigned, and `pcr cycle` never overrides a round
  that another attribute or a text field gives.
  The attributes are kept in `fetch_metadata.json`; if the sample records
  cannot be fetched, `fetch` warns and parses the text fields as before.
  Across the 25,435 runs of the 127 INSDC deposits in the catalogue: 117 runs
  in those three deposits move from rounds 1–39 to their documented cycles;
  4,105 runs that were unassigned get a round from their attributes, 4,080 of
  them in PRJEB76622 and PRJEB61115, whose titles were ambiguous (`CPXCR1`, a
  replicate suffix `_R1`); no run loses its round. Against the ten
  hand-curated round maps in `benchmarks/round_maps/`, no run receives a round
  that differs from the map, and PRJEB49150 goes from 0 to 9 of 9 runs
  assigned.
- **`detect` reports the library constant every round shares, not the earliest
  round's whole flank.** The flank is called on the earliest round from the read
  edge inward, so it took in whatever that round's reads carried outside the
  library constant. Where that sequence changes between runs, the call did not
  fit the other runs: on PRJNA809588, whose runs each carry their own inline
  tags, `extract` kept under 1% of the reads of seven of its ten runs, while
  `detect` and `extract` both reported success. `detect` now aligns every round
  on the flank's inner core (the 12 nt next to the random region) and keeps an
  outer position only if every round agrees on its base with at least 55%
  support, the floor the boundary search already uses for "constant". The inner
  boundary cannot move, so the random-region length is unchanged; position
  consistency, persistence and variants are measured where each round actually
  carries the constant. The same rule drops a heterogeneous base ahead of the
  construct when there is only one round (PRJNA1395820). Found by the
  benchmark's per-run extraction yield; it changes the Tier-1 calls for
  PRJNA809588, PRJEB62495 and PRJNA1395820. A read in which the core occurs
  more than once — a copy in the outer sequence, or a copy inside the random
  region of a later round — is anchored where the round's single-copy reads
  agree, or left out when they do not; a trimmed constant that would occur
  twice in the flank is not used. Found in code review, before any benchmark
  used it: anchoring on the outermost copy placed the constant on an outer
  repeat, and anchoring on the innermost copy let a copy of the core inside the
  random region shorten the constant to the bare core, after which `extract`
  cut inside the random region.
- **`extract` warns when an input's yield collapses.** An input FASTQ keeping
  less than half of its reads, or less than half the yield of the best input,
  is named in a warning with its counts, instead of passing silently.

### Changed

- **PRJEB70964 moved from the adapter-control arm to the specificity arm**
  (arms are now 7 recovery / 8 specificity / 6 adapter-control). It had been
  filed as an adapter-collision control on a spot-check recording 81-nt reads
  with both constants; the all-read check shows every read is the N35 region
  alone. Its outcome (no call) is unchanged. The adapter-control description in
  the benchmark README, `metrics.py` and the library-report example notebook is
  corrected to what those deposits actually carry.

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

[Unreleased]: https://github.com/marcorotanegroni/selexprep/compare/v0.4.5...HEAD
[0.4.5]: https://github.com/marcorotanegroni/selexprep/compare/v0.4.1...v0.4.5
[0.4.1]: https://github.com/marcorotanegroni/selexprep/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/marcorotanegroni/selexprep/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/marcorotanegroni/selexprep/compare/v0.2.1...v0.3.0
[0.2.1]: https://github.com/marcorotanegroni/selexprep/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/marcorotanegroni/selexprep/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/marcorotanegroni/selexprep/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/marcorotanegroni/selexprep/releases/tag/v0.1.0
