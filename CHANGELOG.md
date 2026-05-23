# Changelog

All notable changes to `selexprep` are documented here.

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versioning: [SemVer](https://semver.org/).

## [Unreleased]

### Fixed

**Phase 6b.3a HPC dependency-resolution hotfix (2026-05-23)** — surfaced
during the Phase 6b.4 HPC audit attempt; blocks
``uv sync --extra bench`` from producing a working Snakemake on a
fresh venv (including the DGX5 venv the user ran into).

- **`pyproject.toml` ``bench`` extra pinned.** Previously
  ``bench = ["snakemake >= 7.0"]`` left PuLP transitively unpinned, so
  a fresh ``uv sync --extra bench`` resolved PuLP 3.x; Snakemake 7
  calls ``pulp.list_solvers()`` at import time and PuLP 2.8 removed
  the snake_case alias (only ``listSolvers`` remains), so
  ``snakemake --version`` crashed with ``AttributeError``. Now pinned
  to ``snakemake >= 7.0, < 8`` (ceiling matches the existing comment
  about Python 3.10 compatibility — Snakemake 8+ requires Python
  3.11+) and ``pulp < 2.8``. ``uv lock`` refresh dropped the
  snakemake 9 dual-resolve and its ``snakemake-interface-*`` plugin
  family + sqlalchemy/sqlmodel/greenlet — net reduction in the
  ``bench``-extra footprint.
- **`benchmarks/Snakefile` Tier 1 f-string-wildcard pattern**
  replaced with plain string concatenation. The previous
  ``f"{OUTROOT}/{{accession}}/fetch_metadata.json"`` pattern (f-string
  substituting ``OUTROOT`` while the doubled ``{{accession}}``
  produces a literal Snakemake wildcard) parsed cleanly under
  Snakemake 9 but errors under Snakemake 7 with ``NameError: name
  'accession' is not defined``. Snakemake 7's Snakefile preprocessor
  walks AST nodes that Snakemake 9's doesn't. Refactored four
  occurrences to ``OUTROOT + "/{accession}/..."``; semantically
  identical, robust across versions. Both Snakefiles now dry-run
  cleanly with the pinned Snakemake 7.32.4. **No Tier 1 benchmark
  logic changed.** ``audit.smk`` was already wildcard-free in its
  f-strings so didn't need a patch.

CI gates (ruff / ruff-format / mypy / pytest 504+1 xfailed) and both
Snakefile dry-runs all green after the fix.

### Added

**Phase 6b.3a — Tier 2 corpus-audit scaffolding + Figure B pipeline (2026-05-23)**

Ships the code + Snakefile + tests for the Tier 2 corpus audit. This
is the second figure in the paper's two-figure story:

- **Figure A** (Tier 1, shipped in 6b.1+6b.2): primer-recovery
  validation on N=11 curated source-backed accessions — selexprep
  recovers paper-reported primers exactly / partially / safely
  refuses.
- **Figure B** (Tier 2, *pipeline* shipped here, *results* in 6b.4):
  corpus characterization on N≈30 sampled INSDC accessions, no
  per-row ground truth — fetchability, `LibraryReport.status`
  distribution, `extraction_mode` + `required_action` distributions,
  inference-stage safe-failure rate.

Framing: Tier 2 is an **audit / corpus characterization**, not a
validation in the strict sense — there is no per-row ground truth.
This is selexprep's unique claim: comparator tools (AptaPLEX,
EasyDIVER+) cannot produce this kind of audit by construction —
they require known primers as input, so they cannot be pointed at a
random catalog accession and asked "what does the LibraryReport
say?".

**Scaffolding vs results split.** 6b.3a is deliberately scoped to
the pipeline only: code + Snakefile + tests, reviewable + runnable.
The real-data HPC run that produces measured `audit_metrics.json` +
`audit_accessions.tsv` + `figure_b.{pdf,png}` is the **6b.4
follow-up commit** with no code changes — just the Snakefile output
artifacts after a manual HPC run + curator review. CI does NOT
execute the Snakefile.

**Methodological correction folded in (Codex peer-review + user,
5 rounds of revisions).** `inference_safe_failure_rate` is computed
**only** among rows with a LibraryReport (denominator =
`n_with_library_report`), NEVER mixed with fetch failures. Mixing
the two would inflate the rate with ENA / network / regional-
restriction problems — conflating "selexprep refused" (a feature)
with "the dataset was unreachable" (an external problem). Each
Figure B panel labels its denominator in the subtitle so a reviewer
never has to guess the normalization.

**New modules:**

- **`selexprep.benchmark.corpus_audit`** —
  - `sample_corpus(n, *, sources=None, exclude=(), seed=42)` —
    deterministic uniform-random sample from
    `filter_catalog(insdc_only=True)`. Excluded accessions never
    appear; output is sorted (draw order is implementation detail).
  - `write_accessions_tsv(accessions, path)` — emits the 2-column
    TSV (`accession\tnotes`) that `selexprep run` consumes.
  - `accessions_sha256(accessions)` — sha256 over the sorted
    accession list, used as the reproducibility fingerprint.
  - `CorpusAuditReport` dataclass with **explicit denominator
    partition**: reproducibility envelope
    (`catalog_version`, `sample_seed`, `sample_accessions_sha256`),
    fetchability metrics (denominator `n_sampled`), inference
    metrics (denominator `n_with_library_report`), QC metrics
    (denominator `n_with_qc_run`), per-accession traceability with
    `is_in_ground_truth` annotation.
  - `aggregate_audit_from_run_outputs(run_summary_tsv,
    ground_truth_tsv, *, catalog_version, sample_seed,
    sample_accessions_sha256)` — parses `run_summary.tsv` (Phase
    6a batch driver output), optionally overlaps with
    `ground_truth.tsv` for the overlap count, stamps the
    reproducibility envelope. The reproducibility kwargs are
    required (not defaulted) so callers cannot forget to thread
    them through.
  - `write_audit_json(report, path)` — deterministic sorted-keys
    JSON. Same discipline as `benchmark.metrics.write_metrics_json`.
  - `main(argv)` — CLI with `sample` and `aggregate` subcommands.
    `sample` writes both the accessions TSV and a sidecar
    `*.manifest.json` with the reproducibility envelope; `aggregate`
    reads the sidecar so the Snakefile doesn't have to thread
    catalog version / seed / sha through shell args.

- **`selexprep.benchmark.figure_b`** —
  - `plot_figure_b(audit_json, outdir) → (pdf, png)` — 2×2 panel
    mirroring Figure A's layout:
    - **Panel A** · Fetch outcomes · denominator `n_sampled`.
    - **Panel B** · `LibraryReport.status` (HIGH/MEDIUM/LOW/UNABLE
      _TO_INFER) · denominator `n_with_library_report`.
    - **Panel C** · `extraction_mode` distribution · denominator
      `n_with_library_report`.
    - **Panel D** · `required_action` distribution + inference
      safe-failure rate overlay annotation · denominator
      `n_with_library_report`.
    Each panel labels its denominator in the subtitle (per-panel
    transparency on the normalization). Title carries `n_sampled` /
    `n_fetchable` / `n_with_library_report` + catalog version +
    seed. PNG byte-output not deterministic across matplotlib
    versions (same discipline as Figure A); the audit JSON IS
    deterministic and is the source of truth.
  - `main(argv)` — CLI entry for the Snakefile `figure_b` rule.

**New benchmark infrastructure:**

- **`benchmarks/audit.smk`** — separate Snakefile (NOT merged with
  Tier 1's `Snakefile` — different outdir, different DAG shape,
  different cohort). DAG:
  - `rule sample_corpus` → `audit_accessions.tsv` +
    `audit_accessions.manifest.json` sidecar.
  - `rule run_corpus` → `run_summary.tsv` via
    `selexprep run --resume` (reuses the Phase 6a batch driver —
    no new runner code).
  - `rule aggregate_audit` → `audit_metrics.json` (reads sidecar
    manifest for the reproducibility envelope).
  - `rule figure_b` → `figure_b.{pdf,png}`.
  Configurable via `--config n_sample=30 seed=42`. Excludes
  `ground_truth.tsv` accessions by construction so Tier 1 rows
  don't double-count in Tier 2.

- **`benchmarks/README.md`** — new "Tier 2: corpus audit
  **pipeline**" section (deliberate wording: 6b.3a ships the
  pipeline, NOT the results). Documents what Tier 2 measures, what
  it does NOT claim ("descriptive distributional metrics, no
  per-row ground truth"), the explicit per-panel denominators, and
  the reproduction recipe (`snakemake -s audit.smk --cores N`).
  Calls out 6b.4 as the follow-up commit that lands actual
  `audit_metrics.json` + `audit_accessions.tsv` + `figure_b.{pdf,png}`
  after the HPC run.

**Codex post-implementation peer-review fixes (2026-05-23):**

One blocking + three non-blocking findings from the Codex peer-review
pass, all applied before the commit:

- **[blocking] `benchmarks/audit.smk` paths now anchored to
  `workflow.basedir`.** Originally `GROUND_TRUTH = "ground_truth.tsv"`
  and `OUTROOT = "audit_results"` were CWD-relative, so
  `snakemake -s benchmarks/audit.smk` from the repo root (the plan's
  verification command) died with `MissingInputException` because
  Snakemake looked for `ground_truth.tsv` in the repo root. Fixed by
  anchoring both to `Path(workflow.basedir).resolve()` —
  `workflow.basedir` is the standard Snakemake idiom for "the
  directory of the main Snakefile", stable across snakemake 6+.
  Both `cd benchmarks/ && snakemake -s audit.smk` and
  `snakemake -s benchmarks/audit.smk` from repo root now dry-run
  cleanly.
- **[non-blocking] `corpus_audit aggregate` now refuses when the
  reproducibility envelope is empty.** Previously, if neither
  `--sample-manifest` nor `--sample-sha` was provided, the CLI would
  silently emit an audit JSON with empty `sample_accessions_sha256`
  — defeating the whole point of the envelope (a paper-defensible
  audit fingerprinted by the sampled accession list). Now raises
  `SystemExit` with an actionable message; the function signature
  was already strict (positional kwarg, not defaulted), the CLI
  now enforces the same at its boundary. Regression test
  `test_main_aggregate_refuses_when_no_sha_available`.
- **[non-blocking] `corpus_audit aggregate` now warns on sidecar
  `n_sampled` mismatch.** If the sidecar manifest's `n_sampled`
  differs from the run summary's row count, the audit JSON still
  reports the actually-processed count (the source of truth for
  the per-accession breakdown), but a loud WARNING is printed and
  logged. Surfaces either a selexprep run bug (rows lost) or
  operator error (re-ran with a different TSV). Regression test
  `test_main_aggregate_warns_on_n_sampled_mismatch`.
- **[non-blocking] `write_audit_json` now defensively re-sorts
  `per_accession`.** The aggregator already sorts upstream, but
  the writer no longer trusts caller-supplied ordering — direct
  callers (tests, future hand-rolled scripts) cannot accidentally
  emit non-deterministic JSON. The existing
  `test_write_audit_json_deterministic` was extended to pass
  intentionally-unsorted input and assert sorted output.
- **[non-blocking] `benchmarks/README.md` duplicate-row dedupe.**
  `PRJNA558191` appeared twice in the candidate worklist (full row
  + a second "additional candidate discovered" mini-table). The
  second occurrence was pure documentation duplication and was
  removed.
- **[opportunistic, pre-existing] Tier 1 `benchmarks/Snakefile`
  path anchoring.** Tier 1 had the same cwd-relative quirk that
  the Tier 2 fix addressed (`GROUND_TRUTH = "ground_truth.tsv"`,
  `OUTROOT = "results"`), so
  `snakemake -s benchmarks/Snakefile ...` from repo root failed
  while `cd benchmarks && snakemake -s Snakefile ...` worked. Not
  a 6b.3a blocker (pre-existing Tier 1 behavior), but the fix was
  one-line and folded in here for consistency:
  `_BENCHMARKS_DIR = Path(workflow.basedir).resolve()` anchors
  `GROUND_TRUTH` + `OUTROOT` + `_GROUND_TRUTH_DIR`. Both
  Snakefiles now dry-run cleanly from either invocation point.
  No benchmark logic changed.
  Note: Ruff parses Snakefile rule syntax as plain Python and
  emits false syntax errors on `rule ...:` blocks, so
  `benchmarks/Snakefile` is NOT in the ruff scope. The standard
  `ruff check src/ tests/` is the contract.

**Tests (+~15, expected 494+ passed + 1 xfailed total):**

- `tests/test_benchmark_corpus_audit.py` (~12 tests):
  - Deterministic sampling (same seed → same accessions).
  - Different seed → different sample.
  - `exclude=(...)` drops accessions from the pool.
  - INSDC-only filter excludes `figshare:` / `zenodo:` / `utexas:` prefixes.
  - `sources` substring filter narrows the pool.
  - Empty pool returns `[]` (no crash).
  - `accessions_sha256` is order-independent.
  - `write_accessions_tsv` emits 2-column sorted TSV.
  - Aggregator's `inference_safe_failure_rate` excludes fetch failures
    (the regression test for the Codex+user methodological correction).
  - Aggregator's `n_fetchable` partition matches the locked plan.
  - QC flags histogram counts only OK rows with a flags_raised value.
  - Ground-truth overlap annotation per row + total overlap count.
  - Deterministic JSON output (sorted keys + stable accession order).
  - CLI: `sample` writes TSV + sidecar with envelope.
  - CLI: `aggregate` picks up sidecar values.

- `tests/test_benchmark_figure_b.py` (~5 tests, mirrors
  `test_benchmark_figure_a.py`):
  - PDF + PNG emission with non-trivial size.
  - Empty audit payload handled.
  - Parametric status-bucket smoke (HIGH/MEDIUM/LOW/UNABLE_TO_INFER).
  - Safe-failure rate overlay render path.
  - Forward-compatibility for an unexpected status label.

**Phase 6b.1 — primer-inference benchmark scaffolding + initial
ground-truth curation (2026-05-22)**

Ships the code + Snakefile + initial ground-truth row for the
selexprep primer-inference benchmark. The benchmark tests
selexprep's **unique claim**:

> Given only public/local HT-SELEX reads plus accession metadata,
> selexprep infers primer / constant regions and N-region length,
> reports confidence, and fails safely when inference is ambiguous.

Existing tools (AptaPLEX, EasyDIVER+, FASTAptameR) require the user
to supply the primers — they cannot benchmark this claim. Phase
6b.2 expands the ground-truth set to 5–10 verified rows; 6b.3 runs
the Snakefile on HPC and commits Figure A with measured numbers.

**Scope decision (deviates from locked plan lines 365–366).** The
locked plan listed count-correlation Pearson+Spearman and AptaPLEX /
EasyDIVER+ as comparators. A pre-commit Codex design review pivoted
the scope: comparator-tool head-to-head reduces to a trimming-code
sanity check given identical primers, not a meaningful scientific
comparison. Count correlation moves to optional Phase 6c as a
**self-consistency check** (inferred-primer counts vs
`--override-primer` counts of selexprep itself — no external tool).
`CountCorrelationReport` + `compute_count_correlation` stay in
`metrics.py` as the Phase 6c entry point (so the union+zero-fill /
no-SciPy methodology isn't re-implemented later); `aggregate_metrics`
does NOT call them in 6b.1.

**Codex pre-implementation review** raised seven additional design
amendments, all folded in BEFORE any code landed:

1. **`ground_truth.tsv` is verified-only.** Unverified candidates
   live in `benchmarks/README.md`'s candidate worklist; they do NOT
   pollute the ground-truth file. The metrics aggregator filters
   `verified=true` and emits a stderr warning per skip. The output
   JSON records both `n_verified` and `n_unverified` so Figure A
   labels honestly.
2. **Snakefile fetch/detect restricted to fetchable accessions.**
   Figshare / Zenodo / processed-data sources are excluded; v0.1
   `selexprep fetch` only handles ENA/SRA/DDBJ.
3. **Curated round-map override columns**: `round_map_source` ∈
   {`auto`,`curated`} + `round_map_path` (relative to
   `ground_truth.tsv`). Curated rows get `--allow-manual-review` on
   fetch so FASTQs download even when ENA metadata is too sparse
   for the auto round-parser; detect then uses the curator's
   `rounds.tsv`. Critical: without this, the benchmark would lose
   otherwise-known-primer datasets to ENA metadata sparsity.
4. **No SciPy.** Pearson via `pandas.Series.corr(method="pearson")`.
   Spearman computed manually as Pearson-of-ranks
   (`obs.rank().corr(ref.rank(), method="pearson")`) —
   mathematically equivalent. `pandas.Series.corr(method="spearman")`
   lazily imports SciPy, which would defeat the no-SciPy rule.
5. **Count correlation on union + zero-fill.** When 6c populates
   it: the Pearson is computed on the **union** of observed +
   reference sequences with zero-fill, NOT the intersection — which
   would bias agreement upward by ignoring sequences only one side
   emitted. A top-K Pearson is labeled as `top_k_pearson`
   (secondary diagnostic), never as the primary `pearson` field.
6. **Snakefile FASTQ enumeration uses `Path.glob`, not shell-glob.**
   `round_unknown/` only exists when fetch handled NONE-confidence
   runs; a bash glob like `round_unknown/*.fastq.gz` would pass the
   literal string to `selexprep detect` on auto-only accessions.
   `Path.glob` returns an empty iterator on non-matching patterns.
7. **Pyoverdine library_kind verified.** PRJNA932049 ENA project
   description explicitly says "2'-FY-RNA" — catalog title was
   correct; Codex's independent DNA flag did not match the ENA
   record. Verified before encoding `library_kind`.

**New modules:**

- **`selexprep.benchmark.equivalence`** — `primer_equivalent(observed,
  truth, *, allow_revcomp, allow_ut, strip_barcodes) →
  EquivalenceResult`. Implements the four locked-plan equivalence
  rules (line 364: revcomp + U-T + barcode-strip + IUPAC reject)
  plus PARTIAL_5P / PARTIAL_3P accounting. Promoted to strict-mypy
  alongside `library.report` (typed contract; metric aggregator
  depends on the EquivalenceKind enum).
- **`selexprep.benchmark.metrics`** — per-side primer recovery
  (`compute_primer_recovery`), pair-level recovery cross-tabbed
  against `LibraryReport.status`
  (`compute_pair_recovery_by_status`), safe-failure rate
  (`compute_safe_failure_rate` — the unique distinguishing
  capability vs known-primer pipelines), N-length recovery
  (`compute_n_length_recovery`), honest-accounting distributions
  for `extraction_mode` and `required_action`
  (`compute_extraction_mode_distribution`,
  `compute_required_action_distribution`), top-level
  `aggregate_metrics`, deterministic `write_metrics_json`, and
  `main()` for the Snakefile `compute_metrics` rule.
  `CountCorrelationReport` + `compute_count_correlation` remain as
  the Phase 6c entry point (union+zero-fill, pandas-only
  correlations) but are NOT populated by `aggregate_metrics`.
- **`selexprep.benchmark.figure_a`** — `plot_figure_a(metrics_json,
  outdir) → (pdf_path, png_path)`. **Four-panel** matplotlib Agg
  figure:
  - Panel A · pair recovery by `LibraryReport.status` (headline)
  - Panel B · per-side recovery breakdown (5'/3' EXACT/partial/miss)
  - Panel C · N-length recovery (±tolerance buckets)
  - Panel D · `extraction_mode` + `required_action` distributions
    (honest accounting)
  Title carries N verified + safe-failure rate.

**New benchmark infrastructure:**

- **`benchmarks/Snakefile`** — fetch → detect per verified
  accession → metrics → figure_a. Honors
  `round_map_source=curated` (passes `--allow-manual-review` to
  fetch and routes the curator's `rounds.tsv` to detect),
  enumerates FASTQs via `Path.glob` (shell-glob safety),
  resolves curated paths relative to `ground_truth.tsv` (Snakemake
  can be launched from any cwd).
- **`benchmarks/ground_truth.tsv`** — **11 verified rows**, each
  source-backed and audit-verified. Modality coverage: 2'-F-Py /
  2'-FY RNA (×4) + DNA (×5) + T7-tx RNA (×2) × protein (×7) /
  cell (×2) / small molecule (×2):

  1. **PRJNA315881** — Hoinka 2015 NAR (PMID 25870409, DOI gkv308)
     "Large scale analysis of the mutational landscape in HT-SELEX
     improves aptamer discovery" / IL-10RA / N40 / 2'-F-Py RNA.
     Primer source: PMC4499121 main M&M.
  2. **PRJEB22637** — Cibiel 2014 PLOS ONE (PMID 24489826) /
     ACE4 cell-SELEX (selection vs ETBR cells, aptamer later
     identified as anti-Annexin A2) / N50 / 2'-F-Py RNA.
     Primer source: PMC3906106 main M&M.
  3. **PRJEB28411** — Pleiko 2019 Sci Rep (PMID 31148584) /
     differential binding cell-SELEX vs ccRCC RCC-MF cells /
     N40 / DNA. Primer source: PMC6544647 main M&M.
  4. **PRJDB9110** — Ozaki 2020 NAR / RaptRanker Data1 (PMID 32537639)
     / transglutaminase 2 (TG2) / N30 / T7-tx RNA.
     Primer source: PMC7641312 main M&M.
  5. **PRJDB9111** — Ozaki 2020 NAR / RaptRanker Data2 (PMID 32537639)
     / integrin αVβ3 / N40 / T7-tx RNA.
     Primer source: PMC7641312 main M&M.
  6. **PRJEB70964** — Bouvier-Müller 2024 NAR (PMID 38917326) /
     α-synuclein fibrillar polymorphs / N35 / 2'-F-Py RNA.
     **Intentional edge case**: 5' constant is revcomp of TruSeq R1
     adapter prefix, exercising the adapter-trap path.
     Primer source: PMC11317169 main M&M.
  7. **PRJNA728693** — Sanford 2021 Chemical Science (PMID 34659704)
     / RE-SELEX for kanamycin A structure-switching aptamers /
     N40 / DNA. Edge case: EcoRI site embedded in the 5' constant.
     Primer source: RSC ESI PDF Fig. 2 (extracted via pdfplumber).
  8. **PRJNA935703** — Anisuzzaman 2024 Front Chem (PMID 39148668)
     / pyoverdine PYO-Pf5 / 2'-FY RNA. **Intentional edge case**:
     split-random-region library (N10 + internal constant + N35);
     N=55 is the total length between outer constants.
     Primer source: PMC11324436 main M&M.
  9. **PRJNA975735** — Halder 2023 Sci Rep (PMID 37666993) /
     SARS-CoV-2 viral oligopeptide/RBD aptamer SELEX / N40 / DNA.
     Primer source: PMC10477244 main M&M.
  10. **PRJNA990511** — Hu 2024 Sci Rep (PMID 38374125) / ASFV p30
      MB-SELEX / N40 / DNA. Primer source: MOESM1 supp docx
      Table S3 (extracted via python-docx).
  11. **PRJNA883192** — Ali 2022 Sci Rep (PMID 36577785) /
      eosinophil peroxidase (EPX) aptamer / N40 / DNA. Primer
      source: MOESM1 supp PDF page 2 Table S1, **verified via
      PyMuPDF raster render + visual inspection** (pdfplumber +
      pypdf silently drop the table because it's rendered as a
      raster image inside the PDF). DNA_L =
      ATGCCATCCTACCAAC-N40-GAGCTCTGAACTGG; FP1 = 5' library
      constant; RP1 = revcomp of 3' library constant. The row
      was initially demoted to `candidates.tsv` during the audit
      pass (when only text-extraction tools were tried) and
      promoted back to `ground_truth.tsv` after Codex pointed
      out the visual-rendering verification route.

  Curation toolkit: ENA project XML XREF_LINK → PubMed ID →
  NCBI PMC ID converter → PMC full-text (when in M&M body) →
  Europe PMC `/<PMCID>/supplementaryFiles` ZIP endpoint → `pypdf`
  + `pdfplumber` for supp PDF text + table extraction →
  `python-docx` for supp `.docx` table extraction → Springer CDN
  direct URL for MOESM files not in PMC → **`PyMuPDF` (`fitz`)
  raster render for tables that text extraction silently drops**
  (the Ali 2022 case demonstrated this failure mode: pdfplumber
  returned only the table title from the PDF's text stream,
  missing the actual table body which lives in a raster image).

- **`benchmarks/candidates.tsv`** — **NEW**. Structured record for
  9 documented candidates that did not pass curation. Schema:
  `accession`, `paper_doi`, `paper_pmid`, `source_attempted`,
  `primer_5p/3p_truth` (empty for unrecovered),
  `n_length_truth`, `status` ∈ {`blocked`, `rejected`},
  `rejection_reason`. Two confirmed **rejections** (Blocker-SELEX
  virtual screening + Anti-EGFR uses existing aptamers — both
  outside SELEX-with-random-N scope; backed by supp-PDF inspection
  via pypdf). Seven **blocked** (Dao Cell Press, PRIESSTESS NAR/OUP,
  SELMAP Sci Rep, Penzar bioRxiv, DL-SELEX, Iowa State pyoverdine,
  Camorani iScience TNBC).

- **Audit fixes:**
  - **PRJNA315881 (Hoinka)** row had wrong DOI
    (`10.1093/nar/gkw010`) + wrong PMID (`26773059`) — those
    identify Chen WH 2016 NAR "Integration of multi-omics data of
    a genome-reduced bacterium" (DOI gkw004), not Hoinka 2015.
    Fixed to `10.1093/nar/gkv308` + PMID `25870409` ("Large scale
    analysis of the mutational landscape in HT-SELEX"), which is
    the actual paper the primer sequences were extracted from.
  - **PRJNA883192 (Ali)** was initially demoted to candidates.tsv
    during the audit pass when pdfplumber returned an empty
    Table S1 from the MOESM1 supp PDF. Codex pointed out that the
    table is a raster image inside the PDF; visual inspection of
    PyMuPDF-rendered page 2 confirmed DNA_L = ATGCCATCCTACCAAC-N40-
    GAGCTCTGAACTGG. Row promoted back to ground_truth.tsv. Lesson
    captured in `benchmarks/README.md`: failed text extraction
    does NOT mean data is absent — try raster render before
    blocking a candidate.
- **`benchmarks/README.md`** — "what this benchmark tests / what
  it deliberately does not test" framing, schema docs,
  reproduction recipe, scope-pivot deviation note, candidate
  worklist for 6b.2.

**Tests (+46, 478 + 1 xfailed total):**

- `tests/test_benchmark_equivalence.py` (16 tests): every
  `EquivalenceKind` outcome + observed=None + IUPAC truth +
  barcode longest-prefix.
- `tests/test_benchmark_metrics.py` (23 tests): primer recovery
  aggregation + pair recovery by status (4 tests) + safe-failure
  rate (3 tests) + N-length tolerance + extraction_mode and
  required_action distributions + deterministic JSON
  serialization + Codex-amendment regressions (verified-row
  filter with stderr warning; count_correlation union+zero-fill;
  top_k_pearson labeled secondary; aggregate_metrics does NOT
  populate count_correlation; metrics JSON schema includes the
  pivot fields).
- `tests/test_benchmark_figure_a.py` (7 tests, includes 5
  parametric status-bucket smokes): 4-panel PDF + PNG emission,
  empty-data path, safe-failure rate in title.

**`pyproject.toml`:**

- New `[project.optional-dependencies] bench` extra with
  `snakemake >= 7.0` (Python 3.10 compatible; no SciPy).
- Strict-mypy override added for `selexprep.benchmark.equivalence`.

CI: 478 passed + 1 xfailed (was 432 + 1; +46 new tests). All four
pre-existing gates stay green (ruff, ruff format, mypy, pytest
3.10/3.11/3.12).

### Fixed

**Phase 6a CI hotfix (2026-05-22)** — `test_run_help_lists_resume_and_stop_on_error`
passed locally on macOS but failed on all three Python versions on
Ubuntu CI. Root cause: Typer/Rich renders the help table at the
terminal's column width, and CI's non-TTY default of 80 columns is
narrow enough that the Rich renderer silently truncates `--resume`,
`--stop-on-error`, and `--backend` option names out of the rendered
text. macOS local runs detect a wider effective width (Rich + CliRunner
interaction) and show all options. The test grepped the rendered help
text, so it was width-dependent by construction.

Fix: rewrite the test to introspect the Click command's `params` list
(`flag_names = {opt for p in run_cmd.params for opt in p.opts}`)
instead of grepping rendered output. Renamed to
`test_run_command_registers_resume_and_stop_on_error_options`. Kept
a smoke that `--help` exits cleanly so the rendering path isn't
unobserved. The contract being tested ("these options exist") is
captured exactly; the rendered help width remains
environment-dependent (matplotlib-style informational, per the same
discipline as Phase 5 plots).

**Phase 6a Codex pass 1 (2026-05-22)** — two blocking semantic bugs +
two non-blocking polish items. Calibration verdict: N/A (Phase 6a
introduces no calibration constants).

Blocking:

- **Outer except in `run_batch` mislabeled every unexpected exception
  as `EXTRACT_FAILED` (run/runner.py).** If `build_fetch_plan` or
  `download_srr` raised an `HTTPError`, `RequestException`, or
  `ValueError`, the exception escaped the fetch block to the outer
  safety net and `run_summary.tsv` recorded `status=EXTRACT_FAILED`
  even though the extract stage never ran. Fix: stage-classify
  expected fetch-stage exceptions inside `_process_one_accession`
  (new `try/except` around the `run_fetch` call → `FETCH_FAILED`),
  and rename the outer-net status to a new literal
  `UNEXPECTED_FAILURE` so a truly-escaped exception is no longer
  misreported as extract. The `RunStatus` literal gains
  `UNEXPECTED_FAILURE`; nothing produces `EXTRACT_FAILED` from the
  outer net any more. Resumed reads of `fetch_metadata.json` now also
  guard `read_fetch_metadata_json` against malformed audit-trail
  files (`ValueError` / `KeyError` / `OSError`) and downgrade to
  `FETCH_FAILED` instead of leaking.
- **Fetch resume oracle accepted corrupt `.fastq.gz` files
  (fetch/runner.py).** Both the per-run skip-already-present check
  (`run_fetch`) and `check_fetch_inventory` used bare `Path.exists()`.
  A `SIGKILL`'d download leaves a truncated `.fastq.gz` on disk that
  exists but cannot be decompressed; `--resume` would skip it and
  downstream `detect`/`extract` would hit gzip errors deep in the
  pipeline. Fix: both call sites now use
  `selexprep.fetch.download.validate_fastq_gz` (existing helper at
  fetch/download.py:96 — checks exists, ≥1024 bytes, and `gzip -t`
  decompresses cleanly). Mirrors the same check
  `download_srr_ena_direct` performs at fetch/download.py:288.

Non-blocking:

- **`_count_yaml_flags` used `line.startswith("- name:")` which
  miscounts under `safe_dump(sort_keys=True)` (run/runner.py).** When
  a flag has non-empty `evidence`, the YAML emission orders dict keys
  alphabetically inside each list entry (`evidence` < `name` <
  `severity`), so the first line of each entry can be
  `- evidence:` rather than `- name:`. Resumed `run_summary.tsv`
  would undercount flags. Fix: parse the YAML via `yaml.safe_load`
  and count list entries that contain a `name` key.
- **`fastq_filenames_for_run` synthesized
  `{SRR}.fastq.gz`/`{SRR}_{1,2}.fastq.gz` instead of deriving from
  URLs (fetch/plan.py).** Right for ENA's canonical naming, but
  brittle if a future ENA URL convention or Zenodo mirror used a
  different scheme. `download_srr_ena_direct` writes
  `output_dir / Path(url_path).name` (fetch/download.py:351), so the
  URLs are the source of truth. Fix: derive basenames from
  `run.fastq_urls`; keep the SRR-based synthesis as a defensive
  fallback for the rare case where ENA returns no `fastq_ftp` (which
  also means the run isn't downloadable, so the fallback path is
  effectively dead).

Six new regression tests:

- `test_run_batch_fetch_http_error_records_FETCH_FAILED_not_extract_failed`
  — `requests.HTTPError` from `build_fetch_plan` ⇒ `FETCH_FAILED`
  (not `EXTRACT_FAILED`).
- `test_run_batch_corrupt_fastq_is_redownloaded_by_resume_oracle` —
  garbled `.fastq.gz` on disk triggers re-download under `--resume`;
  the result passes `validate_fastq_gz`.
- `test_run_batch_count_yaml_flags_uses_yaml_parser` — flags YAML
  with `evidence`-first entry ordering counts to 2, not 0.
- `test_fastq_filenames_for_run_derives_from_urls` — hypothetical
  alternate URL naming (`foo_subdir/SRR_X_unusual.fastq.gz`) is
  honored verbatim.
- `test_fastq_filenames_for_run_fallback_when_urls_empty` — defensive
  synthesis path still works when `fastq_urls=[]`.
- `test_run_batch_unexpected_error_uses_unexpected_failure_status` —
  outer-net safety catch uses `UNEXPECTED_FAILURE`, not
  `EXTRACT_FAILED`.

Test fixtures updated: `_stub_download_writes_files` and `_write_fastq`
now emit 800 records of 60-base pseudo-random sequences so the gzipped
outputs exceed `validate_fastq_gz`'s 1024-byte floor. The previous
small fixtures (~38 bytes raw) would now be rejected as too-small to
be valid downloads.

CI: 432 passed + 1 xfailed (was 426; +6 new tests). All four
pre-existing CI gates stay green.

### Added

**Phase 6a — close v0.1 CLI surface: `selexprep fetch` + `selexprep run`
(2026-05-22)**

Wires the two remaining `_not_implemented` CLI verbs from the locked
plan (lines 181–189), closing the v0.1 single-dataset + batch CLI
surface. With Phase 6a landed, the public CLI exercises end-to-end
(`inspect → fetch → detect → extract → count → qc`); Phase 6b
(benchmark Snakefile + Figure A) can now drive the comparison against
AptaPLEX + EasyDIVER+ through the public CLI rather than through
library imports.

- **NEW: `selexprep.fetch.plan`** — `build_fetch_plan(accession,
  timeout_s) → FetchPlan` hits ENA filereport with extended fields
  (`sample_title`, `library_name`, `experiment_title`,
  `sample_accession` on top of the inspect set) so the 5-level
  cascade in `fetch.metadata.parse_round` runs on real metadata
  instead of bare-bones inspect fields. Emits a deterministic
  `fetch_metadata.json` audit trail. Single source of truth for
  "what does this accession look like and how do its runs map to
  rounds?", reused by both `fetch` and `run` CLIs.
- **NEW: `selexprep.fetch.runner`** — `run_fetch(accession, outdir,
  *, backend, allow_manual_review, dry_run, timeout_s) → FetchResult`
  orchestrator. Emits `rounds.tsv` (HIGH/MEDIUM-confidence
  contract that downstream `detect`/`extract` consume; sorted by
  filename for determinism) + per-round FASTQs under `round_NN/` +
  optional `manual_review.tsv` + `fetch_metadata.json`. Cardinal
  rule (per `fetch/metadata.py` line 14): never guess a round
  assignment — NONE-confidence runs refused by default; opt-in via
  `--allow-manual-review` routes them to `round_unknown/` + a
  separate `manual_review.tsv` and **never** lets them enter
  `rounds.tsv`. If every run is NONE-confidence, refuses fail-fast
  before any download. Shared `query_ena_filereport()` helper
  factored out of `inspect.py` so the HTTP call shape is defined
  once.
- **NEW: `selexprep.run.runner`** — `run_batch(accessions_tsv,
  outdir, *, resume, stop_on_error, backend, allow_manual_review,
  timeout_s) → RunReport` batch driver. Per-accession pipeline
  (fetch → detect → extract → count → qc) with file-inventory
  resume oracles per stage (not sentinel-flag based):
  - **fetch**: `fetch_metadata.json` present AND every expected
    FASTQ from `FetchPlan.runs` on disk AND `rounds.tsv` present.
    Mirrors `download_bioproject._missing_srrs` discipline at
    `fetch/download.py:683`.
  - **detect**: `library_report.json` present.
  - **extract**: `selexprep_manifest.json` present.
  - **count**: per-round `round_NN/counts.parquet` (granular —
    only missing rounds re-run).
  - **qc**: `qc/flags.yaml` present.
  Paired-end FASTQs are grouped by ENA naming convention
  (`<srr>_1.fastq.gz` = R1, `<srr>_2.fastq.gz` = R2); R1-only
  sequences feed `compute_library_report`'s primary stream, R2
  sequences are passed as `paired_mate_streams`, and both are
  threaded through `run_extract` as `paired_r2_inputs`. Per-accession
  errors default to log-and-continue with `status=FAILED_<stage>`
  in `run_summary.tsv`; `--stop-on-error` flips to fail-fast.
- **WIRED: `selexprep fetch <accession> --outdir OUT
  [--backend ena|auto|kingfisher|sra] [--allow-manual-review]
  [--dry-run] [--timeout-s N]`** — replaces the
  `_not_implemented` stub. **CLI default `--backend ena`** (paper-
  grade reproducibility; fail-fast if ENA can't serve; never
  silently fall back to GPL tools); `--backend auto` is the
  explicit opt-in for the convenience fallback chain. Library-level
  `download_srr(backend="auto")` default unchanged (Python API
  ergonomic).
- **WIRED: `selexprep run accessions.tsv --outdir OUT [--resume]
  [--stop-on-error] [--backend ena|auto|...] [--allow-manual-review]
  [--timeout-s N]`** — replaces the `_not_implemented` stub. Emits
  a deterministic `run_summary.tsv` (sorted by accession) with
  per-accession status: `OK` /
  `SKIPPED_READ_MERGING_RECOMMENDED` /
  `FETCH_REFUSED` / `FETCH_FAILED` / `DETECT_FAILED` /
  `EXTRACT_REFUSED` / `EXTRACT_FAILED` / `COUNT_FAILED` /
  `QC_FAILED`. Split-primer guard (`required_action ==
  READ_MERGING_RECOMMENDED`) skips count + qc and records the
  status rather than producing misleading half-insert counts
  (locked plan line 325).

**Tests added (+42, 427 + 1 xfailed total):**

- `tests/test_fetch_plan.py` (12 — mocked-ENA paired-vs-single,
  round confidence propagation, deterministic JSON, frozen
  dataclass, extended-field assertion, empty-response ValueError,
  timeout pass-through).
- `tests/test_fetch_cli.py` (15 — orchestrator happy paths,
  refusal paths (all-NONE, mixed-without-flag, allow-manual-review
  rounds.tsv cleanliness), resume oracle (`check_fetch_inventory`
  + skip-already-present), CLI dry-run, invalid backend,
  refusal exit code, propagated ValueError, all four documented
  backends accepted).
- `tests/test_run_runner.py` (12 — duplicate-accession refusal,
  missing-column refusal, fetch refusal status, paired-end
  R1/R2 threading through detect+extract, split-primer skip
  (no count/qc), fetch resume oracle re-fetches missing FASTQ,
  per-stage resume (no work on re-run with all sentinels),
  log-and-continue default, `--stop-on-error` halts loop,
  detect-stage exception → `DETECT_FAILED`, summary TSV sorted
  by accession, manual-review separation keeps `rounds.tsv`
  clean).
- `tests/test_cli.py` (+3 — fetch dry-run smoke, run missing-
  accession-column rejection, run --help lists `--resume` +
  `--stop-on-error`). Replaced the obsolete
  `test_fetch_stub_exits_with_code_2`.

**Calibration:** N/A — Phase 6a introduces no new heuristic thresholds.
The `CALIBRATION-TODO` inventory stays at 19. The fetch refusal
threshold (NONE-confidence) is the locked-plan cardinal-rule binary,
not a tunable.

### Fixed

**Phase 5 Codex pass 1 (2026-05-21)** — three blocking semantic bugs +
three non-blocking polish items. Calibration verdict: all six v0.1
qc constants ✅ confirmed against published HT-SELEX conventions
(RAPID-SELEX, AptaPLEX, Hoinka 2015, EasyDIVER+, FASTAptamer).

Blocking:

- **`rarefy()` returned the original pool unchanged when `depth >=
  total`, reintroducing the depth-confounding the flag was meant to
  avoid (qc/flags.py).** A round with 2k reads was being compared
  against a round with 100k reads rarefied to 10k — exactly the trap
  the locked plan flagged at line 351. Fix: clamp to
  `effective_depth = min(RAREFACTION_DEPTH, min_total_reads_per_round)`.
  Rarefy all rounds to that common depth. Evidence dict now includes
  both `configured_depth` and `effective_depth` for audit. When a
  round is empty (`min_total == 0`), return None (the `low_total_reads`
  flag handles that separately).
- **`check_adapter_contamination_high` mixed measurement universes
  (qc/flags.py).** Numerator was `LR.known_adapter_hits` (Phase 2
  inference on earliest-round subsampled reads); denominator was
  `sum(counts_by_round)` (post-extraction reads across all rounds).
  Apples-to-oranges → fraction could be inflated or deflated by an
  arbitrary factor. Fix: switched the denominator to
  `trim_reports_by_round[earliest_round]["n_in"] *
  lr.read_fraction_used_for_inference`, matching the Phase 2
  inference universe. When `trim_reports.json` is unavailable, surface
  the flag as `severity="info"` with `evidence={"reason":
  "denominator_unavailable"}` rather than guessing. Threaded a
  `trim_reports_by_round` kwarg through `compute_all_flags`; the qc
  runner already aggregates the JSON via `_load_trim_reports_by_round`
  and now passes it down.
- **`selexprep count` silently parsed FASTQ as FASTA (cli.py).** The
  CLI called `count_fasta` unconditionally; a `.fastq.gz` input would
  be misparsed (every other line treated as a sequence header). Fix:
  hard-reject FASTQ extensions by default with a clear error pointing
  the user at `selexprep extract` first OR at the new
  `--from-pretrimmed-fastq` opt-in flag. The v0.1 contract is "`count`
  accepts only extracted FASTA from `extract`"; users with externally
  pre-trimmed FASTQ (e.g., from AptaPLEX, EasyDIVER+, external
  cutadapt) can opt in via `--from-pretrimmed-fastq`, which routes
  through a new `count_fastq_pretrimmed()` library function and emits
  a loud "cannot verify trimming state" `logger.warning`. The legacy
  Phase-1 `count_round` (inline cutadapt trimming) stays a library-
  only API.

Non-blocking:

- **`flags.yaml` was promised to contribute to `output_sha256` but
  didn't (qc/plots.py + manifest.py).** Phase 4's manifest is sealed
  by `extract` before `qc` runs, so the QC artifacts are outside the
  extract manifest lifecycle by construction. Added `.yaml` to
  `_HASHABLE_SUFFIXES` (so future `selexprep qc-amend` could append
  the hash) and rewrote the qc/plots.py docstring to clarify the
  current lifecycle and the optional hash path.
- **`_parse_round_label` accepted negative integers (cli.py).** `R-1`
  produced `round_-1/`. Fix: explicit `r < 0` rejection with
  `typer.BadParameter`.
- **YAML emission lacked a regression test for nested-dict evidence
  (tests/test_flags.py).** `safe_dump(sort_keys=True)` does sort
  nested dict keys recursively, but the determinism guarantee for
  nested LIST-of-DICT evidence (as in
  `adapter_contamination_high`'s `adapters_above_threshold` field)
  wasn't exercised. Added
  `test_write_flags_yaml_deterministic_with_nested_evidence` —
  builds the same flag twice with intentionally-different insertion
  orders and asserts byte-identical output.

### Carry-forward to v0.2 (documented now, not yet implemented)

- **`check_adapter_contamination_high` assumes `min(trim_reports_by_round)`
  is the same round Phase 2 used for primer inference.** Holds for the
  documented v0.1 CLI flow (same FASTQ set passed to detect + extract).
  Could break for partial extract runs, manually-stitched
  ``trim_reports.json``, or the v0.2 batch driver. v0.2 fix: add
  ``LibraryReport.earliest_inference_round`` (set in
  ``compute_library_report``) and key the denominator off that. Schema
  bump is disproportionate for v0.1 since the bug is unreachable through
  the documented CLI flow. Inline comment in
  ``qc/flags.py:check_adapter_contamination_high`` flags the assumption
  for the v0.2 implementor.
- **`_iter_fastq_sequences` only validates record completeness, not
  per-line format conformance** (no ``@`` header check, no ``+``
  separator check, no ``len(seq) == len(qual)``). Acceptable for the
  ``--from-pretrimmed-fastq`` power-user opt-in; if this path is
  promoted to a first-class advertised workflow in v0.2, add the
  validators here matching the truncation ValueError policy.

### Added

**`selexprep.count.counter.count_fastq_pretrimmed`** — library function
for the `--from-pretrimmed-fastq` CLI opt-in. Parses FASTQ records,
extracts the sequence line, builds a Counter, emits the standard
counts.parquet schema (`sequence`, `reads`, `rank`, `rpm`). Same
truncation policy as `extract/strand.py:reorient_fastq_gz` — raises
`ValueError` on incomplete records rather than silently producing
partial output. Used for the external-tool interop story (AptaPLEX +
EasyDIVER+ produce pre-trimmed FASTQ that doesn't need selexprep's
extraction).

Ten new regression tests:
- `test_diversity_increase_clamps_depth_to_min_total`
- `test_diversity_increase_records_effective_depth_in_evidence`
- `test_adapter_contamination_denominator_unavailable`
- `test_adapter_contamination_honors_read_fraction_used_for_inference`
- `test_count_rejects_negative_round_label`
- `test_count_rejects_fastq_input` (default reject path)
- `test_count_with_from_pretrimmed_fastq_succeeds`
- `test_count_from_pretrimmed_fastq_logs_unverified_warning`
- `test_count_from_pretrimmed_fastq_flag_rejects_fasta_input`
- `test_write_flags_yaml_deterministic_with_nested_evidence`

Plus three library-level regressions for `count_fastq_pretrimmed`:
- `test_count_fastq_pretrimmed_basic`
- `test_count_fastq_pretrimmed_raises_on_truncated_record`
- `test_count_fastq_pretrimmed_handles_uncompressed`

Two existing tests updated to reflect the new
`check_adapter_contamination_high` signature (now takes
`trim_reports_by_round`, not `counts_by_round`).

CI: 385 passed + 1 xfailed (was 372; +13 new tests). All four
pre-existing CI gates stay green.

**Phase 4 Codex pass 1 (2026-05-21)** — one blocking semantic bug +
three non-blocking polish items. Calibration verdict: N/A (Phase 4
introduces no calibration constants).

Blocking:

- **`--override-primer-{5p,3p}` was applied AFTER the refusal check
  (runner.py).** Locked plan line 311 explicitly says
  `--override-primer-*` should bypass the UNABLE_TO_INFER /
  UNABLE_TO_EXTRACT refusal — but the runner was checking refusal
  first, so an override could never reach the LR. Worse, the refusal
  message suggested using override as the fix, creating a circular
  error path. Refactored to:
  1. Apply override first (clone LR via `model_copy(update=...)`).
  2. When baseline `extraction_mode == "UNABLE_TO_EXTRACT"`, infer a
     new mode from which primer sides are now defined: both →
     `BOTH_PRIMERS_SINGLE_READ` + `full_insert_recovered=True` +
     `required_action=NONE`; 5p-only → `FIVE_PRIME_ONLY`; 3p-only →
     `THREE_PRIME_ONLY` (`full_insert_recovered=False`).
  3. Promote `status` from `UNABLE_TO_INFER` → `MEDIUM` (manual
     override = medium confidence by convention; hand-editing the LR
     JSON is the path to claim HIGH). Clear `failure_reason`.
  4. THEN run the refusal check. With override applied, UNABLE state
     is naturally cleared and the dispatcher proceeds. Refusal still
     fires if no override was given AND baseline is UNABLE.
  Refusal message also reworded to clarify the override mechanism
  (dropped the stale `(Phase 4)` parenthetical).

Non-blocking:

- **Override didn't check `KNOWN_ADAPTERS` (runner.py).** An override
  matching a known sequencing adapter prefix is now allowed (explicit
  escape hatch) but surfaces a `logger.warning` — the adapter-trap
  story still holds for the auto-inference path, and the override
  path gets a foot-gun guard without blocking deliberate use.
- **`HTTPError` in `inspect_accession` propagates without
  differentiating 404 vs 503 (fetch/inspect.py).** OK for the current
  CLI flow; flagged for the future `selexprep run` batch driver.
  Carry-forward only — not fixed in this round.
- **`parameters` / `runtime_seconds_per_stage` not normalized in
  deterministic JSON (manifest.py).** Currently deterministic because
  the CLI always builds the dict in the same order, but fragile if the
  library API gets a non-CLI caller. Carry-forward.

### Changed

**Architecture hygiene (Codex pass 1 NB follow-up).** Migrated the
`_matches_known_adapter` helper from `library/detect.py` (where it was
underscore-private) to `library/adapters.py` (where `KNOWN_ADAPTERS`,
`KNOWN_ADAPTERS_RC`, and `count_adapter_hits` already live). Renamed
to public `matches_known_adapter_prefix()` — the `_prefix` suffix
makes the semantics explicit (it's a first-`k`-bp match, not arbitrary
substring). Also moved the `ADAPTER_PROBE_K = 13` constant to
`library/adapters.py` so the shared default of `count_adapter_hits`
and `matches_known_adapter_prefix` has a single source of truth.
Cross-package private imports eliminated. No call-site behavior
change; the function is byte-identical.

### Tests (Phase 4 follow-up)

Five new regression tests in `tests/test_extract_override.py`:
- `test_override_with_both_primers_promotes_unable_to_extract`
- `test_override_5p_only_promotes_to_five_prime_only`
- `test_override_3p_only_promotes_to_three_prime_only` (symmetry
  mirror of the 5p-only test)
- `test_override_logs_warning_when_matches_known_adapter`
- `test_no_override_on_unable_still_refuses_with_clear_message`

CI: 372 passed + 1 xfailed (was 367; +5 new tests). All four
pre-existing CI gates stay green. `_make_library_report` in
`test_extract_override.py` extended with kwargs to support multiple
test scenarios (matches the helper pattern in `test_extract_runner.py`).

**Phase 3 Codex pass 1 (2026-05-21)** — three blocking bugs +
three non-blocking polish items, all caught in the second-pass code
review of the extraction pipeline. Calibration verdict for the one
Phase 3 constant: `STRAND_REPORT_PER_READ = False` ✅ CONFIRMED.

Blocking:

- **Multi-FASTQ-per-round overwrite (runner.py).** Each per-mode
  trim loop wrote every input FASTQ in a round to the SAME target
  path — second iteration overwrote first, dropping all but the last
  input's reads. Added `_trim_round_single_end` + `_trim_round_paired_split`
  helpers that trim each input to a per-input temp `.part_NNN.fasta.gz`,
  then concatenate deterministically via `open_gzip_text_deterministic`
  into the final per-round target. Single-input case stays on the
  fast path (no temps).
- **Sample-sheet paired-end demux didn't rebuild `paired_r2_inputs`
  (runner.py).** After demux, only R1-shaped files were collected;
  `paired_r2_inputs` stayed at whatever the caller passed (typically
  None). `PAIRED_END_SPLIT_PRIMERS + --sample-sheet` would refuse with
  "requires --paired-r2". Now collects both `_1.fastq.gz` (R1) and
  `_2.fastq.gz` (R2) from the demuxed layout, builds `round_inputs`
  path-aware from parent dirs (basenames collide across rounds after
  demux), and populates `paired_r2_inputs` when R2 files exist.
- **`output_sha256` collapsed per-round files with same basename
  (manifest.py).** `compute_sha256s` used `p.name` as the dict key, so
  `round_00/extracted.fasta.gz` and `round_01/extracted.fasta.gz`
  collided — last-write-wins. Added `root` kwarg; when given, keys are
  `path.relative_to(root).as_posix()`. Builder passes `output_root=outdir`
  so per-round outputs stay distinct in the manifest.

Non-blocking:

- **Orphaned temp files on trim failure (trim.py).** All four trim
  entry points now wrap cutadapt + repack in `try/finally` that unlinks
  the intermediate `.fa` and `.cutadapt.json` files, even on cutadapt
  raise or repack raise.
- **Truncated FASTQ silently produced partial output (strand.py).**
  `reorient_fastq_gz` now raises `ValueError` on a truncated record
  instead of logging a warning and breaking the loop (the caller has
  no way to distinguish "graceful early termination" from
  "incomplete extraction"). Aligns with the project's no-silent-miscalls
  discipline.
- **R2 basename collision (cli.py).** `extract` now refuses early
  with an explicit list of duplicate basenames when `--fastq` +
  `--paired-r2` collectively contain repeated names (which would
  collapse in the basename-keyed round-map lookup).

Follow-up to the manifest pass (Codex spotted in a second review pass
of the same commit, folded in here):

- **`input_sha256` basename-collision in sample-sheet mode (manifest.py +
  runner.py).** The `root` kwarg now applies to inputs too:
  `build_manifest_from_extract_result` accepts `input_root=...`, and
  the runner passes `input_root = outdir / "demux"` when sample-sheet
  mode is active (demuxed files share basenames across rounds —
  `srr_1.fastq.gz` in every `round_NN/` folder). Normal CLI flow keeps
  `input_root=None` because the CLI's basename-collision guard
  guarantees uniqueness upfront.
- **Latent: `input_sha256` was silently empty for ALL CLI flows
  (manifest.py).** `_HASHABLE_SUFFIXES` only included
  `{.fasta, .fa, .tsv, .json}`, so FASTQ inputs got filtered out and
  the manifest's input-side hashes were missing entirely — contradicting
  the locked plan line 168 ("input_sha256 (FASTQ files)"). Added
  `.fastq` + `.fq` to the suffix set. Caught while writing the regression
  test above.

Six new regression tests:
- `test_run_extract_multi_fastq_same_round_aggregates`
- `test_run_extract_sample_sheet_paired_end_demux_rebuilds_r2_inputs`
- `test_run_extract_sample_sheet_input_sha256_distinct_per_round`
- `test_compute_sha256s_distinct_keys_per_round_with_root`
- `test_compute_sha256s_without_root_falls_back_to_basename`
- `test_reorient_fastq_gz_raises_on_truncated_record`

CI: 367 passed + 1 xfailed (was 361; +6 new tests). All four pre-existing
CI gates stay green.

**Phase 2 Codex pass 1 semantic bugs (2026-05-20)** — surfaced during
the calibration review and folded into the same commit:

- **`match_rate_*` and `position_consistency_*` were aliased.** The
  orchestrator was setting `position_consistency_5p = match_rate_5p`
  literally, double-counting the same evidence in the composite
  confidence formula (which weights them as 4 distinct signals).
  Added `_substring_match_rate()` (substring-anywhere with Hamming ≤
  1), and now `match_rate_*` comes from substring scan while
  `position_consistency_*` comes from `_position_consistency()`
  (flank-anchored ± tolerance). They diverge by construction.
- **`adapter_clean` ignored the drop flag.** Previously the signal was
  `1.0 if (primer_5p_seq is not None or primer_3p_seq is not None)
  else 0.0` — true whenever ANY primer survived, hiding the adapter
  trap. Now tracks per-side `adapter_drop_{5p,3p}` flags at the drop
  site; the signal is `0.0` if either side was dropped as an adapter
  match.
- **Paired-split 3p signals were measured against R1.** In
  `PAIRED_END_SPLIT_PRIMERS` mode, `match_rate_3p`,
  `position_consistency_3p`, and `variants_3p` were computed by
  looking for the (R2-derived) `primer_3p_seq` at R1's 3' end — where
  it cannot exist by construction (that's why it's a split). Refactored
  the orchestrator to toggle a 3p-signal context (R1/3'-end normally,
  R2/5'-end-using-revcomp for paired-split). Per-round persistence
  input for the 3p side now uses the R2 stream in split mode.
- **`_persistence_score` returned `0.0` instead of `None` when mean <
  0.1.** Semantically: `None = "not evaluable"` vs `0.0 = "evaluated
  and zero"`. The composite-confidence formula treats both the same
  numerically (None contributes 0 via `if v is None: continue`; 0.0
  contributes 0 via 0.0 × weight), so the user-visible composite is
  unchanged. The semantic alignment matters for the manifest's
  audit trail.

In `PAIRED_END_SPLIT_PRIMERS` mode the `n_length_*` fields are now set
to `None` / `{}` / `0.0` rather than computed from R1 alone — the full
insert spans R1+R2 and cannot be measured from either read.

Three new regression tests in `tests/test_report.py`:
- `test_match_rate_distinct_from_position_consistency`
- `test_adapter_clean_flag_demotes_confidence_when_primer_dropped`
- `test_paired_split_match_rate_3p_reflects_r2_not_r1`

Behavior-based classification tests stay green throughout (361 passed
+ 1 xfailed; +3 vs previous baseline) — the LR's classification
fields (`extraction_mode`, `required_action`, `status`) are unchanged
for clean inputs; only the numeric signal fields change.

### Calibration

**Phase 2 Codex peer-review, pass 1 (2026-05-20).** Eight constants in
`library/detect.py` were reviewed; six confirmed at locked-plan
defaults, four revised:

| Constant | Before | After | Rationale (Codex) |
|---|---|---|---|
| `POSITION_CONSISTENCY_TOLERANCE` | 2 nt | **3 nt** | AptaPLEX default mismatch tolerance is 3; small public-data offset noise should not over-penalize an otherwise stable flank. |
| `STATUS_HIGH_CUTOFF` | 0.80 | **0.85** | "HIGH" should mean paper-grade high-confidence, harder to reach via additive secondary signals before benchmark calibration. |
| `COMPOSITE_WEIGHTS` (with round map) | match 0.20/0.20 · pos 0.10/0.10 · persistence 0.20 · n_len 0.10 · adapter_clean 0.10 | **match 0.15/0.15 · pos 0.15/0.15 · persistence 0.25 · n_len 0.10 · adapter_clean 0.05** | Position consistency deserves parity with raw match rate; cross-round persistence is the unique SELEX-specific signal (Hoinka et al. 2015, AptaTRACE / AptaTools); adapter_clean is already enforced upstream as a blacklist, so it should be a small confidence bonus, not a driver. |
| `COMPOSITE_WEIGHTS_NO_ROUND_MAP` | match 0.30/0.30 · pos 0.15/0.15 · persistence 0.00 · n_len 0.05 · adapter_clean 0.05 | **match 0.225/0.225 · pos 0.225/0.225 · persistence 0.00 · n_len 0.05 · adapter_clean 0.05** | With persistence absent and status already capped at MEDIUM, within-round evidence (match + position) should carry equal parity weight; n_len and adapter_clean remain supporting signals, not drivers. |

Confirmed as-is: `PRIMER_FOUND_MATCH_RATE_THRESHOLD = 0.70`,
`N_LENGTH_CONFIDENT_FRACTION = 0.80`, `UNABLE_TO_EXTRACT_MATCH_RATE =
0.40`, `STATUS_MEDIUM_CUTOFF = 0.60`, `STATUS_LOW_CUTOFF = 0.30`,
`ORIENTATION_REVERSED_FORWARD_MAX = 0.05`,
`ORIENTATION_REVERSED_REVERSE_MIN = 0.95`.

Tests stayed green across the threshold flip (358 + 1 xfailed) because
they assert on `extraction_mode` / `required_action` / `status`, never
on the threshold numbers themselves. Marker convention shifted:
`CALIBRATION-TODO` → `CALIBRATION-REVIEWED (Codex 2026-05-20, pass 1)`
for the 8 reviewed entries. Eleven `CALIBRATION-TODO` markers still
pending (Phase 5 qc flags, adapter blacklist composition, strand
report granularity).

**Phase 6 benchmark inputs (structural flags Codex raised, read-only).**
These are NOT actionable in Phase 2 calibration; they inform the Phase
6 benchmark dataset selection so empirical ground truth covers the
known edge cases:

1. **Inline-barcoded / messy round-map cases**: earliest-round
   consensus can fail when demux/round-map assumptions are wrong. The
   `--round-map` requirement (locked plan line 289) mitigates the
   common case, but the benchmark should include at least one messy-
   barcode dataset to stress the failure mode.
2. **One-sided primer detection vs partial homology**: orientation
   classification based solely on `primer_5p` vs `revcomp(primer_3p)`
   can be ambiguous if only one primer is inferred, or if the primers
   share local homology. Benchmark should include a dataset with
   homologous flanking regions.
3. **First-13-nt adapter blacklist** is appropriately conservative,
   but rare real primers sharing an Illumina-like prefix would be
   dropped. Acceptable trade-off; just documented for the methods
   section.

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

### Phase 5 — QC plots + flags + count CLI (2026-05-20)

Closes the v0.1 single-dataset CLI surface. The workflow is now
feature-complete end-to-end: ``detect`` -> ``extract`` -> ``count`` ->
``qc``. ``selexprep run`` (batch driver) and ``selexprep fetch``
(accession download) are deferred to Phase 6 / v0.2.

- **NEW: `selexprep.qc.diversity`** — `rarefy()` (deterministic
  multivariate hypergeometric subsampling via numpy default RNG;
  seeded for reproducibility), `shannon_entropy()` (base-2),
  `unique_count()`, `top_n_coverage()`. Used by both `flags.py` and
  `plots.py`.
- **NEW: `selexprep.qc.flags`** — eight depth-aware suspicion flags
  per locked plan lines 350-358:
  `unexpected_rarefied_diversity_increase` (rarefied uniques per
  round; not raw counts),
  `low_primer_match` (threshold imported from
  `library.detect.UNABLE_TO_EXTRACT_MATCH_RATE` — single source of
  truth),
  `n_length_variation_across_rounds`,
  `strand_mix` (from Phase 3 strand_report.tsv),
  `low_total_reads`,
  `adapter_contamination_high`,
  `extraction_mode_changed_across_rounds` (v0.1 inert — placeholder
  for Phase 6 `selexprep run` batch driver),
  `requires_read_merging_for_full_insert` (informational).
  `compute_all_flags()` aggregator + `write_flags_yaml()` deterministic
  emitter (sorted by flag name; float-rounded for cross-platform
  stability).
- **NEW: `selexprep.qc.plots`** — four per-dataset matplotlib plots
  (PNG, 150 DPI, Agg backend, tight bbox): read_retention.png,
  primer_match_per_round.png, n_length_distribution.png,
  per_round_panel.png. **Plots are informational only** — matplotlib
  PNG output is not byte-deterministic across versions; they do not
  contribute to `output_sha256`.
- **NEW: `selexprep.qc.runner`** — `run_qc(manifest_path, ...)`
  orchestrator. Auto-discovers `round_*/counts.parquet` under the
  manifest's directory; reads `trim_reports.json` for read-retention
  plot data; optionally reads `strand_report.tsv` for the strand-mix
  flag. Returns `QcResult` with the list of flags raised, the
  flags.yaml path, and the four plot paths.
- **NEW: `selexprep.count.counter.count_fasta`** — FASTA-aware
  per-round counter (the Phase 3 extract pipeline emits FASTA, not
  FASTQ). Reuses `_counter_to_parquet` so the output schema matches
  `count_round`: `sequence`, `reads`, `rank`, `rpm`.
- **WIRED: `selexprep count <extracted-fasta> --round R0 --outdir OUT`**
  — accepts `R0`, `r0`, `round_0`, or just `0` for the round label;
  writes to `OUT/round_NN/counts.parquet`.
- **WIRED: `selexprep qc <manifest> [--counts-dir DIR] [--outdir OUT]`**
  — prints a one-line summary plus per-flag severity to stdout;
  emits `flags.yaml` and the four PNG plots.
- **Tests added (+56, 358 + 1 xfailed total):**
  `tests/test_diversity.py` (20 — rarefy determinism + edge cases,
  Shannon entropy, top-N coverage monotonicity, depth-aware sanity
  check), `tests/test_flags.py` (21 — positive + negative case per
  flag + aggregator + YAML determinism), `tests/test_plots.py` (5 —
  PNG smoke), `tests/test_qc_runner.py` (6 — end-to-end manifest ->
  flags.yaml + 4 PNGs with realistic synthetic data), `tests/test_cli.py`
  (+4 — count + qc smoke).
- **CALIBRATION-TODO inventory: 19** (was 12). New tunables in
  `qc/flags.py`: rarefaction depth, max modal lengths, strand-mix
  max reverse fraction, low-total-reads minimum, adapter-contamination
  max fraction; in `qc/plots.py`: top-N coverage N. The match-rate
  threshold is **imported** from `library.detect` (not redeclared)
  so a single Codex pass tunes both QC and classifier.

### Phase 4 — Manifest + inspect + extract override/rebuild (2026-05-19)

Closes the v0.1 CLI surface (except QC, Phase 5). Every `extract` run
now emits a `selexprep_manifest.json` — the reproducibility anchor that
Phase 5 `qc` and future v0.2 AnnData export will consume.

- **NEW: `selexprep.manifest`** — `SelexprepManifestV1` pydantic model
  (frozen, extra=forbid) with the locked schema (plan lines 162-175):
  `manifest_version`, dep versions (selexprep / python / cutadapt /
  dnaio / pyarrow), provenance (accession / bioproject_id / runs),
  `input_sha256` + `output_sha256` (FASTA/TSV/JSON only — Parquet
  hashes intentionally absent per locked plan line 28), nested
  `LibraryReport` + denormalized scan fields, CLI argv capture in
  `parameters`, runtime/flags/sampling_seed. Helpers: `compute_sha256s`,
  `write_manifest_json` / `read_manifest_json` with deterministic JSON
  (same numeric-int-key + alphabetical-sha256-keys discipline as
  `library/report.py`), `build_manifest_from_extract_result`.
- **NEW: `selexprep.fetch.inspect`** — `inspect_accession()` hits ENA
  Portal filereport REST (`https://www.ebi.ac.uk/ena/portal/api/filereport`);
  parses run/study metadata into `InspectReport` + `RunFileInfo`
  dataclasses; tolerant of missing fields. Reports
  `library_strategy` / `library_source` **verbatim from SRA** — NOT a
  DNA/RNA classification (locked plan line 332 explicit on this;
  classification deferred to v0.2's library-type-classifier).
- **EXTENDED: `selexprep.extract.runner`** — `run_extract()` now accepts
  `override_primer_{5p,3p}` (cloned via `LibraryReport.model_copy`),
  plus provenance kwargs (`accession`, `bioproject_id`, `runs`,
  `parameters`) for the auto-emitted manifest. Override without
  `--rebuild` routes outputs to `<outdir>/overridden/` (preserves
  baseline). Override + `--rebuild` overwrites baseline AND emits
  `extract_diff.tsv` comparing baseline vs override per-round read
  counts. The diff TSV is read from the baseline `selexprep_manifest.json`
  + `trim_reports.json` BEFORE overwrite; gracefully degrades if either
  baseline artifact is missing/malformed.
- **WIRED: `selexprep inspect <accession>`** — full CLI. Prints a
  human-readable metadata summary; `--outdir` also writes a sorted-keys
  `inspect.json`. `--timeout-s` controls the HTTP timeout (default 30s).
- **EXTENDED: `selexprep extract`** — `--override-primer-{5p,3p}` now
  works (lifts the Phase-3 informative error); CLI argv is captured
  into the emitted manifest's `parameters` field.
- **Tests added (+28, 302 + 1 xfailed total):** `tests/test_manifest.py`
  (12 — schema fields + frozen + extra-forbid + deterministic JSON +
  int-key sort + sha256 helper FASTA/TSV/JSON-only behavior),
  `tests/test_inspect.py` (8 — mocked ENA REST: single run / paired
  semicolon lists / multi-run study / missing-field tolerance / empty
  response → ValueError / HTTP error propagation / timeout pass-through
  / JSON sort-keys),  `tests/test_extract_override.py` (7 — override 5p
  + 3p without rebuild → subtree, override + rebuild → diff TSV + in-
  place overwrite, rebuild alone → no diff, manifest emission, override
  primer recorded in manifest), `tests/test_cli.py` (+2 — inspect smoke
  with mocked REST + override smoke). The old "Phase 4 error" test
  replaced with a real override-works smoke test.
- **CALIBRATION-TODO inventory: 12** (unchanged — Phase 4 is
  serialization + REST + I/O wiring, no new heuristic thresholds).

### Phase 3 — extract: paired-end + strand orientation (2026-05-19)

Turns the Phase 2 `LibraryReport` contract into actual extracted FASTAs.
Cutadapt is invoked as a subprocess (per locked plan); dnaio is available
for paired I/O. **No read merging in v0.1** — paired-end split-primer
mode emits two separate files and flags `READ_MERGING_RECOMMENDED`.

- **NEW: `selexprep.extract.trim`** — cutadapt subprocess wrapper with
  four public entry points (`trim_single_end_linked` for
  BOTH_PRIMERS_SINGLE_READ, `trim_single_end_5p` for FIVE_PRIME_ONLY,
  `trim_single_end_3p` for THREE_PRIME_ONLY, `trim_paired_split` for
  PAIRED_END_SPLIT_PRIMERS). Each returns a `TrimReport` carrying the
  exact cutadapt argv + read counts (Phase 4 manifest precursor).
  Cutadapt writes uncompressed FASTA; this module re-gzips with
  `_io.open_gzip_text_deterministic` (mtime=0 header) so `output_sha256`
  is bit-identical across reruns.
- **NEW: `selexprep.extract.strand`** — strand-orientation handler.
  `detect_strand_distribution()` counts forward/reverse/ambiguous reads;
  `reorient_fastq_gz()` reverse-complements every record (sequence +
  reversed quality string) for `LibraryReport.orientation == "REVERSE"`;
  `write_strand_report()` emits a sorted TSV for the QC trail.
- **NEW: `selexprep.extract.runner`** — `run_extract()` orchestrator.
  Refuses if `LibraryReport.status == "UNABLE_TO_INFER"` or
  `extraction_mode == "UNABLE_TO_EXTRACT"` (no silent miscalls).
  Optional sample-sheet pre-step demuxes multiplexed input; strand
  pre-step rewrites all reads when orientation is `REVERSE` and emits
  `strand_report.tsv` for `MIXED` or `REVERSE`. Per-mode trim dispatch
  writes per-round outputs to `<outdir>/round_NN/<filename>.fasta.gz` +
  `trim_reports.json` (manifest precursor).
- **Output filename contract** (locked plan lines 321-326):
  `extracted.fasta.gz` (full insert), `partial_5p_extracted.fasta.gz` /
  `partial_3p_extracted.fasta.gz` (one-sided), `partial_5p_extracted_R1.fasta.gz`
  + `partial_3p_extracted_R2.fasta.gz` (paired split). Filenames signal
  to downstream ML pipelines whether a full insert was recovered;
  joining R1+R2 by read ID alone is biologically wrong, so
  `joined_counts.tsv` is **not** emitted in v0.1.
- **WIRED: `selexprep extract`** — full CLI. Accepts `--library-report`,
  `--round-map`, `--sample-sheet`, `--paired-r2`, `--rebuild`.
  `--override-primer-{5p,3p}` emit a Phase-4 informative error (full
  diff TSV lands in Phase 4). `--rebuild` toggles the no-clobber guard.
- **Tests added (+37, 274 + 1 xfailed total):** `tests/test_strand.py`
  (13 — distribution + revcomp + deterministic gzip + TSV sort),
  `tests/test_trim.py` (7 — per extraction_mode + determinism + temp
  cleanup, skips if cutadapt absent), `tests/test_extract_runner.py`
  (13 — happy path per mode + UNABLE refusal + no-clobber + rebuild +
  strand-report emission + trim_reports JSON + multi-round),
  `tests/test_cli.py` (+4 — missing-round-map, override Phase-4 error,
  UNABLE refusal, end-to-end smoke).
- **CALIBRATION-TODO inventory: 12** (was 11). New tunable:
  `STRAND_REPORT_PER_READ = False` in `extract/strand.py`. Strand
  classification thresholds (`ORIENTATION_REVERSED_FORWARD_MAX`,
  `ORIENTATION_REVERSED_REVERSE_MIN`) stay in Phase 2's
  `library/detect.py` — no duplication.

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
