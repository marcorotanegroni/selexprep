# selexprep benchmark (Phase 6b)

This directory implements selexprep's two-tier benchmark:

- **Tier 1 — Figure A** (`Snakefile` + `ground_truth.tsv`): curated
  primer-recovery validation against paper-reported primers on N=11
  source-verified accessions.
- **Tier 2 — Figure B** (`audit.smk` + the bundled catalog): a corpus
  audit / utility characterization across N≈30 sampled INSDC
  accessions, no per-row ground truth. Distributional metrics only.

Both tiers share the metric + figure entry points under
`src/selexprep/benchmark/`.

## What this benchmark tests

> **Given only public/local HT-SELEX reads plus accession metadata,
> can `selexprep` infer the primer / constant regions and N-region
> length, then report confidence and fail safely when inference is
> ambiguous?**

This is the **unique claim** of selexprep. Existing tools (AptaPLEX,
EasyDIVER+, FASTAptameR) require the user to supply the primers as
input — they cannot answer this question by construction.

## What this benchmark deliberately does NOT test

- **Comparator-tool count agreement.** AptaPLEX / EasyDIVER+ need the
  paper-reported primers as input, so a head-to-head on count tables
  reduces to "does selexprep's trimming subprocess produce the same
  output as another tool given the same primers?" — a trimming-code
  sanity check, not a meaningful scientific comparison.
- **Downstream count correlation.** Same reasoning. The honest version
  of this question is self-consistency: do inferred-primer counts
  match `--override-primer`-driven counts of the same dataset? That's
  Phase 6c (see *Roadmap* below).

## Scope split: 6b.1 → 6b.2 → 6b.3a → 6b.4 → 6c

| Phase | Deliverable | Status |
|---|---|---|
| **6b.1** | Tier 1 schema + metrics + Figure A + `Snakefile` | ✅ shipped |
| **6b.2** | Curate verified rows from public papers via WebFetch (+ENA XML PubMed links + Europe PMC search → PMC full-text → M&M extraction) | ✅ 11 rows shipped |
| **6b.3a** | Tier 2 corpus-audit **pipeline**: `benchmark.corpus_audit` + `benchmark.figure_b` + `audit.smk` + tests | ✅ scaffolding shipped (this section) |
| **6b.3** | Re-run **Tier 1** benchmark on HPC; commit `results/figure_a.{pdf,png}` + measured `metrics.json` | data work |
| **6b.4** | Run **Tier 2** audit on HPC; commit `audit_results/audit_accessions.tsv` + `audit_metrics.json` + `figure_b.{pdf,png}` (no code changes) | data work |
| **6c (optional)** | Self-consistency check: inferred-primer counts vs `--override-primer` counts (Pearson on union+zero-fill) | post-v0.1 |

## Headline metrics

The Snakefile aggregates these per run; Figure A renders the highlights.

- **Primer recovery** (per-side equivalence): exact, revcomp,
  U-T-normalized, barcode-stripped, partial 5'/3', mismatch, IUPAC
  unsupported.
- **Pair recovery by status**: of HIGH-confidence calls, what fraction
  recovered the pair exactly? Same for MEDIUM / LOW / UNABLE_TO_INFER.
  This is the Codex-honest headline panel.
- **Safe-failure rate**: rows where selexprep refused
  (`status=UNABLE_TO_INFER`, `extraction_mode=UNABLE_TO_EXTRACT`, or
  `required_action=MANUAL_PRIMERS_REQUIRED`). The unique
  distinguishing metric vs known-primer pipelines.
- **N-length recovery within tolerance**.
- **Honest accounting**: `extraction_mode` + `required_action`
  distributions (locked plan line 369).

## Codex pre-implementation amendments folded in

A Codex peer-review pass before any code landed corrected seven
design choices. They are referenced inline in the source files:

1. **Verified-only ground truth.** `ground_truth.tsv` contains
   `verified=true` rows only. Unverified candidates live in the
   *Candidate worklist* below and stay out of the metrics until a
   curator confirms their primer sequences.
2. **Accession-fetchable only.** The Snakefile invokes
   `selexprep fetch`, which v0.1 supports for ENA / SRA / DDBJ.
   Figshare / Zenodo / processed-data sources are deferred to v0.2.
3. **Curated round-map override.** Rows where ENA metadata is too
   sparse for auto round-parsing can set
   `round_map_source=curated` and `round_map_path=<path>`. The
   Snakefile then passes `--allow-manual-review` to fetch and
   routes the curator's TSV to `selexprep detect`.
4. **No SciPy.** Pearson via `pandas.Series.corr(method="pearson")`;
   Spearman via Pearson-of-ranks (`.rank().corr(...)`).
5. **Union + zero-fill count correlation** (kept as a Phase 6c
   entry point in `metrics.py`; not called from
   `aggregate_metrics`). Intersection-only correlation would
   inflate apparent agreement.
6. **Figure A: 4-panel primer-inference focus.** Pair recovery
   x status, per-side breakdown, N-length recovery, and
   `extraction_mode` + `required_action` distribution. Count
   correlation moves to Phase 6c.
7. **Library-kind verification.** Each row's `library_kind` is
   encoded from the paper's Materials & Methods, not the catalog
   title alone (pyoverdine PRJNA932049 was specifically cross-
   checked because Codex flagged a possible DNA-vs-RNA ambiguity;
   the ENA project XML confirms 2'-FY RNA).

### Scope pivot vs locked plan (2026-05-22)

The locked plan (`~/.claude/plans/unified-seeking-treehouse.md` line
365–366) listed count correlation Pearson+Spearman AND AptaPLEX /
EasyDIVER+ as comparators. Both are deprioritized here:

- **Comparators are not benchmarked.** Both require known primers as
  input → cannot test selexprep's unique claim.
- **Count correlation moves to Phase 6c** as a self-consistency
  check (inferred vs `--override-primer` runs of selexprep itself),
  not a comparator-tool oracle. The `CountCorrelationReport`
  dataclass + `compute_count_correlation` function remain in
  `metrics.py` for the 6c entry point; `aggregate_metrics` does not
  call them in 6b.1.

This deviation is conscious and serves the headline claim more
directly. The Application Note framing becomes:

> selexprep recovers published HT-SELEX primer / constant regions
> from accession-derived reads and fails safely on ambiguous
> datasets, removing the manual primer-curation step required by
> existing tools.

## Ground-truth schema

```
accession            ENA/SRA/DDBJ accession that selexprep fetch supports.
library_kind         RNA | DNA | 2'-FY RNA | ssDNA | ...
target_kind          protein | cell | small molecule | DNA (TF) | ...
primer_5p_truth      Paper-reported 5' constant primer (DNA letters).
primer_3p_truth      Paper-reported 3' constant primer (DNA letters).
n_length_truth       Paper-reported random N-region length (integer; 0 if unmeasurable).
paper_doi            10.xxxx/yyyy
paper_pmid           PubMed ID (integer string).
round_map_source     "auto" (default) | "curated".
round_map_path       Relative path (from this file) to a curator-validated rounds.tsv.
verified             "true" iff a curator has confirmed primer_5p/3p_truth against the paper.
notes                Curation citation + any caveats.
```

The metric aggregator filters `verified=true` rows;
`verified=false` rows would emit a stderr warning per skip. To keep
`ground_truth.tsv` clean (Codex amendment 1), unverified rows are
NOT stored here — candidates live in the worklist below until
promoted.

## Tier 2: corpus audit **pipeline**

> **Note (6b.3a):** what ships here is the **pipeline**, not the
> results. The Snakefile + module are reviewable + runnable, but CI
> does not execute the Snakefile (real-data fetch is heavy and not a
> CI workload). The committed `audit_metrics.json` +
> `audit_accessions.tsv` + `figure_b.{pdf,png}` will land in the
> **6b.4 follow-up commit** with no code changes — just the run
> output + curator review.

### What Tier 2 measures

Across N≈30 INSDC accessions sampled deterministically from the
bundled catalog, the audit reports the **distributions** of:

- **Fetch outcomes** — every `RunStatus` value (`OK`,
  `FETCH_FAILED`, `FETCH_REFUSED`, `EXTRACT_REFUSED`, etc.). Answers
  "can public data be obtained?".
- **LibraryReport status** — HIGH / MEDIUM / LOW / UNABLE_TO_INFER
  for rows where `detect` ran successfully.
- **Extraction mode** — `BOTH_PRIMERS_SINGLE_READ` /
  `FIVE_PRIME_ONLY` / `THREE_PRIME_ONLY` /
  `PAIRED_END_SPLIT_PRIMERS` / `UNABLE_TO_EXTRACT`.
- **Required action** — `NONE` / `MANUAL_PRIMERS_REQUIRED` /
  `READ_MERGING_RECOMMENDED`.
- **Inference safe-failure rate** — the headline number. The fraction
  of rows-with-a-LibraryReport where selexprep refused (UNABLE_TO_INFER
  OR UNABLE_TO_EXTRACT OR MANUAL_PRIMERS_REQUIRED).

### What Tier 2 does NOT claim

This is **descriptive distributional metrics, no per-row ground
truth**. The audit cannot say a given accession's primer call is
"correct" or "incorrect" because Tier 2 accessions have no
paper-reported primers attached (that's Tier 1's domain). The
methods text in the paper says so explicitly.

### Per-panel denominators (Figure B)

Each panel labels its denominator in the subtitle so a reviewer
never has to guess the normalization. Two distinct denominators live
in the same figure:

| Panel | Question | Denominator |
|---|---|---|
| **A · Fetch outcomes** | Can public data be obtained? | `n_sampled` |
| **B · Inference confidence** | Given detect ran, what did the report say? | `n_with_library_report` |
| **C · Extraction mode** | Honest accounting of inferred biology | `n_with_library_report` |
| **D · Required action** | Workflow guidance + safe-failure rate overlay | `n_with_library_report` |

The methodological correction (Codex peer-review + user pass):
`inference_safe_failure_rate` is computed **only** among rows with a
LibraryReport. Mixing fetch failures into this rate would inflate the
metric with ENA / network / regional-restriction problems —
conflating "selexprep refused" (a feature) with "the dataset was
unreachable" (an external problem).

### Reproducing the Tier 2 audit

```bash
snakemake -s audit.smk --cores 4
```

Outputs land in `audit_results/`:

- `audit_accessions.tsv` — the sampled accession list (first-class
  artifact; committed in 6b.4 so the figure is reproducible across
  future catalog refreshes).
- `audit_accessions.manifest.json` — sidecar with the reproducibility
  envelope (`catalog_version`, `sample_seed`, `sample_accessions_sha256`).
- `run_summary.tsv` — `selexprep run --resume` per-accession status.
- `audit_metrics.json` — distributional metrics (sorted-keys JSON;
  the source of truth — matplotlib PNG bytes are not deterministic).
- `figure_b.{pdf,png}` — the four-panel Figure B.

Sampling determinism: `sample_corpus(n, seed)` deterministically
samples from `filter_catalog(insdc_only=True)`. A future `selexprep
catalog refresh` will shift row indices and the same seed will draw
different accessions — which is why `audit_accessions.tsv` ships as
a first-class committed artifact in 6b.4: reviewers reproduce the
figure by `selexprep run audit_accessions.tsv` regardless of any
catalog drift.

## Reproducing the benchmark

```bash
# from this directory, with verified rows in ground_truth.tsv:
snakemake -s Snakefile --cores 4
```

Dry-run to inspect the DAG:

```bash
snakemake -s Snakefile --cores 1 --dry-run
```

Snakemake is in the optional `bench` extra:

```bash
uv sync --extra bench
# or
pip install -e ".[bench]"
```

## Currently verified (11 rows, modality-diverse)

See `ground_truth.tsv` for the full table. Each row was source-verified
during a post-curation audit (DOI/PMID confirmed against PubMed esummary;
primer sequences extracted from main text M&M, PMC-bundled supplement
PDFs via pdfplumber, supplement docx via python-docx, RSC ESI PDFs, or
**visual inspection of raster-rendered supplement tables via
PyMuPDF/fitz** when text extraction silently drops the table content;
no row enters this table unless the exact primer sequence is
source-backed).

| Accession | Library | Target | Paper (PMID) | N | Source location |
|---|---|---|---|---|---|
| `PRJNA315881` | 2'-F-Py RNA | protein (IL-10RA) | Hoinka 2015 NAR (25870409) | 40 | PMC4499121 main M&M |
| `PRJEB22637` | 2'-F-Py RNA | cell (Annexin A2 / ETBR) | Cibiel 2014 PLOS ONE (24489826) | 50 | PMC3906106 main M&M |
| `PRJEB28411` | DNA | cell (ccRCC RCC-MF) | Pleiko 2019 Sci Rep (31148584) | 40 | PMC6544647 main M&M |
| `PRJDB9110` | RNA (T7-tx) | protein (TG2) | Ozaki 2020 NAR — RaptRanker Data1 (32537639) | 30 | PMC7641312 main M&M |
| `PRJDB9111` | RNA (T7-tx) | protein (integrin αVβ3) | Ozaki 2020 NAR — RaptRanker Data2 (32537639) | 40 | PMC7641312 main M&M |
| `PRJEB70964` | 2'-F-Py RNA | protein (α-syn fibrils) | NAR 2024 footprints (38917326) | 35 | PMC11317169 main M&M |
| `PRJNA728693` | DNA | small molecule (kanamycin A) | Sanford 2021 Chem Sci — RE-SELEX (34659704) | 40 | RSC ESI PDF Fig. 2 (pdfplumber) |
| `PRJNA883192` | DNA | protein (EPX) | Ali 2022 Sci Rep eosinophil peroxidase (36577785) | 40 | MOESM1 supp PDF page 2 Table S1 (PyMuPDF raster render) |
| `PRJNA935703` | 2'-FY RNA | small molecule (pyoverdine PYO-Pf5) | Anisuzzaman 2024 Front Chem (39148668) | 55\* | PMC11324436 main M&M |
| `PRJNA975735` | DNA | protein / viral peptide (SARS-CoV-2 RBD) | Halder 2023 Sci Rep (37666993) | 40 | PMC10477244 main M&M |
| `PRJNA990511` | DNA | protein (ASFV p30) | Hu 2024 Sci Rep (38374125) | 40 | MOESM1 supp docx Table S3 (python-docx) |

\* PRJNA935703 has a split-random-region library (N10 + internal
constant + N35). The N=55 is the total length between the OUTER
constants, which is what selexprep would measure.

The α-synuclein row's 5' constant matches the revcomp of the TruSeq
R1 adapter prefix that selexprep blacklists — it intentionally
exercises the adapter-trap path. The pyoverdine row tests an
unusual split-N library. Both are valuable edge cases for the
benchmark.

## Candidate worklist — see `candidates.tsv`

Unverified or rejected candidates are tracked in
`benchmarks/candidates.tsv` (NOT in `ground_truth.tsv` per Codex
amendment 1). Each row records: accession, paper_doi, paper_pmid,
the source attempted, the (unrecovered) primer fields, status
(`blocked` / `rejected`), and a structured rejection reason.

Summary of 9 documented candidates:

| Accession | Status | Reason |
|---|---|---|
| `PRJNA321551` | blocked | Dao 2016 *Cell Systems* AptaTRACE — primers in supplemental experimental procedures; Cell Press paywall on supp PDF |
| `PRJEB47428` | blocked | PRIESSTESS NAR — primers in Supp Table S1; OUP blocks bot access to supp |
| `PRJEB9897` | blocked | SELMAP Sci Rep — neither main text nor supp PDF (inspected) prints the adapter/key/barcode sequences |
| `PRJNA1054605` | blocked | Penzar PADIT-seq bioRxiv — license bars PMC archiving; bioRxiv blocks bot PDF download |
| `PRJNA1071201` | **rejected** | Blocker-SELEX — virtual screening of a 1024-sequence pool, NOT a random-N library; outside scope |
| `PRJEB89545` | blocked | DL-SELEX 2025 — recent deposit, no publication cross-link yet |
| `PRJNA932049` | blocked | Iowa State pyoverdine — GEO marks citation as missing; author contact required (a separate pyoverdine SELEX, PRJNA935703, IS verified — see ground_truth.tsv) |
| `PRJNA1081432` | **rejected** | Anti-EGFR Nat Precision Oncol — supp PDF Table 2 confirms paper uses existing aptamers (EGFRapt/MinE07, C36, E3, Waz extensions), not a SELEX library |
| `PRJNA558191` | blocked | Camorani 2020 iScience TNBC cell-SELEX — main text Fig. 1 caption confirms N40 + 21+23 nt constants; exact sequences in Transparent Methods supp not retrievable via current methods (try PyMuPDF raster render if the supp turns out to be a raster table) |

Curation lesson learned during audit: text-extraction tools (`pypdf`,
`pdfplumber`) silently drop table contents when the table is rendered
as a raster image inside the PDF. A failed text extraction does NOT
mean the data is absent — try `PyMuPDF` (`fitz`) to render the page
as an image and visually inspect. Ali 2022 Sci Rep (PRJNA883192) was
initially demoted to this list when pdfplumber returned an empty
Table S1; the row was promoted back to `ground_truth.tsv` after
visual inspection of the raster-rendered page 2.

To promote a blocked candidate: open the paper's supplementary PDF
(usually downloadable from the journal page in a browser session),
extract the exact primer-sequence-containing sentence, add a row
to `ground_truth.tsv` with `verified=true` and the citation in
`notes`, and remove the row from `candidates.tsv`. The metrics
aggregator picks it up on the next Snakefile run with no code
changes needed.

To find new candidates not in this list: scan the bundled catalog
(`src/selexprep/catalog/data/bioprojects.csv`, 273 rows) for
HT-SELEX studies, look up the ENA XML for `XREF_LINK` → `PUBMED`,
then check Europe PMC for PMC full-text availability:

```bash
curl -s "https://www.ebi.ac.uk/ena/browser/api/xml/<accession>" | grep -A1 PUBMED
curl -sL "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=<pmid>&format=json&tool=...&email=..."
```

`pypdf` + Europe PMC's `/<PMCID>/supplementaryFiles` ZIP endpoint
is a fruitful path for journals that bundle their supps with the
PMC record (BMC, Sci Rep, Nat Precision Oncol). NAR and Cell Press
do not — for those, manual browser-session downloads are required.
