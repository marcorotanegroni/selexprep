# selexprep benchmark

This directory holds selexprep's benchmark: two tiers plus a standalone
catalog-completeness audit. Each tier emits deterministic, sorted-key JSON —
`metrics.json` (Tier 1) and `audit_metrics.json` (Tier 2) — which the paper
presents as a per-deposit **scorecard table** and a corpus-audit **summary table**.

- **Tier 1 — primer recovery** (`Snakefile` + `ground_truth.tsv`): per-deposit
  recovery of paper-reported primers from accession-derived reads, on 11
  source-verified deposits. Produces the scorecard (paper **Table 1**).
- **Tier 2 — corpus audit** (`audit.smk` + the bundled catalog): a descriptive
  utility audit over a deterministic sample of audit-eligible INSDC accessions.
  Distributional metrics only — no per-row ground truth. Paper **supplement**.
- **Catalog-completeness audit** (`catalog_completeness_audit.py`): diffs
  `bioprojects.csv` against ENA's `library_strategy="SELEX"` set.

Tier 1 + Tier 2 share metric/table entry points under `src/selexprep/benchmark/`.

## What this benchmark tests

> **Given only public/local HT-SELEX reads plus accession metadata, can
> `selexprep` infer the primer / constant regions and N-region length, report
> its confidence, and fail safely when inference is ambiguous?**

That is selexprep's unique claim. Tools like AptaPLEX, EasyDIVER+, and
FASTAptameR require the user to supply primers as input — they cannot answer
this question by construction.

## What it deliberately does NOT test

- **Comparator-tool count agreement.** AptaPLEX / EasyDIVER+ need the
  paper-reported primers as input, so a head-to-head on count tables reduces to
  "does selexprep's trimming subprocess match another tool given the same
  primers?" — a trimming sanity check, not a scientific comparison.
- **Downstream count correlation.** The honest version is self-consistency:
  do inferred-primer counts match `--override-primer`-driven counts of the same
  dataset? The `CountCorrelationReport` + `compute_count_correlation` scaffolding
  stays in `metrics.py` as a future entry point; `aggregate_metrics` does not
  call it.

## Tier 1 — primer-recovery scorecard

**Headline:** *selexprep recovered the random region at exactly the
paper-reported length on every recovery deposit where a single-read extraction
is possible (6/6, 0 out of tolerance), recovered the paper-reported primer
strings exactly on 4 of 7 and with informative partial recovery on the other 3,
and made zero false-positive primer calls on the 4 primer-absent /
adapter-collision controls.*

Two distinct measurements are reported, because they answer different questions.
The **random-region boundary** is what the tool must get right to trim
correctly; the **exact primer string** additionally requires that the constant
called from the reads coincides with the constant as written in the paper, which
can differ when the deposit carries extra constant technical sequence.

The source-verified deposits, by arm (the table below is the per-row scorecard;
`snakemake -s Snakefile` regenerates it as `metrics.json`). Rows marked
*pending* were added to balance the arms and have not yet been through a
benchmark run; their scores are blank until they have.

| Accession | Chemistry | Target | Arm | `status` | 5′ / 3′ | *N* obs / truth | Note |
|---|---|---|---|---|---|---|---|
| PRJDB9110 | RNA (T7) | TG2 (protein) | recovery | HIGH | EXACT / EXACT | 30 / 30 | counted |
| PRJDB9111 | RNA (T7) | αVβ3 (protein) | recovery | HIGH | EXACT / EXACT | 40 / 40 | counted |
| PRJDB19098 | RNA (T7) | FGF-9 (protein) | recovery | HIGH | EXACT / EXACT | 35 / 35 | counted; truth from patent WO2020204151 |
| PRJNA615076 | DNA | *E. faecalis* (cell) | recovery | HIGH | EXACT / EXACT | 40 / 40 | counted |
| PRJEB62495 | DNA | *Anaplasma* (cell) | recovery | HIGH | EXACT / MISMATCH | 40 / 40 | counted (partial); read-level 3′ boundary divergence |
| PRJNA1395820 | DNA | Co²⁺ (small molecule) | recovery | MEDIUM | PARTIAL / EXACT | — | counted (partial); paired split primers → merge recommended |
| PRJNA809588 | 2′-F-Py RNA | islet (cell) | recovery | MEDIUM | PARTIAL / PARTIAL | 40 / 40 | counted (partial) |
| PRJEB22637 | 2′-F-Py RNA | cell (Annexin A2) | specificity | UNABLE_TO_INFER | null / null | — | correct refusal (N-region-only reads) |
| PRJEB28411 | DNA | cell (ccRCC) | specificity | UNABLE_TO_INFER | null / null | — | correct refusal |
| PRJNA990511 | DNA | protein (ASFV p30) | specificity | UNABLE_TO_INFER | null / null | — | correct refusal (single-round) |
| PRJEB47428 | RNA | RNA-binding proteins | specificity | *pending* | — | — | reads are N40 exactly (92 runs) |
| PRJEB49150 | DNA | BEN-domain TFs | specificity | *pending* | — | — | reads are N40 exactly (9 runs) |
| PRJEB14550 | DNA | HOXB13 / FLI1 | specificity | *pending* | — | — | reads are N40 exactly (8 runs) |
| PRJNA360902 | DNA | *Ciona* TF DBDs | specificity | *pending* | — | — | reads are N20 exactly (14 runs) |
| PRJEB70964 | 2′-F-Py RNA | protein (α-syn) | adapter-control | UNABLE_TO_INFER | null / null | — | correct refusal (5′ const = revcomp TruSeq R1) |
| PRJNA678231 | n/a (ncRNA-Seq) | n/a (*B. mori*) | adapter-control | *pending* | — | — | 51 nt reads over ~20–30 nt inserts |
| PRJDB7022 | n/a (miRNA-Seq) | n/a (*D. melanogaster*) | adapter-control | *pending* | — | — | 51 nt reads over ~20–30 nt inserts |
| PRJNA746278 | n/a (ncRNA-Seq) | n/a (*H. sapiens*) | adapter-control | *pending* | — | — | 59 nt reads, SMARTer smRNA kit |
| PRJDB2183 | n/a (miRNA-Seq) | n/a (*A. thaliana*) | adapter-control | *pending* | — | — | 69 nt reads over ~20–30 nt inserts |
| PRJEB50674 | n/a (miRNA-Seq) | n/a (*T. vaginalis*) | adapter-control | *pending* | — | — | 50 nt reads over ~20–30 nt inserts |
| PRJNA591605 | n/a (ncRNA-Seq) | n/a (*M. musculus*) | adapter-control | *pending* | — | — | 75 nt reads over ~20–30 nt inserts |

**How to read it.** The *recovery* arm asks whether selexprep recovers the
paper primer from raw reads; the *specificity* and *adapter-control* arms are
negative controls where the correct behavior is **no call** (constants are
absent, or collide with a known adapter). A "partial" or "mismatch" usually
reflects the **deposit** (reads carrying an extended / heterogeneous constant
region), not a tool error — `detect` reports what is physically in the reads.
The *N* column is the mode of the extracted random-region length against the
paper-reported length; it is blank where a single-read extraction is not
attempted (PRJNA1395820 is paired-end with split primers, so `detect` asks for
read merging rather than forcing an R1-only call).

PRJNA883192 was **withdrawn from the scored set** (`verified=false`) because its
3′ constant could only be resolved from the deposited reads themselves: scoring
an inferred call against a read-derived truth is circular, and keeping it would
have inflated the recovery denominator with a case the benchmark cannot honestly
adjudicate. It stays in `ground_truth.tsv` with the full reasoning.

### How the two control arms were selected

Both arms were built from archive metadata alone, before any inference was run,
so membership cannot have been conditioned on how selexprep happened to behave.

**Specificity (pre-trimmed).** Every INSDC deposit in the curated catalog that
states a randomised-region length was screened by comparing that length against
the archive-reported read length (`base_count / read_count` per run, halved for
paired layouts — no FASTQ is downloaded). A deposit whose reads *are* the
randomised region carries no library constant, so the correct behaviour is to
make no primer call. Thirteen deposits passed; four were added, chosen for
chemistry and target diversity and for having identical read length in every
run. Deposits sharing a publication with a row already in the arm were
rejected as near-duplicates, as were deposits with mixed randomised-region
lengths by design, where the read-length argument does not hold cleanly.

**Adapter control.** Two different kinds of negative control sit here.
PRJEB70964 is the hard case — a genuine SELEX deposit whose 5′ constant is the
reverse complement of TruSeq R1, so a correct tool must not mistake one for the
other. The other six are non-SELEX small-RNA libraries (`miRNA-Seq` /
`ncRNA-Seq`, six organisms, at least three library-prep kits) where the read
runs 20–50 nt past the insert into 3′ sequencing adapter: a perfectly conserved
block sitting exactly where a library constant would sit. Any primer call on
them is a fabrication. Amplicon libraries were deliberately **not** used: a 16S
amplicon has real constant primers flanking a variable region, so calling them
would be correct behaviour, and scoring it as a false positive would punish the
right answer.

The randomised-region length for three of the four new specificity rows comes
from a publication that cites the accession; for PRJNA360902 no citing
publication was found, so its length is submitter-stated in the SRA record —
external to the reads, but archive-sourced rather than paper-sourced, and the
row says so.

**Excluded deposits** live in `excluded_datasets.tsv` with an evidence-based,
pre-inference reason (e.g. nonstandard read architecture, multiplexing without
documented barcodes). They are removed *before* analysis and documented — not
silently dropped — so the recovery denominator is honest. The pre-detect
screening decision for every candidate is recorded in `screening_log.tsv`.

### Ground-truth schema

```
accession            ENA/SRA/DDBJ accession that `selexprep fetch` supports.
library_kind         RNA | DNA | 2'-F-Py RNA | ...
target_kind          protein | cell | small molecule | ...
primer_5p_truth      Paper-reported 5' constant (DNA letters).
primer_3p_truth      Paper-reported 3' constant (DNA letters).
score_3p             "false" if the 3' truth is read-resolved (excluded from the
                     paper-grounded 3' tally); blank/"true" otherwise.
n_length_truth       Paper-reported random N-region length (0 if unmeasurable).
paper_doi / paper_pmid
round_map_source     "auto" (default) | "curated".
round_map_path       Relative path to a curator-validated rounds.tsv.
verified             "true" iff a curator confirmed the primers against the paper.
read_state           raw_standard (recovery) | pre_trimmed (specificity) | adapter_control.
notes                Curation citation + caveats.
```

Only `verified=true` rows are scored; the metric aggregator filters the rest
and warns per skip.

### Curated round-maps + paired-end policy

Two deposits (`PRJEB28411`, `PRJNA883192`) have no round structure parseable
from ENA metadata; they set `round_map_source=curated` with a hand-supplied TSV
under `round_maps/`. `rule fetch` then passes `--allow-manual-review` (so
all-unassigned runs download into `round_unknown/`) and `detect` consumes the
curated map. Round numbers there are **inferred from sample aliases** and used
only to enable recovery inference — not as per-round biological claims (primers
are constant across rounds, so recovery is robust to the exact numbering).

`rule fetch` writes separate primary/R2 manifests (`fastqs.manifest`,
`fastqs.r2.manifest`). R1 / single-end inputs go positionally to `detect`; R2
mates go via `--paired-r2`, so paired split-primer datasets are not forced into
R1-only partial recovery. Read merging is a v0.2 item.

## Tier 2 — corpus audit

Over a deterministic sample of audit-eligible INSDC accessions (only
`ELIGIBLE_HT_SELEX_ROUNDS` rows from the audit-generated
`audit_results/eligibility.tsv`), the audit reports the
**distributions** of fetch outcomes, `LibraryReport.status`, `extraction_mode`,
`required_action`, and the inference **safe-failure rate**. It is descriptive
(no per-row ground truth): it cannot call any accession's primer "correct" —
that is Tier 1's job.

Discipline: the safe-failure rate is computed **only** among rows that produced
a `LibraryReport`, so unreachable-data (ENA / network) failures don't inflate
"selexprep refused" (a feature) with "the dataset was unreachable" (an external
problem). Each summary denominator is labeled explicitly.

Shipped artifacts live in `audit_results/` with a reproducibility envelope —
`audit_accessions.tsv` + `audit_accessions.manifest.json`
(`catalog_version`, `sample_seed`, `sample_accessions_sha256`) — so a reviewer
reproduces the result via `selexprep run audit_accessions.tsv` regardless of any
later catalog drift.

## Catalog-completeness audit

A standalone, reusable script that hits ENA at the data-type level
(`library_strategy="SELEX"`) and diffs the result against `bioprojects.csv` +
`bioprojects_excluded.csv` — orthogonal to selexprep's text-pattern discovery.

```bash
uv run python -m benchmarks.catalog_completeness_audit
```

Current snapshot:

| Metric | Value |
|---|---|
| ENA `library_strategy="SELEX"` studies | 103 |
| In `bioprojects.csv` (auditable) | 82 (79.6%) |
| In `bioprojects_excluded.csv` (documented exclusions) | 21 (20.4%) |
| Unaccounted for | **0** |

Exclusions split into submission-metadata mis-labels + gSELEX/genomic-fragment
variants; per-row reason strings live in
`selexprep.catalog.rebuild.MANUAL_EXCLUSIONS`.

## Project metadata table

`build_project_metadata.py` emits a static, browsable table of *experiment
characteristics* — one row per catalog deposit — joined so each source owns one
slice of the truth and they never drift:

```bash
uv run python -m benchmarks.build_project_metadata        # hits OpenAlex; --no-network caps tier at ABSTRACT
uv run python -m benchmarks.build_project_metadata --results-dir out/   # + per-round trajectory
```

Outputs: `project_metadata.csv` (flat) + `project_metadata.json` (same rows plus
a per-round `rounds` trajectory). Passing `--results-dir` (a `selexprep run`
output tree) populates each deposit's `rounds` with `{n_reads, n_unique,
singleton_frac}` per cycle, recomputed from `<acc>/round_*/counts.parquet`;
deposits with no count run keep `rounds: null`. The trajectory is nested, so it
lives in the JSON only — the CSV stays flat. A scalar **`round_structure`**
column summarises it in both files and disambiguates the overloaded null:
`multi` / `mono` (counted), `unassigned` (run refused — no round parsable from
metadata), `not_fetchable` (discovery-only deposit), or empty (INSDC not yet
counted). Every row also carries two honesty signals so a reader always knows
how a value was obtained:

- **`curation_level`** — `verified` (the 11 benchmark deposits, primer-checked
  against the paper) / `extracted` (target, `study_type`, format derived from the
  ENA/DDBJ or figshare/Zenodo **title/abstract only — review-grade, not
  full-paper-verified**) / `none`.
- **`metadata_tier`** — `RECORD_ONLY` (no linked paper) / `ABSTRACT` / `FULL_TEXT`
  (open access). Paper DOIs for figshare/Zenodo deposits are resolved from their
  host APIs (the catalog rarely stores them); OA status from OpenAlex.

`study_type` ∈ {`aptamer_selection`, `tf_ht_selex`, `method_or_other`,
`not_selex`}. The curated source columns live in `project_annotations.tsv`
(verified 11) and `catalog_annotations.tsv` (everything else); the OpenAlex /
host-API lookups are cached in git-ignored `.oa_cache.json` /
`.discovery_doi_cache.json`.

**Canonical vs trajectory snapshot.** `project_metadata.{csv,json}` is the
**canonical, sources-reproducible** table: `build_project_metadata` (no
`--results-dir`) regenerates it byte-for-byte from the committed catalog +
annotations, so `rounds` is `null` and `round_structure` only resolves
`not_fetchable`. `project_metadata.trajectories.json` is a **dated run-product
snapshot** (committed for reference, not regenerable from sources alone): the
per-round trajectory + full `round_structure` from a `selexprep run` over the
fetchable INSDC subset. The shipped snapshot is from the 2026-06-11 run (filtered
to aptamer + method INSDC; the giant TF/GHT studies were dropped) — **18
trajectories** (17 multi, 1 mono); the rest of the attempted deposits lacked a
parsable round structure (`unassigned`) or an inferable primer, or their FASTQs
were unavailable from ENA at run time. Re-running `build_project_metadata
--results-dir <selexprep run outdir>` refreshes it.

## Reproducing

Snakemake is in the optional `bench` extra:

```bash
uv sync --extra bench            # or: pip install -e ".[bench]"

snakemake -s Snakefile --cores 4         # Tier 1 → metrics.json + table_1.md
snakemake -s Snakefile --cores 1 --dry-run
snakemake -s audit.smk  --cores 4        # Tier 2 → audit_metrics.json + table_audit.md
```

CI does not execute the Snakefiles (real-data fetch is heavy). Tier 1's
`metrics.json` + `table_1.md` are **regenerated on demand** (`snakemake -s Snakefile`) —
they are reproducible outputs of the committed code + `ground_truth.tsv`, not
committed artifacts. Tier 2 ships its result-of-record under `audit_results/`
(with the reproducibility envelope above), because its deterministic sample is
pinned to a catalog snapshot and so is not trivially regenerable.

## Curation methodology + candidate worklist

Unverified or rejected candidates are tracked in `candidates.tsv` (never in
`ground_truth.tsv`): accession, DOI/PMID, source attempted, status
(`blocked` / `rejected`), and a structured reason. To promote one, extract the
exact primer-sequence sentence from the paper, add a `verified=true` row to
`ground_truth.tsv` with the citation, and remove it from `candidates.tsv`; the
aggregator picks it up on the next run.

**Curation lesson:** text-extraction tools (`pypdf`, `pdfplumber`) silently drop
table contents when a supplement renders the table as a raster image. A failed
text extraction does **not** mean the data is absent — render the page with
PyMuPDF (`fitz`) and inspect it visually. One verified row (PRJNA883192) was
recovered exactly this way after `pdfplumber` returned an empty Table S1.
