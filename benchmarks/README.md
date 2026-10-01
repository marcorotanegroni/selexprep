# selexprep benchmark

This directory holds selexprep's benchmark: two tiers plus a standalone
catalog-completeness audit. Each tier emits deterministic, sorted-key JSON —
`metrics.json` (Tier 1) and `audit_metrics.json` (Tier 2) — which the paper
presents as a per-deposit **scorecard table** and a corpus-audit **summary table**.

- **Tier 1 — primer recovery** (`Snakefile` + `ground_truth.tsv`): per-deposit
  recovery of paper-reported primers from accession-derived reads, on 21
  source-verified deposits in three arms (7 recovery, 8 specificity,
  6 adapter-control). Produces the scorecard (paper **Table 1**).
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

**Two results, kept apart.** The benchmark was specified and first scored with
v0.4.1. Its per-run extraction yield then exposed a defect in `detect` (see
[Extraction yield](#extraction-yield-per-run)), which was fixed; the fix changes
three recovery calls. The v0.4.1 result is the pre-specified one. The post-fix
result is measured on the same deposits that revealed the defect, so it is
in-sample: the rule and its constants were committed before this benchmark was
rerun and were not tuned afterwards, but they were designed with these deposits
in view. Evidence that the rule generalises has to come from data that did not
motivate it. Two later changes to how reads are anchored on the constant's core
(a106c3d and b8ded35, found in code review, not on these deposits) were rerun on
the full benchmark: calls, metrics and per-run yields are identical to the
post-fix run.

**Release 0.4.5 (tag `v0.4.5`, rerun from scratch 2026-09-30, jobs 147833 and
147895).** Every arm's metrics and every per-run yield are identical to the
post-fix run. `detect` and `extract` did not change; `fetch` did (rounds read from
the archive's sample attributes), and it now downloads 4 runs of PRJEB47428 and
1 of PRJDB7022 that it previously left without a round, which changes no call.
In the first job two downloads failed for reasons on the archive's side (an md5
mismatch for SRR18110626, an empty FASTQ link for ERR2764569, both available
again within hours); those two deposits were rerun on the same tag and every
run came through. The sensitivity analysis with the deposit excluded after
inference (below) is measured with this version: PRJNA1244358 is no longer
refused, its 3′ constant is called exactly and its 5′ call is the inner part of
the published constant (PARTIAL_3P), with the random region at its 16 nt. With
it counted in, the recovery arm is 5 exact and 3 partial of 8, and no call is
wrong.

*v0.4.1 (pre-specified):* on 7 recovery deposits selexprep reproduced both
paper-reported primer strings exactly on 4 and partially on 3; in every one of
the 6 where a single-read extraction is possible it recovered the random region
at exactly the paper-reported length, 0 out of tolerance; it made zero
false-positive primer calls on 14 negative controls (8 pre-trimmed, 6 adapter).
Applied to every round, the called constants kept at least 78% of the reads of
every run of five deposits and of 11 of the 12 runs of PRJEB62495, but under 1%
of the reads in 7 of the 10 runs of PRJNA809588.

*After the fix (in-sample):* both strings exact on 5 of 7 (per flank 6 of 7 at
5′ and 6 of 7 at 3′), the random-region length unchanged (exact in all 6
measurable deposits), the controls unchanged (0 false calls on 14), and every
run of every recovery deposit keeps at least 63% of its reads — PRJNA809588 now
84–99% in all ten runs.

Two measurements are reported because they answer different questions. The
**random-region boundary** is what the tool must get right to trim correctly;
the **exact primer string** additionally requires that the constant called from
the reads coincides with the constant as written in the paper. The boundary was
right in every measurable deposit in both versions. The two partials left after
the fix are differences between the reads and the paper, not boundary errors:
in PRJEB62495 one base of the published 3′ constant differs between runs, so the
3′ call now stops before it; in PRJNA809588 the 5′ call still extends past the
published constant at the outer edge (on reconstructed reads, by the one G that
ends every run's inline tag).

The source-verified deposits, by arm (the table below is the per-row scorecard
after the fix; `snakemake -s Snakefile` regenerates it as `metrics.json`):

| Accession | Chemistry | Target | Arm | `status` | 5′ / 3′ | *N* obs / truth | Note |
|---|---|---|---|---|---|---|---|
| PRJDB9110 | RNA (T7) | TG2 (protein) | recovery | HIGH | EXACT / EXACT | 30 / 30 | counted |
| PRJDB9111 | RNA (T7) | αVβ3 (protein) | recovery | HIGH | EXACT / EXACT | 40 / 40 | counted |
| PRJDB19098 | RNA (T7) | FGF-9 (protein) | recovery | HIGH | EXACT / EXACT | 35 / 35 | counted; truth from patent WO2020204151 |
| PRJNA615076 | DNA | *E. faecalis* (cell) | recovery | HIGH | EXACT / EXACT | 40 / 40 | counted |
| PRJEB62495 | DNA | *Anaplasma* (cell) | recovery | HIGH | EXACT / PARTIAL | 40 / 40 | counted (partial); v0.4.1: EXACT / MISMATCH. The 3′ base that differs between runs is left out |
| PRJNA1395820 | DNA | Co²⁺ (small molecule) | recovery | MEDIUM | EXACT / EXACT | — | counted; v0.4.1: PARTIAL / EXACT. The heterogeneous leading base is left out; paired split primers → merge recommended |
| PRJNA809588 | 2′-F-Py RNA | islet (cell) | recovery | HIGH | PARTIAL / EXACT | 40 / 40 | counted (partial); v0.4.1: MEDIUM, PARTIAL / PARTIAL. The run-specific tags are left out |
| PRJEB22637 | 2′-F-Py RNA | cell (Annexin A2) | specificity | UNABLE_TO_INFER | null / null | — | correct refusal (N-region-only reads) |
| PRJEB28411 | DNA | cell (ccRCC) | specificity | UNABLE_TO_INFER | null / null | — | correct refusal |
| PRJNA990511 | DNA | protein (ASFV p30) | specificity | UNABLE_TO_INFER | null / null | — | correct refusal (single-round) |
| PRJEB47428 | RNA | RNA-binding proteins (HTR-SELEX) | specificity | UNABLE_TO_INFER | null / null | — | correct refusal; reads are N40 exactly (92 runs) |
| PRJEB49150 | DNA | BEN-domain TFs | specificity | UNABLE_TO_INFER | null / null | — | correct refusal; submitter states reads carry only the randomised region |
| PRJEB14550 | DNA | HOXB13 / FLI1 | specificity | UNABLE_TO_INFER | null / null | — | correct refusal; reads are N40 exactly (8 runs) |
| PRJNA360902 | DNA | *Ciona* TF DBDs | specificity | UNABLE_TO_INFER | null / null | — | correct refusal; reads are N20 exactly (14 runs) |
| PRJEB70964 | 2′-F-Py RNA | protein (α-syn) | specificity | UNABLE_TO_INFER | null / null | — | correct refusal; reads are N35 exactly (reclassified from adapter-control, see below) |
| PRJNA678231 | n/a (ncRNA-Seq) | n/a (*B. mori*) | adapter-control | UNABLE_TO_INFER | null / null | — | correct refusal; 51 nt reads, TruSeq R1 in 97% at variable positions |
| PRJDB7022 | n/a (miRNA-Seq) | n/a (*D. melanogaster*) | adapter-control | UNABLE_TO_INFER | null / null | — | correct refusal **after the zero-length-insert guard**; 83% of reads are one identical 51-mer |
| PRJNA746278 | n/a (ncRNA-Seq) | n/a (*H. sapiens*) | adapter-control | UNABLE_TO_INFER | null / null | — | correct refusal; 35–76 nt reads ending in poly-A, SMARTer smRNA kit |
| PRJDB2183 | n/a (miRNA-Seq) | n/a (*A. thaliana*) | adapter-control | UNABLE_TO_INFER | null / null | — | correct refusal; 69 nt reads ending in poly-A |
| PRJEB50674 | n/a (miRNA-Seq) | n/a (*T. vaginalis*) | adapter-control | UNABLE_TO_INFER | null / null | — | correct refusal; 50 nt reads, no technical sequence at a fixed position |
| PRJNA591605 | n/a (ncRNA-Seq) | n/a (*M. musculus*) | adapter-control | UNABLE_TO_INFER | null / null | — | correct refusal; 75 nt reads, no technical sequence at a fixed position |

**How to read it.** The *recovery* arm asks whether selexprep recovers the
paper primer from raw reads; the *specificity* and *adapter-control* arms are
negative controls where the correct behavior is **no call** (the reads carry
no library constant, or are not SELEX at all). A "partial" or "mismatch" is a
difference between what the reads carry and the constant as published, not a
wrong boundary; whether the call then extracts every round is measured
separately, under [Extraction yield](#extraction-yield-per-run).
The *N* column is the mode of the extracted random-region length against the
paper-reported length; it is blank where a single-read extraction is not
attempted (PRJNA1395820 is paired-end with split primers, so `detect` asks for
read merging rather than forcing an R1-only call).

PRJNA883192 was **withdrawn from the scored set** (`verified=false`) because its
3′ constant could only be resolved from the deposited reads themselves: scoring
an inferred call against a read-derived truth is circular, and keeping it would
have inflated the recovery denominator with a case the benchmark cannot honestly
adjudicate. It stays in `ground_truth.tsv` with the full reasoning.

### Extraction yield per run

`detect` infers the constants from the earliest round, and the scorecard above
scores that call. The workflow then runs `extract` on every round of the
recovery deposits, the way a user would, and `results/extraction_yield.tsv`
records reads in and reads kept for each input FASTQ, with the modal extracted
length of each round (current file: after the fix).

| Accession | FASTQs | Reads kept per FASTQ, v0.4.1 | After the fix | Modal length / truth |
|---|---|---|---|---|
| PRJDB19098 | 6 | 91.9–93.9% | 91.9–93.9% | 35 / 35 |
| PRJDB9110 | 9 | 78.4–92.0% | 78.4–92.0% | 30 / 30 |
| PRJDB9111 | 4 | 80.9–84.4% | 80.9–84.4% | 40 / 40 |
| PRJNA615076 | 11 | 89.7–95.5% | 89.7–95.5% | 40 / 40 |
| PRJNA1395820 | 1 | 99.5% | 99.6% | partial R1 (paired split primers, no length) |
| PRJEB62495 | 12 | 79.1–93.8%; ERR11470033 36.7% | 83.1–98.7%; ERR11470033 63.1% | 40 / 40 |
| PRJNA809588 | 10 | 89.9–98.5% in 3 runs; 0.03–1.0% in 7 | 83.8–99.2% | 40 / 40 |

**In v0.4.1 PRJNA809588 lost seven of its ten runs.** Every run carries its own
inline tags outside both library constants, differing in sequence and length
between runs (`read_state_evidence.tsv`). `detect` read the earliest run,
SRR18110617, and its call included that run's tags; `extract` applied the call
unchanged to every round, and cutadapt's linked-adapter mode discards any read
in which either called flank is not found. The three runs whose tags match or
nearly match the earliest run's kept 90–99% of their reads; the other seven kept
1% or less, and 31% of the deposit's 332M reads survived overall, while
`extract` exited normally. This was a limitation of `selexprep`, not of the
deposit: the library constant is the part of the flank that every run shares.

`detect` now keeps only that part (`_shared_library_constant`; see the
CHANGELOG), and `extract` warns when an input keeps less than half the yield of
the best one. After the fix all ten runs keep 84–99% of their reads, 91% of the
332M overall. These numbers come from the deposit that motivated the fix and
are in-sample.

In PRJEB62495, ERR11470033 (round 12) kept 36.7% of its reads in v0.4.1 and
63.1% after the fix, which no longer calls the part of the 3′ end that varies
between runs; it remains the run with the lowest positional conservation at its
3′ end (`read_state_full.tsv`).

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

**Adapter control.** Six non-SELEX small-RNA libraries (`miRNA-Seq` /
`ncRNA-Seq`, six organisms, at least three library-prep kits), chosen because
the archive-reported read length exceeds a small-RNA insert, so the reads
should run into sequencing adapter. Any primer call on them is a fabrication.
The all-read check (`read_state_evidence.tsv`) shows what they actually carry:
TruSeq R1 in 97% of PRJNA678231 reads but at variable positions, a poly-A tail
at the read end in PRJDB2183 and PRJNA746278, one dominant sequence in
PRJDB7022, and no technical sequence at a fixed position in PRJEB50674 and
PRJNA591605. The arm therefore tests that `detect` makes no call on non-SELEX
data; it does not test an adapter sitting exactly where a library constant
would sit. Amplicon libraries were deliberately **not** used: a 16S
amplicon has real constant primers flanking a variable region, so calling them
would be correct behaviour, and scoring it as a false positive would punish the
right answer.

The randomised-region length for three of the four new specificity rows comes
from a publication that cites the accession; for PRJNA360902 no citing
publication was found, so its length is submitter-stated in the SRA record —
external to the reads, but archive-sourced rather than paper-sourced, and the
row says so.

**PRJEB70964 moved from the adapter arm to the specificity arm.** It was
filed as an adapter-collision control — a SELEX library whose 5′ constant ends
in the reverse complement of TruSeq R1 — on a read-state spot-check that
recorded 81-nt reads carrying both constants. A check of every read
(`verify_read_state.py`: 155M reads in the 17 downloadable runs, and ENA's
per-run base/read counts for all 27) shows the reads are the N35 region alone,
with neither constant present. Its refusal was always correct, but for the
reason the specificity arm tests, so it now sits there; no outcome changes. The
adapter-collision case itself was never exercised by this deposit. The same
check replaced every spot-check record in `read_state_evidence.tsv` with
all-read measurements.

**The adapter arm earned its place immediately.** On its first run PRJDB7022
produced a false positive: `detect` returned the same 51 nt string as both
constants and an inferred random region of **zero** nucleotides, while
reporting `full_insert_recovered=True` and `required_action=NONE` — a state in
which `run` would have proceeded to `count` and written a table of zero-length
sequences without raising anything. The cause is that 4.7M of the deposit's
5.7M reads are one identical 51-mer, so positional conservation is ~100% at
every position and there is no constant/random boundary for the inward walk to
stop at. `detect` now refuses when the inferred random-region length is 0
(`test_zero_length_insert_guard.py`), and the deposit scores a correct refusal.
The guard keys on exactly 0 and not on falsiness, because paired split-primer
deposits report `None` there and must keep asking for read merging.

**Excluded deposits** live in `excluded_datasets.tsv` with an evidence-based,
pre-inference reason (e.g. nonstandard read architecture, multiplexing without
documented barcodes). They are removed *before* analysis and documented — not
silently dropped — so the recovery denominator is honest. The pre-detect
screening decision for every candidate is recorded in `screening_log.tsv`.

**One candidate was excluded after inference, and is reported here for that
reason.** PRJNA1244358 (HT-SELEX of ASCL2 and HES2 variants, Nucleic Acids Res
2025, doi:10.1093/nar/gkaf831) has both library constants in its publication
(5′ `GTTCAGAGTTCTACAGTCCGACCTAA`, 3′ `TTAGGACTCGGACCTGGACTAGG`, N16), and an
anchored presence test on one run (SRR32924589) found both. On 2026-08-10 it
was added to the recovery arm on a development branch, and `detect` (v0.4.0)
was run on its nine runs: the 3′ constant was recovered exactly, no 5′
constant was called, and the deposit was refused. A decision rule was written
only after that result: the deposit was excluded as outside the model
`selexprep` assumes, one selection trajectory analysed across rounds. Its
archive records show a shared 16-mer starting library (`16mer-SELEX-R0`) and
eight parallel round-1 selections against different protein variants
(`ASCL2_WT-SELEX-R1`, `HES2_K5R-SELEX-R1`, …). The modal 5′ offsets measured
afterwards (7, 0 and 6 nt in SRR32924589, SRR32924596 and SRR32924594) are
consistent with sample-specific sequence before the constant, but do not by
themselves show the design. Because the exclusion was decided after
`detect`'s result was known, it is not recorded in `excluded_datasets.tsv`; it
is disclosed here, and the deposit is evaluated with the released version as a
sensitivity analysis alongside the recovery arm (Table S1).

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

## Round-assignment coverage

Descriptive and in-sample, as fixed in `round_validation/PROTOCOL.md`
(Amendment 1): how often `fetch` assigns a round, never whether the assignment
is right. Run on the tag `v0.4.5` (job 147832) over every run of the 127 INSDC
deposits in the catalogue, from archive metadata only; results in
`round_validation/results/`.

- 25,435 runs: 17,693 receive a round (17,169 HIGH, 524 MEDIUM), 7,742 are left
  unassigned. 55 deposits have every run assigned, 55 none, 17 some.
- Curated round maps, on the 4 catalogue deposits that have one: of 31 runs, 9
  receive the same round, none a different one, and 22 none. Six other maps
  belong to Tier-1 deposits outside the catalogue and are listed apart.
- The run's `summary.json` reads `dirty: true`: the only tracked files that
  differed from the tag were `results/metrics.json` and
  `results/extraction_yield.tsv`, moved aside for the Tier-1 rerun; no code
  differed. The comparison with the curated maps was redone afterwards from the
  same `runs.tsv` (`--agreement-only`), to label the six maps outside the
  catalogue correctly; the commit that redid it is recorded in `summary.json`.

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
