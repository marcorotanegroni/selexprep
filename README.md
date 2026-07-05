# selexprep

[![Tests](https://github.com/marcorotanegroni/selexprep/actions/workflows/tests.yml/badge.svg)](https://github.com/marcorotanegroni/selexprep/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/marcorotanegroni/selexprep/blob/main/LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)

**Accession-first preprocessing for public HT-SELEX, with primer auto-inference
and safe failure modes.**

Give `selexprep` a public accession (or your own FASTQs) and it **infers the
primer/constant regions from the reads**, strips them, and emits clean per-round
random-region FASTAs + count tables — with an explicit confidence level and a
reproducible manifest. When it can't infer the library with confidence, it
**says so and refuses to guess** rather than emitting fabricated primers.

It's the step *before* the existing SELEX toolchain. Tools like FASTAptameR,
AptaSUITE, RaptGen, and AptaTrans assume you already have trimmed reads and
already know your primers. `selexprep` produces exactly those inputs, starting
from a raw public deposit.

```
  accession ──▶ fetch ──▶ detect ────▶ extract ───▶ count ────▶ qc
 (or local             (infer        (trim per     (unique     (flags +
  FASTQs)               primers)      round)        seqs)       4 plots)
```

`run` chains the whole thing for you; the individual verbs (`fetch`, `detect`,
`extract`, `count`, `qc`) are there for when you need to drive one stage by hand.

> **Status — beta.** The core accession/local-FASTQ preprocessing pipeline is
> feature-complete and CI-tested (Python 3.10 / 3.11 / 3.12). The CLI commands and
> the primary output schemas (`library_report.json`, `selexprep_manifest.json`,
> `counts.parquet`) are treated as **stable within the 0.x series**; experimental
> features are listed separately. See [`STABILITY.md`](STABILITY.md) for the full
> stable-vs-experimental surface.

---

## Install

```bash
pip install selexprep        # pulls in cutadapt, the only external dependency
selexprep --version
```

Python 3.10+. The only external tool is `cutadapt`, installed automatically.
Working from a clone instead? `uv pip install -e .` from the repo root.

## Quick start

**One accession, no primers supplied.** The whole pipeline is **one command**.
Write a TSV listing one or more accessions, then point `run` at it. It fetches
each dataset, infers the primers from the read content, then extracts, counts,
and QCs — with nothing about the library supplied by you:

```bash
printf 'accession\nPRJNA615076\n' > accessions.tsv
selexprep run accessions.tsv --outdir out
```

That's the happy path. Everything for each accession lands under
`out/<accession>/`:

| Output | What it is |
|---|---|
| `library_report.json` | inferred primers + `status` (`HIGH` / `MEDIUM` / `LOW` / `UNABLE_TO_INFER`) |
| `round_NN/extracted.fasta.gz` | primer-stripped random region, one folder per round |
| `round_NN/counts.parquet` | unique sequences per round (reads, rank, RPM) |
| `qc/flags.yaml` + 4 PNGs | depth-aware QC flags + diagnostic plots |
| `selexprep_manifest.json` | reproducibility manifest (output hashes, dependency versions, CLI argv) |

A corpus-level `out/run_summary.tsv` spans every accession. If a batch is
interrupted, re-run the **same** command with `--resume` and finished datasets
are skipped.

> **Watch it run, offline, in under a minute.**
> [`examples/01_offline_toy_pipeline.ipynb`](https://github.com/marcorotanegroni/selexprep/blob/main/examples/01_offline_toy_pipeline.ipynb)
> executes the full `detect → extract → count → qc` path on a tiny synthetic
> library — no network, no accession, no HPC. It's the fastest way to *see* what
> each stage produces before pointing the tool at real data.

## Common situations

The situations people actually hit. Find yours.

### "I don't have an accession yet"

Browse the bundled catalog of public SELEX deposits, then feed the hits to `run`:

```bash
selexprep catalog list --target IL-10RA --insdc-only
selexprep catalog show PRJEB12345
```

INSDC accessions are fetchable; figshare/zenodo rows are discovery-only pointers
(flagged in the `fetchable` column) that `run` can't retrieve yet —
`--insdc-only` filters the list to just the retrievable ones.

### "I only have my own local FASTQs"

`run` is for accessions; for files already on disk, drive the stages yourself.
This is also how you inspect an intermediate or re-extract with corrected
primers:

```bash
# detect needs a round map — cross-round persistence is the key primer signal
printf 'file\tround_number\nround_00.fastq.gz\t0\nround_01.fastq.gz\t1\nround_02.fastq.gz\t2\n' > rounds.tsv

selexprep detect round_00.fastq.gz round_01.fastq.gz round_02.fastq.gz \
    --round-map rounds.tsv --outdir ./out

selexprep extract round_00.fastq.gz round_01.fastq.gz round_02.fastq.gz \
    --library-report ./out/library_report.json --round-map rounds.tsv --outdir ./out

for r in 0 1 2; do
    selexprep count ./out/round_$(printf '%02d' $r)/extracted.fasta.gz --round R$r --outdir ./out
done

selexprep qc ./out/selexprep_manifest.json
```

The runnable, annotated version of exactly this path is
[`examples/01_offline_toy_pipeline.ipynb`](https://github.com/marcorotanegroni/selexprep/blob/main/examples/01_offline_toy_pipeline.ipynb);
every flag and the `--round-map` contract are in the [CLI reference](https://github.com/marcorotanegroni/selexprep/blob/main/docs/cli.md).

### "My reads are paired-end (R1 + R2)"

Common for modern Illumina SELEX — and **handled for you from an accession**:
`run` detects the paired layout at fetch time and threads R2 through with no
extra flags. By hand on local files, list the **R1** files in the round map and
the matching **R2** files via `--paired-r2`, in the same order, one R2 per R1:

```bash
printf 'file\tround_number\nr0_R1.fastq.gz\t0\nr1_R1.fastq.gz\t1\n' > rounds.tsv

selexprep detect r0_R1.fastq.gz r1_R1.fastq.gz \
    --paired-r2 r0_R2.fastq.gz r1_R2.fastq.gz \
    --round-map rounds.tsv --outdir ./out

selexprep extract r0_R1.fastq.gz r1_R1.fastq.gz \
    --paired-r2 r0_R2.fastq.gz r1_R2.fastq.gz \
    --library-report ./out/library_report.json --round-map rounds.tsv --outdir ./out
```

If the 5' primer sits on R1 and the 3' on R2, `detect` reports
`extraction_mode: PAIRED_END_SPLIT_PRIMERS` and `extract` writes the two trimmed
sides separately (`partial_5p_extracted_R1…` + `partial_3p_extracted_R2…`). It
stops there on purpose: stitching the sides into one full-length insert is **read
merging, not yet implemented**, so `count`/`qc` are skipped for split-primer rounds
rather than counting half-inserts.

### "I only have the final enriched pool (a single round)"

The common case — HT-SELEX is costly, so many experiments sequence only the
final pool. `selexprep` still works; pass a one-row round map:

```bash
printf 'file\tround_number\nfinal_pool.fastq.gz\t0\n' > rounds.tsv
selexprep detect final_pool.fastq.gz --round-map rounds.tsv --outdir ./out
```

The one difference: **confidence is capped at `MEDIUM`**, because cross-round
persistence (the strongest SELEX-specific primer signal) needs ≥2 rounds.
`detect` warns and falls back to within-round signals (primer match rate, flank
position, low-entropy region, adapter blacklist). Verify the inferred primers
before trusting extraction — or, if you designed the library and know the
primers, pass them with `--override-primer-5p/3p` and the confidence cap no
longer applies (you're not inferring anything).

### "detect returned UNABLE_TO_INFER"

By design, `extract` then **refuses rather than guessing**. Supply the primers
explicitly:

```bash
selexprep extract round_*.fastq.gz \
    --library-report ./out/library_report.json --round-map rounds.tsv \
    --override-primer-5p GGTAATACGACTCACTATAGGG \
    --override-primer-3p CCATGCATGCATGCATGCAT \
    --rebuild --outdir ./out
# --rebuild emits extract_diff.tsv comparing baseline vs override per round.
```

You can also edit `library_report.json` by hand and re-run without overrides.

> **Multiplexed deposit** (several samples or targets pooled in one FASTQ)? v0.1
> doesn't auto-split them — pass a `--sample-sheet` to `extract` in place of
> `--round-map`. Auto-demultiplexing is not yet implemented.

## CLI reference

The core pipeline is `inspect → fetch → detect → extract → count → qc`; `run`
batches it and `catalog` browses the discovery catalog. Every command supports
`--help`; full reference in [`docs/cli.md`](https://github.com/marcorotanegroni/selexprep/blob/main/docs/cli.md).

| Command | What it does |
|---|---|
| `selexprep inspect <ACC>` | ENA metadata preview — round count, `library_strategy` (SRA verbatim), file sizes + MD5s. No download. |
| `selexprep fetch <ACC> --outdir OUT` | Download FASTQ + auto-populate `rounds.tsv`. Unassigned runs go to `round_unknown/` only with `--allow-manual-review`. |
| `selexprep detect <fastq...> --round-map rounds.tsv --outdir OUT` | Auto-infer primers + library structure → `library_report.json`. |
| `selexprep extract <fastq...> --library-report LR.json --round-map rounds.tsv --outdir OUT` | Cutadapt trim + strand reorient + per-round FASTA + manifest. Accepts `--paired-r2`, `--sample-sheet`, `--override-primer-{5p,3p}`, `--rebuild`. |
| `selexprep count <extracted.fasta.gz> --round R0 --outdir OUT` | FASTA → `counts.parquet` (sequence, reads, rank, RPM). |
| `selexprep qc <manifest> --outdir OUT` | Depth-aware suspicion flags (YAML) + 4 PNG plots. |
| `selexprep run <accessions.tsv> --outdir OUT [--resume]` | Batch driver across many accessions; emits `run_summary.tsv` + per-accession outputs. |
| `selexprep catalog list \| show \| version \| refresh` | Browse / refresh the bundled discovery catalog (240 entries; `refresh` hits live ENA). |

## Full output layout

After a complete single-dataset run:

```
out/
├── library_report.json            # primer inference + extraction_mode
├── selexprep_manifest.json        # reproducibility anchor
├── trim_reports.json              # cutadapt argv + n_in/n_out per round
├── strand_report.tsv              # only if orientation ∈ {MIXED, REVERSE}
├── extract_diff.tsv               # only with --rebuild + --override-primer-*
├── round_00/
│   ├── extracted.fasta.gz         # full-insert recovery (BOTH_PRIMERS_SINGLE_READ)
│   │   # or partial_5p_extracted.fasta.gz (FIVE_PRIME_ONLY)
│   │   # or partial_3p_extracted.fasta.gz (THREE_PRIME_ONLY)
│   │   # or partial_5p_extracted_R1.fasta.gz + partial_3p_extracted_R2.fasta.gz (PAIRED_END_SPLIT_PRIMERS)
│   └── counts.parquet             # unique sequences
├── round_NN/...
└── qc/
    ├── flags.yaml                 # depth-aware suspicion flags
    ├── read_retention.png         # input vs extracted per round
    ├── primer_match_per_round.png
    ├── n_length_distribution.png  # per-round N-region histograms (faceted)
    └── per_round_panel.png        # unique seqs + Shannon entropy + top-100 coverage
```

Output filenames deliberately distinguish full inserts (`extracted.fasta.gz`)
from one-sided partials (`partial_5p_…` / `partial_3p_…`) so downstream ML
pipelines never accidentally mix them.

---

## Why it exists (the gap)

No maintained, pip-installable Python tool takes a public HT-SELEX accession or
local FASTQ and produces trimmed per-round FASTAs + count tables + a provenance
manifest **without requiring the user to supply primers**.

| Tool | Accession fetch? | Primer auto-inference? | Safe failure mode? |
|---|---|---|---|
| AptaSUITE (CLI + GUI) | ✗ | ✗ — user supplies | ✗ |
| FASTAptameR | ✗ — consumes already-trimmed | ✗ | ✗ |
| EasyDIVER+ (2025) | ✗ | ✗ — paired-end with known primers | ✗ |
| ht-selex-demo | ✗ | ✗ — Illumina adapters only (a trap `selexprep` guards against via `known_adapter_hits`) | ✗ |
| nf-core | ✗ (no SELEX pipelines exist) | — | — |
| **`selexprep`** | ✓ — INSDC (figshare/Zenodo deferred) | ✓ — cross-round persistence + adapter blacklist | ✓ — explicit `LibraryReport.status` ∈ {HIGH, MEDIUM, LOW, UNABLE_TO_INFER} |

The core inference idea: **a true primer appears at a similar rate across all
rounds; late-round-enriched aptamer motifs do not.** `detect` exploits this
cross-round persistence to separate constants from binders, and reports its
confidence in a typed `LibraryReport` (pydantic, strict-mypy) with explicit
`extraction_mode`, `read_source`, `required_action`, `full_insert_recovered`,
and `status` — no silent miscalls.

**Benchmark headline:** *on a curated, paper-grounded set, `selexprep` recovers
SELEX primers directly from raw reads — with none supplied — and fails safe (no
fabricated primers) where inference is ambiguous.* See
[Benchmarks](#benchmarks).

## What v0.1 ships

`catalog` (240 public SELEX entries: 127 INSDC + 113 figshare/Zenodo
passthrough; documented exclusions in `bioprojects_excluded.csv`) · `inspect`
(ENA metadata preview) · `fetch` (accession download with relaxed
partial-parseability contract) · `detect` (primer inference) · `extract`
(cutadapt-driven trimming + paired-end + strand handling) · `count` (per-round
unique sequences) · `qc` (suspicion flags + 4 PNG plots) · `run` (batch driver
with `--resume`).

Catalog discovery completeness is **auditable** via
`benchmarks/catalog_completeness_audit.py` (it re-queries live ENA): as of the
v0.2.1 snapshot (2026-07-03) every ENA-typed-SELEX study is accounted for —
either included, or documented in `bioprojects_excluded.csv` (mis-labeled and
gSELEX/genomic-fragment variants). Because ENA grows over time, the audit is the
source of truth rather than a frozen percentage.

**Not yet implemented** (deferred to a future release): read merging for
paired-end full-insert recovery · multiplex auto-detection (currently needs a
user-supplied sample sheet) · figshare/Zenodo fetch backends ·
`SELEXPREP_CATALOG_PATH` env var for user-supplied catalogs · AnnData export ·
BibTeX auto-citation · library-type classification.

## What v0.2 adds

**Curated metadata layer** (`selexprep.catalog.load_metadata` /
`load_metadata_records`) — each of the 240 catalog deposits now ships with curated
experimental metadata: `study_type`, `target`, `target_class`, `chemistry`,
`n_random`, `n_rounds`, `selection_format`, `counter_selection` (240 × 8 = 1920
cells, 1491 populated — up from ~4 experimental-field cells in v0.1). Built by **two
independent LLM extractions** (Claude + Codex/GPT), reconciled, with **every value
source-cited** (evidence quote + DOI/URL + location); where the two extractions
genuinely disagreed, both are kept. Bundled as `curated_metadata.json` (canonical,
provenance-rich) + `curated_metadata.csv` (flat view) under `catalog/data/`, accessed
via the API above. Spot-checked against the hand-curated benchmark
ground truth (overlapping fields on 11 verified deposits) — no residual disagreements
after adjudication. The extraction contract, both raw arms, and the reconciliation
live under `benchmarks/dual_extraction/` in the source repository.

## What v0.1 does *not* do

By design — these are handled by mature tools that consume `selexprep`'s outputs:

| Step | Use |
|---|---|
| Clustering | FASTAptameR |
| Motif discovery | MEME · RaptGen-UI |
| Binding-affinity prediction | RaptGen · DeepSELEX · AptaTrans |
| 3D structure | ViennaRNA · RNAfold |
| Aptamer design | MAWS · RNAtranslator |
| Read merging (paired-end full-insert recovery) | bbmerge · vsearch · pear (planned hook) |

## Determinism + reproducibility

All `.fasta.gz`, `.tsv`, and `.json` outputs are **byte-deterministic** given
the same `--sampling-seed`: gzip headers use `mtime=0`, JSON keys are sorted
(int keys in numeric order), TSV rows are sorted. Two runs with identical inputs
+ seed produce identical SHA256s; the manifest's `output_sha256` is the audit
anchor.

Parquet hashes are **advisory** (not guaranteed across pyarrow versions — the
manifest pins `pyarrow_version` instead). PNG plots are **informational**
(matplotlib output isn't byte-deterministic) and don't contribute to
`output_sha256`.

## Benchmarks

Two-tier benchmark under [`benchmarks/`](https://github.com/marcorotanegroni/selexprep/blob/main/benchmarks/README.md):

- **Tier 1 — primer-recovery scorecard** (`Snakefile` + `ground_truth.tsv`):
  curated primer-recovery validation against paper-reported primers (N=11
  source-verified accessions, modality-diverse), reported as a per-deposit table.
- **Tier 2 — corpus audit** (`audit.smk`): corpus utility audit over a random sample
  of audit-eligible INSDC catalog rows. Shipped artifacts under
  `benchmarks/audit_results/`; the reproducibility envelope (catalog version +
  sample sha256 + seed) is committed.
- **Catalog completeness audit** (`catalog_completeness_audit.py`): re-runnable
  script diffing `bioprojects.csv` against ENA's `library_strategy="SELEX"` set,
  so coverage can be re-verified against a current, dated snapshot rather than a
  frozen figure (ENA grows over time).

## Calibration status

Tests assert on **behavior**, never on threshold values (e.g.
`assert report.status == "HIGH"`, never `assert HIGH_CUTOFF == 0.85`), so tuning
the numbers is safe under the suite. `CALIBRATION-REVIEWED` markers in
`library/detect.py` document vetted v0.1 thresholds with rationale;
`CALIBRATION-TODO` markers (in `qc/flags.py`, `library/adapters.py`,
`extract/strand.py`) flag what's still pending.

```bash
grep -rn "CALIBRATION-TODO" src/      # what's left to tune
grep -rn "CALIBRATION-REVIEWED" src/  # what's already vetted with rationale
```

## Development

```bash
uv run pytest                       # full suite + 1 reserved xfail
uv run ruff check src/ tests/
uv run ruff format --check src/ tests/
uv run mypy src/                     # strict on library.report
```

CI matrix: Python 3.10 / 3.11 / 3.12.

## License

MIT. See [`LICENSE`](https://github.com/marcorotanegroni/selexprep/blob/main/LICENSE). Full development log in [`CHANGELOG.md`](https://github.com/marcorotanegroni/selexprep/blob/main/CHANGELOG.md).
