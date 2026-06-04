# selexprep

[![Tests](https://github.com/marcorotanegroni/selexprep/actions/workflows/tests.yml/badge.svg)](https://github.com/marcorotanegroni/selexprep/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)

**Accession-first preprocessing for public HT-SELEX, with primer auto-inference and safe failure modes.**

> *`selexprep` fills the missing preprocessing layer for public datasets by starting from accessions and empirically inferring primer/constant regions, extracting random regions, and emitting confidence-aware, reproducible count tables and manifests.*
>

## Status

**v0.1 RC.** Full single-dataset + batch workflow is feature-complete
and tested (601 passing tests + 1 strict-xfail reserved for v0.2 read
merging). Tier 2 corpus audit shipped against a v0.1.6 catalog
snapshot; catalog discovery completeness measured at 100% of
ENA-typed-SELEX deposits in v0.1.7. The Phase 6 primer-recovery
benchmark (Figure A) is complete — paper-grounded recovery on a
curated multi-chemistry set, alongside safe-failure and specificity
arms. Not yet on PyPI.

**v0.1 ships:** `catalog` (250 public SELEX entries: 125 INSDC +
125 figshare/zenodo passthrough; 21 documented exclusions in
`bioprojects_excluded.csv`) · `inspect` (ENA metadata preview) ·
`fetch` (accession download with relaxed partial-parseability
contract) · `detect` (primer inference) · `extract` (cutadapt-driven
trimming + paired-end + strand handling) · `count` (per-round unique
sequences) · `qc` (suspicion flags + 4 PNG plots) · `run` (batch
driver across many accessions with `--resume`).

**v0.2 deferrals:** read merging for paired-end full-insert recovery ·
multiplex auto-detection (v0.1 needs a user-supplied sample sheet) ·
figshare/zenodo fetch backends · `SELEXPREP_CATALOG_PATH` env var
for user-supplied catalogs · AnnData export · BibTeX auto-citation ·
library-type classification.

## Why it exists (the gap)

No maintained, pip-installable Python tool takes a public HT-SELEX
accession or local FASTQ and produces trimmed per-round FASTAs +
count tables + a provenance manifest **without requiring the user to
supply primers**.

| Tool | Accession fetch? | Primer auto-inference? | Safe failure mode? |
|---|---|---|---|
| AptaSUITE (CLI + GUI) | ✗ | ✗ — user supplies | ✗ |
| FASTAptameR | ✗ — consumes already-trimmed | ✗ | ✗ |
| EasyDIVER+ (2025) | ✗ | ✗ — paired-end with known primers | ✗ |
| ht-selex-demo | ✗ | ✗ — Illumina adapters only (a trap `selexprep` guards against via `known_adapter_hits`) | ✗ |
| nf-core | ✗ (no SELEX pipelines exist) | — | — |
| **`selexprep`** | ✓ (v0.2; catalog ships v0.1) | ✓ — cross-round persistence + adapter blacklist | ✓ — explicit `LibraryReport.status` ∈ {HIGH, MEDIUM, LOW, UNABLE_TO_INFER} |

**Benchmark headline (Phase 6):** *`selexprep` approaches
known-primer pipelines on datasets where primer inference is
high-confidence, while explicitly failing safe on ambiguous ones.*

## What v0.1 does

> *`selexprep` converts public or local HT-SELEX reads into primer-stripped
> random-region FASTA/FASTQ files, per-round count tables, QC reports,
> and provenance manifests, with explicit primer-inference confidence
> and safe failure modes.*

Concretely:

1. **Discovers** public SELEX BioProjects from a bundled catalog (250
   entries: text-pattern queries + `library_strategy="SELEX"` positive
   query — measured at 100% coverage of ENA-typed-SELEX deposits, with
   21 documented exclusions for mis-labeled and gSELEX/genomic-fragment
   variants; `selexprep catalog list / show`).
2. **Infers** the SELEX library structure — primer pair, N-region
   length, orientation, paired-end layout — directly from read content,
   using **cross-round persistence** (a true primer appears at a
   similar rate across all rounds; late-round-enriched aptamer motifs
   do not).
3. **Reports** inference confidence in a typed `LibraryReport`
   (pydantic, strict-mypy) with explicit `extraction_mode`,
   `read_source`, `required_action`, `full_insert_recovered`, and
   `status`. No silent miscalls — ambiguous datasets surface as
   `UNABLE_TO_INFER`, and `extract` refuses without an explicit
   override.
4. **Extracts** the random region per round via cutadapt (subprocess;
   CLI is the stable contract). Output filenames distinguish full
   inserts (`extracted.fasta.gz`) from one-sided partials
   (`partial_5p_extracted.fasta.gz` / `partial_3p_extracted.fasta.gz`)
   to protect downstream ML pipelines from accidentally mixing them.
5. **Counts** unique sequences per round (raw reads + RPM + rank,
   parquet output) for downstream consumption by FASTAptameR /
   RaptGen-UI / AptaTrans.
6. **QC**: depth-aware suspicion flags (rarefied diversity, not raw
   counts; primer match per round; N-length variation; strand mix;
   adapter contamination; …) plus four per-dataset PNG plots.
7. **Manifests** every run in a versioned pydantic
   `SelexprepManifestV1` with SHA256s of FASTA/TSV/JSON outputs,
   dependency-version pins, CLI argv capture, and the full nested
   `LibraryReport` — making reruns byte-identical and audits trivial.

## Quick start

```bash
# Dev install (PyPI release pending)
uv pip install -e .

# 0. Discover available public SELEX datasets
selexprep catalog list --target IL-10RA --insdc-only
selexprep catalog show PRJEB12345

# 1. Preview an accession's metadata without downloading
selexprep inspect SRR12647619

# 2. Detect: auto-infer primers + library structure from local FASTQs
#    (requires a round map — cross-round persistence is the key signal)
cat > rounds.tsv <<EOF
file	round_number
round_00.fastq.gz	0
round_01.fastq.gz	1
round_02.fastq.gz	2
EOF
selexprep detect round_00.fastq.gz round_01.fastq.gz round_02.fastq.gz \
    --round-map rounds.tsv --outdir ./out

# 3. Extract: trim primers, emit per-round FASTAs + manifest
selexprep extract round_00.fastq.gz round_01.fastq.gz round_02.fastq.gz \
    --library-report ./out/library_report.json \
    --round-map rounds.tsv --outdir ./out

# 4. Count: per-round unique sequences -> parquet
for r in 0 1 2; do
    selexprep count ./out/round_$(printf '%02d' $r)/extracted.fasta.gz \
        --round R$r --outdir ./out
done

# 5. QC: depth-aware suspicion flags + 4 PNG plots
selexprep qc ./out/selexprep_manifest.json
```

When primer inference is ambiguous, `extract` refuses with a pointer
to override:

```bash
# Either edit library_report.json by hand…
# …or pass overrides at the CLI:
selexprep extract round_*.fastq.gz \
    --library-report ./out/library_report.json \
    --round-map rounds.tsv \
    --override-primer-5p GGTAATACGACTCACTATAGGG \
    --override-primer-3p CCATGCATGCATGCATGCAT \
    --rebuild --outdir ./out
# Emits extract_diff.tsv comparing baseline vs override per round.
```

### Only have the final pool / a single round?

That's the common case — HT-SELEX is costly, so many experiments
sequence only the final enriched pool. selexprep still works: pass a
one-row round map and run the same `detect` → `extract` → `count` →
`qc` flow.

```bash
printf 'file\tround_number\nfinal_pool.fastq.gz\t0\n' > rounds.tsv
selexprep detect final_pool.fastq.gz --round-map rounds.tsv --outdir ./out
```

The one difference: **confidence is capped at `MEDIUM`**, because
cross-round persistence (the strongest SELEX-specific primer signal)
needs ≥2 rounds. `detect` logs a warning saying so, and inference falls
back to within-round signals only (primer match rate, flank position,
low-entropy region, adapter blacklist). Verify the inferred primers
before trusting extraction — or, if you designed the library and
already know the primers, pass `--override-primer-5p/3p`: extraction
then uses your sequences directly and the inference confidence cap
no longer applies (you're not inferring anything).

## CLI surface

| Command | Status | What it does |
|---|---|---|
| `selexprep catalog list \| show \| version \| refresh` | ✅ v0.1 | Browse / refresh the bundled discovery catalog (250 entries; refresh hits live ENA). |
| `selexprep inspect <ACC>` | ✅ v0.1 | ENA filereport REST preview — round count, `library_strategy` (SRA verbatim, not classified), file sizes + MD5s. No download. |
| `selexprep fetch <ACC> --outdir OUT [--allow-manual-review]` | ✅ v0.1 | Download FASTQ + auto-populate round map. Partial-parseability is warn-and-skip (Phase 6b.5d); unassigned runs go to `round_unknown/` only with `--allow-manual-review`. |
| `selexprep detect <fastq...> --round-map rounds.tsv --outdir OUT` | ✅ v0.1 | Auto-infer primers + library structure → `library_report.json`. |
| `selexprep extract <fastq...> --library-report LR.json --round-map rounds.tsv --outdir OUT [--sample-sheet samples.tsv] [--paired-r2 ...] [--override-primer-{5p,3p} ...] [--rebuild]` | ✅ v0.1 | Cutadapt-driven trim + strand reorient + per-round FASTA + manifest. |
| `selexprep count <extracted.fasta.gz> --round R0 --outdir OUT` | ✅ v0.1 | FASTA → counts.parquet (sequence, reads, rank, RPM). |
| `selexprep qc <manifest> [--counts-dir DIR] [--outdir OUT]` | ✅ v0.1 | Depth-aware suspicion flags (YAML) + 4 PNG plots. |
| `selexprep run <accessions.tsv> --outdir OUT --resume` | ✅ v0.1 | Batch driver across many accessions; emits `run_summary.tsv` + per-accession outputs. Drives both Tier 1 (Figure A) and Tier 2 (audit) pipelines. |

## Output layout (after a full single-dataset run)

```
out/
├── library_report.json            # primer inference + extraction_mode (Phase 2)
├── selexprep_manifest.json        # reproducibility anchor (Phase 4)
├── trim_reports.json              # cutadapt argv + n_in/n_out per round (Phase 3)
├── strand_report.tsv              # only if orientation ∈ {MIXED, REVERSE}
├── extract_diff.tsv               # only with --rebuild + --override-primer-*
├── round_00/
│   ├── extracted.fasta.gz         # full-insert recovery (BOTH_PRIMERS_SINGLE_READ)
│   │   # or partial_5p_extracted.fasta.gz (FIVE_PRIME_ONLY)
│   │   # or partial_3p_extracted.fasta.gz (THREE_PRIME_ONLY)
│   │   # or partial_5p_extracted_R1.fasta.gz + partial_3p_extracted_R2.fasta.gz (PAIRED_END_SPLIT_PRIMERS)
│   └── counts.parquet             # unique sequences (Phase 5)
├── round_NN/...
└── qc/
    ├── flags.yaml                 # depth-aware suspicion flags
    ├── read_retention.png         # input vs extracted per round
    ├── primer_match_per_round.png
    ├── n_length_distribution.png  # per-round N-region histograms (faceted)
    └── per_round_panel.png        # unique seqs + Shannon entropy + top-100 coverage
```

## Determinism + reproducibility

All `.fasta.gz`, `.tsv`, and `.json` outputs are **byte-deterministic**
given the same `--sampling-seed`: gzip headers are written with
`mtime=0`, JSON keys are sorted with int-keyed dicts in numeric (not
lexical) order, and TSV rows are sorted. Two runs with identical
inputs + seed produce identical SHA256s; the manifest's
`output_sha256` field is the audit anchor.

Parquet hashes are **advisory** (not guaranteed across pyarrow
versions; the manifest pins `pyarrow_version` instead).

PNG plots are **informational** (matplotlib output is not
byte-deterministic across versions). They do not contribute to
`output_sha256`.

## What v0.1 does *not* do

By design — these are handled by mature existing tools that consume
`selexprep`'s outputs:

| Step | Use |
|---|---|
| Clustering | FASTAptameR |
| Motif discovery | MEME · RaptGen-UI |
| Binding-affinity prediction | RaptGen · DeepSELEX · AptaTrans |
| 3D structure | ViennaRNA · RNAfold |
| Aptamer design | MAWS · RNAtranslator |
| Read merging (paired-end full-insert recovery) | bbmerge · vsearch · pear (v0.2 hook) |

## Calibration status

Tests assert on **behavior**, never on threshold values (e.g.
`assert report.status == "HIGH"` for high-match-rate inputs, never
`assert HIGH_CUTOFF == 0.85`). Tuning the numbers is therefore safe
under the existing test suite.

**Phase 2 (LibraryReport inference)** — `CALIBRATION-REVIEWED` markers
in `library/detect.py` document the v0.1 values + rationale for each
threshold (`STATUS_HIGH_CUTOFF`, `POSITION_CONSISTENCY_TOLERANCE`, the
two `COMPOSITE_WEIGHTS` regimes, etc.). See `CHANGELOG.md` for the
diff history.

**Phase 5 (QC suspicion flags) + adapter blacklist composition** —
still pending. `CALIBRATION-TODO` markers in `qc/flags.py`,
`library/adapters.py` (TruSeq + Nextera vs full Illumina set), and
`extract/strand.py`. Inventory:

```bash
grep -rn "CALIBRATION-TODO" src/      # what's left to tune
grep -rn "CALIBRATION-REVIEWED" src/  # what's already vetted with rationale
```

Final calibration tuning will use Phase 6 benchmark recovery numbers
(15+ known-primer datasets) as empirical ground truth.

## Architecture

```
                     ┌──────────────────┐
                     │  catalog         │
                     │  (250 entries;   │
                     │   v0.1.7)        │
                     └────────┬─────────┘
                              │
                              ▼
   ┌─────────────┐    ┌──────────────┐    ┌───────────────┐
   │  inspect    │    │   fetch      │    │  detect       │
   │  (ENA REST  │    │  (ENA REST   │    │ (LibraryReport│
   │   preview)  │    │   download)  │───►│  schema)      │
   └─────────────┘    └──────────────┘    └───────┬───────┘
                                                  │
                              ┌───────────────────┘
                              ▼
                     ┌────────────────────┐     ┌──────────────────────┐
                     │  extract           │     │ manifest             │
                     │  (cutadapt + PE +  │────►│ (SelexprepManifestV1 │
                     │   strand)          │     │  reproducibility)    │
                     └────────┬───────────┘     └──────────────────────┘
                              │
                              ▼
                     ┌────────────────────┐
                     │  count             │
                     │  (FASTA → parquet) │
                     └────────┬───────────┘
                              │
                              ▼
                     ┌────────────────────┐
                     │  qc                │
                     │  (8 flags + 4 PNG) │
                     └────────────────────┘
```

## Development

```bash
# Pre-commit gates (run all four before pushing)
uv run pytest                       # 587 + 1 xfailed
uv run ruff check src/ tests/
uv run ruff format --check src/ tests/
uv run mypy src/                     # strict on library.report
```

CI matrix: Python 3.10 / 3.11 / 3.12.

## Benchmarks

Two-tier benchmark under `benchmarks/`:

- **Tier 1 — Figure A** (`Snakefile` + `ground_truth.tsv`): curated
  primer-recovery validation against paper-reported primers (N=11
  source-verified accessions, modality-diverse).
- **Tier 2 — Figure B** (`audit.smk`): a corpus utility audit over a
  random sample of audit-eligible INSDC catalog rows. Shipped audit
  artifacts live in `benchmarks/audit_results/`; reproducibility
  envelope (catalog version + sample sha256 + seed) is committed.
- **Catalog completeness audit** (`catalog_completeness_audit.py`):
  re-runnable one-off script that diffs `bioprojects.csv` against
  ENA's `library_strategy="SELEX"` set. Current snapshot:
  100% discovery / 79.6% auditable (82/103) with 21 documented
  exclusions for mis-labels + gSELEX variants.

See [`benchmarks/README.md`](benchmarks/README.md) for the methodology
and curation notes.

## License

MIT. See [`LICENSE`](LICENSE).

## Changelog

Full development log: [`CHANGELOG.md`](CHANGELOG.md).
