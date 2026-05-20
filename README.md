# selexprep

**Accession-first preprocessing for public HT-SELEX, with primer auto-inference and safe failure modes.**

> *`selexprep` fills the missing preprocessing layer for public datasets by starting from accessions and empirically inferring primer/constant regions, extracting random regions, and emitting confidence-aware, reproducible count tables and manifests.*
>

## Status

**v0.1 RC.** The single-dataset workflow is feature-complete and tested
(358 passing tests + 1 strict-xfail reserved for v0.2 read merging).
Calibration constants are documented v0.1 placeholders pending a Codex
peer-review pass + Phase 6 benchmark recovery numbers. Not yet on PyPI.

**v0.1 ships:** `catalog` (273 public SELEX BioProjects) · `detect`
(primer inference) · `extract` (cutadapt-driven trimming + paired-end
+ strand handling) · `count` (per-round unique sequences) · `qc`
(suspicion flags + 4 PNG plots) · `inspect` (ENA metadata preview).

**v0.2 / Phase 6 deferrals:** `fetch` (accession download) · `run`
(batch driver across many accessions) · read merging for paired-end
full-insert recovery · AnnData export · BibTeX auto-citation · library-
type classification.

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

**Benchmark headline (Codex-honest, Phase 6):** *`selexprep` approaches
known-primer pipelines on datasets where primer inference is
high-confidence, while explicitly failing safe on ambiguous ones.*

## What v0.1 does (Codex-frozen core claim)

> *`selexprep` converts public or local HT-SELEX reads into primer-stripped
> random-region FASTA/FASTQ files, per-round count tables, QC reports,
> and provenance manifests, with explicit primer-inference confidence
> and safe failure modes.*

Concretely:

1. **Discovers** public SELEX BioProjects from a bundled catalog (273
   BPs refreshed against broad ENA queries; `selexprep catalog list /
   show`).
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

## CLI surface

| Command | Status | What it does |
|---|---|---|
| `selexprep catalog list \| show \| version \| refresh` | ✅ v0.1 | Browse the bundled discovery catalog (273 BPs). |
| `selexprep inspect <ACC>` | ✅ v0.1 | ENA filereport REST preview — round count, `library_strategy` (SRA verbatim, not classified), file sizes + MD5s. No download. |
| `selexprep detect <fastq...> --round-map rounds.tsv --outdir OUT` | ✅ v0.1 | Auto-infer primers + library structure → `library_report.json`. |
| `selexprep extract <fastq...> --library-report LR.json --round-map rounds.tsv --outdir OUT [--sample-sheet samples.tsv] [--paired-r2 ...] [--override-primer-{5p,3p} ...] [--rebuild]` | ✅ v0.1 | Cutadapt-driven trim + strand reorient + per-round FASTA + manifest. |
| `selexprep count <extracted.fasta.gz> --round R0 --outdir OUT` | ✅ v0.1 | FASTA → counts.parquet (sequence, reads, rank, RPM). |
| `selexprep qc <manifest> [--counts-dir DIR] [--outdir OUT]` | ✅ v0.1 | Depth-aware suspicion flags (YAML) + 4 PNG plots. |
| `selexprep fetch <ACC> --outdir OUT` | ⬜ v0.2 | Download FASTQ + auto-populate round map. |
| `selexprep run <accessions.tsv> --outdir OUT --resume` | ⬜ v0.2 | Batch driver across many accessions; emits corpus-level plot. |

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

## Calibration status (v0.1 placeholder values)

Nineteen heuristic thresholds across `library/detect.py`,
`qc/flags.py`, and `extract/strand.py` are documented v0.1
placeholders. Each carries a `# CALIBRATION-TODO` comment citing the
locked plan line it came from (or "no locked default" otherwise). The
full inventory:

```bash
grep -rn "CALIBRATION-TODO" src/
```

Tests assert on **behavior**, never on threshold values (e.g.
`assert report.status == "HIGH"` for high-match-rate inputs, never
`assert HIGH_CUTOFF == 0.80`). When Codex peer-review tunes the
numbers — together with Phase 6 benchmark recovery numbers — the test
suite stays green by construction.

## Architecture

```
                     ┌──────────────────┐
                     │  catalog (Phase  │
                     │  1.5 — 273 BPs)  │
                     └────────┬─────────┘
                              │
                              ▼
   ┌─────────────┐    ┌──────────────┐    ┌───────────────┐
   │  inspect    │    │   fetch      │    │  detect       │
   │  (ENA REST  │    │  (v0.2 —     │    │ (LibraryReport│
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
uv run pytest                       # 358 + 1 xfailed
uv run ruff check src/ tests/
uv run ruff format --check src/ tests/
uv run mypy src/                     # strict on library.report
```

CI matrix: Python 3.10 / 3.11 / 3.12.

## License

MIT. See [`LICENSE`](LICENSE).

## Citation

Publication pending — Bioinformatics Advances Application Note for
v0.1 (tool); NAR Database Issue for v0.2 (catalog +
LibraryReport-annotated corpus). DOI placeholders to follow.

## Acknowledgments

Built with Claude Code (Anthropic) under the four-round
Codex-frozen implementation plan at
`~/.claude/plans/unified-seeking-treehouse.md`. The full development
log is in [`CHANGELOG.md`](CHANGELOG.md).
