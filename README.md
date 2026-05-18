# selexprep

**Accession-first preprocessing for public HT-SELEX datasets with primer auto-inference.**

> `selexprep` is not another aptamer-analysis suite. It is the missing accession-to-clean-library preprocessing layer, with primer inference, uncertainty reporting, random-region extraction, and reproducible outputs.

## Status

**Pre-release.** v0.1 in active development. Not yet on PyPI.

## What it does

Given an ENA / SRA / DDBJ accession (or local FASTQ files), `selexprep`:

1. Downloads FASTQ from the public archive
2. Auto-infers the SELEX library structure: primer pair, N-region length, orientation, paired-end layout — directly from read content, without supplementary PDFs
3. Reports inference confidence with an explicit `LibraryReport` (`extraction_mode` + `read_source` + `required_action` + `full_insert_recovered`)
4. Extracts the randomized region per round, with safe failure modes when inference is ambiguous
5. Emits trimmed per-round FASTAs + count tables + QC plots + a provenance manifest consumable by downstream tools (FASTAptameR / RaptGen-UI / APTANI2)

## What it does *not* do

By design — these are handled by mature existing tools that consume `selexprep`'s output:

| Step | Use |
|---|---|
| Clustering | FASTAptameR |
| Motif discovery | MEME / RaptGen-UI |
| Binding-affinity prediction | RaptGen / DeepSELEX / AptaTrans |
| 3D structure | ViennaRNA / RNAfold |
| Aptamer design | MAWS / RNAtranslator |

## Quick start

*Coming soon — Phase 0 scaffold only.*

```bash
pip install selexprep  # not yet published
selexprep inspect SRR12647619
```

## Roadmap

v0.1 features (in development):
- Accession or local FASTQ input
- Primer inference with `LibraryReport` (confidence + `extraction_mode`)
- Safe manual override (`--library-report` or `--override-primer-*`)
- Strand/orientation detection + auto-reorient
- Random-region extraction (R1, R2, or split — no merging in v0.1)
- Per-round unique-sequence counts
- QC plots + suspicion flags (depth-aware)
- Reproducibility manifest (deterministic sha256s on FASTA/TSV/JSON)

v0.2 deferrals: read merging, AnnData output, BibTeX auto-citation, HTML report, library-type classification, corpus-level analytics.

## License

MIT. See [`LICENSE`](LICENSE).

## Citation

Application Note in preparation for *Bioinformatics Advances*.
