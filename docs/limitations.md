# Limitations

v0.1 is intentionally narrow. The following are NOT supported and are deferred to v0.2 or are explicitly out of scope:

- **Read merging.** `PAIRED_END_SPLIT_PRIMERS` mode emits separate R1/R2 partial extractions; joined-count reconstruction requires bbmerge / vsearch / pear / flash and is v0.2.
- **IUPAC ambiguous bases in primers.** Reported as unsupported in the `LibraryReport`; excluded from exact-recovery counts.
- **Barcode inference.** Demultiplexing is sample-sheet driven only.
- **DNA vs RNA classification from FASTQ.** RNA SELEX is sequenced as cDNA, so reads are A/C/G/T regardless. Assay type is read from SRA metadata when available; classification is v0.2.
- **Clustering, motif discovery, binding-affinity prediction.** Handled by FASTAptameR, MEME / RaptGen-UI, RaptGen / DeepSELEX / AptaTrans respectively — `selexprep` outputs feed into these.
- **Parquet bit-identical reproducibility across `pyarrow` versions.** Only FASTA/TSV/JSON hashes are guaranteed deterministic; Parquet hashes are version-pinned in the manifest as advisory.
