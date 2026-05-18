# v0.2 roadmap

Deferred from v0.1 (paper mentions, doesn't ship):

- **AnnData output** — `.h5ad` export for ML-pipeline consumers
- **BibTeX auto-citation** — BioProject → PubMed E-utils → `citations.bib`
- **HTML report per dataset** — single self-contained `report.html` with all QC plots + flags + LibraryReport
- **Library-type classification** — metadata-driven DNA/RNA hint surface
- **Optional read merging** — bbmerge / vsearch / pear / flash subprocess wrappers; completes the `PAIRED_END_SPLIT_PRIMERS` path with joined counts
- **Corpus-level analytics** — cross-dataset queries, primer-reuse detection, last-round-uniques histogram with statistical summaries
- **Barcode inference** — automatic barcode discovery when no sample sheet is provided
