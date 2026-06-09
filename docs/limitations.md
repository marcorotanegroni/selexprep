# Limitations

v0.1 is intentionally narrow. The following are NOT supported and are deferred to v0.2 or are explicitly out of scope:

- **Read merging.** `PAIRED_END_SPLIT_PRIMERS` mode emits separate R1/R2 partial extractions; joined-count reconstruction requires bbmerge / vsearch / pear / flash and is v0.2.
- **IUPAC ambiguous bases in primers.** Reported as unsupported in the `LibraryReport`; excluded from exact-recovery counts.
- **Barcode inference.** Demultiplexing is sample-sheet driven only.
- **DNA vs RNA classification from FASTQ.** RNA SELEX is sequenced as cDNA, so reads are A/C/G/T regardless. Assay type is read from SRA metadata when available; classification is v0.2.
- **Clustering, motif discovery, binding-affinity prediction.** Handled by FASTAptameR, MEME / RaptGen-UI, RaptGen / DeepSELEX / AptaTrans respectively — `selexprep` outputs feed into these.
- **Parquet bit-identical reproducibility across `pyarrow` versions.** Only FASTA/TSV/JSON hashes are guaranteed deterministic; Parquet hashes are version-pinned in the manifest as advisory.

## Benchmark (Tier 1 scorecard) — how to read the recovery numbers

The primer-recovery benchmark scores recovery of the **paper-reported** primer from **real deposited reads**. A few honest caveats for interpreting it:

- **Small, non-random set.** The recovery arm is a hand-curated set of deposits whose library is paper-documented (7 in v0.1), not a random sample. Report exact / partial / miss counts separately — never a single headline rate or an "X/N" that implies a representative sample.
- **Limited lab independence.** Some rows share a lab or library family (e.g. the two RaptRanker RNA deposits), so the number of *independent* constructs is smaller than the row count.
- **Reads can differ from the paper template.** Deposited reads sometimes carry an extended or heterogeneous constant region — read-through into adapter/cloning sequence, or several co-dominant variants — that does not equal the clean paper primer. `detect` reports what is physically in the reads (surfacing `MEDIUM` / partial / mismatch) rather than fabricating the paper sequence, so a "partial" or "mismatch" often describes the **deposit**, not a tool error.
- **Read-resolved truths are not scored as paper recovery.** Where the only available truth for one side was derived from the deposited reads themselves, that side is flagged (`score_3p = false`) and excluded from the paper-grounded tally — scoring an inferred call against a read-derived truth would be circular. The pair is scored on the independent (paper-grounded) side alone.
- **Specificity / refusal arms are the complement.** On pre-trimmed (N-region-only) and adapter-colliding deposits the correct behavior is **no primer call**; those arms report false-positive calls against a target of 0, not recovery.
