# Parameters of `detect` and `extract` (selexprep 0.4.5)

Values read from the installed package by `benchmarks/parameters/parameter_table.py`.

| Stage | Parameter | Value | What it controls |
|---|---|---|---|
| detect: input | `DEFAULT_MIN_SEQS_FOR_DETECTION` | 500 | Fewest sequences in the earliest round for detect to infer anything; below it, it refuses. |
| detect: input | `DEFAULT_TOP_N` | None (all) | Sequences used for inference: None means every unique sequence, no subsampling. |
| detect: flank boundary | `DEFAULT_MIN_LEN` | 14 | Shortest flank called as a library constant (nt). |
| detect: flank boundary | `DEFAULT_MAX_LEN` | 60 | Longest flank examined from each read edge (nt). |
| detect: flank boundary | `DEFAULT_CONFIDENCE` | 0.75 | Positional support for a 'strong' position of the flank. |
| detect: flank boundary | `BOUNDARY_SUPPORT_FLOOR` | 0.55 | Support below which a position is no longer constant; the boundary stops at the first run of such positions. |
| detect: flank boundary | `BOUNDARY_DROP_RUN` | 2 | Consecutive positions that must fall below a threshold to end the flank, so one noisy base does not. |
| detect: flank boundary | `BOUNDARY_STRONG_POSITION_FRACTION` | 0.8 | Fraction of the called flank's positions that must be strong. |
| detect: flank boundary | `BOUNDARY_HIGH_SUPPORT_BASELINE` | 0.9 | Drop test: median support of the preceding positions for the test to apply. |
| detect: flank boundary | `BOUNDARY_HIGH_SUPPORT_DROP` | 0.12 | Drop test: fall in support from that median that marks the boundary. |
| detect: flank boundary | `BOUNDARY_HIGH_SUPPORT_POST_MAX` | 0.85 | Drop test: support the following positions must stay below. A sequence that makes up this much of every round, the earliest included, can pass for constant (Known limits). |
| detect: core rescue | `CORE_MIN_SUPPORT` | 0.9 | Support of the high-support core used when the whole flank's match rate fails. |
| detect: core rescue | `CORE_MIN_LEN` | 12 | Shortest core that can be rescued (nt). |
| detect: core rescue | `CORE_STABLE_RUN` | 3 | Consecutive well-supported positions that end the outer-edge scan. |
| detect: shared constant | `CONSTANT_CORE_LEN` | 12 | Inner bases of the flank, next to the random region, on which every round is aligned; never trimmed (nt). |
| detect: shared constant | `CONSTANT_SEARCH_SLACK` | 12 | How far beyond the earliest round's flank the core may sit in another round (nt). |
| detect: shared constant | `CONSTANT_MIN_ANCHORED_FRACTION` | 0.5 | Fraction of a round's reads that must carry the core for the round to vote. |
| detect: shared constant | `CONSTANT_ANCHOR_AGREEMENT` | 0.8 | Agreement among single-copy reads needed to anchor reads that carry the core twice. |
| detect: classification | `PRIMER_FOUND_MATCH_RATE_THRESHOLD` | 0.7 | Fraction of reads matching a flank for that side to count as found. |
| detect: classification | `UNABLE_TO_EXTRACT_MATCH_RATE` | 0.4 | Match rate below which, on both sides, extraction is refused. |
| detect: classification | `N_LENGTH_CONFIDENT_FRACTION` | 0.8 | Share of reads at the modal random-region length for a one-sided extraction to be allowed. |
| detect: classification | `POSITION_CONSISTENCY_TOLERANCE` | 3 | Offset of a flank, in nt, still counted as the same position. |
| detect: status | `STATUS_HIGH_CUTOFF` | 0.85 | Composite confidence for HIGH. Rates support for the inferred primers, not a verified boundary. |
| detect: status | `STATUS_MEDIUM_CUTOFF` | 0.6 | Composite confidence for MEDIUM; a single round caps the status at MEDIUM. |
| detect: status | `STATUS_LOW_CUTOFF` | 0.3 | Composite confidence for LOW; below it, UNABLE_TO_INFER. |
| detect: status | `COMPOSITE_WEIGHTS` | match_5p 0.15, match_3p 0.15, pos_5p 0.15, pos_3p 0.15, persistence 0.25, n_len 0.1, adapter_clean 0.05 | Weights of the signals in the composite confidence, with a round map. |
| detect: status | `COMPOSITE_WEIGHTS_NO_ROUND_MAP` | match_5p 0.225, match_3p 0.225, pos_5p 0.225, pos_3p 0.225, persistence 0.0, n_len 0.05, adapter_clean 0.05 | The same weights without a round map (no cross-round persistence). |
| detect: orientation | `ORIENTATION_REVERSED_FORWARD_MAX` | 0.05 | Share of reverse reads below which the deposit is FORWARD. |
| detect: orientation | `ORIENTATION_REVERSED_REVERSE_MIN` | 0.95 | Share of reverse reads above which it is REVERSE (between the two: MIXED). |
| detect: report | `VARIANTS_TOP_K` | 3 | Flank variants listed in the report, per side. |
| detect: adapters | `ADAPTER_PROBE_K` | 13 | Leading bases of a known adapter that mark a call or a read as adapter (nt). |
| detect: adapters | `KNOWN_ADAPTERS` | NEXTERA, TRUSEQ_R1 | Sequencing adapters excluded from primer candidates (names listed). |
| extract: yield warning | `LOW_YIELD_RATIO` | 0.5 | An input keeping less than this fraction of the best input's yield is reported. |
| extract: yield warning | `LOW_YIELD_FLOOR` | 0.5 | An input keeping less than this fraction of its own reads is reported. |
| extract: cutadapt | `--error-rate` | 0.1 | Maximum error rate in a primer match (cutadapt default). |
| extract: cutadapt | `--overlap` | 3 | Shortest primer overlap at a read end (cutadapt default). |
