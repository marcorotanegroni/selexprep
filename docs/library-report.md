# LibraryReport reference

`detect` emits a `LibraryReport` (JSON) — the heart of `selexprep`. It records
what primers were found, how confident the call is, and what to do next. The key
design choice: **biology, workflow, and file layout are separate axes**, read
together. The full pydantic schema is in the [API reference](api.md); worked
real examples are in [`02_library_report_interpretation.ipynb`](examples.md).

## Fields

### Primers
- `primer_5p`, `primer_3p` — inferred constants (`null` if not recoverable).
  Each is the part of the flank that **every round shares**: the flank is
  called on the earliest round, then its outer edge is trimmed back to the
  sequence all rounds agree on, so a run-specific inline tag, a base that
  differs between rounds or a heterogeneous base ahead of the construct is
  left out. The boundary with the random region never moves.
- `variants_5p`, `variants_3p` — top candidates with read counts.

### Biology — `extraction_mode`
- `BOTH_PRIMERS_SINGLE_READ` — both flanks in one read (full insert).
- `FIVE_PRIME_ONLY` / `THREE_PRIME_ONLY` — one flank; N bounded on one side.
- `PAIRED_END_SPLIT_PRIMERS` — 5' on R1, 3' on R2, no overlap (needs merging).
- `UNABLE_TO_EXTRACT` — neither flank usable.

`full_insert_recovered` (bool) is `True` only for `BOTH_PRIMERS_SINGLE_READ`.

### Workflow — `required_action`
- `NONE` — proceed.
- `MANUAL_PRIMERS_REQUIRED` — supply `--override-primer-*` (or hand-edit the report).
- `READ_MERGING_RECOMMENDED` — paired split; full-insert recovery needs a merger (v0.2).

### File layout — `read_source`
`R1` / `R2` / `R1_AND_R2` / `INTERLEAVED` / `UNKNOWN`.

### Adapters — `known_adapter_hits`
Sequencing adapters seen (TruSeq, Nextera, …), recorded as counts and
**excluded from primer candidates** — never silently called as a primer, never
used to discard reads.

### Quality + N-region
- `match_rate_5p/3p`, `position_consistency_5p/3p` — how cleanly each constant
  sits at its flank, measured where it actually is (a run-specific tag outside
  the constant shifts it away from the read edge).
- `n_length_mode`, `n_length_distribution`, `n_length_confidence`.

### Confidence + reproducibility
- `confidence` (0–1) and `status`: `HIGH` / `MEDIUM` / `LOW` / `UNABLE_TO_INFER`.
- `failure_reason` — set when inference is refused.
- `read_fraction_used_for_inference`, `sampling_seed`.

## How detection maps to mode + action

| Detection result | `extraction_mode` | `full_insert_recovered` | `required_action` |
|---|---|---|---|
| 5' + 3' both > 0.7, same read | `BOTH_PRIMERS_SINGLE_READ` | `True` | `NONE` |
| 5' > 0.7, N-length conf > 0.8 | `FIVE_PRIME_ONLY` | `False` | `NONE` |
| 5' > 0.7, N-length conf ≤ 0.8 | `UNABLE_TO_EXTRACT` | `False` | `MANUAL_PRIMERS_REQUIRED` |
| 3' > 0.7, N-length conf > 0.8 | `THREE_PRIME_ONLY` | `False` | `NONE` |
| 3' > 0.7, N-length conf ≤ 0.8 | `UNABLE_TO_EXTRACT` | `False` | `MANUAL_PRIMERS_REQUIRED` |
| 5' on R1, 3' on R2, no overlap | `PAIRED_END_SPLIT_PRIMERS` | `False` | `READ_MERGING_RECOMMENDED` |
| both match rates < 0.4 | `UNABLE_TO_EXTRACT` | `False` | `MANUAL_PRIMERS_REQUIRED` (`status = UNABLE_TO_INFER`) |

## Safe failure

When `status == UNABLE_TO_INFER` or `extraction_mode == UNABLE_TO_EXTRACT`,
downstream `extract` refuses without an explicit `--override-primer-*` or a
hand-edited report, rather than cutting where it cannot place the boundary.
Refusal covers the cases `detect` recognises; the limits below are cases it
does not.

`extract` also warns when an input FASTQ keeps less than half of its reads, or
less than half the yield of the best input (`trim_reports.json` holds every input's reads in and out): the
primers fit that run's reads badly, for example because it carries a different
construct.

!!! note
    Single-round deposits cap `status` at `MEDIUM` (cross-round persistence is
    unavailable). Below ~500 sequences in the earliest round, `detect` refuses
    outright.

## Known limits

`detect` finds the boundary from how conserved each position is across the
reads of the earliest round. Sequence that is conserved in the random region
itself can therefore be taken for constant, and the call then reaches a few
bases into the random region. Seen in the semi-synthetic benchmark
(`benchmarks/synthetic/`):

- **One sequence dominating every round, the earliest included.** With one
  sequence making up 70% of every round, 1 run in 24 placed the boundary inside
  the random region; at 80%, 17 in 24, all with `status` `HIGH`. From 90% on,
  every run was refused. What protects the call is an earliest round in which
  no sequence dominates.
- **Only a late round provided, with a motif next to the constant.** The motif
  enters the call (`status` `MEDIUM`, since one round caps it).

A high extraction yield does not show that the boundary is right: a call a few
bases into the random region still trims almost every read. When the deposit
has an early, diverse round, provide it; when it has none, check the called
constants against the publication.
