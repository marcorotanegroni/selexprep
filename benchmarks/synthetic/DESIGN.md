# Semi-synthetic benchmark — design (pre-registration)

This file fixes the design **before** any simulated library is generated or
scored. It is committed first. The donor screen, the generator, the evaluator
and the runs come after it, and the commits they ran at are recorded with every
result. What may change afterwards, and how, is set out in
[Freezing and changes](#freezing-and-changes).

Companion files, fixed with this one: `constants.tsv` (constant families),
`donors.tsv` (donors known today, completed by the donor screen) and
`configurations.tsv` (one row per configuration, with the families it applies
to).

## Why this benchmark exists

The real Tier-1 benchmark (`../README.md`) has seven recovery deposits and,
since the per-run-tag fix of 2026-09-27, is a development set: its post-fix
numbers are in-sample. It also cannot know, read by read, where the random
region is. This benchmark asks:

1. **Does `detect` place the constant/random boundary correctly and, when it
   cannot, does `extract` refuse instead of cutting in the wrong place?** With
   the true boundary known in every read.
2. **Do the rule that keeps only the constant every round shares, and its
   anchoring on the inner core, hold beyond the deposits and the code review
   that motivated them?**
3. **Where does the current version stop working?** Short constants, sequencing
   error, rounds without a diverse pool, converged pools, sequence next to the
   boundary: the operating envelope, measured rather than asserted.

**Scope of the conclusions.** The benchmark can show that a failure mechanism
exists and whether it is handled, under the conditions simulated. It cannot
estimate how often any mechanism occurs among public deposits: the conditions
are chosen by us, and the biological material comes from a few donors. Many
simulated runs on few donors are not many independent experiments. It does not
check called constants against publications (Tier-1 does), and does not cover
paired-end layouts or demultiplexing inside a run.

## Simulated libraries

Each simulated deposit has one or more rounds; each round is a pool of
single-end reads built as

```
[5' outer][5' constant][N][3' constant][3' outer][read-through]
```

**N, the random region, is taken from real pre-trimmed deposits** (the donors),
so it keeps their base composition, biases, sequencing errors, conservation and
round-to-round dynamics. From one documented selection trajectory per donor,
three rounds are used: the earliest, the middle (index ⌊(n−1)/2⌋ of the n
rounds with reads available) and the latest. The trajectory is the one with the
most rounds with reads available; ties go to the trajectory with the lowest run
accession. The first available round is not necessarily a naive pool, and the
design does not assume it is.

**Which donor reads are used.** Only reads of length N ± 2 nt (N as stated
outside the reads) and without independent evidence of truncation; the rest are
excluded before sampling, and the excluded fraction is recorded per round. The
retained reads may differ in composition from the whole pool; that fraction is
what the benchmark represents. 10,000 reads per round are sampled without
replacement from the retained reads (seeded); a round with fewer than 10,000
retained reads is replaced by the next round inward and the substitution is
recorded.

**Constants are published library constants** (`constants.tsv`, from
`ground_truth.tsv`, in read orientation). Shorter constants keep the inner part
next to N.

**Construction order.** For every read: assemble the true molecule (outer
sequences, constants, N); apply the configuration's edits to the molecule
(round-specific base, adapter replacement, core copy, motif, dominant clone,
biological N±1); truncate (F); add sequencing indels (I), then substitutions
(E); reverse-complement last (O). Substitutions occur at rate `r · (0.5 + i/L)`
at read position `i` of a read of length `L` (mean `r`, rising toward the 3′
end), the new base drawn uniformly from the other three; indels at the stated
per-base rate, split equally between insertions and deletions. Every
coordinate is carried through these steps.

### Truth, per read

The read name carries the coordinates `[start, end)` of N **in the observed
read**, after every step, and whether the read is recoverable by two-sided,
5′-only and 3′-only extraction (the constants that extraction needs are present
and not truncated). Evaluation compares extracted sequences with that observed
span, errors included: `selexprep` is asked to cut in the right place, not to
correct sequencing errors. `extract` keeps read names, so each extracted
sequence is traced to its truth.

## Development and test sets

Donors and constant families are split between the two sets, with the overlaps
declared below. Seeds 1–3 are used in development and 101–103 in test.

| | Donors (random regions) | Constant families |
|---|---|---|
| **dev** | PRJNA360902 (DNA, N20, GATA4_5_6 rounds 1/4/7); PRJEB22637 (2′-F RNA, N50, starting library / round 7 / round 15) | PRJDB9110 (40/19 nt, T7), PRJNA615076 (23/23), PRJEB62495 (20/20), PRJNA990511 (20/20) |
| **test** | PRJEB70964 (2′-F RNA, N35, S1 rounds 1/8/15); PRJEB47428 (RNA, N40, ELAVL1 cycles 1/2/4); PRJEB49150 (DNA, N40, NACC1 rounds 1/2/3); up to three more from the donor screen | PRJEB22637 (23/24), PRJEB28411 (18/18), PRJEB70964 (23/23) |

Declared overlaps and limits:

- PRJEB22637 is a dev donor and a test constant family; its random regions and
  its constants are different sequences.
- The 3′ constants of PRJDB9110 and PRJNA990511 share their first 13 nt, so both
  are dev families; the test families share no such prefix.
- PRJEB70964 is a test donor and a test constant family: some test runs rebuild
  that library from its own reads and its own published constants.
- The only documented naive pool (PRJEB22637's starting library) is in dev. The
  test set measures generalisation to trajectories whose first available round
  is already selected.
- The dev families come from Tier-1 recovery deposits; the test families from
  deposits whose constants were never used to develop `detect`. All current
  donors are Tier-1 specificity deposits, which did not shape the
  shared-constant rule.
- The length thresholds and the donor rules below were written after examining
  the current donors; their eligibility is in-sample. Freezing protects only the
  candidates the screen adds.

**Replication.** Three kinds, reported separately: seeds (technical variation:
sampling and simulated error), configurations (conditions and interactions),
donors and families (biological diversity). Donors and families are crossed,
so two cells can share a donor or a family; cells are not independent
experiments, and neither are runs. Results are descriptive and stratified;
no statistic treats runs or cells as independent.

## Donor screen (to add test donors)

Run once, after this file is committed, on data and metadata only; never on
`selexprep`'s behaviour.

1. **Candidates:** INSDC deposits in the bundled catalogue whose curated
   metadata state the randomised-region length N, other than the current
   donors.
2. **Length screen:** the archive's mean read length (`base_count /
   read_count`, halved for paired layouts) is within N ± 2 nt.
3. **Evidence that the reads are N**, from independent sources taken together:
   documentation of the preprocessing (submitter or publication), the published
   architecture, and an all-read check (`../verify_read_state.py`) in which the
   modal read length is N and at least 70% of reads lie within N ± 2. Sufficient
   evidence admits the donor even when N carries conserved blocks.
4. **Technical residue:** an exact fragment of at least 8 nt from the part of a
   published constant next to N, or from the start of a known adapter (TruSeq
   R1/R2, Nextera, small-RNA 3′ adapter), at the read edge where a residue would
   sit, in at least 5% of the reads of a chosen round, **triggers an
   investigation**, not an exclusion. The investigation records the diversity
   of the reads that carry it (distinct sequences / carrying reads), whether it
   is present in the earliest round, and the documentation; diversity is a clue,
   not a test. The donor is excluded as technical only when the evidence taken
   together shows a residue; the reason is recorded.
5. **Insufficient evidence** on completeness or boundaries makes the candidate
   **inconclusive**, recorded with its reason and not used. Conservation,
   diversity and convergence of N are recorded as donor descriptors and are
   never a reason for exclusion.
6. **Trajectory:** a documented biological trajectory with at least three
   rounds, each chosen round with at least 10,000 retained reads.
7. **Independence from dev:** no shared BioSample or selection with a dev donor.
8. **Selection:** at most six test donors in total. Priority to a chemistry not
   yet represented, then an N length not yet represented, then a different
   laboratory or publication; remaining ties go to the lowest accession.

Every candidate is written to `donor_screening.tsv` with each criterion's result
and the decision; admitted donors are added to `donors.tsv`. If the screen
brings the test donors to five or six, test uses seeds 101–102 only.

## Configurations

A baseline, one factor changed at a time, and a few targeted interactions.
Baseline: three rounds, 10,000 reads per round, full-length published
constants, substitution rate 1%, no indels, no outer sequence, reads ending at
the 3′ constant, donor N as retained. `configurations.tsv` lists every
configuration with the families it applies to; a configuration that would
reproduce the baseline for a family (for instance length 20 on a 20/20 family)
does not apply to that family.

| ID | Factor | Levels | Exact transformation |
|---|---|---|---|
| L | Constant length | 8, 10, 12, 14, 16, 20 | each constant keeps its inner `k` nt (5′: last `k`; 3′: first `k`); a flank shorter than `k` stays whole |
| E | Substitution rate | 0, 0.5, 2, 5% (baseline 1%) | error model above |
| I | Sequencing indels | 0.1, 0.5% per base | error model above |
| T1 | Tag per round | — | per round, a random 5′ tag and 3′ tag, lengths drawn independently from 4–8 nt |
| T2 | Tag shared by all rounds | — | the same random 6 nt at each end in every round |
| T3 | Heterogeneous leading base | — | one uniformly random base per read before the 5′ constant |
| T4 | Round-specific base | — | in the earliest round, the outermost 3′-constant base replaced by a different base, the same in every read |
| T5 | Adapter read-through | — | after the 3′ constant, the first 20 nt of the TruSeq Read 1 read-through, `AGATCGGAAGAGCACACGTC` |
| A | TruSeq R1 reverse complement in the 5′ constant | outer start; middle; inner end | `GCTCTTCCGATCT` replaces 13 nt of the 5′ constant at positions 0–12; centred; the last 13 |
| K | Core copy | 5′ outer; 3′ outer; inside N at 5′; inside N at 3′ | the constant's inner 12 nt, separated from the constant by 4 random nt, in the outer sequence; or overwriting the first (5′) or last (3′) 12 nt of N in 60% of the reads of the middle and latest rounds |
| R | Rounds provided | earliest only; latest only; middle + latest | subsets of the three rounds |
| C | Dominant clone | 50, 70, 80, 90, 99, 100% | in every provided round, that fraction of reads replaced by one donor read drawn with the seed from the earliest round |
| D | Depth and diversity | 500; 2,000; 1,000 distinct | subsample of the 10,000 reads; 1,000 distinct N resampled with replacement to 10,000 |
| M | Motif at the start of N | 50%; 90% | 6 random nt (fixed per seed) overwrite the first 6 nt of N in that fraction of the reads of the middle and latest rounds |
| F | Missing or truncated flank | no 3′ constant; 30% truncated in the 3′ constant; no 5′ constant | constant removed; read cut at a uniform position inside the 3′ constant |
| B | Biological N±1 | — | 10% of molecules lose one base of N, 10% gain one random base |
| O | Orientation | — | every read reverse-complemented |
| X1 | T1 × substitution rate | 0.5, 2, 5% | |
| X2 | Length × dominant clone | 14 and 16 nt at 70% | |
| X3 | Latest round only × motif | motif in 90% | |
| NC | Negative controls | donor reads alone; one sequence alone | no constants |

52 configurations with the baseline. Development runs every configuration for
every dev donor × applicable dev family × seed; test likewise with test donors
and families; negative controls, which carry no constants, run once per donor ×
seed. From `configurations.tsv`: 1,200 development runs, and 450 test runs per
test donor with three seeds (1,350 with the three current test donors). Every
run is also extracted with a report built from the configuration's truth (the
oracle, below).

## Evaluation

The evaluator reads the generator's truth and `selexprep`'s output files only;
it imports nothing from `selexprep.library` or `selexprep.extract`.

**Three levels of correctness**, kept apart:

1. **The call.** A called 5′ constant is correct when it is a suffix of the true
   sequence before N (outer + constant, read orientation) in the earliest round
   provided; a 3′ constant, when it is a prefix of the true sequence after N.
2. **The extraction, against the declared mode**, per read: two-sided output
   must equal the observed N span; `FIVE_PRIME_ONLY` output must equal the read
   from the start of N to its end; `THREE_PRIME_ONLY` output, the read from its
   start to the end of N. Precision among emitted reads; recall over the reads
   recoverable for that mode.
3. **Full recovery of N**: reads whose emitted sequence equals their observed N
   span, over the reads recoverable by two-sided extraction.

**One outcome per run**, the first that applies:

1. **Error** — the pipeline crashed or timed out (not `extract`'s refusal).
   Reported separately, kept in every denominator.
2. **Negative control** — refused is correct; anything else is **wrong**.
3. **Refused** — `extract` refuses: `status == UNABLE_TO_INFER` or
   `extraction_mode == UNABLE_TO_EXTRACT`.
4. **Wrong** — a side used for extraction has an incorrect call. A `MEDIUM` or
   `LOW` status does not change this.
5. **Complete, correct** — two-sided mode, both calls correct.
6. **Partial, correct and declared** — one-sided mode, the used side correct.

An outcome is never relabelled after the fact: a call that reaches into N is
wrong even where N carries a block conserved in every round, and is then
interpreted within that donor's stratum. Calling such a case unidentifiable
from the reads is a separate claim that needs its own argument.

Zero denominators are reported as not applicable; a run with no emitted read
has precision not applicable and recall 0.

**Reported per configuration**, stratified by donor and by family: outcome
counts; coverage (runs not refused) and accuracy among accepted runs; per-read
precision and recall (median and range over runs) for the inferred constants
and for the oracle; lowest per-round yield; the called constant relative to the
true constant (exact, includes outer material, shorter, other); outcomes by
`status`. `confidence` is a weighted sum of signals, not a calibrated
probability, and the status table is presented as descriptive.

**The oracle.** For each run, a `LibraryReport` written with the public schema
from the configuration's truth: the true library constants in read
orientation, the extraction mode the truth allows (two-sided, or one-sided when
F removes a flank), orientation FORWARD, status HIGH. `extract` runs on it
directly, without inheriting `detect`'s decisions or relying on primer
overrides (which do not promote a one-sided mode). The difference between the
oracle and the inferred run is the cost of inference. Not applicable to
negative controls.

## Expectations, fixed in advance

For each factor, the **desirable** outcome and the outcome **expected from the
current version** (`main` at b8ded35), kept apart: a limit met exactly as
predicted is reported as a limit. "Uncertain" means we do not know.

| ID | Desirable | Expected from the current version |
|---|---|---|
| Baseline | complete, correct | complete, correct; per-read precision near 1 |
| L ≥ 14 | complete, correct | complete, correct |
| L 8–12 | correct, or refused | refused: a flank shorter than `min_len` (14) flush with the read edge is not called (outer sequence would lengthen the called flank; L adds none). **Risk:** if the first bases of N reach 0.55 support, the call can reach into N |
| E 0–2% | complete, correct | complete, correct |
| E 5% | correct, or partial / refused | uncertain: whole-constant match rates approach the 0.70 threshold; the 3′ side, where errors are highest, may be dropped |
| I | complete, correct | complete, correct; per-read precision slightly lower |
| T1 | complete, correct; constant without the tags | complete, correct; the constant is the shared part, which can include tag bases all rounds happen to share |
| T2 | boundary correct | boundary correct; the call is tag + constant, including outer material |
| T3 | complete, correct; constant without the base | complete, correct; the base is left out |
| T4 | complete, correct; constant without the differing base | complete, correct; the call stops before that base |
| T5 | complete, correct; constant without the adapter | boundary correct; the read-through stays in the 3′ call: only a call that *starts* like an adapter is dropped |
| A outer start | complete, correct | uncertain between partial and refused: a 5′ call whose first 13 nt match an adapter is dropped |
| A middle / inner end | complete, correct | complete, correct: the adapter check does not see it |
| K | complete, correct | complete, correct: reads with two copies are anchored where single-copy reads agree; a core repeated in the outer sequence keeps the flank whole (includes outer material) |
| R earliest only | complete, correct | complete, correct; status at most MEDIUM |
| R latest only; middle + latest | correct, or refused | uncertain; depends on the donor's late rounds |
| C 50% | complete, correct | complete, correct: support in N reaches 0.85 only if the rest of the pool agrees at 70% |
| C 70%, 80% | correct, or refused | uncertain: support in N near the 0.85 of the cliff test |
| C 90–99% | correct, or refused | refused: the call runs through N, the flanks meet and the zero-length guard fires |
| C 100% | refused: the boundaries are not identifiable from conservation | refused, by the zero-length guard |
| D 500 | correct, or refused | uncertain: 500 reads is the detection floor |
| D 2,000; 1,000 distinct | complete, correct | complete, correct |
| M | complete, correct | complete, correct: the calls are made on the earliest round, where the motif is absent |
| F no 3′ / no 5′ | partial, correct and declared | partial, correct and declared |
| F 30% truncated | complete, or partial correct | uncertain: the 3′ match rate lands near 0.70 |
| B | complete, correct; N±1 reads extracted exactly | complete, correct |
| O | complete, correct | complete, correct, constants in read orientation — **except the PRJEB70964 family**, whose reverse-complemented 3′ constant starts with TruSeq R1 (`AGATCGGAAGAGC`) and is dropped: partial or refused |
| X1 | as T1 until errors make the constant unreliable | as E |
| X2 | complete, correct | uncertain |
| X3 | correct, or refused | **wrong:** with only the latest round, a motif conserved next to the constant enters the call |
| NC | refused | refused |

## Freezing and changes

1. This file and its companions are committed before the donor screen.
2. The donor screen runs; its results are committed (`donor_screening.tsv`,
   `donors.tsv`) as a dated amendment.
3. The generator, the harness and the evaluator are written with unit tests and
   committed.
4. **Development**: dev runs may lead to changes in `selexprep`, the generator
   or the evaluator; each is recorded below, validated by rerunning dev, and dev
   results before and after are both kept.
5. **Test**: `selexprep`, the generator and the evaluator are frozen at recorded
   commits; the test set is generated and run **once** and reported as it is.
6. A defect found by test is fixed for users afterwards; the test results stay
   tied to the commit that was evaluated, and the test set becomes development
   material.

Changes after this file is committed are allowed only before the test run, and
only for reasons that do not come from test outcomes.

## What goes into the paper

Test-set results only: one sentence in the Results (outcomes, accuracy among
accepted runs, per-read precision, and where wrong calls fell relative to the
expectations above), with conclusions limited to the conditions simulated; a
supplementary table with every configuration, stratified by donor and family,
and the three kinds of replication stated separately. Development results, and
any change they led to, are reported in the supplement as development.

## Outputs for the package

- A small deterministic regression suite in `tests/` built from dev
  configurations, with the expected outcome of each.
- The outcomes-by-status table, to state in the documentation how often each
  status accompanied a correct call.
- The absolute yield check in `extract` (`LOW_YIELD_FLOOR = 0.5`, fixed a priori
  in a106c3d: below every healthy Tier-1 run, the lowest being 63%) is
  evaluated on dev and test, not tuned.

## Amendments

### Amendment 1 — 2026-09-27, before any read of a candidate donor was examined

**What.** The trajectory rule (step 6 of the donor screen, and "Simulated
libraries") is made precise:

1. **Round documentation.** A run's round is taken from an explicit round or
   cycle label in the archive records (ENA titles, aliases and library names;
   the ArrayExpress SDRF; GEO sample characteristics) or from the depositor's
   publication. `selexprep`'s round parser is only a first pass.
2. **Only trajectories selected against a target.** No-target controls
   (bead-only, no protein) are not donor trajectories.
3. **Input library.** A documented input library (round or cycle 0, no target)
   is round 0 of every trajectory it seeded.
4. **Replicates.** Of several runs of one round, the lowest run accession.
5. **Ties** in the number of rounds go to the lowest run accession among the
   runs not shared with another trajectory.

**Why.** Running the metadata screen showed that `selexprep`'s parser reads
sample identifiers as rounds (PRJEB51212 and PRJEB51473: "RV01"…"RV39" became
rounds 1–39, while ArrayExpress documents cycles 0, 1, 3 and 6 with three
replicates each) and misses round labels in sample aliases, SDRFs and GEO
records; and that "the trajectory with the most rounds" picked no-target
controls (PRJEB61115, PRJEB51473), whose weak enrichment would make the test
easier than intended.

**Disclosure.** Written after seeing the candidates' metadata, before any of
their reads was examined.

**Applied** to all ten deposits that pass the length screen; the review is in
`trajectory_documentation.tsv`, and `screen_donors.py` regenerates
`donor_screening.tsv` from it. Three candidates go on to the read check
(PRJEB51212, PRJEB51473, PRJNA741127). Two are pending the publication
(PRJEB14550, PRJEB25907), whose cycle numbers are known only from the naming of
the samples. Rule 1 also applies to the current donors: PRJEB49150's rounds were
inferred from the trailing digit of its sample titles, so it is pending the same
check; the other four current donors carry explicit round labels and are
unchanged.

**Publication check (same day, before any read was examined).** One standard for
all three: the publication states that the deposited samples of a trajectory are
its successive selection cycles, and the run suffixes are consistent with them.
PRJEB49150 (Genes Dev 2022, doi:10.1101/gad.348993.121): the bound populations
were sequenced after all three selection cycles, three runs per protein with
suffixes 1–3 — it stays a test donor. PRJEB25907 (Genome Res 2020,
doi:10.1101/gr.258848.119): three or four selection cycles, ligands from the
cycles sequenced; aliases end in 1–4 — passes. PRJEB14550 (NAR 2018,
doi:10.1093/nar/gky027): the paper analyses an early cycle of FLI1 (143,389 reads,
the run suffixed 2) against the previous cycle, which is deposited with suffix 1;
that supports suffixes 1 and 2 and not the suffix-4 run chosen as the latest
round, so its trajectory is **inconclusive** — the reasoning is in amendment 2.
Four candidates go on to the read check (`donor_screening.tsv`).

### Amendment 2 — 2026-09-27, before any simulated library was generated

An external review of the written design, its screen and its tables found eleven
points that needed a decision: five in the design, four in the screen, and two
in how the donor pools and the read job are described. None comes from a
simulated result, because nothing has been generated yet. The body above is left
as it was committed and each item below names the sentence it replaces.

**Design.**

1. **Recall denominators** (replaces the recall clauses of levels 2 and 3 of
   "Evaluation"). The numerator must lie inside the denominator. Measured, not
   assumed: with the linked non-anchored adapters `extract` uses, a read cut
   inside the 3′ constant is still trimmed correctly from the residual fragment.
   Of 200 reads cut at a uniform position inside a 20-nt 3′ constant (cutadapt
   5.2), every read that kept at least 3 nt of it came out equal to its true
   random region, none came out wrong, and only those left with 0–2 nt (below
   cutadapt's 3-nt minimum overlap) were lost. Those reads are not
   "recoverable" under the truth definition, so they would have entered the
   numerator of `F_trunc30` and not the denominator, giving more than 100%.
   Therefore: for each mode, recall = reads that are eligible **and** correct /
   eligible reads; and the recovery of the random region overall is reported
   separately over **all reads that contain it**, which is where a recovery
   beyond the eligible set is credited. Precision stays over every emitted read.
   A zero denominator is not applicable, with or without output.
2. **The oracle with round-specific constants** (adds to "The oracle"). T4 is
   the only configuration in which the true constants differ between rounds. The
   oracle keeps one report per run — the constraint `extract` itself has — and
   carries the longest constant that is true in **every** provided round, which
   for T4 is the 3′ constant without the round-specific base. The review
   suggested a per-round oracle instead; that would measure something else,
   because the difference between the oracle and the inferred run would then mix
   the cost of inference with a freedom `selexprep` does not have.
3. **The dominant clone is a share, not a substitution rate** (replaces
   configuration C in the table, and follows `configurations.tsv`, which already
   said "makes up"). The level `p` is the clone's share of the round: its count
   is set to exactly ⌈p·M⌉ by replacing a random subset of the reads that are
   not already that read, so reads that already equal the clone count toward
   `p`. When the round already holds more copies than that, a random subset of
   the surplus copies is replaced by reads drawn, with the seed and with
   replacement, from the round's reads that are not the clone. C100 leaves only
   the clone. Recorded per run: the realised share and the
   pool's own top-read share before the replacement.
4. **The K spacer is fixed** (replaces "separated from the constant by 4 random
   nt" in configuration K). The 4 nt are drawn once per cell and seed and are the
   same in every read and round, so the outer region is a fixed 16 nt. A spacer
   drawn per read would break the consensus 12 nt into the read, below
   `DEFAULT_MIN_LEN` (14), and the case would never reach the repeated-core
   anchoring it exists to exercise.
5. **Indel coordinates** (adds to "Construction order"). An insertion is placed
   before a sampled position and belongs to the block of the position it
   precedes; a deletion removes the sampled position from its block. Both
   conventions are arbitrary at the two boundaries of the random region, and
   so is any placement inside a run of identical bases, where the same indel at
   another position gives the same read. A read is therefore flagged in its
   name, and reported separately, when any placement that gives the same read
   (the indel shifted left and right as far as it stays equivalent, as in
   variant normalisation) touches a boundary of the random region: a convention
   must not move a headline number.

**Screen.**

6. **Residue probes of every length ≥ 8** (replaces "an exact fragment of at
   least 8 nt … at the read edge" in step 4 with its intended meaning). A single
   8-mer tested with `startswith` / `endswith` misses longer residues, which
   shift it away from the edge: `N + AGATCGGAA` does not end with `AGATCGGA`.
   Every length from 8 to the full probe is tested, anchored at the edge, the
   longest match is recorded, and a read counts once per signal. Signals are
   counted per probe family: TruSeq R1 and R2 share their first 13 nt, so they
   are one signal; counted apart, one residue would be split between two labels
   and each part could stay under the 5% trigger.
7. **Which side was checked is recorded** (replaces the single `none` value).
   Constants are loaded per side, so one known side is still used, and each side
   is recorded with the probe and its length, as `unknown`, or as `too short`. A
   documented fragment shorter than 8 nt is still probed when it is at least 6
   nt, where the background expectation (1/4096 ≈ 0.02%) is fifty times below
   the 5% trigger; below 6 nt it is recorded as too short. `ground_truth.tsv`
   covers no candidate, so the candidates' documented fragments are recorded,
   with their source, in `candidate_constants.tsv`.
8. **The documented runs are re-verified** (adds to Amendment 1, rule 1).
   Substituting the runs the documentation names skipped the checks the parser's
   runs had passed. Each chosen run is now re-checked against the archive
   (FASTQ present, mean length within N ± 2, at least 10,000 reads) and the
   result is recorded. Step 7 is split: `step7_biosample` (no BioSample shared
   with a development donor) and `step7_selection` (no shared selection or input
   library), which a shared BioSample does not cover.
9. **A stated length inside a longer molecule is a stated length** (replaces the
   integer rule of step 1). Taking the single integer of the free-text field
   excluded PRJEB9897 ("12-mer random sequence (within a 71-bp dsDNA)"), which
   states N = 12: an interpretation error, not an absent N. All 71 deposits
   excluded as "no single stated N" are re-reviewed by hand into
   `n_random_review.tsv` (accession, text, reviewed N or none, reason) and the
   screen reads that table. Metadata only; no read of theirs has been examined.

**Consequences for the current candidates.**

- **PRJEB14550 is inconclusive** on its trajectory, by the same standard applied
  to the other two: the publication documents the FLI1 cycle of 143,389 reads
  (the run suffixed 2) and the previous cycle, which supports suffixes 1 and 2,
  not the suffix-4 run chosen as the latest round. Restricted to what is
  documented it has two rounds, below step 6. Recorded, not used.
- **Its reads were already examined**, which the Amendment 1 disclosure did not
  say: PRJEB14550 is in the all-read check of commit f9d8bd7
  (`../read_state_full.tsv`, ERR2270887: 186,388 reads, all 40 nt, no constant
  at either edge). PRJEB28411, excluded at step 6 for having two rounds, is in
  that check too. "Before any read of a candidate donor was examined" holds for
  the four candidates still under evaluation; for PRJEB14550 and PRJEB28411 it
  does not, and both are out for reasons that have nothing to do with their
  reads.
- **The re-review of item 9 changed no candidate.** Five deposits do state one
  randomised length after all — PRJEB9897 (12), PRJNA491840 (4), PRJNA854957
  (45), PRJNA1091700 (6), PRJNA1246497 (40) — and none of them has a single run
  whose mean read length is within N ± 2, so all five stop at step 2. What
  changes is the record: they are now excluded by measured read lengths instead
  of by a misreading of their metadata.
- A candidate that passes step 2 without a row in `trajectory_documentation.tsv`
  is now **pending**, not decided by the parser, as Amendment 1 rule 1 requires.
- The four remaining candidates all pass the re-verification of item 8
  (PRJEB25907 mean 40.00; PRJEB51212 and PRJEB51473 mean 21.00; PRJNA741127
  mean 16.00; every run above 170,000 reads). **PRJEB51212 and PRJEB51473 are
  expected to fail step 3**: 21.00 nt for a 20-nt random region is one extra
  base in every read. The all-read check measures it; if the modal length is 21
  they fail step 3 and are not used. What the extra base is (a residue, or a
  stated length that is wrong) is recorded from the evidence, not assumed.

**Sampling replicates** (replaces "10,000 reads per round are sampled … (seeded)"
in "Which donor reads are used"). One pool per run **and seed** — 1–3 for
development, 101–103 for test — each with its own hash, so a seed varies the
donor sampling as well as the simulated error, as "Replication" states. A single
pool would have made the seeds vary only the error model.

**Status of the read job's output.** `donor_checks.tsv` and the pools it writes
are measurements and provisional material: admission (steps 3–5, the trajectory,
the cap of six and the priority order) is a separate decision, recorded in
`donor_screening.tsv` with its evidence. The job ranks nothing and excludes
nothing by itself, and conservation, diversity and convergence stay descriptors.

**Disclosure.** Every item of this amendment was written after reading the
design, the screen's code and its tables, after the metadata of every candidate,
and after the reads of PRJEB14550 and PRJEB28411 as recorded above; before any
read of the four candidates still under evaluation and before any simulated
library.

### Amendment 3 — 2026-09-28, donor screen results, before any simulated library

**Read job.** `check_and_sample_donors.py` at the Amendment 2 commit, SLURM job
121232: 27 runs (15 of the current donors, 12 of the four candidates), every
read of each. Results in `donor_checks.tsv`.

**Correction made after the job, for one donor.** PRJNA741127's constants had
been recorded as undocumented; the GEO sample records do document them, in the
depositor's trimming command (test samples: `cutadapt -g ^CGC -a
CCTGGAATTCTCGGGTGCCAAGGAACTCCAG$ … -m 16 -M 16`). `candidate_constants.tsv` was
corrected, and a defect found while doing so was fixed: a documented constant
longer than 8 nt had been probed only at its full length instead of from 8 nt.
PRJNA741127 alone was checked again with the corrected probes (no residue
signal); its rows in `donor_checks.tsv` come from that rerun, and its pools are
bit-identical to the job's (same SHA-256), so only the probe columns changed. No
other donor had a documented constant longer than 8 nt taken from
`candidate_constants.tsv`.

**Current donors.** All five pass the mechanical part of step 3 in every chosen
round (modal length N, at least 99.6% of reads within N ± 2) and show no residue
signal. Retained reads, excluded fraction and the pool hashes are now in
`donors.tsv`.

**Candidates** (`donor_screening.tsv`, with the evidence):

- **PRJEB25907 — admitted** (RNA, N40). Every read 40 nt; no residue signal
  (adapters, and the 5′ fragment TCCCAA); at least 171,491 reads per round.
- **PRJNA741127 — admitted** (DNA, N16). Every read 16 nt, as the depositor's
  documented trimming requires; no residue signal (adapters, and the documented
  3′ constant); at least 2.69 million reads per round.
- **PRJEB51212 — inconclusive.** Every read 21 nt for a stated N of 20, in all
  three rounds; no edge position above 0.44 support, so the extra base is not a
  fixed residue, and what it is cannot be told from the reads or the archive.
- **PRJEB51473 — excluded.** Every read 21 nt for a stated N of 20; a position
  among the last ten carries one base in at least 97.7% of the reads of every
  round, the cycle-0 input library included: a fixed base, not selection.

**Test set.** Five test donors: PRJEB70964, PRJEB47428, PRJEB49150, PRJEB25907,
PRJNA741127. Below the cap of six, so the priority order was not needed. As
fixed in "Donor screen", test uses seeds 101–102 only: 300 runs per test donor,
1,500 in all; the seed-103 pools are kept and not used. Declared limit:
PRJEB25907 and PRJEB47428 are both HTR-SELEX RNA with N40, so they add depth on
one assay rather than a new one; PRJNA741127 adds a new N length (16).

### Amendment 4 — 2026-09-28, implementation choices, before any development run

The generator (`simulate.py`), the evaluator (`evaluate.py`) and the harness
(`run_synthetic.py`, `run_synthetic.sh`) are written, with unit tests in
`tests/test_synthetic_benchmark.py`. They were exercised only on random fake
pools; no run has used a donor pool. Where the text above left a detail open,
the code fixes it as follows.

1. **A correct call** is a suffix (5′) or prefix (3′) of the error-free,
   untruncated sequence before or after N in **every** read of the earliest
   round provided, in read orientation. Where that sequence varies between
   reads (T3's leading base), a call that includes the varying part is wrong.
2. **Truncation (F)** cuts at a uniform position that keeps between 1 nt and
   all but 1 nt of the 3′ constant.
3. **The motif (M, X3)** is drawn from the seed alone, so it is the same in
   every configuration run with that seed; the K spacer, the tags (T1, T2) and
   the T4 base are drawn per cell and seed.
4. **N±1 (B)** deletes a base at a uniform position of N, or inserts a uniform
   random base at a uniform position.
5. **Round numbers** passed to `selexprep` are the donor's own (for example 1,
   4, 7), and rounds are passed earliest first.
6. **`selexprep` runs through its command line**, as a user runs it; a detect or
   extract call that exits non-zero or exceeds 30 minutes is an error. The
   oracle report is written with the public schema (`LibraryReport`) and
   extracted with the same command.
7. **Per-run output** is one row of `results.tsv`. The inputs and the extracted
   reads are deleted after evaluation and can be regenerated from the seed; the
   truth, the reports and the logs are kept.
8. **A run whose report says REVERSE** is evaluated in the frame `extract`
   works in, after reverse-complementing the reads.
