# Round-assignment validation — protocol (pre-registration)

This file fixes the question, the populations, the admissible truth and the
outcomes **before** any truth is collected or compared. It is committed first;
the truth table, the comparison script and the results come after it, each in
its own commit. Changes after this commit are recorded, dated, under
[Amendments](#amendments), with the reason.

## Question

When `selexprep fetch` assigns a round to a run, is it the round the depositors
meant? And how often does it assign none?

## Code under test

- **Before:** the round parser of v0.4.1 (unchanged up to 05c21cc; the only
  change to `fetch/plan.py` in between is a warning that does not touch rounds).
- **After:** the parser of fd1a28d, which reads ENA sample attributes, weighs
  every round attribute before deciding, and no longer reads `RV<n>` as a round.
- The final numbers are rerun on the 0.4.5 tag and reported with its commit.

## Disclosures

- **In-sample.** The parser's patterns and its attribute rule were developed on
  this catalogue. The attribute-key vocabulary was widened on 2026-09-28 after a
  census of the attribute keys in these same deposits. Nothing here estimates
  performance on deposits outside the catalogue.
- The before/after catalogue snapshots of 2026-09-28 (round per run, no truth)
  were seen before this protocol was written.

## Populations (three views, reported separately)

| View | Records | Measures | Purpose |
|---|---|---|---|
| **A. Universe** | all 127 INSDC deposits of the bundled catalogue, every run | coverage only | how often a round is assigned at all |
| **B. Fixed before/after set** | the 28 records classified `ELIGIBLE_HT_SELEX_ROUNDS` in the committed audit (`benchmarks/audit_results/eligibility.tsv`, cae12f3), every run | coverage and accuracy, before and after, on the same runs | the effect of the fix |
| **C. New eligible set** | records classified `ELIGIBLE_HT_SELEX_ROUNDS` by the eligibility module with the code after the fix | coverage and accuracy | what the S2 audit will stand on |

Eligibility is computed from the parser's own rounds, so B and C are not a
sample of the parser's failures; A is the only view that sees every run.

**Prediction, to be checked, not assumed:** PRJEB51212 and PRJEB51473 leave the
eligible set (their attributes give cycles 0/1/3/6 across several constructs,
which the eligibility module should classify as a mixed project).

## Truth

For each run of B ∪ C, the round the depositors meant, taken from a source whose
**provenance** differs from the field the parser read for that run. Provenance
is where the value originates, not the site that shows it.

Admissible sources:

1. The publication: text, figures, methods, supplementary tables.
2. An archive record the parser did not read for that run, and that is not
   derived from what it read.
3. A curated round map in `benchmarks/round_maps/`, only where its own
   provenance is documented and differs from the parser's field.

Not independent (examples): an ArrayExpress SDRF against the ENA `selex cycle`
attribute mirrored from it; a GEO sample page repeating the BioSample `round`
attribute; a curated map built from the same attribute or title.

Every truth row records: accession, run, true round (or a category below), the
source (DOI or URL), the exact location (table, figure, field), the provenance
of the value, and whether it is independent of the field the parser read.
A row whose provenance coincides with the parser's field is recorded as
**archive-consistent** and kept out of the accuracy counts.

**Truth categories.** A run is either a numbered round in the depositors'
numbering (an input library counts as round 0 only where the depositors number
it so), or not a selection round (a no-target control, an unrelated library, a
run of another experiment). Assigning a number to a run that is not a selection
round is wrong.

**Numbering.** Correct means the depositors' number. A number that differs from
it by the same offset across a whole trajectory is wrong, and is also counted
as a consistent offset, reported apart.

## Outcomes

**Coverage** (every run of each view): assigned or unassigned. A run is
unassigned when `RoundRecord.is_unassigned` is true.

**Accuracy** (runs with an independent truth only): correct or wrong. An
assigned run without an independent truth is **not verified**; it is counted,
never folded into correct.

**Per record:**

- records with at least one wrong run;
- records with at least one unassigned run (the two can coexist);
- records **all correct**: every run verifiable and correct.

Reported for B before and after, and for C after, as run counts and record
counts, with the denominators stated.

## Order of work

1. This protocol is committed.
2. The truth table (`truth.tsv`) is compiled from the sources without looking
   at the parser's rounds for those runs, and committed.
3. The comparison script is written with tests and committed.
4. The comparison runs, first on the code after the fix, then on the 0.4.5 tag;
   the results are committed with the commits they ran at.

## Reporting

A separate section of Supplementary Data S2, "Round-assignment validation",
keeping the statement that the S2 audit does not test primer correctness. The
paper states the views, the denominators, the in-sample disclosure and the
number of runs not verified.

## Amendments

### Amendment 1 — 2026-09-30, before any truth was compiled or compared

**What changes.** The endpoint "accuracy against a source of independent
provenance" on views B and C is **suspended**. No truth table is compiled and
no accuracy is reported. What is kept is descriptive: coverage, the before/after
comparison of the parser on the catalogue, the agreement with the curated round
maps on the runs they can be compared on, and the specific cases whose
correction is documented. No assignment without verification is counted as
correct.

**Why.** A feasibility review of view B, done before any truth was collected,
showed that an independent reference is hard to find systematically:

- B has 28 records and 679 runs.
- In the 9 DDBJ records the round appears only in the library name, which the
  parser reads; their sample description names the target, not the round.
- In 14 of the other 19 records the sample description repeats the sample
  title or the library name in every run, or is empty; it is a copy of what the
  parser reads, not a second source.
- 22 of the 28 records have no DOI or PMID in the catalogue metadata. This does
  not show that no publication exists; it shows that finding one, and a
  run-level mapping in it, would be a search per record.

The accuracy that could be measured would therefore cover few records, and most
runs would be "not verified". An interpretation reference (a person reading the
depositors' label without the parser's output) was considered: it would measure
a different endpoint, the interpretation of the labels rather than their
correctness, and the difficult cases are already known from development, so it
would not be a blind validation. It is not adopted now.

**What was known when this was written.** The before/after catalogue snapshots
of 2026-09-28 and 2026-09-30 (round per run, no truth); the development cases
(`RV<n>` sample numbers in PRJEB51212, PRJEB51473 and PRJEB38961; titles with a
protein name or replicate suffix read as round 1 in PRJEB76622 and PRJEB61115);
the raw metadata of the 679 runs of B, extracted without the parser's rounds;
and a classification of the catalogue with the code after the fix, run for view
C and not examined.

**What is reported instead**, for the released version (0.4.5), as in-sample
and descriptive, with denominators:

- coverage on every run of the 127 INSDC deposits: runs that receive a round
  and runs left unassigned, without treating an assignment as correct;
- agreement with the curated maps in `benchmarks/round_maps/`, counting only the
  runs that receive a round, with the unassigned runs stated;
- that every assignment keeps its source field, confidence and notes, and that a
  user-supplied round map overrides it.

The paper describes the released version, not how it was reached: the
before/after comparison of the fix stays in `CHANGELOG.md`. The numbers are
produced by a script run on the 0.4.5 tag together with Supplementary Data S2.
None of this is an estimate of accuracy on deposits outside the catalogue.
