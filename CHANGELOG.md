# Changelog

All notable changes to `selexprep` are documented here.

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Versioning: [SemVer](https://semver.org/).

## [Unreleased]

### Fixed

**Phase 6a CI hotfix (2026-05-22)** — `test_run_help_lists_resume_and_stop_on_error`
passed locally on macOS but failed on all three Python versions on
Ubuntu CI. Root cause: Typer/Rich renders the help table at the
terminal's column width, and CI's non-TTY default of 80 columns is
narrow enough that the Rich renderer silently truncates `--resume`,
`--stop-on-error`, and `--backend` option names out of the rendered
text. macOS local runs detect a wider effective width (Rich + CliRunner
interaction) and show all options. The test grepped the rendered help
text, so it was width-dependent by construction.

Fix: rewrite the test to introspect the Click command's `params` list
(`flag_names = {opt for p in run_cmd.params for opt in p.opts}`)
instead of grepping rendered output. Renamed to
`test_run_command_registers_resume_and_stop_on_error_options`. Kept
a smoke that `--help` exits cleanly so the rendering path isn't
unobserved. The contract being tested ("these options exist") is
captured exactly; the rendered help width remains
environment-dependent (matplotlib-style informational, per the same
discipline as Phase 5 plots).

**Phase 6a Codex pass 1 (2026-05-22)** — two blocking semantic bugs +
two non-blocking polish items. Calibration verdict: N/A (Phase 6a
introduces no calibration constants).

Blocking:

- **Outer except in `run_batch` mislabeled every unexpected exception
  as `EXTRACT_FAILED` (run/runner.py).** If `build_fetch_plan` or
  `download_srr` raised an `HTTPError`, `RequestException`, or
  `ValueError`, the exception escaped the fetch block to the outer
  safety net and `run_summary.tsv` recorded `status=EXTRACT_FAILED`
  even though the extract stage never ran. Fix: stage-classify
  expected fetch-stage exceptions inside `_process_one_accession`
  (new `try/except` around the `run_fetch` call → `FETCH_FAILED`),
  and rename the outer-net status to a new literal
  `UNEXPECTED_FAILURE` so a truly-escaped exception is no longer
  misreported as extract. The `RunStatus` literal gains
  `UNEXPECTED_FAILURE`; nothing produces `EXTRACT_FAILED` from the
  outer net any more. Resumed reads of `fetch_metadata.json` now also
  guard `read_fetch_metadata_json` against malformed audit-trail
  files (`ValueError` / `KeyError` / `OSError`) and downgrade to
  `FETCH_FAILED` instead of leaking.
- **Fetch resume oracle accepted corrupt `.fastq.gz` files
  (fetch/runner.py).** Both the per-run skip-already-present check
  (`run_fetch`) and `check_fetch_inventory` used bare `Path.exists()`.
  A `SIGKILL`'d download leaves a truncated `.fastq.gz` on disk that
  exists but cannot be decompressed; `--resume` would skip it and
  downstream `detect`/`extract` would hit gzip errors deep in the
  pipeline. Fix: both call sites now use
  `selexprep.fetch.download.validate_fastq_gz` (existing helper at
  fetch/download.py:96 — checks exists, ≥1024 bytes, and `gzip -t`
  decompresses cleanly). Mirrors the same check
  `download_srr_ena_direct` performs at fetch/download.py:288.

Non-blocking:

- **`_count_yaml_flags` used `line.startswith("- name:")` which
  miscounts under `safe_dump(sort_keys=True)` (run/runner.py).** When
  a flag has non-empty `evidence`, the YAML emission orders dict keys
  alphabetically inside each list entry (`evidence` < `name` <
  `severity`), so the first line of each entry can be
  `- evidence:` rather than `- name:`. Resumed `run_summary.tsv`
  would undercount flags. Fix: parse the YAML via `yaml.safe_load`
  and count list entries that contain a `name` key.
- **`fastq_filenames_for_run` synthesized
  `{SRR}.fastq.gz`/`{SRR}_{1,2}.fastq.gz` instead of deriving from
  URLs (fetch/plan.py).** Right for ENA's canonical naming, but
  brittle if a future ENA URL convention or Zenodo mirror used a
  different scheme. `download_srr_ena_direct` writes
  `output_dir / Path(url_path).name` (fetch/download.py:351), so the
  URLs are the source of truth. Fix: derive basenames from
  `run.fastq_urls`; keep the SRR-based synthesis as a defensive
  fallback for the rare case where ENA returns no `fastq_ftp` (which
  also means the run isn't downloadable, so the fallback path is
  effectively dead).

Six new regression tests:

- `test_run_batch_fetch_http_error_records_FETCH_FAILED_not_extract_failed`
  — `requests.HTTPError` from `build_fetch_plan` ⇒ `FETCH_FAILED`
  (not `EXTRACT_FAILED`).
- `test_run_batch_corrupt_fastq_is_redownloaded_by_resume_oracle` —
  garbled `.fastq.gz` on disk triggers re-download under `--resume`;
  the result passes `validate_fastq_gz`.
- `test_run_batch_count_yaml_flags_uses_yaml_parser` — flags YAML
  with `evidence`-first entry ordering counts to 2, not 0.
- `test_fastq_filenames_for_run_derives_from_urls` — hypothetical
  alternate URL naming (`foo_subdir/SRR_X_unusual.fastq.gz`) is
  honored verbatim.
- `test_fastq_filenames_for_run_fallback_when_urls_empty` — defensive
  synthesis path still works when `fastq_urls=[]`.
- `test_run_batch_unexpected_error_uses_unexpected_failure_status` —
  outer-net safety catch uses `UNEXPECTED_FAILURE`, not
  `EXTRACT_FAILED`.

Test fixtures updated: `_stub_download_writes_files` and `_write_fastq`
now emit 800 records of 60-base pseudo-random sequences so the gzipped
outputs exceed `validate_fastq_gz`'s 1024-byte floor. The previous
small fixtures (~38 bytes raw) would now be rejected as too-small to
be valid downloads.

CI: 432 passed + 1 xfailed (was 426; +6 new tests). All four
pre-existing CI gates stay green.

### Added

**Phase 6a — close v0.1 CLI surface: `selexprep fetch` + `selexprep run`
(2026-05-22)**

Wires the two remaining `_not_implemented` CLI verbs from the locked
plan (lines 181–189), closing the v0.1 single-dataset + batch CLI
surface. With Phase 6a landed, the public CLI exercises end-to-end
(`inspect → fetch → detect → extract → count → qc`); Phase 6b
(benchmark Snakefile + Figure A) can now drive the comparison against
AptaPLEX + EasyDIVER+ through the public CLI rather than through
library imports.

- **NEW: `selexprep.fetch.plan`** — `build_fetch_plan(accession,
  timeout_s) → FetchPlan` hits ENA filereport with extended fields
  (`sample_title`, `library_name`, `experiment_title`,
  `sample_accession` on top of the inspect set) so the 5-level
  cascade in `fetch.metadata.parse_round` runs on real metadata
  instead of bare-bones inspect fields. Emits a deterministic
  `fetch_metadata.json` audit trail. Single source of truth for
  "what does this accession look like and how do its runs map to
  rounds?", reused by both `fetch` and `run` CLIs.
- **NEW: `selexprep.fetch.runner`** — `run_fetch(accession, outdir,
  *, backend, allow_manual_review, dry_run, timeout_s) → FetchResult`
  orchestrator. Emits `rounds.tsv` (HIGH/MEDIUM-confidence
  contract that downstream `detect`/`extract` consume; sorted by
  filename for determinism) + per-round FASTQs under `round_NN/` +
  optional `manual_review.tsv` + `fetch_metadata.json`. Cardinal
  rule (per `fetch/metadata.py` line 14): never guess a round
  assignment — NONE-confidence runs refused by default; opt-in via
  `--allow-manual-review` routes them to `round_unknown/` + a
  separate `manual_review.tsv` and **never** lets them enter
  `rounds.tsv`. If every run is NONE-confidence, refuses fail-fast
  before any download. Shared `query_ena_filereport()` helper
  factored out of `inspect.py` so the HTTP call shape is defined
  once.
- **NEW: `selexprep.run.runner`** — `run_batch(accessions_tsv,
  outdir, *, resume, stop_on_error, backend, allow_manual_review,
  timeout_s) → RunReport` batch driver. Per-accession pipeline
  (fetch → detect → extract → count → qc) with file-inventory
  resume oracles per stage (not sentinel-flag based):
  - **fetch**: `fetch_metadata.json` present AND every expected
    FASTQ from `FetchPlan.runs` on disk AND `rounds.tsv` present.
    Mirrors `download_bioproject._missing_srrs` discipline at
    `fetch/download.py:683`.
  - **detect**: `library_report.json` present.
  - **extract**: `selexprep_manifest.json` present.
  - **count**: per-round `round_NN/counts.parquet` (granular —
    only missing rounds re-run).
  - **qc**: `qc/flags.yaml` present.
  Paired-end FASTQs are grouped by ENA naming convention
  (`<srr>_1.fastq.gz` = R1, `<srr>_2.fastq.gz` = R2); R1-only
  sequences feed `compute_library_report`'s primary stream, R2
  sequences are passed as `paired_mate_streams`, and both are
  threaded through `run_extract` as `paired_r2_inputs`. Per-accession
  errors default to log-and-continue with `status=FAILED_<stage>`
  in `run_summary.tsv`; `--stop-on-error` flips to fail-fast.
- **WIRED: `selexprep fetch <accession> --outdir OUT
  [--backend ena|auto|kingfisher|sra] [--allow-manual-review]
  [--dry-run] [--timeout-s N]`** — replaces the
  `_not_implemented` stub. **CLI default `--backend ena`** (paper-
  grade reproducibility; fail-fast if ENA can't serve; never
  silently fall back to GPL tools); `--backend auto` is the
  explicit opt-in for the convenience fallback chain. Library-level
  `download_srr(backend="auto")` default unchanged (Python API
  ergonomic).
- **WIRED: `selexprep run accessions.tsv --outdir OUT [--resume]
  [--stop-on-error] [--backend ena|auto|...] [--allow-manual-review]
  [--timeout-s N]`** — replaces the `_not_implemented` stub. Emits
  a deterministic `run_summary.tsv` (sorted by accession) with
  per-accession status: `OK` /
  `SKIPPED_READ_MERGING_RECOMMENDED` /
  `FETCH_REFUSED` / `FETCH_FAILED` / `DETECT_FAILED` /
  `EXTRACT_REFUSED` / `EXTRACT_FAILED` / `COUNT_FAILED` /
  `QC_FAILED`. Split-primer guard (`required_action ==
  READ_MERGING_RECOMMENDED`) skips count + qc and records the
  status rather than producing misleading half-insert counts
  (locked plan line 325).

**Tests added (+42, 427 + 1 xfailed total):**

- `tests/test_fetch_plan.py` (12 — mocked-ENA paired-vs-single,
  round confidence propagation, deterministic JSON, frozen
  dataclass, extended-field assertion, empty-response ValueError,
  timeout pass-through).
- `tests/test_fetch_cli.py` (15 — orchestrator happy paths,
  refusal paths (all-NONE, mixed-without-flag, allow-manual-review
  rounds.tsv cleanliness), resume oracle (`check_fetch_inventory`
  + skip-already-present), CLI dry-run, invalid backend,
  refusal exit code, propagated ValueError, all four documented
  backends accepted).
- `tests/test_run_runner.py` (12 — duplicate-accession refusal,
  missing-column refusal, fetch refusal status, paired-end
  R1/R2 threading through detect+extract, split-primer skip
  (no count/qc), fetch resume oracle re-fetches missing FASTQ,
  per-stage resume (no work on re-run with all sentinels),
  log-and-continue default, `--stop-on-error` halts loop,
  detect-stage exception → `DETECT_FAILED`, summary TSV sorted
  by accession, manual-review separation keeps `rounds.tsv`
  clean).
- `tests/test_cli.py` (+3 — fetch dry-run smoke, run missing-
  accession-column rejection, run --help lists `--resume` +
  `--stop-on-error`). Replaced the obsolete
  `test_fetch_stub_exits_with_code_2`.

**Calibration:** N/A — Phase 6a introduces no new heuristic thresholds.
The `CALIBRATION-TODO` inventory stays at 19. The fetch refusal
threshold (NONE-confidence) is the locked-plan cardinal-rule binary,
not a tunable.

### Fixed

**Phase 5 Codex pass 1 (2026-05-21)** — three blocking semantic bugs +
three non-blocking polish items. Calibration verdict: all six v0.1
qc constants ✅ confirmed against published HT-SELEX conventions
(RAPID-SELEX, AptaPLEX, Hoinka 2015, EasyDIVER+, FASTAptamer).

Blocking:

- **`rarefy()` returned the original pool unchanged when `depth >=
  total`, reintroducing the depth-confounding the flag was meant to
  avoid (qc/flags.py).** A round with 2k reads was being compared
  against a round with 100k reads rarefied to 10k — exactly the trap
  the locked plan flagged at line 351. Fix: clamp to
  `effective_depth = min(RAREFACTION_DEPTH, min_total_reads_per_round)`.
  Rarefy all rounds to that common depth. Evidence dict now includes
  both `configured_depth` and `effective_depth` for audit. When a
  round is empty (`min_total == 0`), return None (the `low_total_reads`
  flag handles that separately).
- **`check_adapter_contamination_high` mixed measurement universes
  (qc/flags.py).** Numerator was `LR.known_adapter_hits` (Phase 2
  inference on earliest-round subsampled reads); denominator was
  `sum(counts_by_round)` (post-extraction reads across all rounds).
  Apples-to-oranges → fraction could be inflated or deflated by an
  arbitrary factor. Fix: switched the denominator to
  `trim_reports_by_round[earliest_round]["n_in"] *
  lr.read_fraction_used_for_inference`, matching the Phase 2
  inference universe. When `trim_reports.json` is unavailable, surface
  the flag as `severity="info"` with `evidence={"reason":
  "denominator_unavailable"}` rather than guessing. Threaded a
  `trim_reports_by_round` kwarg through `compute_all_flags`; the qc
  runner already aggregates the JSON via `_load_trim_reports_by_round`
  and now passes it down.
- **`selexprep count` silently parsed FASTQ as FASTA (cli.py).** The
  CLI called `count_fasta` unconditionally; a `.fastq.gz` input would
  be misparsed (every other line treated as a sequence header). Fix:
  hard-reject FASTQ extensions by default with a clear error pointing
  the user at `selexprep extract` first OR at the new
  `--from-pretrimmed-fastq` opt-in flag. The v0.1 contract is "`count`
  accepts only extracted FASTA from `extract`"; users with externally
  pre-trimmed FASTQ (e.g., from AptaPLEX, EasyDIVER+, external
  cutadapt) can opt in via `--from-pretrimmed-fastq`, which routes
  through a new `count_fastq_pretrimmed()` library function and emits
  a loud "cannot verify trimming state" `logger.warning`. The legacy
  Phase-1 `count_round` (inline cutadapt trimming) stays a library-
  only API.

Non-blocking:

- **`flags.yaml` was promised to contribute to `output_sha256` but
  didn't (qc/plots.py + manifest.py).** Phase 4's manifest is sealed
  by `extract` before `qc` runs, so the QC artifacts are outside the
  extract manifest lifecycle by construction. Added `.yaml` to
  `_HASHABLE_SUFFIXES` (so future `selexprep qc-amend` could append
  the hash) and rewrote the qc/plots.py docstring to clarify the
  current lifecycle and the optional hash path.
- **`_parse_round_label` accepted negative integers (cli.py).** `R-1`
  produced `round_-1/`. Fix: explicit `r < 0` rejection with
  `typer.BadParameter`.
- **YAML emission lacked a regression test for nested-dict evidence
  (tests/test_flags.py).** `safe_dump(sort_keys=True)` does sort
  nested dict keys recursively, but the determinism guarantee for
  nested LIST-of-DICT evidence (as in
  `adapter_contamination_high`'s `adapters_above_threshold` field)
  wasn't exercised. Added
  `test_write_flags_yaml_deterministic_with_nested_evidence` —
  builds the same flag twice with intentionally-different insertion
  orders and asserts byte-identical output.

### Carry-forward to v0.2 (documented now, not yet implemented)

- **`check_adapter_contamination_high` assumes `min(trim_reports_by_round)`
  is the same round Phase 2 used for primer inference.** Holds for the
  documented v0.1 CLI flow (same FASTQ set passed to detect + extract).
  Could break for partial extract runs, manually-stitched
  ``trim_reports.json``, or the v0.2 batch driver. v0.2 fix: add
  ``LibraryReport.earliest_inference_round`` (set in
  ``compute_library_report``) and key the denominator off that. Schema
  bump is disproportionate for v0.1 since the bug is unreachable through
  the documented CLI flow. Inline comment in
  ``qc/flags.py:check_adapter_contamination_high`` flags the assumption
  for the v0.2 implementor.
- **`_iter_fastq_sequences` only validates record completeness, not
  per-line format conformance** (no ``@`` header check, no ``+``
  separator check, no ``len(seq) == len(qual)``). Acceptable for the
  ``--from-pretrimmed-fastq`` power-user opt-in; if this path is
  promoted to a first-class advertised workflow in v0.2, add the
  validators here matching the truncation ValueError policy.

### Added

**`selexprep.count.counter.count_fastq_pretrimmed`** — library function
for the `--from-pretrimmed-fastq` CLI opt-in. Parses FASTQ records,
extracts the sequence line, builds a Counter, emits the standard
counts.parquet schema (`sequence`, `reads`, `rank`, `rpm`). Same
truncation policy as `extract/strand.py:reorient_fastq_gz` — raises
`ValueError` on incomplete records rather than silently producing
partial output. Used for the external-tool interop story (AptaPLEX +
EasyDIVER+ produce pre-trimmed FASTQ that doesn't need selexprep's
extraction).

Ten new regression tests:
- `test_diversity_increase_clamps_depth_to_min_total`
- `test_diversity_increase_records_effective_depth_in_evidence`
- `test_adapter_contamination_denominator_unavailable`
- `test_adapter_contamination_honors_read_fraction_used_for_inference`
- `test_count_rejects_negative_round_label`
- `test_count_rejects_fastq_input` (default reject path)
- `test_count_with_from_pretrimmed_fastq_succeeds`
- `test_count_from_pretrimmed_fastq_logs_unverified_warning`
- `test_count_from_pretrimmed_fastq_flag_rejects_fasta_input`
- `test_write_flags_yaml_deterministic_with_nested_evidence`

Plus three library-level regressions for `count_fastq_pretrimmed`:
- `test_count_fastq_pretrimmed_basic`
- `test_count_fastq_pretrimmed_raises_on_truncated_record`
- `test_count_fastq_pretrimmed_handles_uncompressed`

Two existing tests updated to reflect the new
`check_adapter_contamination_high` signature (now takes
`trim_reports_by_round`, not `counts_by_round`).

CI: 385 passed + 1 xfailed (was 372; +13 new tests). All four
pre-existing CI gates stay green.

**Phase 4 Codex pass 1 (2026-05-21)** — one blocking semantic bug +
three non-blocking polish items. Calibration verdict: N/A (Phase 4
introduces no calibration constants).

Blocking:

- **`--override-primer-{5p,3p}` was applied AFTER the refusal check
  (runner.py).** Locked plan line 311 explicitly says
  `--override-primer-*` should bypass the UNABLE_TO_INFER /
  UNABLE_TO_EXTRACT refusal — but the runner was checking refusal
  first, so an override could never reach the LR. Worse, the refusal
  message suggested using override as the fix, creating a circular
  error path. Refactored to:
  1. Apply override first (clone LR via `model_copy(update=...)`).
  2. When baseline `extraction_mode == "UNABLE_TO_EXTRACT"`, infer a
     new mode from which primer sides are now defined: both →
     `BOTH_PRIMERS_SINGLE_READ` + `full_insert_recovered=True` +
     `required_action=NONE`; 5p-only → `FIVE_PRIME_ONLY`; 3p-only →
     `THREE_PRIME_ONLY` (`full_insert_recovered=False`).
  3. Promote `status` from `UNABLE_TO_INFER` → `MEDIUM` (manual
     override = medium confidence by convention; hand-editing the LR
     JSON is the path to claim HIGH). Clear `failure_reason`.
  4. THEN run the refusal check. With override applied, UNABLE state
     is naturally cleared and the dispatcher proceeds. Refusal still
     fires if no override was given AND baseline is UNABLE.
  Refusal message also reworded to clarify the override mechanism
  (dropped the stale `(Phase 4)` parenthetical).

Non-blocking:

- **Override didn't check `KNOWN_ADAPTERS` (runner.py).** An override
  matching a known sequencing adapter prefix is now allowed (explicit
  escape hatch) but surfaces a `logger.warning` — the adapter-trap
  story still holds for the auto-inference path, and the override
  path gets a foot-gun guard without blocking deliberate use.
- **`HTTPError` in `inspect_accession` propagates without
  differentiating 404 vs 503 (fetch/inspect.py).** OK for the current
  CLI flow; flagged for the future `selexprep run` batch driver.
  Carry-forward only — not fixed in this round.
- **`parameters` / `runtime_seconds_per_stage` not normalized in
  deterministic JSON (manifest.py).** Currently deterministic because
  the CLI always builds the dict in the same order, but fragile if the
  library API gets a non-CLI caller. Carry-forward.

### Changed

**Architecture hygiene (Codex pass 1 NB follow-up).** Migrated the
`_matches_known_adapter` helper from `library/detect.py` (where it was
underscore-private) to `library/adapters.py` (where `KNOWN_ADAPTERS`,
`KNOWN_ADAPTERS_RC`, and `count_adapter_hits` already live). Renamed
to public `matches_known_adapter_prefix()` — the `_prefix` suffix
makes the semantics explicit (it's a first-`k`-bp match, not arbitrary
substring). Also moved the `ADAPTER_PROBE_K = 13` constant to
`library/adapters.py` so the shared default of `count_adapter_hits`
and `matches_known_adapter_prefix` has a single source of truth.
Cross-package private imports eliminated. No call-site behavior
change; the function is byte-identical.

### Tests (Phase 4 follow-up)

Five new regression tests in `tests/test_extract_override.py`:
- `test_override_with_both_primers_promotes_unable_to_extract`
- `test_override_5p_only_promotes_to_five_prime_only`
- `test_override_3p_only_promotes_to_three_prime_only` (symmetry
  mirror of the 5p-only test)
- `test_override_logs_warning_when_matches_known_adapter`
- `test_no_override_on_unable_still_refuses_with_clear_message`

CI: 372 passed + 1 xfailed (was 367; +5 new tests). All four
pre-existing CI gates stay green. `_make_library_report` in
`test_extract_override.py` extended with kwargs to support multiple
test scenarios (matches the helper pattern in `test_extract_runner.py`).

**Phase 3 Codex pass 1 (2026-05-21)** — three blocking bugs +
three non-blocking polish items, all caught in the second-pass code
review of the extraction pipeline. Calibration verdict for the one
Phase 3 constant: `STRAND_REPORT_PER_READ = False` ✅ CONFIRMED.

Blocking:

- **Multi-FASTQ-per-round overwrite (runner.py).** Each per-mode
  trim loop wrote every input FASTQ in a round to the SAME target
  path — second iteration overwrote first, dropping all but the last
  input's reads. Added `_trim_round_single_end` + `_trim_round_paired_split`
  helpers that trim each input to a per-input temp `.part_NNN.fasta.gz`,
  then concatenate deterministically via `open_gzip_text_deterministic`
  into the final per-round target. Single-input case stays on the
  fast path (no temps).
- **Sample-sheet paired-end demux didn't rebuild `paired_r2_inputs`
  (runner.py).** After demux, only R1-shaped files were collected;
  `paired_r2_inputs` stayed at whatever the caller passed (typically
  None). `PAIRED_END_SPLIT_PRIMERS + --sample-sheet` would refuse with
  "requires --paired-r2". Now collects both `_1.fastq.gz` (R1) and
  `_2.fastq.gz` (R2) from the demuxed layout, builds `round_inputs`
  path-aware from parent dirs (basenames collide across rounds after
  demux), and populates `paired_r2_inputs` when R2 files exist.
- **`output_sha256` collapsed per-round files with same basename
  (manifest.py).** `compute_sha256s` used `p.name` as the dict key, so
  `round_00/extracted.fasta.gz` and `round_01/extracted.fasta.gz`
  collided — last-write-wins. Added `root` kwarg; when given, keys are
  `path.relative_to(root).as_posix()`. Builder passes `output_root=outdir`
  so per-round outputs stay distinct in the manifest.

Non-blocking:

- **Orphaned temp files on trim failure (trim.py).** All four trim
  entry points now wrap cutadapt + repack in `try/finally` that unlinks
  the intermediate `.fa` and `.cutadapt.json` files, even on cutadapt
  raise or repack raise.
- **Truncated FASTQ silently produced partial output (strand.py).**
  `reorient_fastq_gz` now raises `ValueError` on a truncated record
  instead of logging a warning and breaking the loop (the caller has
  no way to distinguish "graceful early termination" from
  "incomplete extraction"). Aligns with the project's no-silent-miscalls
  discipline.
- **R2 basename collision (cli.py).** `extract` now refuses early
  with an explicit list of duplicate basenames when `--fastq` +
  `--paired-r2` collectively contain repeated names (which would
  collapse in the basename-keyed round-map lookup).

Follow-up to the manifest pass (Codex spotted in a second review pass
of the same commit, folded in here):

- **`input_sha256` basename-collision in sample-sheet mode (manifest.py +
  runner.py).** The `root` kwarg now applies to inputs too:
  `build_manifest_from_extract_result` accepts `input_root=...`, and
  the runner passes `input_root = outdir / "demux"` when sample-sheet
  mode is active (demuxed files share basenames across rounds —
  `srr_1.fastq.gz` in every `round_NN/` folder). Normal CLI flow keeps
  `input_root=None` because the CLI's basename-collision guard
  guarantees uniqueness upfront.
- **Latent: `input_sha256` was silently empty for ALL CLI flows
  (manifest.py).** `_HASHABLE_SUFFIXES` only included
  `{.fasta, .fa, .tsv, .json}`, so FASTQ inputs got filtered out and
  the manifest's input-side hashes were missing entirely — contradicting
  the locked plan line 168 ("input_sha256 (FASTQ files)"). Added
  `.fastq` + `.fq` to the suffix set. Caught while writing the regression
  test above.

Six new regression tests:
- `test_run_extract_multi_fastq_same_round_aggregates`
- `test_run_extract_sample_sheet_paired_end_demux_rebuilds_r2_inputs`
- `test_run_extract_sample_sheet_input_sha256_distinct_per_round`
- `test_compute_sha256s_distinct_keys_per_round_with_root`
- `test_compute_sha256s_without_root_falls_back_to_basename`
- `test_reorient_fastq_gz_raises_on_truncated_record`

CI: 367 passed + 1 xfailed (was 361; +6 new tests). All four pre-existing
CI gates stay green.

**Phase 2 Codex pass 1 semantic bugs (2026-05-20)** — surfaced during
the calibration review and folded into the same commit:

- **`match_rate_*` and `position_consistency_*` were aliased.** The
  orchestrator was setting `position_consistency_5p = match_rate_5p`
  literally, double-counting the same evidence in the composite
  confidence formula (which weights them as 4 distinct signals).
  Added `_substring_match_rate()` (substring-anywhere with Hamming ≤
  1), and now `match_rate_*` comes from substring scan while
  `position_consistency_*` comes from `_position_consistency()`
  (flank-anchored ± tolerance). They diverge by construction.
- **`adapter_clean` ignored the drop flag.** Previously the signal was
  `1.0 if (primer_5p_seq is not None or primer_3p_seq is not None)
  else 0.0` — true whenever ANY primer survived, hiding the adapter
  trap. Now tracks per-side `adapter_drop_{5p,3p}` flags at the drop
  site; the signal is `0.0` if either side was dropped as an adapter
  match.
- **Paired-split 3p signals were measured against R1.** In
  `PAIRED_END_SPLIT_PRIMERS` mode, `match_rate_3p`,
  `position_consistency_3p`, and `variants_3p` were computed by
  looking for the (R2-derived) `primer_3p_seq` at R1's 3' end — where
  it cannot exist by construction (that's why it's a split). Refactored
  the orchestrator to toggle a 3p-signal context (R1/3'-end normally,
  R2/5'-end-using-revcomp for paired-split). Per-round persistence
  input for the 3p side now uses the R2 stream in split mode.
- **`_persistence_score` returned `0.0` instead of `None` when mean <
  0.1.** Semantically: `None = "not evaluable"` vs `0.0 = "evaluated
  and zero"`. The composite-confidence formula treats both the same
  numerically (None contributes 0 via `if v is None: continue`; 0.0
  contributes 0 via 0.0 × weight), so the user-visible composite is
  unchanged. The semantic alignment matters for the manifest's
  audit trail.

In `PAIRED_END_SPLIT_PRIMERS` mode the `n_length_*` fields are now set
to `None` / `{}` / `0.0` rather than computed from R1 alone — the full
insert spans R1+R2 and cannot be measured from either read.

Three new regression tests in `tests/test_report.py`:
- `test_match_rate_distinct_from_position_consistency`
- `test_adapter_clean_flag_demotes_confidence_when_primer_dropped`
- `test_paired_split_match_rate_3p_reflects_r2_not_r1`

Behavior-based classification tests stay green throughout (361 passed
+ 1 xfailed; +3 vs previous baseline) — the LR's classification
fields (`extraction_mode`, `required_action`, `status`) are unchanged
for clean inputs; only the numeric signal fields change.

### Calibration

**Phase 2 Codex peer-review, pass 1 (2026-05-20).** Eight constants in
`library/detect.py` were reviewed; six confirmed at locked-plan
defaults, four revised:

| Constant | Before | After | Rationale (Codex) |
|---|---|---|---|
| `POSITION_CONSISTENCY_TOLERANCE` | 2 nt | **3 nt** | AptaPLEX default mismatch tolerance is 3; small public-data offset noise should not over-penalize an otherwise stable flank. |
| `STATUS_HIGH_CUTOFF` | 0.80 | **0.85** | "HIGH" should mean paper-grade high-confidence, harder to reach via additive secondary signals before benchmark calibration. |
| `COMPOSITE_WEIGHTS` (with round map) | match 0.20/0.20 · pos 0.10/0.10 · persistence 0.20 · n_len 0.10 · adapter_clean 0.10 | **match 0.15/0.15 · pos 0.15/0.15 · persistence 0.25 · n_len 0.10 · adapter_clean 0.05** | Position consistency deserves parity with raw match rate; cross-round persistence is the unique SELEX-specific signal (Hoinka et al. 2015, AptaTRACE / AptaTools); adapter_clean is already enforced upstream as a blacklist, so it should be a small confidence bonus, not a driver. |
| `COMPOSITE_WEIGHTS_NO_ROUND_MAP` | match 0.30/0.30 · pos 0.15/0.15 · persistence 0.00 · n_len 0.05 · adapter_clean 0.05 | **match 0.225/0.225 · pos 0.225/0.225 · persistence 0.00 · n_len 0.05 · adapter_clean 0.05** | With persistence absent and status already capped at MEDIUM, within-round evidence (match + position) should carry equal parity weight; n_len and adapter_clean remain supporting signals, not drivers. |

Confirmed as-is: `PRIMER_FOUND_MATCH_RATE_THRESHOLD = 0.70`,
`N_LENGTH_CONFIDENT_FRACTION = 0.80`, `UNABLE_TO_EXTRACT_MATCH_RATE =
0.40`, `STATUS_MEDIUM_CUTOFF = 0.60`, `STATUS_LOW_CUTOFF = 0.30`,
`ORIENTATION_REVERSED_FORWARD_MAX = 0.05`,
`ORIENTATION_REVERSED_REVERSE_MIN = 0.95`.

Tests stayed green across the threshold flip (358 + 1 xfailed) because
they assert on `extraction_mode` / `required_action` / `status`, never
on the threshold numbers themselves. Marker convention shifted:
`CALIBRATION-TODO` → `CALIBRATION-REVIEWED (Codex 2026-05-20, pass 1)`
for the 8 reviewed entries. Eleven `CALIBRATION-TODO` markers still
pending (Phase 5 qc flags, adapter blacklist composition, strand
report granularity).

**Phase 6 benchmark inputs (structural flags Codex raised, read-only).**
These are NOT actionable in Phase 2 calibration; they inform the Phase
6 benchmark dataset selection so empirical ground truth covers the
known edge cases:

1. **Inline-barcoded / messy round-map cases**: earliest-round
   consensus can fail when demux/round-map assumptions are wrong. The
   `--round-map` requirement (locked plan line 289) mitigates the
   common case, but the benchmark should include at least one messy-
   barcode dataset to stress the failure mode.
2. **One-sided primer detection vs partial homology**: orientation
   classification based solely on `primer_5p` vs `revcomp(primer_3p)`
   can be ambiguous if only one primer is inferred, or if the primers
   share local homology. Benchmark should include a dataset with
   homologous flanking regions.
3. **First-13-nt adapter blacklist** is appropriately conservative,
   but rare real primers sharing an Illumina-like prefix would be
   dropped. Acceptable trade-off; just documented for the methods
   section.

### Added

**Phase 1 — library modules ported (166 tests, all green):**

- `selexprep._common` — shared utilities: `iter_srr_files` (exact-name match, no SRR1234↔SRR12345 collisions), `load_csv`, `parse_round_number`, `setup_logging`.
- `selexprep.count.counter` — FASTQ.gz → parquet sequence counting with anchored linked-adapter cutadapt trimming; paired-end pair-sync; pyfastx fast path with gzip fallback; Shannon-entropy + singleton-fraction pool stats; multi-target BioProject layout support.
- `selexprep.fetch.metadata` — deterministic 5-level cascade for round inference (sample_attributes → sample_title → library_name/experiment_title/design_description → abstract count → manual review). Never guesses; 20% unknowns preferred over silent miscalls.
- `selexprep.fetch.discover` — multi-source SELEX discovery across nine sources (Seed YAML, ENA, NCBI SRA, GEO, UTexas Aptamer DB, Zenodo, Figshare, Crossref, OpenAlex). Optional dependencies (pysradb, Bio.Entrez) soft-fail. Library-type classification is conditional (deferred to v0.2).
- `selexprep.fetch.download` — MIT-compatible-first download dispatcher: **ENA-direct (default)** → kingfisher (optional, GPL-3.0) → sra-toolkit (optional). Pure-`requests` ENA filereport API with Range-resumable streaming + MD5 verification; no external tools required for the default install.
- `selexprep.library.detect` — empirical primer/constant-region inference. **Default: scans every unique sequence in the parquet (no top-N subsampling)** so the long tail of rare unique sequences also confirms the consensus.
- `selexprep.library.audit` — pre/post-extraction structure audits: raw-FASTQ sampling with 3'-aligned positional base frequencies (variable-length-tolerant), trimmed-parquet length distribution + TruSeq R1 contamination probe. Pure dataclasses; no CLI.
- `selexprep.extract.demux` — sample-sheet-driven barcode demultiplexer for pooled multi-round SELEX runs. Validates barcodes for Hamming distance, keeps paired-end R1/R2 in lockstep, never trims R2.
- `selexprep.qc.coverage` — per-BioProject round-coverage classification (all_rounds_public / partial / multiplexed_unrecoverable / unknown). Filter is a configurable callable (no thesis-specific hard-coding).
- `selexprep.qc.consistency` — k-mer Jaccard distance + monotonicity check across rounds. Strictly diagnostic; never reassigns rounds from enrichment signal.
- `selexprep.qc.readiness` — eight-section sequence-level readiness review (`pre / alphabet / lengths / trim_seq / composition / diversity / selection / consistency`). Tag-aware composition/diversity thresholds; tag is a per-call parameter (not a hard-coded BP map).

**CLI:** seven stub subcommands (`inspect`, `fetch`, `detect`, `extract`, `count`, `qc`, `run`) ship in the Typer dispatcher; wiring to the ported library modules lands in Phase 2.

### Notes
- v0.1 packaging note: `discover.py` keeps its nine adapter classes in one file. Splitting into a `selexprep.fetch.sources.*` subpackage is a v0.2 cleanup.
- `kingfisher` (GPL-3.0) was dropped from `pyproject.toml` — it remains a runtime-detected optional subprocess backend so a default `pip install selexprep` stays MIT-only.

### Known v0.1 follow-ups (documented Codex peer-review findings)

These were flagged during the Phase 0/1 peer-review and are **not** blocking Phase 2; capturing them so they don't get lost.

- **`count.counter` still trims raw FASTQs inline.** The Phase 1 port preserved the original `selex_corpus` behavior — `count_round()` can run cutadapt on raw FASTQ inputs. The plan's final shape is *`extract` produces primer-stripped FASTAs and `count` only counts those*. The separation lands when the `extract` step is wired in Phase 2/3; until then, `count.counter` is dual-purpose.
- **`qc.readiness` requires clusters / enrich parquets.** The module is a faithful port and still expects `round_*.clusters.parquet`, `enrich_*.parquet`, `summary.json`, and `cluster_stats.json` — artifacts that v0.1 does *not* produce (clustering / enrichment are out of v0.1 scope). It remains exposed as a library API so the thesis pipeline can use it, but it is **not** wired into the `selexprep qc` CLI verb. The v0.1 `qc` verb will get a thinner, manifest-driven implementation when `extract`/`count` are fully separated.
- **Mocked-HTTP coverage gap.** The nine network adapters in `fetch.discover` and `download_srr_*` paths beyond ENA-direct don't yet have offline mocked tests (only their parsing helpers + dispatcher + SeedAdapter are covered). To be addressed before PyPI release.

### Phase 5 — QC plots + flags + count CLI (2026-05-20)

Closes the v0.1 single-dataset CLI surface. The workflow is now
feature-complete end-to-end: ``detect`` -> ``extract`` -> ``count`` ->
``qc``. ``selexprep run`` (batch driver) and ``selexprep fetch``
(accession download) are deferred to Phase 6 / v0.2.

- **NEW: `selexprep.qc.diversity`** — `rarefy()` (deterministic
  multivariate hypergeometric subsampling via numpy default RNG;
  seeded for reproducibility), `shannon_entropy()` (base-2),
  `unique_count()`, `top_n_coverage()`. Used by both `flags.py` and
  `plots.py`.
- **NEW: `selexprep.qc.flags`** — eight depth-aware suspicion flags
  per locked plan lines 350-358:
  `unexpected_rarefied_diversity_increase` (rarefied uniques per
  round; not raw counts),
  `low_primer_match` (threshold imported from
  `library.detect.UNABLE_TO_EXTRACT_MATCH_RATE` — single source of
  truth),
  `n_length_variation_across_rounds`,
  `strand_mix` (from Phase 3 strand_report.tsv),
  `low_total_reads`,
  `adapter_contamination_high`,
  `extraction_mode_changed_across_rounds` (v0.1 inert — placeholder
  for Phase 6 `selexprep run` batch driver),
  `requires_read_merging_for_full_insert` (informational).
  `compute_all_flags()` aggregator + `write_flags_yaml()` deterministic
  emitter (sorted by flag name; float-rounded for cross-platform
  stability).
- **NEW: `selexprep.qc.plots`** — four per-dataset matplotlib plots
  (PNG, 150 DPI, Agg backend, tight bbox): read_retention.png,
  primer_match_per_round.png, n_length_distribution.png,
  per_round_panel.png. **Plots are informational only** — matplotlib
  PNG output is not byte-deterministic across versions; they do not
  contribute to `output_sha256`.
- **NEW: `selexprep.qc.runner`** — `run_qc(manifest_path, ...)`
  orchestrator. Auto-discovers `round_*/counts.parquet` under the
  manifest's directory; reads `trim_reports.json` for read-retention
  plot data; optionally reads `strand_report.tsv` for the strand-mix
  flag. Returns `QcResult` with the list of flags raised, the
  flags.yaml path, and the four plot paths.
- **NEW: `selexprep.count.counter.count_fasta`** — FASTA-aware
  per-round counter (the Phase 3 extract pipeline emits FASTA, not
  FASTQ). Reuses `_counter_to_parquet` so the output schema matches
  `count_round`: `sequence`, `reads`, `rank`, `rpm`.
- **WIRED: `selexprep count <extracted-fasta> --round R0 --outdir OUT`**
  — accepts `R0`, `r0`, `round_0`, or just `0` for the round label;
  writes to `OUT/round_NN/counts.parquet`.
- **WIRED: `selexprep qc <manifest> [--counts-dir DIR] [--outdir OUT]`**
  — prints a one-line summary plus per-flag severity to stdout;
  emits `flags.yaml` and the four PNG plots.
- **Tests added (+56, 358 + 1 xfailed total):**
  `tests/test_diversity.py` (20 — rarefy determinism + edge cases,
  Shannon entropy, top-N coverage monotonicity, depth-aware sanity
  check), `tests/test_flags.py` (21 — positive + negative case per
  flag + aggregator + YAML determinism), `tests/test_plots.py` (5 —
  PNG smoke), `tests/test_qc_runner.py` (6 — end-to-end manifest ->
  flags.yaml + 4 PNGs with realistic synthetic data), `tests/test_cli.py`
  (+4 — count + qc smoke).
- **CALIBRATION-TODO inventory: 19** (was 12). New tunables in
  `qc/flags.py`: rarefaction depth, max modal lengths, strand-mix
  max reverse fraction, low-total-reads minimum, adapter-contamination
  max fraction; in `qc/plots.py`: top-N coverage N. The match-rate
  threshold is **imported** from `library.detect` (not redeclared)
  so a single Codex pass tunes both QC and classifier.

### Phase 4 — Manifest + inspect + extract override/rebuild (2026-05-19)

Closes the v0.1 CLI surface (except QC, Phase 5). Every `extract` run
now emits a `selexprep_manifest.json` — the reproducibility anchor that
Phase 5 `qc` and future v0.2 AnnData export will consume.

- **NEW: `selexprep.manifest`** — `SelexprepManifestV1` pydantic model
  (frozen, extra=forbid) with the locked schema (plan lines 162-175):
  `manifest_version`, dep versions (selexprep / python / cutadapt /
  dnaio / pyarrow), provenance (accession / bioproject_id / runs),
  `input_sha256` + `output_sha256` (FASTA/TSV/JSON only — Parquet
  hashes intentionally absent per locked plan line 28), nested
  `LibraryReport` + denormalized scan fields, CLI argv capture in
  `parameters`, runtime/flags/sampling_seed. Helpers: `compute_sha256s`,
  `write_manifest_json` / `read_manifest_json` with deterministic JSON
  (same numeric-int-key + alphabetical-sha256-keys discipline as
  `library/report.py`), `build_manifest_from_extract_result`.
- **NEW: `selexprep.fetch.inspect`** — `inspect_accession()` hits ENA
  Portal filereport REST (`https://www.ebi.ac.uk/ena/portal/api/filereport`);
  parses run/study metadata into `InspectReport` + `RunFileInfo`
  dataclasses; tolerant of missing fields. Reports
  `library_strategy` / `library_source` **verbatim from SRA** — NOT a
  DNA/RNA classification (locked plan line 332 explicit on this;
  classification deferred to v0.2's library-type-classifier).
- **EXTENDED: `selexprep.extract.runner`** — `run_extract()` now accepts
  `override_primer_{5p,3p}` (cloned via `LibraryReport.model_copy`),
  plus provenance kwargs (`accession`, `bioproject_id`, `runs`,
  `parameters`) for the auto-emitted manifest. Override without
  `--rebuild` routes outputs to `<outdir>/overridden/` (preserves
  baseline). Override + `--rebuild` overwrites baseline AND emits
  `extract_diff.tsv` comparing baseline vs override per-round read
  counts. The diff TSV is read from the baseline `selexprep_manifest.json`
  + `trim_reports.json` BEFORE overwrite; gracefully degrades if either
  baseline artifact is missing/malformed.
- **WIRED: `selexprep inspect <accession>`** — full CLI. Prints a
  human-readable metadata summary; `--outdir` also writes a sorted-keys
  `inspect.json`. `--timeout-s` controls the HTTP timeout (default 30s).
- **EXTENDED: `selexprep extract`** — `--override-primer-{5p,3p}` now
  works (lifts the Phase-3 informative error); CLI argv is captured
  into the emitted manifest's `parameters` field.
- **Tests added (+28, 302 + 1 xfailed total):** `tests/test_manifest.py`
  (12 — schema fields + frozen + extra-forbid + deterministic JSON +
  int-key sort + sha256 helper FASTA/TSV/JSON-only behavior),
  `tests/test_inspect.py` (8 — mocked ENA REST: single run / paired
  semicolon lists / multi-run study / missing-field tolerance / empty
  response → ValueError / HTTP error propagation / timeout pass-through
  / JSON sort-keys),  `tests/test_extract_override.py` (7 — override 5p
  + 3p without rebuild → subtree, override + rebuild → diff TSV + in-
  place overwrite, rebuild alone → no diff, manifest emission, override
  primer recorded in manifest), `tests/test_cli.py` (+2 — inspect smoke
  with mocked REST + override smoke). The old "Phase 4 error" test
  replaced with a real override-works smoke test.
- **CALIBRATION-TODO inventory: 12** (unchanged — Phase 4 is
  serialization + REST + I/O wiring, no new heuristic thresholds).

### Phase 3 — extract: paired-end + strand orientation (2026-05-19)

Turns the Phase 2 `LibraryReport` contract into actual extracted FASTAs.
Cutadapt is invoked as a subprocess (per locked plan); dnaio is available
for paired I/O. **No read merging in v0.1** — paired-end split-primer
mode emits two separate files and flags `READ_MERGING_RECOMMENDED`.

- **NEW: `selexprep.extract.trim`** — cutadapt subprocess wrapper with
  four public entry points (`trim_single_end_linked` for
  BOTH_PRIMERS_SINGLE_READ, `trim_single_end_5p` for FIVE_PRIME_ONLY,
  `trim_single_end_3p` for THREE_PRIME_ONLY, `trim_paired_split` for
  PAIRED_END_SPLIT_PRIMERS). Each returns a `TrimReport` carrying the
  exact cutadapt argv + read counts (Phase 4 manifest precursor).
  Cutadapt writes uncompressed FASTA; this module re-gzips with
  `_io.open_gzip_text_deterministic` (mtime=0 header) so `output_sha256`
  is bit-identical across reruns.
- **NEW: `selexprep.extract.strand`** — strand-orientation handler.
  `detect_strand_distribution()` counts forward/reverse/ambiguous reads;
  `reorient_fastq_gz()` reverse-complements every record (sequence +
  reversed quality string) for `LibraryReport.orientation == "REVERSE"`;
  `write_strand_report()` emits a sorted TSV for the QC trail.
- **NEW: `selexprep.extract.runner`** — `run_extract()` orchestrator.
  Refuses if `LibraryReport.status == "UNABLE_TO_INFER"` or
  `extraction_mode == "UNABLE_TO_EXTRACT"` (no silent miscalls).
  Optional sample-sheet pre-step demuxes multiplexed input; strand
  pre-step rewrites all reads when orientation is `REVERSE` and emits
  `strand_report.tsv` for `MIXED` or `REVERSE`. Per-mode trim dispatch
  writes per-round outputs to `<outdir>/round_NN/<filename>.fasta.gz` +
  `trim_reports.json` (manifest precursor).
- **Output filename contract** (locked plan lines 321-326):
  `extracted.fasta.gz` (full insert), `partial_5p_extracted.fasta.gz` /
  `partial_3p_extracted.fasta.gz` (one-sided), `partial_5p_extracted_R1.fasta.gz`
  + `partial_3p_extracted_R2.fasta.gz` (paired split). Filenames signal
  to downstream ML pipelines whether a full insert was recovered;
  joining R1+R2 by read ID alone is biologically wrong, so
  `joined_counts.tsv` is **not** emitted in v0.1.
- **WIRED: `selexprep extract`** — full CLI. Accepts `--library-report`,
  `--round-map`, `--sample-sheet`, `--paired-r2`, `--rebuild`.
  `--override-primer-{5p,3p}` emit a Phase-4 informative error (full
  diff TSV lands in Phase 4). `--rebuild` toggles the no-clobber guard.
- **Tests added (+37, 274 + 1 xfailed total):** `tests/test_strand.py`
  (13 — distribution + revcomp + deterministic gzip + TSV sort),
  `tests/test_trim.py` (7 — per extraction_mode + determinism + temp
  cleanup, skips if cutadapt absent), `tests/test_extract_runner.py`
  (13 — happy path per mode + UNABLE refusal + no-clobber + rebuild +
  strand-report emission + trim_reports JSON + multi-round),
  `tests/test_cli.py` (+4 — missing-round-map, override Phase-4 error,
  UNABLE refusal, end-to-end smoke).
- **CALIBRATION-TODO inventory: 12** (was 11). New tunable:
  `STRAND_REPORT_PER_READ = False` in `extract/strand.py`. Strand
  classification thresholds (`ORIENTATION_REVERSED_FORWARD_MAX`,
  `ORIENTATION_REVERSED_REVERSE_MIN`) stay in Phase 2's
  `library/detect.py` — no duplication.

### Phase 2 — LibraryReport schema + cross-round inference (2026-05-19)

Adds the `LibraryReport` pydantic schema and the cross-round primer
inference pipeline that turns Phase 1's single-pool flank detector into
the typed contract every downstream stage (`extract`, `count`, `qc`,
`manifest`) consumes.

- **NEW: `selexprep.library.report`** — `LibraryReport(BaseModel)`
  with the locked schema (plan lines 233-285), `Literal` aliases for the
  five categorical fields (`ExtractionMode`, `ReadSource`,
  `RequiredAction`, `Orientation`, `Status`), the `_classify` pure
  function implementing the locked decision table (plan lines 300-309)
  with the no-round-map status cap (line 289), and deterministic JSON
  I/O (`write_library_report_json` / `read_library_report_json` —
  bit-identical output across reruns, numeric ordering for int-keyed
  dicts). Strict-mypy clean (pydantic plugin enabled).
- **NEW: `selexprep.library.adapters`** — conservative v0.1 blacklist
  (TruSeq R1 + Nextera) with auto-computed reverse complements,
  `reverse_complement()` helper (rejects IUPAC ambiguity in v0.1), and
  `count_adapter_hits()` substring scanner (records hits; does NOT
  filter reads).
- **EXTENDED: `selexprep.library.detect`** — `compute_library_report()`
  orchestrator. Cross-round persistence as `1 - clip(stdev/mean, 0, 1)`,
  position consistency with ±2 nt tolerance, U→T normalization of RNA
  primers, paired-end split detection (R1 5' + R2 5' = revcomp(3'
  primer)), MIXED/FORWARD/REVERSE orientation diagnostic, composite
  confidence via weighted sum (two regimes: with vs without round map).
  Phase 1 functions (`detect_flank`, `detect_primers`,
  `detect_from_parquet`, `earliest_round_parquet`) are unchanged and
  remain the algorithmic primitives the new orchestrator consumes.
- **WIRED: `selexprep detect`** — CLI command parses `--round-map` TSV
  (columns `file<TAB>round_number`), groups FASTQs by round, runs
  `compute_library_report`, writes `library_report.json` to `--outdir`.
  Refuses to run without `--round-map` (cross-round persistence is a
  core inference signal). Single-end only in Phase 2; paired-end via
  `compute_library_report`'s `paired_mate_streams` kwarg awaits CLI
  surface in Phase 3.
- **INFRASTRUCTURE:** `pyproject.toml` declares
  `plugins = ["pydantic.mypy"]` under `[tool.mypy]` (required for
  strict-mypy to resolve `BaseModel` field types — without it every
  field decays to `Any`).
- **Tests added (32 new, 234 + 1 xfailed total):** `tests/test_adapters.py`
  (14 tests covering revcomp + blacklist + substring scan);
  `tests/test_report.py` (19 tests covering every row of the locked
  classification table plus edge cases — status cap, adapter
  demotion, orientation, U→T, deterministic serialization, schema
  immutability, numeric int-key ordering, empty input, sub-floor
  input); `tests/test_cli.py` (3 new tests for `detect` CLI:
  missing-round-map, end-to-end round-trip, FASTQ-not-in-map).

#### Calibration status — placeholder pending Codex peer review

Codex usage was rate-limited 2026-05-19 → 2026-05-26 when Phase 2
shipped, so calibration numbers ship as locked-plan literals (match
rates `> 0.7`, n_length confidence `> 0.8`, UNABLE floor `< 0.4`) plus
placeholders for everything the locked plan does not pin down
(composite weights, status cutoffs, position-consistency tolerance,
adapter list exact composition, persistence formula). Every numeric
placeholder carries a `# CALIBRATION-TODO` comment naming the locked
plan line (where applicable) or "not in locked plan - placeholder
pending Codex" otherwise.

Test discipline: **all tests assert on behavior, never on threshold
values** — e.g. `assert report.extraction_mode == "BOTH_PRIMERS_SINGLE_READ"`,
never `assert HIGH_CONFIDENCE_CUTOFF == 0.80`. So when Codex tuning
lands (or Phase 6 benchmark numbers update the constants), the test
suite stays green by construction.

Recovery list: `grep -rn "CALIBRATION-TODO" src/` returns the full
inventory.

### Phase 1.5.1 — catalog refresh against broad ENA queries

A Codex / sanity-check pass after Phase 1.5 revealed that the initial
bundled catalog (219 rows, sourced from the thesis-specific
`selex_corpus.discover` run) under-counted INSDC studies by ~50%: 94
INSDC accessions vs ~120 unique studies surfaced by broader keyword
queries against ENA. The thesis queries combined keywords with AND
clauses to maximize precision; a generic-tool catalog wants the
broader OR-style net.

- **NEW: `selexprep.catalog.rebuild`** — reproducible refresh script
  that runs 13 broad ENA queries (HT-SELEX, SELEX-seq, SELEX, aptamer,
  Cell-SELEX, RNA aptamer, DNA aptamer, systematic evolution, SELEX
  rounds, …), unions the studies, merges hand-enriched fields
  (`protein_target` / `paper_doi` / `paper_pmid` /
  `n_rounds_declared`) forward from the previous catalog when an
  accession is still upstream, and carries non-INSDC deposits
  (zenodo/figshare/utexas processed-data entries) across refreshes
  unchanged.
- **NEW: `selexprep catalog refresh [--out PATH --no-preserve-enrichment]`**
  — CLI verb that wraps `rebuild_catalog`. Lets users (or CI) refresh
  the catalog on demand without touching the package source.
- **Catalog refreshed in-place** for v0.1.5: 273 bioprojects (148
  ENA-discovered INSDC studies + 125 carried-forward
  Zenodo/Figshare deposits). Snapshot bumped to
  `v0.1.5-snapshot-2026-05-19`. The 4 seed entries (Hoinka IL-10RA,
  Dao CCR7, …) keep their hand-curated enrichment.
- **No curation flags.** Confirmed with the PI: the package never
  ships `include` / `manual_curation_notes` columns. Curation is the
  user's downstream job; the catalog reflects the public archives
  only.

### Phase 1.5 — discovery catalog (new)

The biggest "where do I even start?" UX gap in v0.1 was: a user installs the
package, knows nothing about which public SELEX accessions exist, and is
expected to run a multi-API discovery scan before anything useful happens.
Phase 1.5 fills this gap by shipping a **bundled discovery catalog** as
package data.

- **`selexprep.catalog`** — new subpackage with `load_catalog()`,
  `filter_catalog()` (target / organism / source / min-rounds / INSDC-only),
  and a `catalog_version()` snapshot identifier.
- **`selexprep catalog list/show/version`** — new Typer subapp wired into the
  root CLI. `list` supports the same filters as the Python API; `show
  <accession>` prints the full row including the study abstract.
- **Catalog content (v0.1 snapshot):** 219 public SELEX bioprojects with
  bioproject_id / source / study_title / protein_target / target_organism /
  paper_doi / paper_pmid / n_rounds_declared / abstract. Thesis-specific
  columns (`include`, `manual_curation_notes`,
  `library_type_verification`, `library_type_evidence`, `has_processed_counts`)
  are intentionally stripped so the catalog reflects the public archives,
  not any single researcher's curation.
- **v0.2 plan:** enrich each catalog row with the inferred `LibraryReport`
  (primer pair, N-region length, extraction_mode, confidence) once the full
  Phase 2 pipeline runs end-to-end. The enriched catalog will be deposited
  to Zenodo with a DOI and unlocks the NAR Database Issue paper venue
  alongside the planned Bioinformatics Advances Application Note.

### Improvements landed alongside the Phase 0/1 wrap

- `selexprep._io` — new module with `open_gzip_text_deterministic()` and `sha256_file()`. All `.gz` writes now produce bit-identical bytes across reruns (gzip header `mtime=0`), making manifest SHA256 hashes reproducible.
- `extract.demux` switched to the deterministic gzip writer; new regression test (`test_demux_output_is_deterministic_across_reruns`) compares byte-for-byte.
- `download_srr(backend="auto|ena|kingfisher|sra")` — explicit backend selector. `auto` (default) preserves the ENA-first dispatch; `ena` is the pure-MIT path with no fallback; `kingfisher` and `sra` force a specific backend. Any kingfisher invocation now logs an explicit GPL-3.0 notice.
- `fetch.discover._classify_all` writes a `NOT_ASSESSED_V0_1` sentinel + evidence JSON when the v0.2 classifier is absent. Distinguishes "deferred to v0.2" from an empty-string verdict downstream callers might silently treat as success.
- `qc.readiness.review_bioproject(tag=None)` now WARNs at call time when callers omit the tag (silently defaulting to `"untagged"` is the most likely call-site bug).
- `numpy` added as an explicit core dependency (it was being imported directly while only present transitively via `pandas`).

## [0.1.0] — TBD

First public release. See implementation plan for the locked v0.1 feature set.
