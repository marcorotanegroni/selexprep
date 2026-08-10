# Dual-extraction → published metadata — PLAN

Contract: `extract-v1.4` — SHA256 `c069d2dfb5a9a19d8a7aa26cc0beac68f5a19d0e14f8f0baeb7273de2a87631f`
Two independent arms: **Claude** (`~/Estrazione_claude_batch`, 01–16) + **Codex**
(**canonical = `~/Estrazione_codex_batch`**, one json per batch; `~/Estrazione_codex_batch_2`
+ `~/Estrazione_codex_batch/_archive/` = originals/backups).
Extraction rule: **medium** effort, per-batch (~15 accessions), fresh session each. Never max.
Keep at the end: **two raw extractions + agreement report = EVIDENCE; reconciled table = PRODUCT.**

## Phase 0 — Codex re-run of the weak batches (14/15/16)
Reason: 14/15/16 under-extracted vs Claude (paper-only fields missing though the data exists).
Output to `*.v2.json` (don't clobber originals), medium, one accession at a time.
- [x] 14 → `accessions_14.v2.json`
- [x] 15 → `accessions_15.v2.json`
- [ ] 16 → `accessions_16.v2.json`

## Phase 1 — Consolidate the Codex arm (best per record)
- [x] 09 → keep `codex_batch` (folder2's was shallow, all-ABSTRACT)
- [x] 10 → keep `codex_batch` (complete + marginally richer)
- [x] 11 → best-of-both: 10 deep (`codex_batch`) + 5 figshare (`codex_batch_2`) — orig → `*.incomplete.bak.json`
- [x] 12 → best-of-both — orig → `*.incomplete.bak.json`
- [x] 13 → copied into folder1 (canonical); already deep (FULL=11, chemistry 15/15)
- [x] 14 → best-of(new re-run, old) → `accessions_14.json` (backup `.rerun.bak.json`); re-run OK: paper_src 36 (trails Claude 60 ma solido)
- [x] 15 → best-of(new re-run, old) → `accessions_15.json` (backup `.rerun.bak.json`); stated 60 (≥ Claude 56) MA provenienza debole (paper_src 8 vs 35: Codex ha trovato meno paper)
- [ ] 16 → compare `v2` vs old vs **Claude** → keep best per record (after Phase 0)
- [ ] 16 → consolidate into folder1 (after Phase 0 re-run)
- [ ] Assemble `codex_extractions.json` (238) = `glob ~/Estrazione_codex_batch/accessions_*.json` (canonical, one/batch, no `.bak`) — **needs 16** (now 225/238)

## Phase 2 — Assemble the Claude arm
- [x] 238 extracted; provenance clean (0 real traceability defects; the 2 figshare = theses, correct)
- [ ] Merge `~/Estrazione_claude_batch/accessions_*.json` → `claude_extractions.json` (238)

## Phase 3 — Reconciliation (Claude vs Codex) — PIVOT: curated resource, not a reproducible metric
- [x] `compute_agreement.py` run as a **scaffold** → flagged 747 cells (`out/disagreements.tsv`)
- [x] **Claude did the semantic comparison** of the 747: ~700 = wording the mind collapses → concordant; **~47 genuine** → keep both
- Note: "reproducible" (PI's sense) = the *method* is repeatable by anyone with model access, not an identical number; value to the note = its presence

## Phase 4 — Reconcile → published metadata
Per (accession, field): concordant → `verified`; conflict → **adjudicate from source** → `adjudicated`;
presence-gap → take the deeper arm (usually Claude) → `single_source`; both absent → `not_stated`
(the old CSV is the tertiary safety net for the rare both-absent-but-CSV-has case).
- [x] Built `curated_metadata.json` (238): per cell `status` (concordant/discordant/single_source/not_stated) + value + **inline provenance**; the **47 discordant carry BOTH arms** (value+source+location)
- Counts at that point: 1194 concordant / 238 single_source / 47 discordant / 425 not_stated
  (of 1904 cells, 238 deposits)
- **Superseded.** The catalog then grew to 240 deposits and the 47 disagreements were
  adjudicated one by one against the sources, so the shipped layer
  (`v0.3.1-dual-extraction-adjudicated-2026-08-09`) is 1206 concordant / 236 single_source /
  47 adjudicated / 2 verified / 429 not_stated of 1920 cells, with **no discordant cells left**.
  Adjudicated cells keep both arms plus the rule applied and the reasoning; the per-cell trail is
  in `adjudication_worklist.tsv`.
- [x] Emitted flat `curated_metadata.csv` (value + `_curation` per field; discordant → `claude || codex`); provenance companion = `curated_metadata.json`
- Information gain vs old `bioprojects.csv`: experimental fields **4 → 1479 cells (~370×)**, all source-cited (old had 0 evidence-cited values); `target_organism` artifact dropped

## Phase 4.5 — Catalog cross-analysis + release QC (`src/selexprep/catalog/data/bioprojects.csv`)
The shipped catalog (238/238 same accessions) is API-derived, **NOT a third extraction arm** → OUT of the
agreement metric, and **do NOT backfill its DOIs** into the published table (a DOI with no extraction behind
it is a hollow reference, and injecting an unread source would fake provenance we never used).
- [ ] **`paper_doi` (114/238) = INTERNAL EVALUATION of "did source-grounded dual extraction pay off"** (NOT a dataset field).
  Preview: extraction independently recovered **104/114 = 91%** of the API DOIs, and found a paper for **+98**
  deposits the API had NO DOI for → ~doubles paper coverage (114 → ~202) vs deterministic accession→DOI linking.
  (Full pass: also verify the 104 cite the *same* DOI, not just *a* paper.) → a figure/number for the Application Note.
- [ ] The **10** "catalog-DOI but extraction found no paper" = gaps to optionally **verify/re-extract** (targeted), NOT blind-inject.
- [ ] Carry forward the deterministic, non-extracted fields: `source`, `study_title`, `abstract` (238/238); confirm **no regression**.
- [ ] Drop the wrong `target_organism` ("synthetic construct" artifact, 119/238) — already decided.
- [ ] Changelog for release notes: v2 ADDS chemistry/n_random/counter_selection/study_type/target_class,
  fills protein_target (was 2/238) + real paper DOIs, and FIXES target_organism.

## Phase 5 — Package release
- [x] Wired the curated layer into the package: `catalog/metadata.py` (`load_metadata` / `load_metadata_records` / `metadata_version`), data bundled via wheel force-include, exported from `catalog/__init__`; **v0.1.1 → 0.2.0** (pyproject + `__init__`), CHANGELOG `[0.2.0]`, Dockerfile + `.def` bumped. Tests **617 passed**, ruff + mypy clean.
- [ ] USER: `git add` + commit + push; GitHub Release tag **v0.2.0** → OIDC publish-pypi (same flow as 0.1.1)
- [ ] Bioconda 0.2.0 PR (later): bump `conda-recipe/meta.yaml` version + recompute sha256 from the 0.2.0 sdist
- [ ] Application Note cites the resource + the dual-extraction-with-provenance method + the 91% recall / +98-papers result
- [ ] Supplementary/repo: `benchmarks/dual_extraction/` (contract, both raw arms, reconciliation, PLAN) = the evidence trail
