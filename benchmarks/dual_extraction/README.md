# Dual independent extraction — protocol (v1.4)

Methodology for the per-deposit SELEX metadata that ships with `selexprep` and is
described in the Application Note. Each text-derived field is extracted by **two
independent LLMs** (Claude and Codex/GPT), then reconciled. Pre-registered here so
there is no post-hoc tuning.

## Why two extractors
Single-LLM curation is a credibility liability and an AI-disclosure problem.
Double independent extraction + reconciliation is the systematic-review standard;
it turns "an LLM made it" into a reproducible, source-cited, human-verified
resource with a reportable inter-extractor agreement rate. AI use is disclosed
per OUP policy; the authors remain responsible for the final dataset.

## Independence (hard rules)
- Each arm runs in its **own clean directory outside `~/Documents`** (e.g.
  `~/extract_claude`, `~/extract_codex`), in **fresh agentic sessions** with web access —
  not in this repo dir, so neither `project_metadata.*` nor `CLAUDE.md` is in scope → no
  contamination from the prior curation.
- **Small batches in fresh sessions.** Accessions are processed in batches of ~10–15 per
  fresh session (the size validated by the pilot). A long list in a single session degrades
  into a metadata scrape — the failure mode this protocol exists to avoid. Do NOT feed the
  whole list at once, and do NOT instruct the model to "loop / append / resume", which reads
  like a program spec and makes it write a scraper.
- **Source-grounded, ID-only**: each model is given just the accessions and looks up the
  deposit record + associated publication itself; a value is reported only with a verbatim
  quote + the source (URL/DOI) consulted; if absent after reading → `not_stated` (never invented).
- **Identical contract** ([`PROMPT_v1.md`](PROMPT_v1.md)) for both arms; each looks up its
  sources independently (a difference in the source found surfaces as a disagreement at
  reconciliation).
- `temperature = 0`. Each arm is **blind** to the other. The contaminated main session only
  *designs* the protocol; it does **not** extract.

## Two steps: extract, then normalize (decided 2026-06-16; Claude+Codex unanimous)
1. **Extraction = free-text / verbatim**, no controlled-vocab menu in the prompt.
   The extractor reports what the source says, with an evidence quote + location.
   (A menu in the prompt forces borderline cases identically across models →
   inflates fake agreement and hides normalization.)
2. **Normalization = separate deterministic step** ([`normalize.py`](normalize.py) +
   [`normalization.yaml`](normalization.yaml)): maps verbatim values onto controlled
   vocabularies. Inter-extractor agreement is computed on **normalized** values.
   An off-table label is flagged `UNMAPPED` → a human adjudication decision (add it
   to the table), never silently forced. The table is versioned and citable.

`study_type` is the exception: it is a project-defined classification with no
verbatim value in the paper. The prompt gives the **definitions/criteria** (not a
forced pick) and requires a one-sentence rationale + the supporting quote →
reasoned-but-grounded classification.

## Fields extracted from text
`study_type`, `target`, `target_class`, `chemistry`, `n_random`, `n_rounds`,
`selection_format`, and `counter_selection` (the last as `present`/`absent` + a
free-text `note`). `target_organism` was **dropped**: redundant with
`target` for whole-organism targets, meaningless for ions/small molecules, and the
only field where the prior curation was wrong (it had carried the INSDC "synthetic
construct" library organism, not the target's organism); species, when needed, is
recoverable from `target`. It was likewise removed from `project_metadata.{csv,json}`.
Everything else (`accession`, `study_title`, `paper_*`, computed signals, `rounds`)
is deterministic/derived — **not** double-extracted.

## Run (replicate the validated pilot, in small batches)
The pilot (11 deposits) extracted at full depth — it read the deposits *and* the papers
(down to supplementary methods). The first attempt at the full 238 degraded to a metadata
scrape **because the run instruction had become a program spec** ("process one-at-a-time,
append to `extractions.json`, skip if already present, resumable"): both models implemented
it literally as a script (`append_record.py`, `extract_accessions.py`) and stopped at the
deposit record. The fix is to run the corpus exactly like the pilot — just in batches.

1. In a clean dir per arm (outside `~/Documents`) containing `PROMPT_v1.md` and a batch file
   `accessions_N.txt` (~10–15 accessions from `accessions_full.txt`), open a **fresh agentic
   session with web access** (Claude Code / Codex CLI) and give it the **launch prompt below**
   (identical for both arms). It reads `PROMPT_v1.md` from disk (the `extract-v1.4` contract),
   processes each accession in `accessions_N.txt`, and writes `accessions_N.json`:

   ```
   Read PROMPT_v1.md in this directory and first print its CONTRACT VERSION line to
   confirm you've read it. Then, following its rules and JSON schema exactly, extract
   SELEX metadata for every accession in accessions_N.txt: for each, look up the deposit
   and read its associated paper's full text yourself; fill each field only with a
   verbatim quote + source URL/DOI, otherwise not_stated. Write the results as one JSON
   array to accessions_N.json. Work autonomously, no questions. Stop at the end of the
   batch, working one accession at a time (including the JSON write).
   ```
2. Repeat with a **fresh session per batch** until all 238 are done (~16–24 batches per
   arm). Fresh sessions keep each batch small enough to stay deep; "resumability" is simply
   which batches you have saved. Do **not** add loop / append / resume language.
3. Merge each arm's batch arrays into one file (run inside each arm's dir, over its
   `accessions_*.json` batches):
   ```bash
   # in the Claude arm dir:
   python3 -c "import json,glob; r=[x for f in sorted(glob.glob('accessions_*.json')) for x in json.load(open(f))]; json.dump(r,open('claude_extractions.json','w'),indent=2,ensure_ascii=False)"
   # in the Codex arm dir:
   python3 -c "import json,glob; r=[x for f in sorted(glob.glob('accessions_*.json')) for x in json.load(open(f))]; json.dump(r,open('codex_extractions.json','w'),indent=2,ensure_ascii=False)"
   ```
4. Compute agreement + the adjudication worklist:
   ```bash
   python compute_agreement.py claude_extractions.json codex_extractions.json --outdir out
   ```
   → `out/agreement_summary.tsv` (per-field agreement) + `out/disagreements.tsv`
   (both values + both quotes, for human adjudication).
5. Adjudicate `out/disagreements.tsv` by re-reading the cited locations → final verified
   values. Concordant + adjudicated = `curation_level = verified`.

**Record `PROMPT_v1.md`'s SHA256 in the Methods** so reviewers can audit the exact contract.

**Sanity check before trusting a batch:** for fields normally only in the paper
(`n_random`, `chemistry`, `counter_selection`), a healthy batch has most records citing a
paper `source` with `tier = FULL_TEXT`. If those fields are overwhelmingly `not_stated`
with paper `tier = RECORD_ONLY`, the model scraped instead of reading — discard that batch
and re-run it in a fresh session.

## Adjudication conventions (for the Methods)
Two conventions resolve the recurring borderline / free-text cases, fixed from the pilot:
- **`target`** — when both extractors name the same entity at different granularity
  (full name vs abbreviation, e.g. "African swine fever virus (ASFV) p30" vs "ASFV p30"),
  keep the **more complete form**. These are not real disagreements.
- **`study_type`** — **data-centric**: a deposit with genuine selection-round reads is
  `aptamer_selection` even if the paper frames itself as a method/protocol;
  `method_or_other` is reserved for deposits without their own selection reads
  (simulation, reanalysis, tool/benchmark). This is also baked into `PROMPT_v1.md`.

## Record shape (one per deposit, per extractor)
```json
{
  "accession": "PRJDB9110",
  "extractor": "claude",
  "model": "claude-opus-4-8",
  "prompt_version": "extract-v1.4",
  "source_documents": [
    {"type": "paper", "id": "doi:10.xxxx/xxxx", "tier": "FULL_TEXT"},
    {"type": "deposit", "id": "PRJDB9110"}
  ],
  "fields": {
    "study_type":  {"value": "aptamer_selection", "rationale": "in vitro selection of RNA aptamers against a protein target", "evidence_quote": "we performed HT-SELEX to generate RNA aptamers against ...", "source": "doi:10.xxxx/xxxx", "location": "Abstract"},
    "target":      {"value": "human transglutaminase 2", "evidence_quote": "aptamers against human transglutaminase 2", "source": "doi:10.xxxx/xxxx", "location": "Title"},
    "target_class":{"value": "protein", "evidence_quote": "the protein transglutaminase 2", "source": "doi:10.xxxx/xxxx", "location": "Intro, p.2"},
    "chemistry":   {"value": "RNA", "evidence_quote": "an RNA library", "source": "doi:10.xxxx/xxxx", "location": "Methods"},
    "n_random":    {"value": "40", "evidence_quote": "40-nt randomized region", "source": "doi:10.xxxx/xxxx", "location": "Methods"},
    "n_rounds":    {"value": "8", "evidence_quote": "eight rounds of selection", "source": "doi:10.xxxx/xxxx", "location": "Methods"},
    "selection_format": {"value": "HT-SELEX", "evidence_quote": "high-throughput SELEX", "source": "doi:10.xxxx/xxxx", "location": "Abstract"},
    "counter_selection": {"value": "present", "note": "negative selection against mock beads", "evidence_quote": "a counter-selection against mock beads was performed", "source": "doi:10.xxxx/xxxx", "location": "Methods"}
  }
}
```
A stated field carries `value` + `evidence_quote` + `source` + `location`
(+ `rationale` for `study_type`). `counter_selection.value` is `present`/`absent`
(+ a free-text `note`). An absent field is `{"value": null, "status": "not_stated"}`.
The paper `tier` (`FULL_TEXT` / `ABSTRACT` / `RECORD_ONLY`) records how deeply the
source was read; a paper-only field may be `not_stated` only when the tier explains it.

## Files
- `PROMPT_v1.md` — the identical extraction contract (`extract-v1.4`): a batch of accessions in, a JSON array out
- `accessions_pilot.txt` / `accessions_full.txt` — the 11 benchmark deposits / all 238
- `normalization.yaml` — controlled vocabularies + verbatim→canonical mappings (versioned)
- `normalize.py` — deterministic normalizer used by the comparison scripts
- `compute_agreement.py` — per-field inter-extractor agreement + adjudication worklist
- `compare_to_existing.py` — 3-way comparison vs the existing `project_metadata.json` curation
