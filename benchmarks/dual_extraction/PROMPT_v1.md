# Extraction prompt — `extract-v1.4`

**How to use (this header is NOT part of the prompt):**
- Runs in an **agentic session with web access** (Claude Code / Codex CLI): the model
  must open and read the deposit record and the associated paper itself — a plain
  non-browsing chat will not work, and an automated metadata scrape is forbidden below.
- Paste a **small batch** of accessions (~10–15; see README → "Run") where
  `<<PASTE ACCESSIONS HERE>>` is, in a **fresh** session, and send. The model reads each
  deposit + its paper and returns one JSON array. Small batches in fresh sessions are
  what keep the per-deposit depth — a long list in one session degrades into a scrape.
- Use the **identical** contract for both extractors (Claude and Codex); do not edit it
  between models.

=== COPY FROM HERE ===

CONTRACT VERSION: extract-v1.4

You are a careful, human-style data curator for a systematic curation of public HT-SELEX sequencing deposits. You are given a LIST OF ACCESSIONS at the end of this message. For EACH accession, look it up yourself, read its sources, and extract its metadata.

THIS IS MANUAL CURATION, NOT AN AUTOMATED SCRAPE:
- Do NOT write or run any script, program, regex, notebook, or bulk API harvester to populate the fields. Treat each accession independently and thoroughly, reading the sources yourself the way a human curator would.
- For EACH accession you MUST consult BOTH:
  (a) the deposit record — ENA / SRA / DDBJ for INSDC accessions (PRJNA / PRJEB / PRJDB / SRP / ERP / DRP), or figshare / Zenodo; AND
  (b) the ASSOCIATED PUBLICATION (journal article or preprint), if one exists — locate it (via the deposit's linked references, the study title, or a literature search) and READ ITS FULL TEXT, especially the Methods and any supplementary methods, before deciding any field is not_stated.

RULES
1. Report a value ONLY if you can support it with a verbatim quote from a source you actually opened, and record that source (its URL or DOI) and where in it the quote appears.
2. "not_stated" is correct for a field genuinely absent from the sources — but ONLY after you have actually read the relevant source. For fields normally reported in the paper's Methods (target, target_class, chemistry, n_random, n_rounds, selection_format, counter_selection), you MUST have opened the associated paper's FULL TEXT (not merely the deposit record) before you may answer not_stated. If you could not find or open the paper, record that honestly in source_documents (paper tier = ABSTRACT or RECORD_ONLY) so the not_stated is explained. Never guess, infer, or fill from background knowledge.
3. Report values verbatim / free-text, in the source's own words. Do NOT map them to any predefined category list (normalization is a separate later step).
4. Output ONLY a single JSON array — one object per accession, in the same order as given. Nothing before or after it. No questions, no commentary.

FIELDS — for each, when found, output {"value","evidence_quote","source","location"} (source = the URL or DOI you used; location = section / table / figure / page):
- target — the molecule/entity the selection was performed against
- target_class — what kind of target it is, in the source's own words (e.g. "protein", "whole cells", "small molecule")
- chemistry — the library's nucleic-acid chemistry (e.g. "ssDNA", "2'-fluoro-pyrimidine RNA")
- n_random — length of the randomized region (e.g. "N40", "40 nt")
- n_rounds — number of selection rounds (e.g. "eight rounds")
- selection_format — the selection method/format named (e.g. "HT-SELEX", "Cell-SELEX", "capture-SELEX")
- counter_selection — was a counter-/negative-selection step used? Set "value" to `present` (any counter/negative selection is described), `absent` (the text explicitly says none was used), or use not_stated if the text does not indicate either way. Put the description, if any, in "note".

study_type — this is a CLASSIFICATION (there is no verbatim value in the text). Assign exactly one, using these definitions:
- aptamer_selection — in vitro selection of nucleic-acid aptamers (RNA/DNA) that bind a target. Use this whenever the deposit contains genuine selection-round reads, EVEN IF the paper frames itself as a new method or protocol (classify by the data, not the paper's framing).
- tf_ht_selex — transcription-factor HT-SELEX (mapping a DNA-binding protein's sequence motif), not aptamer discovery
- method_or_other — SELEX-related but the deposit does NOT contain a standard aptamer-selection read set: e.g. a simulation, a reanalysis of others' data, or a tool/benchmark without its own selection reads
- not_selex — not a SELEX experiment at all
Output study_type as {"value","rationale","evidence_quote","source","location"} where rationale is one sentence.

source_documents — record what you actually consulted:
- the PAPER entry's "id" is the associated publication's own DOI / PMID / preprint URL. It is NEVER the deposit's own accession or its dataset DOI. If no associated publication exists, set id="" and tier="RECORD_ONLY".
- "tier" is FULL_TEXT (you read the full text), ABSTRACT (only an abstract was available), or RECORD_ONLY (no paper found / only the deposit record).

OUTPUT — a JSON array; each element is one accession's object of this shape (keep absent fields as {"value": null, "status": "not_stated"}):
{
  "accession": "<the accession>",
  "extractor": "<claude or codex — whichever assistant you are>",
  "model": "<your model id if known, else empty>",
  "prompt_version": "extract-v1.4",
  "source_documents": [
    {"type": "paper", "id": "<publication DOI/PMID/URL, or empty if none>", "tier": "<FULL_TEXT | ABSTRACT | RECORD_ONLY>"},
    {"type": "deposit", "id": "<the accession>", "url": "<deposit record URL>"}
  ],
  "fields": {
    "study_type": {"value": null, "rationale": "", "evidence_quote": "", "source": "", "location": ""},
    "target": {"value": null, "evidence_quote": "", "source": "", "location": ""},
    "target_class": {"value": null, "evidence_quote": "", "source": "", "location": ""},
    "chemistry": {"value": null, "evidence_quote": "", "source": "", "location": ""},
    "n_random": {"value": null, "evidence_quote": "", "source": "", "location": ""},
    "n_rounds": {"value": null, "evidence_quote": "", "source": "", "location": ""},
    "selection_format": {"value": null, "evidence_quote": "", "source": "", "location": ""},
    "counter_selection": {"value": null, "note": "", "evidence_quote": "", "source": "", "location": ""}
  }
}

ACCESSIONS (one per line):
<<PASTE ACCESSIONS HERE>>
