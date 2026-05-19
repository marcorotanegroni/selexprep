# selexprep — Claude Code project context

Accession-first preprocessing for public HT-SELEX datasets with empirical
primer / constant-region inference. v0.1 ships a CLI + bundled discovery
catalog. Target venues: Bioinformatics Advances Application Note (tool) +
NAR Database Issue (catalog, v0.2).

The locked Codex-peer-reviewed implementation plan lives at
`~/.claude/plans/unified-seeking-treehouse.md` (four review passes; treat
as a contract — don't rewrite for additive extensions, just commit +
CHANGELOG entry).

## Commands

| Task                | Command                                  |
| ------------------- | ---------------------------------------- |
| Tests               | `uv run pytest`                          |
| Lint                | `uv run ruff check src/ tests/`          |
| Format              | `uv run ruff format src/ tests/`         |
| Type-check          | `uv run mypy src/`                       |
| Catalog browse      | `uv run selexprep catalog list --target X --insdc-only` |
| Catalog refresh     | `uv run selexprep catalog refresh`       |

CI runs lint + format + mypy + pytest on Python 3.10 / 3.11 / 3.12.
Always run all four locally before committing.

## Layout

```
src/selexprep/
├── _common.py         shared utilities (iter_srr_files, load_csv, …)
├── _io.py             deterministic gzip + sha256 helpers
├── cli.py             Typer dispatcher (root)
├── catalog/           bundled public-SELEX catalog + filters + CLI subapp
├── count/             per-round sequence counting
├── extract/           demux + trim (sample-sheet driven; deterministic gzip)
├── fetch/             accession discovery + download (ENA-first)
├── library/           primer detection + audit (Phase 2: LibraryReport)
└── qc/                round coverage, consistency, readiness, plots
```

## Strict-mypy boundary

Only `selexprep.library.report` is strict-mypy (Phase 2 LibraryReport
schema — type precision matters for confidence calibration). Other
modules use the permissive profile declared in `pyproject.toml`'s
`[tool.mypy] disable_error_code`.

## Critical gotchas (do not relearn)

- **macOS TCC** blocks `~/Documents` on this Mac → repo lives at
  `/Users/marcorotanegroni/selexprep`, NOT under `~/Documents/...`.
- **`detect.detect_from_parquet`** defaults to `top_n=None` (no
  subsampling). Don't reintroduce the 10_000 cap — the long tail of
  rare uniques confirms primer consensus.
- **`download_srr`** defaults to `backend="auto"` → ENA-direct first
  (MIT-licensed path). kingfisher is GPL-3.0 and only invoked as
  opt-in fallback (with a loud GPL notice). Never promote it to default.
- **All `.gz` writes** must go through
  `selexprep._io.open_gzip_text_deterministic` (mtime=0 + suppressed
  FNAME). Plain `gzip.open` breaks SHA256 reproducibility.
- **Curation flags** (`include`, `manual_curation_notes`,
  `library_type_verification`, `library_type_evidence`,
  `has_processed_counts`) are NEVER in the package catalog or in any
  user-facing output. Curation is the user's downstream job; the
  package reflects the public archives.
- **Phase 2 centerpiece** is `LibraryReport` with explicit
  `extraction_mode` × `read_source` × `required_action` matrix
  (see the plan file for the full pydantic schema and the
  cross-round-persistence inference algorithm).

## Workflow

- **Big design decisions** → user cross-checks with Codex
  (independent peer-review). Don't self-approve scope changes. Trigger
  phrase: *"let me run this through Codex"*. Wait for the response
  before applying changes.
- **Per-feature implementations** → use `EnterPlanMode` to get plan
  approval before writing code on anything > a one-file change.
- **Subagents** → `Explore` agent for read-only mapping when the
  surface area is unclear; do edits in the main context.
- **Background polling** (CI, long-running jobs) → use `Bash` with
  `run_in_background: true` + a sha-aware Python poll script at
  `/tmp/wait_for_sha_ci.py`. Don't write nested `python -c` with
  escaped quotes — they break and run forever silently.

## Where context lives

- **Plan (locked contract):** `~/.claude/plans/unified-seeking-treehouse.md`
- **Change log + v0.1 follow-ups:** `CHANGELOG.md` in this repo
- **Project memory** (durable user preferences, project state):
  `~/.claude/projects/-Users-marcorotanegroni-Documents-subtractive-proteomics-pipeline-aptameri-selex-corpus/memory/MEMORY.md`
  (the path is keyed on the original thesis directory, not on this
  repo — index file lists every memory note).
