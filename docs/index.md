# selexprep

**Accession-first preprocessing for public HT-SELEX datasets with primer auto-inference.**

> `selexprep` is not another aptamer-analysis suite. It is the missing accession-to-clean-library preprocessing layer, with primer inference, uncertainty reporting, random-region extraction, and reproducible outputs.

!!! note "v0.1"
    Install with `pip install selexprep`. Source and issues on [GitHub](https://github.com/marcorotanegroni/selexprep).

## Sections

- [Quick start](quickstart.md) — single accession, start to finish
- [Batch tutorial](batch.md) — many accessions with `run`
- [LibraryReport reference](library-report.md) — every field of the inference output explained
- [CLI reference](cli.md) — all eight `selexprep` subcommands
- [API reference](api.md) — `LibraryReport` + manifest schemas
- [Examples](examples.md) — runnable notebooks (toy + cached-real)
- [Limitations](limitations.md) — what v0.1 does NOT do; honest constraints
- [v0.2 roadmap](roadmap.md) — deferred features (read merging, AnnData, HTML report, …)
