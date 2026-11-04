# Examples (notebooks)

Two runnable notebooks live in
[`examples/`](https://github.com/marcorotanegroni/selexprep/tree/main/examples)
(GitHub renders them directly).

### [`01_offline_toy_pipeline.ipynb`](https://github.com/marcorotanegroni/selexprep/blob/main/examples/01_offline_toy_pipeline.ipynb) — offline toy pipeline

Runs the full `detect → extract → count → qc` pipeline on **synthetic,
self-generated** data — no download, fully deterministic. It demonstrates that
the pipeline *runs* end-to-end and reproducibly. Because the constants are
planted, it is **not** evidence about inference accuracy on real reads.

### [`02_library_report_interpretation.ipynb`](https://github.com/marcorotanegroni/selexprep/blob/main/examples/02_library_report_interpretation.ipynb) — interpreting real `LibraryReport`s

Walks through **cached real `detect` outputs** (0.4.7) on four public deposits
spanning the outcome range: a `HIGH` clean recovery, a `HIGH` call whose 3′
constant is shorter than the published one, a `MEDIUM` one-round paired-end
deposit with the constants on different mates, and an `UNABLE_TO_INFER` safe
failure on pre-trimmed reads. No network needed — the reports are bundled under
`examples/data/` and `examples/regenerate_reports.sh` regenerates them.

!!! note
    Neither notebook is a performance benchmark. Recovery performance on
    paper-documented deposits is reported in the **Tier 1 scorecard** (`benchmarks/`).
