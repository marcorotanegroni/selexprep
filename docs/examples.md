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

Walks through **cached real `detect` outputs** on three public deposits spanning
the outcome range: a `HIGH` clean recovery, a `MEDIUM` partial / `FIVE_PRIME_ONLY`
case, and an `UNABLE_TO_INFER` safe failure (adapter collision). No network
needed — the reports are bundled under `examples/data/`.

!!! note
    Neither notebook is a performance benchmark. Recovery performance on
    paper-documented deposits is reported in **Figure A** (`benchmarks/`).
