"""Benchmark scaffolding for primer-recovery metrics + result tables.

Public surface kept minimal — most callers should reach for the named
functions in the submodules directly.
"""

from selexprep.benchmark.corpus_audit import (
    CorpusAuditReport,
    aggregate_audit_from_run_outputs,
    sample_corpus,
)
from selexprep.benchmark.equivalence import EquivalenceResult, primer_equivalent

__all__ = [
    "CorpusAuditReport",
    "EquivalenceResult",
    "aggregate_audit_from_run_outputs",
    "primer_equivalent",
    "sample_corpus",
]
