"""Benchmark scaffolding for Phase 6: primer-recovery metrics + Figure A.

Public surface kept minimal — most callers should reach for the named
functions in the submodules directly.
"""

from selexprep.benchmark.equivalence import EquivalenceResult, primer_equivalent

__all__ = ["EquivalenceResult", "primer_equivalent"]
