"""Batch driver: fetch → detect → extract → count → qc per accession."""

from selexprep.run.runner import RunReport, RunRowReport, run_batch

__all__ = ["RunReport", "RunRowReport", "run_batch"]
