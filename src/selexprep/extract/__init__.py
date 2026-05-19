"""Demux + trim + random-region extraction."""

from selexprep.extract.runner import ExtractResult, run_extract
from selexprep.extract.trim import TrimReport

__all__ = ["ExtractResult", "TrimReport", "run_extract"]
