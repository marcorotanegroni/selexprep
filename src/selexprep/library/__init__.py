"""Primer + library-structure inference (LibraryReport)."""

from selexprep.library.detect import compute_library_report
from selexprep.library.report import (
    ExtractionMode,
    LibraryReport,
    Orientation,
    ReadSource,
    RequiredAction,
    Status,
    read_library_report_json,
    write_library_report_json,
)

__all__ = [
    "ExtractionMode",
    "LibraryReport",
    "Orientation",
    "ReadSource",
    "RequiredAction",
    "Status",
    "compute_library_report",
    "read_library_report_json",
    "write_library_report_json",
]
