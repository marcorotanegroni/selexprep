"""Deterministic post-hoc normalization for dual-extraction metadata.

Maps verbatim extracted values onto controlled vocabularies (``normalization.yaml``)
and applies equivalence rules, so inter-extractor agreement is computed on
normalized values rather than raw wording. Unmapped categorical labels are
returned flagged ``UNMAPPED`` so they surface as a human-adjudication decision
(add them to the table) and are never silently forced into a category.

This module is intentionally dependency-light (PyYAML only) and pure: the same
``normalize`` is used by ``compute_agreement.py`` and by the final dataset build,
so normalization is identical everywhere and auditable.
"""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

import yaml

_HERE = Path(__file__).parent
_DEFAULT_TABLE = _HERE / "normalization.yaml"

CATEGORICAL = {"study_type", "chemistry", "target_class", "selection_format"}
NUMERIC = {"n_random", "n_rounds"}
FREE_TEXT = {"target", "counter_selection"}

_ABSENT = {"", "not_stated", "none", "null", "na", "n/a"}


@lru_cache(maxsize=8)
def _table(path: str) -> dict:
    return yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}


def _clean(s: object) -> str:
    """Normalize for matching: drop parentheticals, unify prime/quote glyphs,
    collapse whitespace, casefold.

    Parentheticals (``"cells (cell line)"`` -> ``"cells"``) and unicode prime /
    apostrophe glyphs (U+2032, U+2019, U+02BC -> ASCII ``'``) are removed so
    wording variants collapse. The verbatim value is preserved upstream in
    ``evidence_quote``.
    """
    t = str(s).replace("′", "'").replace("’", "'").replace("ʼ", "'")  # noqa: RUF001
    t = re.sub(r"\([^)]*\)", " ", t)
    return re.sub(r"\s+", " ", t.strip()).strip(" .;,-").casefold()


def normalize(field: str, value: object, *, table_path: str | None = None) -> tuple[object, str]:
    """Normalize one extracted value.

    Returns ``(normalized_value, flag)`` where ``flag`` is ``""`` normally,
    ``"UNMAPPED"`` for a categorical value absent from the table, or
    ``"UNPARSED"`` for a numeric value with no parseable integer. An absent /
    ``not_stated`` input normalizes to ``(None, "")``.
    """
    if value is None or _clean(value) in _ABSENT:
        return None, ""

    if field in NUMERIC:
        m = re.search(r"\d+", str(value))
        return (int(m.group()), "") if m else (None, "UNPARSED")

    if field in CATEGORICAL:
        tbl = _table(str(table_path or _DEFAULT_TABLE)).get(field, {}) or {}
        key = _clean(value)
        mapping = {_clean(k): v for k, v in (tbl.get("map") or {}).items()}
        if key in mapping:
            return mapping[key], ""
        canonical = {_clean(c): c for c in (tbl.get("canonical") or [])}
        if key in canonical:
            return canonical[key], ""
        return str(value).strip(), "UNMAPPED"

    # free-text: casefold + whitespace only (synonyms resolved at adjudication)
    return _clean(value), ""
