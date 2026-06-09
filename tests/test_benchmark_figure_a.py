"""Tests for ``selexprep.benchmark.figure_a`` — the Tier 1 scorecard table emitter.

``figure_a.py`` now emits a per-deposit Markdown scorecard (Table 1), not a bar chart:
``metrics.json`` (outcomes) joined with ``ground_truth.tsv`` (descriptors).
"""

from __future__ import annotations

import json
from pathlib import Path

from selexprep.benchmark.figure_a import build_scorecard, emit_scorecard
from selexprep.benchmark.metrics import load_ground_truth


def _metrics() -> dict:
    """A two-arm metrics.json: 2 recovery (one with read-resolved 3') + 2 controls."""
    return {
        "recovery_denominator": 2,
        "pair_recovery_by_status": {
            "counts": {"HIGH": {"pair_exact": 1}, "MEDIUM": {"pair_partial": 1}}
        },
        "primer_recovery": {
            "pairs": [
                {
                    "accession": "REC1",
                    "score_3p": True,
                    "status_5p": {"equivalence_kind": "EXACT", "matched": True},
                    "status_3p": {"equivalence_kind": "EXACT", "matched": True},
                },
                {
                    "accession": "REC2",
                    "score_3p": False,
                    "status_5p": {"equivalence_kind": "EXACT", "matched": True},
                    "status_3p": {"equivalence_kind": "MISMATCH", "matched": False},
                },
            ]
        },
        "specificity": {
            "n_evaluated": 1,
            "n_false_positive": 0,
            "per_row": [
                {
                    "accession": "SPEC1",
                    "bucket": "no_false_call",
                    "primer_5p": None,
                    "primer_3p": None,
                }
            ],
        },
        "adapter_control": {
            "n_evaluated": 1,
            "n_no_false_call": 1,
            "per_row": [
                {
                    "accession": "ADPT1",
                    "bucket": "no_false_call",
                    "primer_5p": None,
                    "primer_3p": None,
                }
            ],
        },
    }


_GROUND_TRUTH = (
    "accession\tlibrary_kind\ttarget_kind\tverified\tread_state\tscore_3p\tnotes\n"
    "REC1\tDNA\tprotein\ttrue\traw_standard\t\tfirst recovery\n"
    "REC2\tRNA\tcell\ttrue\traw_standard\tfalse\tread-resolved 3'\n"
    "SPEC1\tDNA\tprotein\ttrue\tpre_trimmed\t\tpre-trimmed\n"
    "ADPT1\t2'-F-Py RNA\tprotein\ttrue\tadapter_control\t\tadapter collision\n"
    "UNVER\tDNA\tprotein\tfalse\traw_standard\t\tunverified — excluded\n"
)


def _write(tmp_path: Path) -> tuple[Path, Path]:
    metrics = tmp_path / "metrics.json"
    metrics.write_text(json.dumps(_metrics()), encoding="utf-8")
    gt = tmp_path / "ground_truth.tsv"
    gt.write_text(_GROUND_TRUTH, encoding="utf-8")
    return metrics, gt


def test_emit_scorecard_writes_table_1(tmp_path: Path) -> None:
    metrics, gt = _write(tmp_path)
    out = emit_scorecard(metrics, gt, tmp_path / "o")
    assert out.name == "table_1.md"
    text = out.read_text(encoding="utf-8")
    assert "Table 1" in text
    assert "| Accession | Chemistry | Target | Arm | 5' | 3' | Note |" in text
    # factual headline (counts only)
    assert "1 exact" in text
    assert "false-positive calls" in text
    # verified deposits present; the unverified row excluded
    for acc in ("REC1", "REC2", "SPEC1", "ADPT1"):
        assert acc in text
    assert "UNVER" not in text
    # control rows refuse; the read-resolved 3' row is scored on 5'
    assert "correct refusal" in text
    assert "3' read-resolved → scored on 5'" in text


def test_scorecard_groups_by_arm(tmp_path: Path) -> None:
    _, gt = _write(tmp_path)
    rows = build_scorecard(_metrics(), load_ground_truth(gt))
    arms = [r["arm"] for r in rows]
    rank = {"recovery": 0, "specificity": 1, "adapter-control": 2}
    assert arms == sorted(arms, key=lambda a: rank[a])
    # recovery rows carry the equivalence outcome; controls carry null/null
    by_acc = {r["accession"]: r for r in rows}
    assert by_acc["REC1"]["five"] == "EXACT"
    assert by_acc["REC2"]["three"] == "—"  # read-resolved → not scored
    assert by_acc["SPEC1"]["five"] == "null"


def test_emit_scorecard_empty_metrics_still_writes_table(tmp_path: Path) -> None:
    _, gt = _write(tmp_path)
    metrics = tmp_path / "empty.json"
    metrics.write_text("{}", encoding="utf-8")
    out = emit_scorecard(metrics, gt, tmp_path / "o")
    text = out.read_text(encoding="utf-8")
    assert "| Accession | Chemistry | Target | Arm | 5' | 3' | Note |" in text
    assert "REC1" in text  # descriptors still come from ground_truth
