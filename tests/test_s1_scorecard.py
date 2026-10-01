"""Supplementary Data S1, generated from the benchmark's data
(benchmarks/supplementary/s1_scorecard.py)."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "benchmarks" / "supplementary" / "s1_scorecard.py"


@pytest.fixture(scope="module")
def s1():
    spec = importlib.util.spec_from_file_location("s1_scorecard", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


EXCLUDED = ("PRJNA728693", "PRJNA935703", "PRJNA975735", "PRJNA315881", "PRJNA1068659")


@pytest.fixture(scope="module")
def text(s1, tmp_path_factory):
    """Legend and notes are manuscript text kept outside the repository; the
    tests use stand-ins, so they never depend on paper/."""
    d = tmp_path_factory.mktemp("paper")
    (d / "legend.md").write_text("**Legend.** Stand-in.")
    rows = [f"{acc}\tclean note for {acc}" for acc in (*EXCLUDED, "PRJNA883192", "PRJNA1244358")]
    (d / "notes.tsv").write_text("accession\tnote\n" + "\n".join(rows) + "\n")
    return s1.build(legend=d / "legend.md", notes_tsv=d / "notes.tsv")


def test_without_legend_and_notes_the_table_is_still_generated(s1, tmp_path):
    bare = s1.build(legend=tmp_path / "missing.md", notes_tsv=tmp_path / "missing.tsv")
    assert "### Recovery arm (7)" in bare and "**Legend.**" not in bare


@pytest.mark.parametrize(
    ("five", "three", "expected"),
    [
        ("EXACT", "EXACT", "both exact"),
        ("EXACT", "MISMATCH", "partial"),
        ("PARTIAL_3P", "PARTIAL_5P", "partial"),
        ("MISMATCH", "MISMATCH", "missed"),
    ],
)
def test_outcome_follows_the_pair_rule_of_the_metrics(s1, five, three, expected):
    assert s1.outcome(five, three) == expected


def test_every_scored_deposit_is_in_the_table_once(text):
    with (ROOT / "benchmarks" / "ground_truth.tsv").open() as fh:
        gt = list(csv.DictReader(fh, delimiter="\t"))
    arms = {"raw_standard", "pre_trimmed", "adapter_control"}
    scored = [r["accession"] for r in gt if r["verified"] == "true" and r["read_state"] in arms]
    assert len(scored) == 21
    for acc in scored:
        assert text.count(f"| {acc} | ") == 1, acc


def test_the_totals_are_those_of_the_metrics(text):
    assert "*v0.4.1, pre-specified:* 4 with both constants exact, 3 partial, 0 missed, of 7" in text
    assert "*0.4.5, after the fix (in-sample):* 5 with both constants exact, 2 partial" in text
    assert "Specificity: 0 calls in 8 deposits. Adapter control: 0 calls in 6 deposits." in text
    assert (
        "with it counted in (0.4.5):* 5 with both constants exact, 3 partial, 0 missed, of 8"
        in text
    )


def test_excluded_and_unscored_deposits_are_listed_apart(text):
    sections = text.split("**Excluded before inference**")
    assert len(sections) == 2
    for acc in EXCLUDED:
        assert f"clean note for {acc}" in sections[1]
    assert "PRJNA1244358" in text.split("**Sensitivity analysis")[1]
    assert "PRJNA883192" in text.split("**Not scored.**")[1]
    # No internal working notes leak into the supplementary.
    assert "v0.1" not in text and "->" not in text
