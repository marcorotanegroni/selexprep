"""Tests for ``benchmarks.build_project_metadata`` — the static SELEX metadata builder.

Covers the merge precedence (verified > extracted > catalog), the two honesty
signals (``curation_level`` / ``metadata_tier``), and an integrity check over the
bundled curated TSVs so a stray hand-edit can't silently ship a bad row.
"""

from __future__ import annotations

import csv
from pathlib import Path

from benchmarks import build_project_metadata as bpm

_STUDY_TYPES = {"aptamer_selection", "tf_ht_selex", "method_or_other", "not_selex"}


def test_metadata_tier_ladder() -> None:
    assert bpm.metadata_tier("", None) == "RECORD_ONLY"
    assert bpm.metadata_tier("10.x/y", None) == "ABSTRACT"  # paper, OA unknown
    assert bpm.metadata_tier("10.x/y", False) == "ABSTRACT"  # paper, not OA
    assert bpm.metadata_tier("10.x/y", True) == "FULL_TEXT"  # paper, OA


def test_build_rows_precedence_and_signals() -> None:
    catalog = {
        "PRJV": {
            "study_title": "verified study",
            "protein_target": "",
            "target_organism": "Homo sapiens",
            "n_rounds_declared": "",
            "paper_doi": "10.cat/verified",
            "paper_pmid": "",
        },
        "PRJE": {
            "study_title": "extracted study",
            "protein_target": "catalog target",
            "target_organism": "",
            "n_rounds_declared": "8",
            "paper_doi": "",
            "paper_pmid": "",
        },
        "PRJN": {
            "study_title": "bare study",
            "protein_target": "",
            "target_organism": "",
            "n_rounds_declared": "",
            "paper_doi": "",
            "paper_pmid": "",
        },
    }
    descriptors = {
        "PRJV": {
            "chemistry": "DNA",
            "target_class": "protein",
            "n_random": "40",
            "paper_doi": "10.gt/verified",
            "paper_pmid": "999",
        },
    }
    project_ann = {"PRJV": {"study_type": "aptamer_selection", "target": "verified target"}}
    catalog_ann = {"PRJE": {"study_type": "tf_ht_selex", "paper_doi": "10.ann/extracted"}}
    oa = {"10.gt/verified": (True, "http://oa"), "10.ann/extracted": (False, None)}

    rows = {
        r["accession"]: r
        for r in bpm.build_rows(
            catalog=catalog,
            descriptors=descriptors,
            project_ann=project_ann,
            catalog_ann=catalog_ann,
            oa=oa,
        )
    }

    # verified: descriptors + project_ann win; ground_truth DOI beats catalog DOI; OA -> FULL_TEXT
    v = rows["PRJV"]
    assert v["curation_level"] == "verified"
    assert v["target"] == "verified target" and v["chemistry"] == "DNA"
    assert v["paper_doi"] == "10.gt/verified" and v["metadata_tier"] == "FULL_TEXT"

    # extracted: annotation DOI overrides empty catalog DOI; not OA -> ABSTRACT; n_rounds from catalog
    e = rows["PRJE"]
    assert e["curation_level"] == "extracted"
    assert e["paper_doi"] == "10.ann/extracted" and e["metadata_tier"] == "ABSTRACT"
    assert e["n_rounds"] == "8" and e["target"] == "catalog target"  # catalog target as fallback

    # none: no annotation, no DOI -> RECORD_ONLY
    n = rows["PRJN"]
    assert n["curation_level"] == "none" and n["metadata_tier"] == "RECORD_ONLY"


def test_load_annotations_reads_discovery_doi(tmp_path: Path) -> None:
    path = tmp_path / "ann.tsv"
    path.write_text(
        "accession\tstudy_type\ttarget\tselection_format\tn_rounds\tpaper_doi\n"
        "zenodo:1\taptamer_selection\tThrombin\t\t\t10.foo/bar\n",
        encoding="utf-8",
    )
    ann = bpm.load_annotations(path)
    assert ann["zenodo:1"]["paper_doi"] == "10.foo/bar"
    assert ann["zenodo:1"]["target"] == "Thrombin"


def test_write_csv_and_json_round_trip(tmp_path: Path) -> None:
    rows = bpm.build_rows(
        catalog={
            "PRJX": {
                "study_title": "t",
                "protein_target": "",
                "target_organism": "",
                "n_rounds_declared": "",
                "paper_doi": "",
                "paper_pmid": "",
            }
        },
        descriptors={},
        project_ann={},
        catalog_ann={},
        oa={},
    )
    csv_path, json_path = tmp_path / "m.csv", tmp_path / "m.json"
    bpm.write_csv(rows, csv_path)
    bpm.write_json(rows, json_path)
    header = csv_path.read_text(encoding="utf-8").splitlines()[0].split(",")
    assert header == bpm.COLUMNS
    assert "false" in csv_path.read_text(encoding="utf-8")  # paper_linked bool rendered
    import json

    assert json.loads(json_path.read_text(encoding="utf-8"))[0]["rounds"] is None


def test_bundled_curated_files_are_valid() -> None:
    """The shipped curated TSVs must keep their column shape and study_type vocabulary."""
    with bpm.DEFAULT_CATALOG_ANNOTATIONS.open(encoding="utf-8", newline="") as fh:
        cat = list(csv.DictReader(fh, delimiter="\t"))
    assert cat, "catalog_annotations.tsv is empty"
    for r in cat:
        assert r["study_type"] in _STUDY_TYPES, (
            f"{r['accession']}: bad study_type {r['study_type']!r}"
        )
        assert r["accession"], "blank accession row"

    with bpm.DEFAULT_ANNOTATIONS.open(encoding="utf-8", newline="") as fh:
        proj = list(csv.DictReader(fh, delimiter="\t"))
    assert all(r["study_type"] == "aptamer_selection" for r in proj), (
        "verified rows are all aptamer selections"
    )


def _write_round(acc_dir: Path, round_name: str, reads: list[int]) -> None:
    import pandas as pd

    d = acc_dir / round_name
    d.mkdir(parents=True)
    pd.DataFrame({"sequence": [f"S{i}" for i in range(len(reads))], "reads": reads}).to_parquet(
        d / "counts.parquet", index=False
    )


def test_load_trajectory_from_run_outputs(tmp_path: Path) -> None:
    """Trajectory recomputed from round_*/counts.parquet (the `selexprep run` layout)."""
    acc_dir = tmp_path / "PRJTEST"
    _write_round(acc_dir, "round_00", [1, 1, 2])  # 3 unique, 4 reads, 2 singletons
    _write_round(acc_dir, "round_01", [10, 5])  # 2 unique, 15 reads, 0 singletons
    _write_round(acc_dir, "round_unknown", [3])  # run skips this dir; loader handles it if present

    traj = bpm.load_trajectory(tmp_path, "PRJTEST")
    assert traj is not None
    assert [r["round"] for r in traj] == [0, 1, "unknown"]  # numeric ascending, label last
    assert traj[0]["n_reads"] == 4 and traj[0]["n_unique"] == 3
    assert traj[0]["singleton_frac"] == 2 / 3
    assert traj[1]["singleton_frac"] == 0.0


def test_load_trajectory_absent_is_none(tmp_path: Path) -> None:
    assert bpm.load_trajectory(tmp_path, "NOT_RUN") is None  # no per-round counts -> rounds: null


def test_load_trajectory_flat_layout(tmp_path: Path) -> None:
    """The standalone counter's flat ``round_NN.counts.parquet`` is also accepted."""
    import pandas as pd

    acc_dir = tmp_path / "PRJFLAT"
    acc_dir.mkdir()
    for label, reads in (("round_00", [1, 1]), ("round_01", [4])):
        pd.DataFrame({"sequence": [f"S{i}" for i in range(len(reads))], "reads": reads}).to_parquet(
            acc_dir / f"{label}.counts.parquet", index=False
        )
    traj = bpm.load_trajectory(tmp_path, "PRJFLAT")
    assert traj is not None
    assert [r["round"] for r in traj] == [0, 1]
    assert traj[0]["n_reads"] == 2 and traj[0]["singleton_frac"] == 1.0


def test_load_trajectory_empty_round(tmp_path: Path) -> None:
    """A present-but-empty parquet yields a zero-read round, not a crash."""
    import pandas as pd

    acc_dir = tmp_path / "PRJEMPTY"
    d = acc_dir / "round_00"
    d.mkdir(parents=True)
    pd.DataFrame({"sequence": [], "reads": []}).to_parquet(d / "counts.parquet", index=False)
    traj = bpm.load_trajectory(tmp_path, "PRJEMPTY")
    assert traj == [{"round": 0, "n_reads": 0, "n_unique": 0, "singleton_frac": 0.0}]


def test_write_json_preserves_attached_trajectory(tmp_path: Path) -> None:
    rows = bpm.build_rows(
        catalog={
            "PRJX": dict.fromkeys(
                (
                    "study_title",
                    "protein_target",
                    "target_organism",
                    "n_rounds_declared",
                    "paper_doi",
                    "paper_pmid",
                ),
                "",
            )
        },
        descriptors={},
        project_ann={},
        catalog_ann={},
        oa={},
    )
    rows[0]["rounds"] = [{"round": 0, "n_reads": 4, "n_unique": 3, "singleton_frac": 0.5}]
    import json

    path = tmp_path / "m.json"
    bpm.write_json(rows, path)
    assert json.loads(path.read_text(encoding="utf-8"))[0]["rounds"][0]["n_reads"] == 4


def test_round_structure_categories() -> None:
    rs = bpm.round_structure
    assert rs("PRJX", [{"round": 0}, {"round": 1}], "OK") == "multi"
    assert rs("PRJX", [{"round": 0}], "OK") == "mono"
    assert rs("PRJX", None, "FETCH_REFUSED") == "unassigned"  # INSDC the run refused
    assert rs("figshare:1", None, "") == "not_fetchable"
    assert rs("zenodo:9", None, "") == "not_fetchable"
    assert rs("PRJX", None, "") == ""  # INSDC not yet counted


def test_load_run_status(tmp_path: Path) -> None:
    (tmp_path / "run_summary.tsv").write_text(
        "accession\tstatus\tnotes\nPRJA\tOK\t\nPRJB\tFETCH_REFUSED\tno round\n", encoding="utf-8"
    )
    assert bpm.load_run_status(tmp_path) == {"PRJA": "OK", "PRJB": "FETCH_REFUSED"}
    assert bpm.load_run_status(tmp_path / "nope") == {}  # absent summary -> empty


def test_build_rows_sets_round_structure_default() -> None:
    """build_rows seeds round_structure (not_fetchable for discovery, '' for INSDC)."""
    catalog = {
        k: dict.fromkeys(
            (
                "study_title",
                "protein_target",
                "target_organism",
                "n_rounds_declared",
                "paper_doi",
                "paper_pmid",
            ),
            "",
        )
        for k in ("PRJX", "zenodo:1")
    }
    rows = {
        r["accession"]: r
        for r in bpm.build_rows(
            catalog=catalog, descriptors={}, project_ann={}, catalog_ann={}, oa={}
        )
    }
    assert rows["PRJX"]["round_structure"] == ""
    assert rows["zenodo:1"]["round_structure"] == "not_fetchable"
