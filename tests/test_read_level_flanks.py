from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any

from benchmarks import read_level_flanks


class _GitResult:
    def __init__(self, stdout: str = "", returncode: int = 0) -> None:
        self.stdout = stdout
        self.returncode = returncode


def _write_fastq(path: Path, seqs: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        for i, seq in enumerate(seqs):
            fh.write(f"@read{i}\n{seq}\n+\n{'I' * len(seq)}\n")


def test_read_level_flanks_can_generate_r2_provenance(tmp_path: Path) -> None:
    """R2-derived truth must be reproducible by the same recipe as R1 truth."""
    accession = "PRJNA883192"
    results_dir = tmp_path / "results"
    acc_dir = results_dir / accession
    r1 = acc_dir / "round_unknown" / "SRR1_1.fastq.gz"
    r2 = acc_dir / "round_unknown" / "SRR1_2.fastq.gz"

    _write_fastq(r1, ["ATGCCATCCTACCAACAAAAA"] * 4)
    r2_prefix = "TCGAGTTCAGAGCTC"
    _write_fastq(
        r2,
        [
            r2_prefix + "AAAAA",
            r2_prefix + "CCCCC",
            r2_prefix + "GGGGG",
            r2_prefix + "TTTTT",
        ],
    )
    (acc_dir / "fastqs.manifest").write_text(f"{r1}\n")
    (acc_dir / "fastqs.r2.manifest").write_text(f"{r2}\n")

    result = read_level_flanks.analyze(
        results_dir,
        accession,
        sample=0,
        flank_len=20,
        threshold=0.6,
        manifest_name="fastqs.r2.manifest",
    )

    assert result is not None
    assert result["fastqs"] == r2.name
    assert result["consensus_5p_flank"] == r2_prefix

    out = tmp_path / "read_level_truth_provenance.tsv"
    read_level_flanks._append_provenance(out, result)
    text = out.read_text()
    assert "R2 mates" in text
    assert "insert 3' = revcomp(this 5' consensus)" in text


def test_git_commit_dirty_flag_is_scoped_to_provenance_script(monkeypatch: Any) -> None:
    """Untracked HPC artifacts must not mark the provenance script as dirty."""

    def fake_run(cmd: list[str], **_: Any) -> _GitResult:
        assert "status" not in cmd
        if cmd[:3] == ["git", "rev-parse", "--short"]:
            return _GitResult(stdout="50d3ca2\n")
        if cmd[:4] == ["git", "diff", "--quiet", "--"]:
            return _GitResult(returncode=0)
        if cmd[:5] == ["git", "diff", "--cached", "--quiet", "--"]:
            return _GitResult(returncode=0)
        raise AssertionError(f"unexpected git command: {cmd}")

    monkeypatch.setattr(read_level_flanks.subprocess, "run", fake_run)

    assert read_level_flanks._git_commit() == "50d3ca2"


def test_git_commit_marks_script_edits_dirty(monkeypatch: Any) -> None:
    def fake_run(cmd: list[str], **_: Any) -> _GitResult:
        if cmd[:3] == ["git", "rev-parse", "--short"]:
            return _GitResult(stdout="50d3ca2\n")
        if cmd[:4] == ["git", "diff", "--quiet", "--"]:
            return _GitResult(returncode=1)
        if cmd[:5] == ["git", "diff", "--cached", "--quiet", "--"]:
            return _GitResult(returncode=0)
        raise AssertionError(f"unexpected git command: {cmd}")

    monkeypatch.setattr(read_level_flanks.subprocess, "run", fake_run)

    assert read_level_flanks._git_commit() == "50d3ca2-dirty"
