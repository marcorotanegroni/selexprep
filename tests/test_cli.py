"""Smoke tests for the Phase 0 Typer CLI scaffold."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

from typer.testing import CliRunner

from selexprep.cli import app
from selexprep.library.adapters import reverse_complement

runner = CliRunner()


# Local copy of _synthetic_pool — keeps test_cli independent of test_detect.
def _synthetic_pool(
    primer_5p: str | None,
    primer_3p: str | None,
    n: int = 1000,
    random_len: int = 30,
) -> list[str]:
    bases = "ACGT"
    out: list[str] = []
    for i in range(n):
        rand = "".join(bases[(i * 7 + j * 13) % 4] for j in range(random_len))
        out.append((primer_5p or "") + rand + (primer_3p or ""))
    return out


def _write_fastq_gz(path: Path, seqs: list[str]) -> None:
    """Write a minimal FASTQ.gz: identical quality line per read."""
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for i, s in enumerate(seqs):
            fh.write(f"@read_{i}\n{s}\n+\n{'I' * len(s)}\n")


def test_version_flag_emits_version_string() -> None:
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert "selexprep" in result.stdout


def test_help_lists_all_subcommands() -> None:
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    # Stub processing verbs + the Phase 1.5 catalog subapp
    for cmd in ("inspect", "fetch", "detect", "extract", "count", "qc", "run", "catalog"):
        assert cmd in result.stdout


def test_inspect_emits_summary_when_ena_returns_data(tmp_path: Path) -> None:
    """Phase 4: inspect now hits ENA (mocked) and prints metadata to stdout."""
    from unittest.mock import MagicMock, patch

    import requests as _requests

    fake_resp = MagicMock(spec=_requests.Response)
    fake_resp.status_code = 200
    fake_resp.raise_for_status.return_value = None
    fake_resp.json.return_value = [
        {
            "run_accession": "SRR000000",
            "study_accession": "PRJEB00000",
            "study_title": "Mocked SELEX study",
            "library_strategy": "OTHER",
            "library_source": "OTHER",
            "read_count": "1000",
            "base_count": "75000",
            "fastq_md5": "abc",
            "fastq_bytes": "1000",
        }
    ]
    outdir = tmp_path / "out"
    with patch("selexprep.fetch.inspect.requests.get", return_value=fake_resp):
        result = runner.invoke(app, ["inspect", "SRR000000", "--outdir", str(outdir)])

    assert result.exit_code == 0, result.output
    assert "SRR000000" in result.output
    assert "Mocked SELEX study" in result.output
    assert "library_strategy" in result.output
    assert (outdir / "inspect.json").exists()


def test_inspect_exits_with_code_2_when_accession_unknown() -> None:
    """ValueError on empty ENA response -> CLI exits 2."""
    from unittest.mock import MagicMock, patch

    import requests as _requests

    fake_resp = MagicMock(spec=_requests.Response)
    fake_resp.status_code = 200
    fake_resp.raise_for_status.return_value = None
    fake_resp.json.return_value = []  # empty -> inspect_accession raises
    with patch("selexprep.fetch.inspect.requests.get", return_value=fake_resp):
        result = runner.invoke(app, ["inspect", "SRR_DOES_NOT_EXIST"])

    assert result.exit_code == 2


# ===========================================================================
# Phase 5 — `count` + `qc` CLI commands
# ===========================================================================


def _write_fasta_gz(path: Path, seqs: list[str]) -> None:
    import gzip as _gzip

    with _gzip.open(path, "wt", encoding="utf-8") as fh:
        for i, s in enumerate(seqs):
            fh.write(f">read_{i}\n{s}\n")


def test_count_produces_counts_parquet(tmp_path: Path) -> None:
    """selexprep count <extracted.fasta.gz> --round R0 --outdir OUT
    writes OUT/round_00/counts.parquet."""
    fa = tmp_path / "extracted.fasta.gz"
    _write_fasta_gz(fa, ["AAAA", "AAAA", "CCCC", "GGGG"])
    outdir = tmp_path / "out"

    result = runner.invoke(
        app,
        ["count", str(fa), "--round", "R0", "--outdir", str(outdir)],
    )
    assert result.exit_code == 0, result.output
    parquet = outdir / "round_00" / "counts.parquet"
    assert parquet.exists()


def test_count_accepts_integer_round_label(tmp_path: Path) -> None:
    fa = tmp_path / "extracted.fasta.gz"
    _write_fasta_gz(fa, ["AAAA", "CCCC"])
    outdir = tmp_path / "out"

    result = runner.invoke(app, ["count", str(fa), "--round", "1", "--outdir", str(outdir)])
    assert result.exit_code == 0, result.output
    assert (outdir / "round_01" / "counts.parquet").exists()


def test_count_rejects_invalid_round_label(tmp_path: Path) -> None:
    fa = tmp_path / "extracted.fasta.gz"
    _write_fasta_gz(fa, ["AAAA"])

    result = runner.invoke(
        app,
        ["count", str(fa), "--round", "notanumber", "--outdir", str(tmp_path / "out")],
    )
    assert result.exit_code != 0


def test_count_rejects_negative_round_label(tmp_path: Path) -> None:
    """Phase 5 Codex pass 1 NB: negative integers like ``R-1`` would
    create ``round_-1/`` directories; CLI must reject."""
    fa = tmp_path / "extracted.fasta.gz"
    _write_fasta_gz(fa, ["AAAA"])

    result = runner.invoke(
        app,
        ["count", str(fa), "--round", "R-1", "--outdir", str(tmp_path / "out")],
    )
    assert result.exit_code != 0


def test_count_rejects_fastq_input(tmp_path: Path) -> None:
    """Phase 5 Codex pass 1 BLOCKING: selexprep count accepts only the
    extracted FASTA from `selexprep extract` (primer-stripped). FASTQ
    inputs would silently mis-parse, so the CLI rejects them by default
    with a clear error pointing the user at `selexprep extract` OR the
    `--from-pretrimmed-fastq` opt-in flag."""
    import gzip as _gzip

    fq = tmp_path / "raw.fastq.gz"
    with _gzip.open(fq, "wt", encoding="utf-8") as fh:
        fh.write("@r0\nACGT\n+\nIIII\n")

    result = runner.invoke(
        app,
        ["count", str(fq), "--round", "R0", "--outdir", str(tmp_path / "out")],
    )
    assert result.exit_code == 2
    output = result.stderr or result.output
    assert "FASTQ" in output
    assert "selexprep extract" in output
    # Error must also mention the explicit opt-in escape hatch.
    assert "--from-pretrimmed-fastq" in output


def test_count_with_from_pretrimmed_fastq_succeeds(tmp_path: Path) -> None:
    """Opt-in path: --from-pretrimmed-fastq routes FASTQ input through
    count_fastq_pretrimmed and produces a counts.parquet."""
    import gzip as _gzip

    fq = tmp_path / "pretrimmed.fastq.gz"
    # Three reads, two unique sequences.
    with _gzip.open(fq, "wt", encoding="utf-8") as fh:
        fh.write("@r0\nACGTACGT\n+\nIIIIIIII\n")
        fh.write("@r1\nACGTACGT\n+\nIIIIIIII\n")
        fh.write("@r2\nGCATGCAT\n+\nIIIIIIII\n")

    outdir = tmp_path / "out"
    result = runner.invoke(
        app,
        [
            "count",
            str(fq),
            "--round",
            "R0",
            "--outdir",
            str(outdir),
            "--from-pretrimmed-fastq",
        ],
    )
    assert result.exit_code == 0, result.output
    parquet = outdir / "round_00" / "counts.parquet"
    assert parquet.exists()

    import pandas as pd

    df = pd.read_parquet(parquet)
    assert set(df["sequence"]) == {"ACGTACGT", "GCATGCAT"}
    # ACGTACGT has 2 reads, GCATGCAT has 1.
    counts_by_seq = dict(zip(df["sequence"], df["reads"], strict=True))
    assert counts_by_seq["ACGTACGT"] == 2
    assert counts_by_seq["GCATGCAT"] == 1


def test_count_from_pretrimmed_fastq_logs_unverified_warning(
    tmp_path: Path,
) -> None:
    """Opt-in path emits a loud warning that selexprep cannot verify the
    trimming state."""
    import gzip as _gzip
    import logging

    fq = tmp_path / "pretrimmed.fastq.gz"
    with _gzip.open(fq, "wt", encoding="utf-8") as fh:
        fh.write("@r0\nACGT\n+\nIIII\n")

    # CliRunner uses a separate logging context; we route the runner's
    # stderr through a handler to capture the warning.
    handler_records: list[logging.LogRecord] = []

    class _Capture(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            handler_records.append(record)

    cli_logger = logging.getLogger("selexprep.cli")
    handler = _Capture(level=logging.WARNING)
    cli_logger.addHandler(handler)
    try:
        result = runner.invoke(
            app,
            [
                "count",
                str(fq),
                "--round",
                "R0",
                "--outdir",
                str(tmp_path / "out"),
                "--from-pretrimmed-fastq",
            ],
        )
    finally:
        cli_logger.removeHandler(handler)

    assert result.exit_code == 0, result.output
    pretrimmed_warnings = [
        rec for rec in handler_records if "pre-trimmed FASTQ" in rec.getMessage()
    ]
    assert pretrimmed_warnings, "expected a 'cannot verify trimming' warning"


def test_count_from_pretrimmed_fastq_flag_rejects_fasta_input(tmp_path: Path) -> None:
    """Symmetry: --from-pretrimmed-fastq passed with a FASTA input is a
    user error (flag/extension mismatch). Rejects with exit 2."""
    fa = tmp_path / "extracted.fasta.gz"
    _write_fasta_gz(fa, ["AAAA", "CCCC"])

    result = runner.invoke(
        app,
        [
            "count",
            str(fa),
            "--round",
            "R0",
            "--outdir",
            str(tmp_path / "out"),
            "--from-pretrimmed-fastq",
        ],
    )
    assert result.exit_code == 2
    output = result.stderr or result.output
    assert "not FASTQ" in output


def test_qc_emits_flags_yaml_and_plots(tmp_path: Path) -> None:
    """End-to-end smoke: synthetic outdir -> selexprep qc -> flags.yaml + 4 PNGs."""
    import json as _json

    import pandas as _pd

    from selexprep.library.report import LibraryReport
    from selexprep.manifest import build_manifest_from_extract_result, write_manifest_json

    outdir = tmp_path / "ds"
    outdir.mkdir()
    for r in range(2):
        round_dir = outdir / f"round_{r:02d}"
        round_dir.mkdir()
        _pd.DataFrame(
            {
                "sequence": [f"s{i}" for i in range(20)],
                "reads": [(20 - i) for i in range(20)],
                "rank": list(range(1, 21)),
                "rpm": [(20 - i) * 1000.0 for i in range(20)],
            }
        ).to_parquet(round_dir / "counts.parquet", index=False)

    (outdir / "trim_reports.json").write_text(
        _json.dumps(
            [
                {
                    "cutadapt_cmd": ["cutadapt", "..."],
                    "n_in": 1000,
                    "n_out": 950,
                    "return_code": 0,
                    "output_paths": [str(outdir / f"round_{r:02d}" / "extracted.fasta.gz")],
                }
                for r in range(2)
            ]
        ),
        encoding="utf-8",
    )

    lr = LibraryReport(
        primer_5p="GGTAATACGACTCACTATAGGG",
        primer_3p="CCATGCATGCATGCATGCAT",
        variants_5p=[],
        variants_3p=[],
        known_adapter_hits={},
        extraction_mode="BOTH_PRIMERS_SINGLE_READ",
        full_insert_recovered=True,
        read_source="R1",
        required_action="NONE",
        orientation="FORWARD",
        n_length_mode=30,
        n_length_distribution={30: 100},
        n_length_confidence=1.0,
        match_rate_5p=0.95,
        match_rate_3p=0.92,
        position_consistency_5p=0.95,
        position_consistency_3p=0.92,
        read_fraction_used_for_inference=1.0,
        sampling_seed=42,
        confidence=0.85,
        status="HIGH",
        failure_reason=None,
    )
    manifest = build_manifest_from_extract_result(
        library_report=lr,
        input_paths=[],
        output_paths=[],
        accession=None,
        bioproject_id=None,
        runs=[],
        parameters={},
    )
    manifest_path = outdir / "selexprep_manifest.json"
    write_manifest_json(manifest, manifest_path)

    result = runner.invoke(app, ["qc", str(manifest_path)])
    assert result.exit_code == 0, result.output
    assert (outdir / "qc" / "flags.yaml").exists()
    assert (outdir / "qc" / "read_retention.png").exists()
    assert (outdir / "qc" / "primer_match_per_round.png").exists()
    assert (outdir / "qc" / "n_length_distribution.png").exists()
    assert (outdir / "qc" / "per_round_panel.png").exists()


def test_no_args_shows_help() -> None:
    result = runner.invoke(app, [])
    assert "Usage:" in result.stdout


# ===========================================================================
# Phase 6a — `fetch` + `run` CLI smokes (wired)
# ===========================================================================


def test_fetch_dry_run_smoke(tmp_path: Path) -> None:
    """End-to-end CLI smoke: --dry-run emits metadata + rounds.tsv with no download."""
    from unittest.mock import MagicMock, patch

    import requests as _requests

    fake_response = MagicMock(spec=_requests.Response)
    fake_response.status_code = 200
    fake_response.raise_for_status.return_value = None
    fake_response.json.return_value = [
        {
            "run_accession": "SRR_DRY",
            "study_accession": "PRJX",
            "study_title": "T",
            "library_strategy": "OTHER",
            "library_source": "OTHER",
            "library_name": "",
            "experiment_title": "",
            "sample_title": "Round 0",
            "sample_accession": "SAM_X",
            "read_count": "10",
            "base_count": "100",
            "fastq_md5": "aaa",
            "fastq_bytes": "100",
            "fastq_ftp": "ftp.ena.example/SRR_DRY.fastq.gz",
        }
    ]

    with patch("selexprep.fetch.inspect.requests.get", return_value=fake_response):
        result = runner.invoke(app, ["fetch", "PRJ", "--outdir", str(tmp_path), "--dry-run"])
    assert result.exit_code == 0, result.output
    assert (tmp_path / "fetch_metadata.json").exists()


def test_run_missing_accession_column_exits_2(tmp_path: Path) -> None:
    tsv = tmp_path / "bad.tsv"
    tsv.write_text("not_accession\nPRJ_A\n", encoding="utf-8")
    result = runner.invoke(app, ["run", str(tsv), "--outdir", str(tmp_path / "out")])
    assert result.exit_code == 2
    assert "accession" in result.output


def test_run_exits_zero_when_all_rows_fail_but_summary_written(tmp_path: Path) -> None:
    """Phase 6b.4 HPC audit fix: per-accession failures are first-class data
    captured in ``run_summary.tsv``. ``selexprep run`` is a batch driver —
    a non-zero exit when every row safely failed would conflate "the runner
    did its job and recorded refusals" (a normal operational outcome on
    a noisy public corpus — the whole point of the Tier 2 audit) with
    "the runner itself crashed" (already handled separately via exit 2).
    The audit Snakefile's ``rule run_corpus`` runs under ``set -e``; a
    non-zero exit here would abort the pipeline before aggregate_audit
    could even read the summary.
    """
    from unittest.mock import patch

    from selexprep.run.runner import RunReport, RunRowReport

    tsv = tmp_path / "accs.tsv"
    tsv.write_text("accession\tnotes\nPRJ_A\t\nPRJ_B\t\n", encoding="utf-8")
    out = tmp_path / "out"

    # Stub run_batch so we get a deterministic all-failed report without
    # touching the network. The summary TSV path is set; the runner
    # otherwise behaves exactly as it would on a real all-FETCH_REFUSED
    # sample.
    summary_path = out / "run_summary.tsv"
    out.mkdir(parents=True, exist_ok=True)
    summary_path.write_text("accession\tstatus\n", encoding="utf-8")
    fake_report = RunReport(
        accessions_tsv=tsv,
        outdir=out,
        rows=[
            RunRowReport(
                accession="PRJ_A",
                status="FETCH_REFUSED",
                last_stage_completed="fetch",
                notes="all NONE-confidence rounds",
            ),
            RunRowReport(
                accession="PRJ_B",
                status="FETCH_REFUSED",
                last_stage_completed="fetch",
                notes="all NONE-confidence rounds",
            ),
        ],
        summary_tsv=summary_path,
    )

    with patch("selexprep.run.run_batch", return_value=fake_report):
        result = runner.invoke(app, ["run", str(tsv), "--outdir", str(out)])

    assert result.exit_code == 0, result.output
    # The summary path is still announced to the user.
    assert "run_summary.tsv" in result.output
    # The per-row failure surface is still emitted for human review.
    assert "FETCH_REFUSED" in result.output


def test_run_exits_2_on_empty_accessions_tsv(tmp_path: Path) -> None:
    """An empty input TSV (zero parseable accessions) is operator error,
    not an audit outcome. Distinct from the all-rows-failed case above.
    """
    tsv = tmp_path / "empty.tsv"
    tsv.write_text("accession\tnotes\n", encoding="utf-8")
    result = runner.invoke(app, ["run", str(tsv), "--outdir", str(tmp_path / "out")])
    assert result.exit_code == 2
    assert "zero rows" in result.output or "zero rows" in (result.stderr or "")


def test_run_command_registers_resume_and_stop_on_error_options() -> None:
    """``selexprep run`` exposes ``--resume`` and ``--stop-on-error``.

    Asserts via Click command introspection rather than grepping the
    rendered ``--help`` text — Rich help tables get truncated in CI's
    80-column non-TTY (Ubuntu runner default), which silently elides
    option names. The contract we care about is "the option exists";
    the rendered help is matplotlib-style informational and varies by
    terminal width.
    """
    import typer.main

    click_app = typer.main.get_command(app)
    run_cmd = click_app.commands["run"]  # type: ignore[attr-defined]
    flag_names: set[str] = set()
    for param in run_cmd.params:
        for opt in getattr(param, "opts", []):
            flag_names.add(opt)
    assert "--resume" in flag_names
    assert "--stop-on-error" in flag_names

    # Smoke: --help still exits cleanly (the rendering itself works,
    # we just don't grep its text).
    result = runner.invoke(app, ["run", "--help"])
    assert result.exit_code == 0


# ===========================================================================
# Phase 2 — `detect` CLI command (Wired)
# ===========================================================================


def test_detect_without_round_map_exits_with_code_2(tmp_path: Path) -> None:
    fq = tmp_path / "round0.fastq.gz"
    _write_fastq_gz(fq, _synthetic_pool("GGTAATACGACTCACTATAGGG", "CCATGCATGCATGCATGCAT"))
    result = runner.invoke(app, ["detect", str(fq), "--outdir", str(tmp_path / "out")])
    assert result.exit_code == 2
    assert "round-map" in result.stderr or "round-map" in result.output


def test_detect_with_round_map_emits_library_report(tmp_path: Path) -> None:
    primer_5p = "GGTAATACGACTCACTATAGGG"
    primer_3p = "CCATGCATGCATGCATGCAT"
    seqs = _synthetic_pool(primer_5p, primer_3p, n=1000)

    fq = tmp_path / "round0.fastq.gz"
    _write_fastq_gz(fq, seqs)

    round_map = tmp_path / "rounds.tsv"
    round_map.write_text("file\tround_number\nround0.fastq.gz\t0\n", encoding="utf-8")

    outdir = tmp_path / "out"
    result = runner.invoke(
        app, ["detect", str(fq), "--round-map", str(round_map), "--outdir", str(outdir)]
    )
    assert result.exit_code == 0, result.output

    out_path = outdir / "library_report.json"
    assert out_path.exists()
    payload = json.loads(out_path.read_text())
    assert payload["primer_5p"] is not None
    assert payload["primer_3p"] is not None
    assert payload["extraction_mode"] == "BOTH_PRIMERS_SINGLE_READ"


def test_detect_with_paired_r2_threads_split_primer_stream(tmp_path: Path) -> None:
    primer_5p = "GGTAATACGACTCACTATAGGG"
    # PRJNA883192 regression: its 3' library constant is only 14 nt.
    # The paired R2 stream must still recover it instead of falling back
    # to a spurious R1 3' suffix.
    primer_3p = "GAGCTCTGAACTGG"
    r1 = tmp_path / "round0_1.fastq.gz"
    r2 = tmp_path / "round0_2.fastq.gz"
    _write_fastq_gz(r1, _synthetic_pool(primer_5p, primer_3p=None, n=1000, random_len=80))
    _write_fastq_gz(
        r2,
        _synthetic_pool(reverse_complement(primer_3p), primer_3p=None, n=1000, random_len=80),
    )

    round_map = tmp_path / "rounds.tsv"
    round_map.write_text(
        "file\tround_number\nround0_1.fastq.gz\t0\nround0_2.fastq.gz\t0\n",
        encoding="utf-8",
    )

    outdir = tmp_path / "out"
    result = runner.invoke(
        app,
        [
            "detect",
            str(r1),
            "--paired-r2",
            str(r2),
            "--round-map",
            str(round_map),
            "--outdir",
            str(outdir),
        ],
    )
    assert result.exit_code == 0, result.output

    payload = json.loads((outdir / "library_report.json").read_text())
    assert payload["read_source"] == "R1_AND_R2"
    assert payload["extraction_mode"] == "PAIRED_END_SPLIT_PRIMERS"
    assert payload["primer_5p"] == primer_5p
    assert payload["primer_3p"] == primer_3p


def test_detect_fastq_not_in_round_map_exits_with_code_2(tmp_path: Path) -> None:
    fq = tmp_path / "round0.fastq.gz"
    _write_fastq_gz(fq, _synthetic_pool("GGTAATACGACTCACTATAGGG", "CCATGCATGCATGCATGCAT"))

    round_map = tmp_path / "rounds.tsv"
    # Map references a different filename than the one we pass.
    round_map.write_text("file\tround_number\nelsewhere.fastq.gz\t0\n", encoding="utf-8")

    outdir = tmp_path / "out"
    result = runner.invoke(
        app, ["detect", str(fq), "--round-map", str(round_map), "--outdir", str(outdir)]
    )
    assert result.exit_code == 2


# ===========================================================================
# Phase 3 - `extract` CLI command
# ===========================================================================


def _write_library_report_json(path: Path, **overrides: object) -> None:
    """Hand-write a minimal LibraryReport JSON for CLI tests."""
    import json as _json

    base = {
        "primer_5p": "GGTAATACGACTCACTATAGGG",
        "primer_3p": "CCATGCATGCATGCATGCAT",
        "variants_5p": [],
        "variants_3p": [],
        "known_adapter_hits": {},
        "extraction_mode": "BOTH_PRIMERS_SINGLE_READ",
        "full_insert_recovered": True,
        "read_source": "R1",
        "required_action": "NONE",
        "orientation": "FORWARD",
        "n_length_mode": 30,
        "n_length_distribution": {"30": 100},
        "n_length_confidence": 1.0,
        "match_rate_5p": 0.95,
        "match_rate_3p": 0.95,
        "position_consistency_5p": 0.95,
        "position_consistency_3p": 0.95,
        "read_fraction_used_for_inference": 1.0,
        "sampling_seed": 42,
        "confidence": 0.85,
        "status": "HIGH",
        "failure_reason": None,
    }
    base.update(overrides)
    path.write_text(_json.dumps(base, indent=2), encoding="utf-8")


def test_extract_without_round_map_or_sample_sheet_errors(tmp_path: Path) -> None:
    fq = tmp_path / "r0.fastq.gz"
    _write_fastq_gz(fq, _synthetic_pool("GGTAATACGACTCACTATAGGG", "CCATGCATGCATGCATGCAT"))
    lr_path = tmp_path / "lr.json"
    _write_library_report_json(lr_path)

    result = runner.invoke(
        app,
        [
            "extract",
            str(fq),
            "--library-report",
            str(lr_path),
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code == 2


def test_extract_override_primer_writes_overridden_subtree(tmp_path: Path) -> None:
    """Phase 4: --override-primer-5p without --rebuild writes to outdir/overridden/."""
    import shutil as _shutil

    if _shutil.which("cutadapt") is None:
        return  # cutadapt not available; this test needs a real trim run

    primer_5p = "GGTAATACGACTCACTATAGGG"
    primer_3p = "CCATGCATGCATGCATGCAT"
    fq = tmp_path / "r0.fastq.gz"
    _write_fastq_gz(fq, _synthetic_pool(primer_5p, primer_3p, n=500))

    lr_path = tmp_path / "lr.json"
    _write_library_report_json(lr_path)
    rm = tmp_path / "rounds.tsv"
    rm.write_text(f"file\tround_number\n{fq.name}\t0\n", encoding="utf-8")
    outdir = tmp_path / "out"

    # Override the 5' primer (still using the same FASTQ; cutadapt will
    # discard untrimmed reads since the new "primer" is not present).
    result = runner.invoke(
        app,
        [
            "extract",
            str(fq),
            "--library-report",
            str(lr_path),
            "--round-map",
            str(rm),
            "--override-primer-5p",
            "AAAAAAAAAAAAAAA",  # 15 nt junk primer to force divergence
            "--outdir",
            str(outdir),
        ],
    )
    assert result.exit_code == 0, result.output
    # Override outputs land in the overridden/ subtree (baseline outdir is untouched).
    assert (outdir / "overridden" / "selexprep_manifest.json").exists()
    assert (outdir / "overridden" / "round_00" / "extracted.fasta.gz").exists()


def test_extract_refuses_unable_library_report(tmp_path: Path) -> None:
    fq = tmp_path / "r0.fastq.gz"
    _write_fastq_gz(fq, _synthetic_pool("GGTAATACGACTCACTATAGGG", "CCATGCATGCATGCATGCAT"))
    lr_path = tmp_path / "lr.json"
    _write_library_report_json(
        lr_path,
        status="UNABLE_TO_INFER",
        extraction_mode="UNABLE_TO_EXTRACT",
        failure_reason="synthetic test failure",
    )
    rm = tmp_path / "rounds.tsv"
    rm.write_text(f"file\tround_number\n{fq.name}\t0\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "extract",
            str(fq),
            "--library-report",
            str(lr_path),
            "--round-map",
            str(rm),
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code == 2


def test_extract_end_to_end_both_primers(tmp_path: Path) -> None:
    """Smoke test: real cutadapt run via the CLI -> deterministic FASTA on disk."""
    import shutil as _shutil

    if _shutil.which("cutadapt") is None:
        return  # silently skip — covered by @pytest.mark.skipif in trim/runner tests

    fq = tmp_path / "r0.fastq.gz"
    _write_fastq_gz(fq, _synthetic_pool("GGTAATACGACTCACTATAGGG", "CCATGCATGCATGCATGCAT"))
    lr_path = tmp_path / "lr.json"
    _write_library_report_json(lr_path)
    rm = tmp_path / "rounds.tsv"
    rm.write_text(f"file\tround_number\n{fq.name}\t0\n", encoding="utf-8")
    outdir = tmp_path / "out"

    result = runner.invoke(
        app,
        [
            "extract",
            str(fq),
            "--library-report",
            str(lr_path),
            "--round-map",
            str(rm),
            "--outdir",
            str(outdir),
        ],
    )
    assert result.exit_code == 0, result.output
    assert (outdir / "round_00" / "extracted.fasta.gz").exists()
    assert (outdir / "trim_reports.json").exists()
