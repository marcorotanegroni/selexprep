"""Generator, evaluator and harness of the semi-synthetic benchmark
(benchmarks/synthetic/DESIGN.md; amendments 2 and 3).

The donor pools are fake here (random reads written to a temporary directory
under the real donors' names), so these tests pin the mechanics — coordinates,
flags, edits, metrics, outcomes — not any benchmark result.
"""

from __future__ import annotations

import gzip
import importlib.util
import json
import math
import random
import shutil
import sys
from pathlib import Path

import pytest

SYN = Path(__file__).resolve().parent.parent / "benchmarks" / "synthetic"
DONOR = "PRJNA360902"  # dev donor, N20
FAMILY = "PRJEB62495"  # dev family, 20/20


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SYN / f"{name}.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


sim = _load("simulate")
ev = _load("evaluate")


def _fake_pools(root: Path, size: int = 400) -> Path:
    donors = sim.load_donors()
    d = donors[DONOR]
    for role, (rnd, run) in d["rounds"].items():
        for seed in sim.DEV_SEEDS:
            rng = random.Random(f"{role}{seed}")
            reads = ["".join(rng.choice("ACGT") for _ in range(d["n"])) for _ in range(size)]
            path = sim.pool_path(root, DONOR, role, rnd, run, seed)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("\n".join(reads) + "\n")
    return root


@pytest.fixture(scope="module")
def pools(tmp_path_factory):
    return _fake_pools(tmp_path_factory.mktemp("pools"))


def _fastq(path: Path) -> list[tuple[str, str]]:
    with gzip.open(path, "rt") as fh:
        lines = fh.read().splitlines()
    return [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines), 4)]


def _build(config: str, pools: Path, out: Path, family: str | None = FAMILY, seed: int = 1):
    cell = sim.Cell(config, DONOR, None if config.startswith("NC") else family, seed)
    truth = sim.generate(cell, pools, out)
    reads = {r["role"]: _fastq(out / r["fastq"]) for r in truth["rounds"]}
    return truth, reads


def _pool(pools: Path, role: str, seed: int = 1) -> list[str]:
    rnd, run = sim.load_donors()[DONOR]["rounds"][role]
    return sim.load_pool(sim.pool_path(pools, DONOR, role, rnd, run, seed))


# ---------------------------------------------------------------------------
# Generator
# ---------------------------------------------------------------------------


def test_run_counts_match_the_design():
    """DESIGN.md: 1,200 development runs; amendment 3: 300 per test donor, 1,500."""
    assert len(sim.cells("dev")) == 1200
    assert len(sim.cells("test")) == 1500


def test_every_configuration_has_a_specification():
    listed = {r["config_id"] for r in sim._tsv("configurations.tsv")}
    assert listed == set(sim.SPECS)
    assert len(listed) == 52


def test_error_free_reads_carry_n_exactly_where_the_name_says(pools, tmp_path):
    truth, reads = _build("E0", pools, tmp_path)
    c5, c3 = sim.load_constants()[FAMILY]
    for role, pool in (("earliest", _pool(pools, "earliest")), ("latest", _pool(pools, "latest"))):
        for (name, seq), n in zip(reads[role], pool, strict=True):
            r = ev.parse_name(name)
            assert seq[r.start : r.end] == n
            assert seq[: r.start] == c5 and seq[r.end :] == c3
            assert (r.two, r.five, r.three, r.flagged) == (True, True, True, False)
    assert truth["oracle"] == {
        "mode": "BOTH_PRIMERS_SINGLE_READ",
        "primer_5p": c5,
        "primer_3p": c3,
    }


def test_reverse_orientation_reverses_reads_truth_and_oracle(pools, tmp_path):
    # O keeps the 1% baseline errors; compare the error-free parts through the truth.
    truth, reads = _build("O", pools, tmp_path / "o")
    c5, c3 = sim.load_constants()[FAMILY]
    assert truth["oracle"]["primer_5p"] == sim.revcomp(c3)
    assert truth["oracle"]["primer_3p"] == sim.revcomp(c5)
    assert set(truth["rounds"][0]["pre"]) == {sim.revcomp(c3)}
    pool = _pool(pools, "earliest")
    errors = 0
    for (name, seq), n in zip(reads["earliest"], pool, strict=True):
        r = ev.parse_name(name)
        errors += seq[r.start : r.end] != sim.revcomp(n)
        assert r.end - r.start == len(n)
    assert errors < 0.5 * len(pool)  # 20 nt at ~1% error


def test_truncation_cuts_inside_the_three_prime_constant(pools, tmp_path):
    _, reads = _build("F_trunc30", pools, tmp_path)
    c3 = sim.load_constants()[FAMILY][1]
    cut = 0
    for name, seq in reads["earliest"]:
        r = ev.parse_name(name)
        tail = len(seq) - r.end
        if not r.two:
            cut += 1
            assert r.five and not r.three
            assert 1 <= tail <= len(c3) - 1
        else:
            assert tail == len(c3)
    assert 0.2 < cut / len(reads["earliest"]) < 0.4


@pytest.mark.parametrize(
    ("config", "mode", "flags"),
    [
        ("F_no3", "FIVE_PRIME_ONLY", (False, True, False)),
        ("F_no5", "THREE_PRIME_ONLY", (False, False, True)),
    ],
)
def test_a_missing_flank_is_declared_in_the_truth(pools, tmp_path, config, mode, flags):
    truth, reads = _build(config, pools, tmp_path)
    assert truth["oracle"]["mode"] == mode
    r = ev.parse_name(reads["earliest"][0][0])
    assert (r.two, r.five, r.three) == flags


def test_t4_changes_one_base_in_the_earliest_round_and_the_oracle_leaves_it_out(pools, tmp_path):
    truth, _ = _build("T4", pools, tmp_path)
    c3 = sim.load_constants()[FAMILY][1]
    earliest, latest = truth["rounds"][0], truth["rounds"][-1]
    (post_early,) = earliest["post"]
    assert post_early[:-1] == c3[:-1] and post_early[-1] != c3[-1]
    assert set(latest["post"]) == {c3}
    assert truth["oracle"]["primer_3p"] == c3[:-1]


def test_the_clone_share_is_exact_and_a_surplus_is_reduced():
    rng = random.Random(0)
    pool = ["A" * 20] * 800 + ["C" * 20] * 200
    shared, before = sim._set_clone_share(pool, "A" * 20, 0.5, rng)
    assert before == 0.8
    assert shared.count("A" * 20) == 500
    raised, _ = sim._set_clone_share(pool, "A" * 20, 0.9, rng)
    assert raised.count("A" * 20) == math.ceil(0.9 * 1000)


def test_dominant_clone_is_applied_to_every_provided_round(pools, tmp_path):
    truth, _ = _build("C70", pools, tmp_path)
    clone = truth["record"]["clone"]
    assert clone in _pool(pools, "earliest")
    for role in ("earliest", "middle", "latest"):
        assert truth["record"]["clone_share"][role] == pytest.approx(0.7)


def test_the_k_spacer_is_fixed_across_reads_and_rounds(pools, tmp_path):
    truth, _ = _build("K_5outer", pools, tmp_path)
    c5 = sim.load_constants()[FAMILY][0]
    spacer = truth["record"]["spacer"]
    for rnd in truth["rounds"]:
        assert set(rnd["pre"]) == {c5[-sim.CORE_LEN :] + spacer + c5}


def test_core_copy_in_n_only_in_later_rounds(pools, tmp_path):
    _, reads = _build("K_N5", pools, tmp_path)
    core = sim.load_constants()[FAMILY][0][-sim.CORE_LEN :]

    def share(role):
        hits = 0
        for name, seq in reads[role]:
            r = ev.parse_name(name)
            hits += seq[r.start : r.start + sim.CORE_LEN] == core
        return hits / len(reads[role])

    assert share("earliest") < 0.05
    assert 0.45 < share("latest") < 0.7


def test_the_motif_is_fixed_per_seed_across_configurations(pools, tmp_path):
    a, _ = _build("M50", pools, tmp_path / "a")
    b, _ = _build("M90", pools, tmp_path / "b")
    c, _ = _build("M90", pools, tmp_path / "c", seed=2)
    assert a["record"]["motif"] == b["record"]["motif"] != c["record"]["motif"]


def test_generation_is_deterministic(pools, tmp_path):
    _build("I0.5", pools, tmp_path / "a")
    _build("I0.5", pools, tmp_path / "b")
    for path in (tmp_path / "a").glob("*.fastq.gz"):
        assert path.read_bytes() == (tmp_path / "b" / path.name).read_bytes()


def test_an_insertion_or_deletion_at_a_boundary_of_n_is_flagged():
    P, NN, Q = sim.PRE, sim.N, sim.POST
    seq = list("CCGA" + "ACGT" + "TTGC")
    blocks = [P] * 4 + [NN] * 4 + [Q] * 4
    # Inserting before the first base of N: the base could belong to either side.
    assert sim.insertion_ambiguous(seq, blocks, 4, "G")
    # Deleting an A of the run "A|A" that straddles the 5' boundary.
    assert sim.deletion_ambiguous(seq, blocks, 3)
    assert sim.deletion_ambiguous(seq, blocks, 4)
    # Deleting a T of "T|TT" at the 3' boundary, from the constant's side.
    assert sim.deletion_ambiguous(seq, blocks, 9)
    # Inserting a T two bases into the 3' constant slides to the boundary.
    assert sim.insertion_ambiguous(seq, blocks, 10, "T")
    # Well inside N, or well inside a constant: not ambiguous.
    assert not sim.deletion_ambiguous(seq, blocks, 5)
    assert not sim.insertion_ambiguous(seq, blocks, 6, "A")
    assert not sim.deletion_ambiguous(seq, blocks, 1)


def test_indels_keep_n_contiguous_and_tracked():
    rng = random.Random(3)
    for _ in range(300):
        seq = list("ACGTACGTAC" + "".join(rng.choice("ACGT") for _ in range(20)) + "GGTACCATGA")
        blocks = [sim.PRE] * 10 + [sim.N] * 20 + [sim.POST] * 10
        seq, blocks, _ = sim._indels(seq, blocks, 0.05, rng)
        assert len(seq) == len(blocks)
        positions = [i for i, b in enumerate(blocks) if b == sim.N]
        assert positions == list(range(positions[0], positions[-1] + 1))
        assert blocks == sorted(blocks)


# ---------------------------------------------------------------------------
# Evaluator
# ---------------------------------------------------------------------------


def _write_outputs(run_dir: Path, sub: str, round_reads: dict[int, list[tuple[str, str]]]):
    for rnd, items in round_reads.items():
        path = run_dir / sub / f"round_{rnd:02d}" / "extracted.fasta.gz"
        path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(path, "wt") as fh:
            for name, seq in items:
                fh.write(f">{name}\n{seq}\n")


def _report(run_dir: Path, **fields) -> None:
    base = {
        "primer_5p": None,
        "primer_3p": None,
        "extraction_mode": "BOTH_PRIMERS_SINGLE_READ",
        "status": "HIGH",
        "confidence": 0.9,
        "orientation": "FORWARD",
    }
    base.update(fields)
    (run_dir / "detect").mkdir(parents=True, exist_ok=True)
    (run_dir / "detect" / "library_report.json").write_text(json.dumps(base))
    (run_dir / "run.json").write_text("{}")


def _perfect(run_dir: Path, truth: dict, sub: str = "extract", mode: str = ev.TWO_SIDED):
    reads = ev.read_inputs(run_dir, truth)
    by_round: dict[int, list] = {}
    for name, r in reads.items():
        expected = r.expected(mode)
        if r.eligible(mode):
            by_round.setdefault(r.round, []).append((name, expected))
    _write_outputs(run_dir, sub, by_round)


def test_a_perfect_extraction_is_complete_and_correct(pools, tmp_path):
    truth, _ = _build("BASE", pools, tmp_path)
    c5, c3 = sim.load_constants()[FAMILY]
    _report(tmp_path, primer_5p=c5, primer_3p=c3)
    _perfect(tmp_path, truth)
    _perfect(tmp_path, truth, "oracle/extract")
    row = ev.evaluate(tmp_path)
    assert row["outcome"] == "complete, correct"
    assert row["call_5p_category"] == row["call_3p_category"] == "exact"
    assert row["inferred_precision"] == row["inferred_recall"] == row["inferred_n_recovery"] == 1
    assert row["oracle_recall"] == 1


def test_a_call_reaching_into_n_is_wrong(pools, tmp_path):
    truth, _ = _build("BASE", pools, tmp_path)
    c5, c3 = sim.load_constants()[FAMILY]
    _report(tmp_path, primer_5p=c5, primer_3p="A" + c3)
    _perfect(tmp_path, truth)
    row = ev.evaluate(tmp_path)
    assert row["outcome"] == "wrong"
    assert row["call_3p_category"] == "other"


def test_outer_material_is_correct_and_labelled(pools, tmp_path):
    truth, _ = _build("K_5outer", pools, tmp_path)
    (pre,) = truth["rounds"][0]["pre"]
    c3 = sim.load_constants()[FAMILY][1]
    _report(tmp_path, primer_5p=pre, primer_3p=c3)
    row = ev.evaluate(tmp_path)
    assert row["call_5p_correct"] and row["call_5p_category"] == "includes outer material"


def test_a_leading_random_base_in_the_call_is_wrong(pools, tmp_path):
    truth, _ = _build("T3", pools, tmp_path)
    c5 = sim.load_constants()[FAMILY][0]
    assert len(truth["rounds"][0]["pre"]) == 4
    assert ev.call_correct(c5, "5p", truth, reverse=False)
    assert not ev.call_correct("A" + c5, "5p", truth, reverse=False)


def test_a_refusal_is_refused_with_no_inferred_metrics(pools, tmp_path):
    _build("C100", pools, tmp_path)
    _report(tmp_path, status="UNABLE_TO_INFER", extraction_mode="UNABLE_TO_EXTRACT")
    row = ev.evaluate(tmp_path)
    assert row["outcome"] == "refused"
    assert "inferred_recall" not in row


def test_recall_never_counts_reads_outside_its_denominator(pools, tmp_path):
    """Amendment 2, item 1: in F_trunc30 extract recovers truncated reads too."""
    truth, _ = _build("F_trunc30", pools, tmp_path)
    c5, c3 = sim.load_constants()[FAMILY]
    _report(tmp_path, primer_5p=c5, primer_3p=c3)
    reads = ev.read_inputs(tmp_path, truth)
    by_round: dict[int, list] = {}
    for name, r in reads.items():  # every read, eligible or not, emitted correctly
        by_round.setdefault(r.round, []).append((name, r.seq[r.start : r.end]))
    _write_outputs(tmp_path, "extract", by_round)
    row = ev.evaluate(tmp_path)
    assert row["inferred_recall"] == 1
    assert row["inferred_eligible"] < len(reads)
    assert row["inferred_n_recovery"] == 1


def test_a_zero_denominator_is_not_applicable(pools, tmp_path):
    truth, _ = _build("F_no5", pools, tmp_path)
    reads = ev.read_inputs(tmp_path, truth)
    metrics = ev.per_read("FIVE_PRIME_ONLY", reads, {})
    assert metrics["recall"] is None and metrics["precision"] is None


def test_flagged_reads_are_kept_apart(pools, tmp_path):
    truth, _ = _build("E0", pools, tmp_path)
    reads = ev.read_inputs(tmp_path, truth)
    name = next(iter(reads))
    flagged = {
        name.replace("_i0", "_i1"): reads[name].__class__(
            **{**reads[name].__dict__, "flagged": True}
        )
    }
    metrics = ev.per_read(ev.TWO_SIDED, flagged, {name.replace("_i0", "_i1"): "WRONG"})
    assert metrics["flagged_reads"] == 1 and metrics["flagged_emitted"] == 1
    assert metrics["precision"] is None


@pytest.mark.parametrize("refused", [True, False])
def test_negative_controls(pools, tmp_path, refused):
    _build("NC_donor", pools, tmp_path)
    if refused:
        _report(tmp_path, status="UNABLE_TO_INFER", extraction_mode="UNABLE_TO_EXTRACT")
    else:
        _report(
            tmp_path,
            primer_5p="ACGTACGTACGTACGT",
            primer_3p=None,
            extraction_mode="FIVE_PRIME_ONLY",
        )
    assert ev.evaluate(tmp_path)["outcome"] == (
        "negative control, refused" if refused else "negative control, wrong"
    )


def test_a_pipeline_error_is_an_error(pools, tmp_path):
    _build("BASE", pools, tmp_path)
    (tmp_path / "run.json").write_text(json.dumps({"error": "detect exited 1"}))
    assert ev.evaluate(tmp_path)["outcome"] == "error"


def test_a_reversed_read_is_evaluated_in_the_frame_extract_uses():
    read = ev.Read(1, "AAAACCGGTTTT", 4, 8, True, True, False, False)
    rev = read.reversed()
    assert rev.seq == "AAAACCGGTTTT"[::-1].translate(str.maketrans("ACGT", "TGCA"))
    assert rev.expected(ev.TWO_SIDED) == sim.revcomp("CCGG")
    assert (rev.five, rev.three) == (False, True)


# ---------------------------------------------------------------------------
# Harness, end to end (needs cutadapt)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    shutil.which("cutadapt") is None and not (Path(sys.executable).parent / "cutadapt").exists(),
    reason="cutadapt not available",
)
def test_one_baseline_run_end_to_end(tmp_path):
    pools = _fake_pools(tmp_path / "pools", size=1500)
    harness = _load("run_synthetic")
    cell = sim.Cell("BASE", DONOR, FAMILY, 1)
    row = harness.run_cell(cell, pools, tmp_path / "out", keep=False)
    assert row["outcome"] == "complete, correct", row
    assert row["oracle_recall"] > 0.9 and row["inferred_precision"] > 0.99
    run_dir = tmp_path / "out" / "runs" / cell.cell_id
    assert not list(run_dir.glob("*.fastq.gz"))  # bulky files removed
    assert (run_dir / "result.json").is_file()
    # Resuming returns the stored result without running again.
    assert harness.run_cell(cell, pools, tmp_path / "out", keep=False) == json.loads(
        (run_dir / "result.json").read_text()
    )
