"""Residue probes of the semi-synthetic donor check (amendment 2, items 6-7).

The first version of the check built one 8-mer per adapter and tested it with
``startswith`` / ``endswith``. A residue longer than 8 nt shifts that 8-mer away
from the read edge, so it was missed: these tests pin the fix.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = (
    Path(__file__).resolve().parent.parent
    / "benchmarks"
    / "synthetic"
    / "check_and_sample_donors.py"
)
TRUSEQ_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
N30 = "ACGTTGCAAGGCTTACGGATCCATGCAAGT"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("check_and_sample_donors", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def adapters(mod, side="end"):
    probes = mod.EdgeProbes(side)
    for name, seq in mod.ADAPTERS.items():
        probes.add(name, seq)
    return probes


@pytest.mark.parametrize("k", [8, 9, 13, 20, 34])
def test_an_adapter_residue_is_found_at_every_length_from_eight(mod, k):
    probes = adapters(mod)
    found = probes.match(N30 + TRUSEQ_R1[:k])
    assert found is not None
    _, length = found
    assert length == k


def test_a_residue_shorter_than_the_floor_is_not_a_signal(mod):
    assert adapters(mod).match(N30 + TRUSEQ_R1[:7]) is None


def test_a_read_with_no_residue_is_not_a_signal(mod):
    assert adapters(mod).match(N30) is None


def test_the_two_truseq_adapters_are_one_family(mod):
    """R1 and R2 share 13 nt: a residue counts toward one TruSeq signal at any length."""
    probes = adapters(mod)
    assert probes.match(N30 + TRUSEQ_R1[:13]) == ("TruSeq", 13)
    assert probes.match(N30 + TRUSEQ_R1[:14]) == ("TruSeq", 14)


def test_the_longest_match_is_reported_once(mod):
    """A read carrying 20 nt of adapter is one signal, not thirteen."""
    probes = adapters(mod)
    label, length = probes.match(N30 + TRUSEQ_R1[:20])
    assert (label, length) == ("TruSeq", 20)


def test_a_five_prime_constant_is_probed_at_the_read_start(mod):
    constant = "GGCAGGTTCTAAGGCTAGCA"
    probes = mod.EdgeProbes("start")
    probes.add("constant_5p", constant)
    assert probes.match(constant[-11:] + N30) == ("constant_5p", 11)
    assert probes.match(constant[-7:] + N30) is None


def test_a_documented_fragment_is_probed_down_to_six_nucleotides(mod):
    probes = mod.EdgeProbes("start")
    probes.add("constant_5p", "TCCCAA", kmin=6)
    assert probes.match("TCCCAA" + N30) == ("constant_5p", 6)
    assert probes.described == ["constant_5p=TCCCAA(>=6nt)"]


def test_a_fragment_below_six_nucleotides_is_not_probed(mod):
    probes = mod.EdgeProbes("end")
    probes.add("constant_3p", "GCGC", kmin=4)
    assert probes.described == ["constant_3p=GCGC(too short, <6nt)"]
    assert probes.match(N30 + "GCGC") is None


def test_short_and_long_probes_coexist_on_one_edge(mod):
    """The index key must be as short as the shortest probe, or it hides it."""
    probes = adapters(mod)
    probes.add("constant_3p", "TTGAGC", kmin=6)
    assert probes.match(N30 + "TTGAGC") == ("constant_3p", 6)
    assert probes.match(N30 + TRUSEQ_R1[:20]) == ("TruSeq", 20)


def test_constants_are_loaded_per_side(mod):
    """A deposit with one documented side is probed on that side (item 7)."""
    loaded = mod.documented_constants()
    assert loaded["PRJEB25907"]["5p"] == ("TCCCAA", 6)
    # Documented but too short to probe: loaded, and recorded as such.
    assert loaded["PRJEB25907"]["3p"] == ("GCGC", 4)
    for accession in ("PRJEB51212", "PRJEB51473", "PRJNA741127"):
        assert loaded.get(accession, {}) == {}
    five, kmin = loaded["PRJEB70964"]["5p"]
    assert kmin == mod.RESIDUE_K and len(five) > mod.RESIDUE_K


def fastq(path: Path, reads: list[str]) -> Path:
    import gzip

    with gzip.open(path, "wt") as fh:
        for i, seq in enumerate(reads):
            fh.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    return path


@pytest.fixture
def run(mod, tmp_path, monkeypatch):
    """12,000 reads of 20 nt, plus 1,200 longer reads carrying 11 nt of TruSeq."""
    import random

    rng = random.Random(7)
    reads = ["".join(rng.choice("ACGT") for _ in range(20)) for _ in range(12_000)]
    reads += [r + rng.choice("ACGT") * 6 + TRUSEQ_R1[:11] for r in reads[:1_200]]
    monkeypatch.setattr(
        mod, "download", lambda accession, cache: fastq(tmp_path / "r.fastq.gz", reads)
    )
    return reads


def test_one_pool_per_seed_with_its_own_sample(mod, run, tmp_path):
    """Amendment 2: a seed varies the donor sampling, not only the error model."""
    item = dict(donor="PRJTEST", set="test", n=20, role="earliest", round="1", run="ERR1")
    row = mod.check_and_sample(item, {}, tmp_path, tmp_path / "pools")
    assert row["sample_seeds"] == "101;102;103"
    digests = dict(pair.split("=") for pair in row["sample_sha256"].split(";"))
    assert len(set(digests.values())) == 3
    pools = sorted((tmp_path / "pools" / "PRJTEST").glob("*.txt"))
    assert [p.name for p in pools] == [
        f"earliest_round1_ERR1_seed{seed}.txt" for seed in (101, 102, 103)
    ]
    assert all(len(p.read_text().splitlines()) == mod.SAMPLE_SIZE for p in pools)


def test_the_mechanical_part_of_step_three_and_the_residue_signal(mod, run, tmp_path):
    item = dict(donor="PRJTEST", set="dev", n=20, role="earliest", round="1", run="ERR1")
    row = mod.check_and_sample(item, {}, tmp_path, tmp_path / "pools")
    assert row["sample_seeds"] == "1;2;3"
    assert row["modal_length"] == 20
    assert row["step3_reads"] == "pass"  # 12,000 of 13,200 reads are 20 nt
    assert row["within_n2_fraction"] == "0.9091"
    assert row["residue_signals"] == "TruSeq=0.0910"
    assert row["residue_max_length"] == "TruSeq=11"
    assert row["probes_end"].startswith("TruSeq_R1=AGATCGGA")
    assert row["probes_end"].endswith("constant_3p=unknown")
    assert row["probes_start"] == "constant_5p=unknown"


def test_step_three_fails_when_the_modal_length_is_not_n(mod, run, tmp_path):
    """The 21-nt case: one base too many in every read is a fail, not a pass."""
    item = dict(donor="PRJTEST", set="dev", n=21, role="earliest", round="1", run="ERR1")
    row = mod.check_and_sample(item, {}, tmp_path, tmp_path / "pools")
    assert row["step3_reads"].startswith("fail: modal length 20 != N 21")


def test_one_residue_split_across_lengths_is_one_signal(mod, tmp_path, monkeypatch):
    """3% of reads with 10 nt of TruSeq (shared by R1 and R2) and 3% with 20 nt of
    R1 are one 6% residue: counted apart, neither half would reach the 5% trigger."""
    import random

    rng = random.Random(11)
    reads = ["".join(rng.choice("ACGT") for _ in range(20)) for _ in range(9_400)]
    reads += [r + TRUSEQ_R1[:10] for r in reads[:300]]
    reads += [r + TRUSEQ_R1[:20] for r in reads[300:600]]
    monkeypatch.setattr(
        mod, "download", lambda accession, cache: fastq(tmp_path / "s.fastq.gz", reads)
    )
    item = dict(donor="PRJTEST", set="dev", n=20, role="earliest", round="1", run="ERR2")
    row = mod.check_and_sample(item, {}, tmp_path, tmp_path / "pools")
    assert row["residue_signals"] == "TruSeq=0.0600"
    assert row["residue_max_length"] == "TruSeq=20"


def test_a_five_prime_probe_is_described_by_the_end_it_matches(mod):
    """The start-side probe matches the constant's last bases, so it is shown by them."""
    probes = mod.EdgeProbes("start")
    probes.add("constant_5p", "GGCAGGTTCTAAGGCTAGCA")
    assert probes.described == ["constant_5p=...GGCTAGCA(>=8nt)"]
