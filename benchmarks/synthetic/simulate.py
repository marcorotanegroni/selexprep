"""Generator of the semi-synthetic benchmark (DESIGN.md, "Simulated libraries",
"Configurations" and amendments 2-3).

One *cell* is a configuration x donor x constant family x seed. For a cell this
module builds, per provided round, a single-end FASTQ.gz of
``[5' outer][5' constant][N][3' constant][3' outer][read-through]`` reads whose
random region N comes from the donor's sampled pool, and a ``truth.json``.

Every read name carries its truth, so an extracted sequence can be traced back:

    R{round}_{index}_s{start}_e{end}_r{two}{five}{three}_i{flag}

``[start, end)`` is N in the observed read, after every step; ``two``, ``five``
and ``three`` say whether the read is recoverable by two-sided, 5'-only and
3'-only extraction (the constants that extraction needs are present and not
truncated), in read orientation; ``flag`` is 1 when an indel of the read has a
placement giving the same read that touches a boundary of N (amendment 2,
item 5).

Construction order, per read: assemble the molecule; apply the configuration's
edits; truncate (F); indels (I), then substitutions (E); reverse-complement last
(O). Nothing here imports ``selexprep``.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import random
from dataclasses import dataclass, field, replace
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
ROLES = ("earliest", "middle", "latest")
DEV_SEEDS = (1, 2, 3)
# Amendment 3: five test donors, so test uses seeds 101-102 only.
TEST_SEEDS = (101, 102)
CORE_LEN = 12
TRUSEQ_R1_RC_13 = "GCTCTTCCGATCT"
READ_THROUGH = "AGATCGGAAGAGCACACGTC"
PRE, N, POST = 0, 1, 2


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


# ---------------------------------------------------------------------------
# Configurations
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Spec:
    """What a configuration changes relative to the baseline."""

    rounds: tuple[str, ...] = ROLES
    constant_len: int | None = None
    substitution: float = 0.01
    indel: float = 0.0
    tag: str | None = None  # "per_round" | "shared"
    lead_base: bool = False
    t4: bool = False
    read_through: bool = False
    adapter: str | None = None  # "outer_start" | "middle" | "inner_end"
    core_copy: str | None = None  # "5outer" | "3outer" | "N5" | "N3"
    clone: float | None = None
    depth: int | str | None = None  # 500 | 2000 | "1000distinct"
    motif: float | None = None
    flank: str | None = None  # "no3" | "no5" | "trunc30"
    biological: bool = False
    reverse: bool = False
    negative: str | None = None  # "donor" | "clone"


BASE = Spec()
SPECS: dict[str, Spec] = {
    "BASE": BASE,
    **{f"L{k:02d}": replace(BASE, constant_len=k) for k in (8, 10, 12, 14, 16, 20)},
    "E0": replace(BASE, substitution=0.0),
    "E0.5": replace(BASE, substitution=0.005),
    "E2": replace(BASE, substitution=0.02),
    "E5": replace(BASE, substitution=0.05),
    "I0.1": replace(BASE, indel=0.001),
    "I0.5": replace(BASE, indel=0.005),
    "T1": replace(BASE, tag="per_round"),
    "T2": replace(BASE, tag="shared"),
    "T3": replace(BASE, lead_base=True),
    "T4": replace(BASE, t4=True),
    "T5": replace(BASE, read_through=True),
    "A_outer_start": replace(BASE, adapter="outer_start"),
    "A_middle": replace(BASE, adapter="middle"),
    "A_inner_end": replace(BASE, adapter="inner_end"),
    "K_5outer": replace(BASE, core_copy="5outer"),
    "K_3outer": replace(BASE, core_copy="3outer"),
    "K_N5": replace(BASE, core_copy="N5"),
    "K_N3": replace(BASE, core_copy="N3"),
    "R_earliest": replace(BASE, rounds=("earliest",)),
    "R_latest": replace(BASE, rounds=("latest",)),
    "R_middle_latest": replace(BASE, rounds=("middle", "latest")),
    **{f"C{p}": replace(BASE, clone=p / 100) for p in (50, 70, 80, 90, 99, 100)},
    "D500": replace(BASE, depth=500),
    "D2000": replace(BASE, depth=2000),
    "D1000distinct": replace(BASE, depth="1000distinct"),
    "M50": replace(BASE, motif=0.5),
    "M90": replace(BASE, motif=0.9),
    "F_no3": replace(BASE, flank="no3"),
    "F_trunc30": replace(BASE, flank="trunc30"),
    "F_no5": replace(BASE, flank="no5"),
    "B": replace(BASE, biological=True),
    "O": replace(BASE, reverse=True),
    "X1_E0.5": replace(BASE, tag="per_round", substitution=0.005),
    "X1_E2": replace(BASE, tag="per_round", substitution=0.02),
    "X1_E5": replace(BASE, tag="per_round", substitution=0.05),
    "X2_L14_C70": replace(BASE, constant_len=14, clone=0.7),
    "X2_L16_C70": replace(BASE, constant_len=16, clone=0.7),
    "X3": replace(BASE, rounds=("latest",), motif=0.9),
    "NC_donor": replace(BASE, negative="donor"),
    "NC_clone": replace(BASE, negative="clone"),
}


@dataclass(frozen=True)
class Cell:
    config: str
    donor: str
    family: str | None  # None for negative controls
    seed: int

    @property
    def cell_id(self) -> str:
        return f"{self.config}__{self.donor}__{self.family or 'none'}__s{self.seed}"


def _tsv(name: str) -> list[dict[str, str]]:
    with (HERE / name).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_constants() -> dict[str, tuple[str, str]]:
    return {r["family"]: (r["constant_5p"], r["constant_3p"]) for r in _tsv("constants.tsv")}


def load_donors() -> dict[str, dict]:
    """Donor -> set, N and its three rounds (role -> round number, run)."""
    donors: dict[str, dict] = {}
    for r in _tsv("donors.tsv"):
        d = donors.setdefault(r["donor"], {"set": r["set"], "n": int(r["n_stated"]), "rounds": {}})
        d["rounds"][r["role"]] = (int(r["round"]), r["run_accession"])
    return donors


def cells(which: str) -> list[Cell]:
    """Every run of the development or the test set, from ``configurations.tsv``."""
    column = "dev_families" if which == "dev" else "test_families"
    seeds = DEV_SEEDS if which == "dev" else TEST_SEEDS
    donors = sorted(d for d, v in load_donors().items() if v["set"] == which)
    out = []
    for row in _tsv("configurations.tsv"):
        spec = SPECS[row["config_id"]]
        families = [None] if spec.negative else row[column].split(",")
        for donor in donors:
            for family in families:
                for seed in seeds:
                    out.append(Cell(row["config_id"], donor, family, seed))
    return out


# ---------------------------------------------------------------------------
# Building one cell
# ---------------------------------------------------------------------------


def _rng(*parts: object) -> random.Random:
    key = "|".join(str(p) for p in parts)
    return random.Random(int(hashlib.sha256(key.encode()).hexdigest()[:16], 16))


def _random_seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice(BASES) for _ in range(n))


def _other_base(rng: random.Random, base: str) -> str:
    return rng.choice([b for b in BASES if b != base])


def pool_path(pools: Path, donor: str, role: str, round_number: int, run: str, seed: int) -> Path:
    return pools / donor / f"{role}_round{round_number}_{run}_seed{seed}.txt"


def load_pool(path: Path) -> list[str]:
    return [line for line in path.read_text().splitlines() if line]


@dataclass
class Library:
    """The per-cell choices every read of the cell shares."""

    c5: str
    c3: str
    outer5: dict[str, str] = field(default_factory=dict)  # role -> outer 5' sequence
    outer3: dict[str, str] = field(default_factory=dict)
    c3_by_role: dict[str, str] = field(default_factory=dict)
    record: dict = field(default_factory=dict)


def _constants(spec: Spec, family: str | None, rng: random.Random) -> Library:
    if spec.negative:
        return Library(
            "", "", {r: "" for r in ROLES}, {r: "" for r in ROLES}, {r: "" for r in ROLES}
        )
    c5, c3 = load_constants()[family]  # type: ignore[index]
    if spec.constant_len is not None:
        k = spec.constant_len
        c5, c3 = c5[-k:], c3[:k]
    if spec.adapter is not None:
        a = TRUSEQ_R1_RC_13
        start = {"outer_start": 0, "middle": (len(c5) - len(a)) // 2, "inner_end": len(c5) - len(a)}
        i = start[spec.adapter]
        c5 = c5[:i] + a + c5[i + len(a) :]
    lib = Library(c5, c3)
    for role in ROLES:
        lib.outer5[role], lib.outer3[role], lib.c3_by_role[role] = "", "", c3
    if spec.t4:
        changed = c3[:-1] + _other_base(rng, c3[-1])
        lib.c3_by_role["earliest"] = changed
        lib.record["t4_earliest_3p"] = changed
    if spec.tag == "shared":
        t5, t3 = _random_seq(rng, 6), _random_seq(rng, 6)
        for role in ROLES:
            lib.outer5[role], lib.outer3[role] = t5, t3
        lib.record["tags"] = {"5p": t5, "3p": t3}
    elif spec.tag == "per_round":
        lib.record["tags"] = {}
        for role in ROLES:
            t5 = _random_seq(rng, rng.randint(4, 8))
            t3 = _random_seq(rng, rng.randint(4, 8))
            lib.outer5[role], lib.outer3[role] = t5, t3
            lib.record["tags"][role] = {"5p": t5, "3p": t3}
    if spec.core_copy in ("5outer", "3outer"):
        # Amendment 2, item 4: the 4-nt spacer is fixed per cell and seed.
        spacer = _random_seq(rng, 4)
        lib.record["spacer"] = spacer
        for role in ROLES:
            if spec.core_copy == "5outer":
                lib.outer5[role] = c5[-CORE_LEN:] + spacer
            else:
                lib.outer3[role] = spacer + c3[:CORE_LEN]
    if spec.read_through:
        for role in ROLES:
            lib.outer3[role] += READ_THROUGH
    return lib


def _set_clone_share(
    ns: list[str], clone: str, share: float, rng: random.Random
) -> tuple[list[str], float]:
    """Amendment 2, item 3: the clone's count becomes exactly ceil(share * M).

    Reads that already equal the clone count toward it; a surplus is replaced by
    reads drawn with replacement from the round's other reads.
    """
    target = math.ceil(share * len(ns))
    ns = list(ns)
    idx_clone = [i for i, s in enumerate(ns) if s == clone]
    idx_other = [i for i, s in enumerate(ns) if s != clone]
    before = len(idx_clone) / len(ns)
    if len(idx_clone) < target:
        for i in rng.sample(idx_other, target - len(idx_clone)):
            ns[i] = clone
    elif len(idx_clone) > target:
        others = [ns[i] for i in idx_other]
        for i in rng.sample(idx_clone, len(idx_clone) - target):
            ns[i] = rng.choice(others)
    return ns, before


def _edit_random_regions(
    spec: Spec, role: str, ns: list[str], lib: Library, clone: str | None, rng: random.Random
) -> list[str]:
    """Edits to N, in the order of "Construction order": core copy, motif, clone, N±1."""
    ns = list(ns)
    late = role in ("middle", "latest")
    if spec.core_copy in ("N5", "N3") and late:
        core = lib.c5[-CORE_LEN:] if spec.core_copy == "N5" else lib.c3[:CORE_LEN]
        for i in rng.sample(range(len(ns)), round(0.6 * len(ns))):
            s = ns[i]
            ns[i] = core + s[CORE_LEN:] if spec.core_copy == "N5" else s[:-CORE_LEN] + core
    if spec.motif is not None and late:
        motif = lib.record["motif"]
        for i in rng.sample(range(len(ns)), round(spec.motif * len(ns))):
            ns[i] = motif + ns[i][len(motif) :]
    if clone is not None and spec.clone is not None:
        ns, before = _set_clone_share(ns, clone, spec.clone, rng)
        lib.record.setdefault("clone_share_before", {})[role] = round(before, 6)
    if spec.biological:
        order = list(range(len(ns)))
        rng.shuffle(order)
        tenth = len(ns) // 10
        for i in order[:tenth]:
            j = rng.randrange(len(ns[i]))
            ns[i] = ns[i][:j] + ns[i][j + 1 :]
        for i in order[tenth : 2 * tenth]:
            j = rng.randrange(len(ns[i]) + 1)
            ns[i] = ns[i][:j] + rng.choice(BASES) + ns[i][j:]
    return ns


def _depth(spec: Spec, ns: list[str], rng: random.Random) -> list[str]:
    if spec.depth is None:
        return ns
    if spec.depth == "1000distinct":
        distinct = sorted(set(ns))
        chosen = rng.sample(distinct, min(1000, len(distinct)))
        return [rng.choice(chosen) for _ in range(len(ns))]
    order = list(range(len(ns)))
    rng.shuffle(order)
    return [ns[i] for i in order[: int(spec.depth)]]


def _boundary_touching(seq: list[str], blocks: list[int], lo: int, hi: int) -> bool:
    """True when a gap in [lo, hi] lies between N and another block."""
    for g in range(max(lo, 1), min(hi, len(seq) - 1) + 1):
        a, b = blocks[g - 1], blocks[g]
        if a != b and N in (a, b):
            return True
    return False


def insertion_ambiguous(seq: list[str], blocks: list[int], i: int, base: str) -> bool:
    """Inserting ``base`` before position ``i`` gives the same read at every gap of
    the run of ``base`` around that gap; ambiguous when one of them borders N."""
    lo = i
    while lo > 0 and seq[lo - 1] == base:
        lo -= 1
    hi = i
    while hi < len(seq) and seq[hi] == base:
        hi += 1
    return _boundary_touching(seq, blocks, lo, hi)


def deletion_ambiguous(seq: list[str], blocks: list[int], i: int) -> bool:
    """Deleting position ``i`` gives the same read as deleting any base of its
    run; ambiguous when the run holds N and another block."""
    a = b = i
    while a > 0 and seq[a - 1] == seq[i]:
        a -= 1
    while b + 1 < len(seq) and seq[b + 1] == seq[i]:
        b += 1
    run = set(blocks[a : b + 1])
    return N in run and len(run) > 1


def _indels(
    seq: list[str], blocks: list[int], rate: float, rng: random.Random
) -> tuple[list[str], list[int], bool]:
    """Indels at ``rate`` per base, half insertions and half deletions.

    An insertion goes before the sampled position and joins that position's
    block; a deletion removes the sampled position. Events are applied from the
    3' end, so earlier coordinates are unaffected. The read is flagged when a
    placement of an indel that gives the same read touches a boundary of N.
    """
    events = []
    for i in range(len(seq)):
        u = rng.random()
        if u < rate:
            events.append((i, "ins" if u < rate / 2 else "del"))
    flagged = False
    for i, kind in reversed(events):
        if kind == "ins":
            x = rng.choice(BASES)
            flagged |= insertion_ambiguous(seq, blocks, i, x)
            seq.insert(i, x)
            blocks.insert(i, blocks[i])
        else:
            flagged |= deletion_ambiguous(seq, blocks, i)
            del seq[i]
            del blocks[i]
    return seq, blocks, flagged


def _substitutions(seq: list[str], rate: float, rng: random.Random) -> list[str]:
    """Rate ``rate * (0.5 + i / L)`` at position i: mean ``rate``, rising to 3'."""
    length = len(seq)
    for i in range(length):
        if rng.random() < rate * (0.5 + i / length):
            seq[i] = _other_base(rng, seq[i])
    return seq


@dataclass
class ReadTruth:
    name: str
    seq: str
    pre: str  # error-free, untruncated sequence before N, read orientation
    post: str


def _build_read(
    spec: Spec,
    role: str,
    round_number: int,
    index: int,
    n: str,
    lib: Library,
    rng: random.Random,
) -> ReadTruth:
    pre = lib.outer5[role] + ("" if spec.flank == "no5" else lib.c5)
    if spec.lead_base:
        pre = rng.choice(BASES) + pre
    c3 = "" if spec.flank == "no3" else lib.c3_by_role[role]
    post = c3 + lib.outer3[role] if c3 else ""
    has5 = spec.negative is None and spec.flank != "no5"
    has3 = spec.negative is None and spec.flank != "no3"
    complete3 = has3
    body = pre + n + post
    blocks = [PRE] * len(pre) + [N] * len(n) + [POST] * len(post)
    if spec.flank == "trunc30" and rng.random() < 0.3:
        cut = len(pre) + len(n) + rng.randint(1, len(c3) - 1)
        body, blocks = body[:cut], blocks[:cut]
        complete3 = False
    seq = list(body)
    flagged = False
    if spec.indel:
        seq, blocks, flagged = _indels(seq, blocks, spec.indel, rng)
    seq = _substitutions(seq, spec.substitution, rng)
    read = "".join(seq)
    ok5, ok3 = has5, complete3
    if spec.reverse:
        read = revcomp(read)
        blocks = blocks[::-1]
        pre, post = revcomp(post), revcomp(pre)
        ok5, ok3 = ok3, ok5
    positions = [i for i, b in enumerate(blocks) if b == N]
    start, end = (positions[0], positions[-1] + 1) if positions else (0, 0)
    flags = f"{int(ok5 and ok3)}{int(ok5)}{int(ok3)}"
    name = f"R{round_number}_{index:05d}_s{start}_e{end}_r{flags}_i{int(flagged)}"
    return ReadTruth(name, read, pre, post)


def _oracle(spec: Spec, lib: Library) -> dict:
    """The report the truth allows (DESIGN.md, "The oracle"; amendment 2, item 2)."""
    if spec.negative:
        return {"mode": None, "primer_5p": None, "primer_3p": None}
    c5 = None if spec.flank == "no5" else lib.c5
    c3 = None if spec.flank == "no3" else (lib.c3[:-1] if spec.t4 else lib.c3)
    if spec.reverse:
        c5, c3 = (revcomp(c3) if c3 else None), (revcomp(c5) if c5 else None)
    mode = (
        "BOTH_PRIMERS_SINGLE_READ" if c5 and c3 else "FIVE_PRIME_ONLY" if c5 else "THREE_PRIME_ONLY"
    )
    return {"mode": mode, "primer_5p": c5, "primer_3p": c3}


def _write_fastq(path: Path, reads: list[ReadTruth]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with (
        path.open("wb") as raw,
        gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as gz,
    ):
        for r in reads:
            gz.write(f"@{r.name}\n{r.seq}\n+\n{'I' * len(r.seq)}\n".encode())


def generate(cell: Cell, pools: Path, outdir: Path) -> dict:
    """Write the cell's FASTQs, ``rounds.tsv`` and ``truth.json`` under ``outdir``."""
    spec = SPECS[cell.config]
    donor = load_donors()[cell.donor]
    rng = _rng(cell.config, cell.donor, cell.family, cell.seed)
    lib = _constants(spec, cell.family, rng)
    if spec.motif is not None:
        # "6 random nt (fixed per seed)".
        lib.record["motif"] = _random_seq(_rng("motif", cell.seed), 6)

    earliest_round, earliest_run = donor["rounds"]["earliest"]
    earliest_pool = load_pool(
        pool_path(pools, cell.donor, "earliest", earliest_round, earliest_run, cell.seed)
    )
    clone = None
    if spec.clone is not None or spec.negative == "clone":
        clone = rng.choice(earliest_pool)
        lib.record["clone"] = clone

    rounds_out = []
    round_map = []
    for role in spec.rounds:
        round_number, run = donor["rounds"][role]
        ns = load_pool(pool_path(pools, cell.donor, role, round_number, run, cell.seed))
        if spec.negative == "clone":
            ns = [clone] * len(ns)  # type: ignore[list-item]
        ns = _edit_random_regions(spec, role, ns, lib, clone, rng)
        ns = _depth(spec, ns, rng)
        if clone is not None and spec.clone is not None:
            share = sum(1 for s in ns if s == clone) / len(ns)
            lib.record.setdefault("clone_share", {})[role] = round(share, 6)
        reads = [_build_read(spec, role, round_number, i, s, lib, rng) for i, s in enumerate(ns)]
        fastq = outdir / f"round_{round_number:02d}.fastq.gz"
        _write_fastq(fastq, reads)
        round_map.append((fastq.name, round_number))
        pre_counts: dict[str, int] = {}
        post_counts: dict[str, int] = {}
        for r in reads:
            pre_counts[r.pre] = pre_counts.get(r.pre, 0) + 1
            post_counts[r.post] = post_counts.get(r.post, 0) + 1
        rounds_out.append(
            {
                "role": role,
                "round": round_number,
                "fastq": fastq.name,
                "reads": len(reads),
                "pre": pre_counts,
                "post": post_counts,
            }
        )

    with (outdir / "rounds.tsv").open("w") as fh:
        fh.write("file\tround_number\n")
        for name, number in round_map:
            fh.write(f"{name}\t{number}\n")
    truth = {
        "cell": cell.cell_id,
        "config": cell.config,
        "donor": cell.donor,
        "family": cell.family,
        "seed": cell.seed,
        "negative_control": spec.negative is not None,
        "rounds": rounds_out,
        "oracle": _oracle(spec, lib),
        "record": lib.record,
    }
    (outdir / "truth.json").write_text(json.dumps(truth, indent=2, sort_keys=True) + "\n")
    return truth
