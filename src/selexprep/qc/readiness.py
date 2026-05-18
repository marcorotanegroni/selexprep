"""Full sequence-level readiness review of processed SELEX rounds.

Scans EVERY sequence in EVERY round parquet (no sampling). For each
BioProject, runs eight diagnostic sections that each yield a status from
``PASS / WARN / FAIL / INFO``:

- ``pre``         — required artifacts exist; loads each round parquet once
- ``alphabet``    — every sequence is pure ACGT
- ``lengths``     — per-round mode_L matches expected random-region length;
                    mode-L coverage; no short primer-dimer outliers; mode_L
                    consistent across rounds
- ``trim_seq``    — sequences do NOT start with primer_5p, do NOT end with
                    primer_3p, and TruSeq R1 (AGATCGGAAGAGC) <1% of reads
- ``composition`` — per-position top-base on mode-L reads of the earliest
                    round; tag-aware thresholds (undoped ≤50% with entropy
                    fallback; doped ≤90%)
- ``diversity``   — per-position Shannon entropy; singleton fraction; top
                    6-mer dominance (undoped: WARN if top 6-mer >10%)
- ``selection``   — final enrich parquet: max_log2FC ≥ 2.0; >100 post-noise
                    sequences; positive log2FC fraction; cross-round top-5
                    trace
- ``consistency`` — n_unique / n_reads ratio across rounds; seed-vs-summary
                    primer alignment

Library tags drive composition/diversity thresholds. The selexprep package
does NOT hardcode a BP→tag dict (that's thesis-specific data); callers
supply the tag per BioProject when invoking ``review_bioproject``.

Supported tag values (passed through from the original selex_corpus rule
set, kept stable so thesis-side callers can reuse them as-is):
    standard / undoped / doped_3pct / doped_8pct /
    variable_length_undoped / multiplexed_origin / aptasim_excluded /
    no_r0_post_selection / untagged
"""

from __future__ import annotations

import json
import math
import re
from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

TRUSEQ_R1_PREFIX = "AGATCGGAAGAGC"
TRUSEQ_MAX_FRAC_READS = 0.01
MODE_L_TOLERANCE = 2
UNDOPED_POS_TOPBASE_MAX = 0.50  # under 50% is library bias / sequencing
DOPED_POS_TOPBASE_MAX = 0.90
UNDOPED_MIN_ENTROPY_BITS = 1.7  # pristine random has H≈2.0
MIN_MAX_LOG2FC = 2.0
MIN_POST_NOISE_SEQS = 100
PRIMER_MATCH_PREFIX_LEN = 15
TOP_KMER_K = 6
TOP_KMER_MAX_FRAC_UNDOPED = 0.10
ALPHABET = set("ACGT")
PRIMER_RESIDUE_MAX_FRAC = 0.02
ALPHABET_BAD_FRAC_FAIL = 0.005
ALPHABET_BAD_FRAC_WARN = 0.001
PRIMER_MIN_LEN = 10
VALID_RRL_SOURCES_LITERAL = {"post_rebuild_verified", "paper"}
NO_PRIMER_3P_TAGS = {"multiplexed_origin"}  # 5p-only by design

PASS, WARN, FAIL, INFO = "PASS", "WARN", "FAIL", "INFO"

SECTIONS = (
    "pre",
    "alphabet",
    "lengths",
    "trim_seq",
    "composition",
    "diversity",
    "selection",
    "consistency",
)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class CheckResult:
    section: str
    status: str  # PASS | WARN | FAIL | INFO
    detail: str = ""


@dataclass
class BPReport:
    bp_id: str
    tag: str
    checks: list[CheckResult] = field(default_factory=list)
    stats: dict[str, Any] = field(default_factory=dict)

    def add(self, section: str, status: str, detail: str = "") -> None:
        self.checks.append(CheckResult(section, status, detail))

    @property
    def worst(self) -> str:
        s = {c.status for c in self.checks}
        if FAIL in s:
            return FAIL
        if WARN in s:
            return WARN
        return PASS

    def to_dict(self) -> dict[str, Any]:
        return {
            "bp_id": self.bp_id,
            "tag": self.tag,
            "worst": self.worst,
            "checks": [asdict(c) for c in self.checks],
            "stats": self.stats,
        }


# ---------------------------------------------------------------------------
# Loaders / helpers
# ---------------------------------------------------------------------------


def _columns(tbl: pq.lib.Table) -> tuple[str, str]:
    cols = tbl.column_names
    seq_col = next((c for c in ("sequence", "seq", "aptamer") if c in cols), None)
    count_col = next((c for c in ("count", "reads", "n", "abundance") if c in cols), None)
    if seq_col is None or count_col is None:
        raise ValueError(f"missing seq/count column; got cols={cols}")
    return seq_col, count_col


def load_parquet_full(path: Path) -> tuple[list[str], list[int]]:
    """Read the entire parquet — no sampling."""
    tbl = pq.read_table(path)
    seq_col, count_col = _columns(tbl)
    return tbl.column(seq_col).to_pylist(), tbl.column(count_col).to_pylist()


def round_parquets(bp_dir: Path, kind: str) -> list[Path]:
    return sorted(bp_dir.glob(f"round_*.{kind}.parquet"))


def round_index(p: Path) -> int:
    m = re.search(r"round_(\d+)", p.name)
    return int(m.group(1)) if m else -1


def final_enrich_parquet(bp_dir: Path) -> Path | None:
    """Pick the widest-span enrich parquet (final round_0 → round_last comparison)."""
    enrichs = sorted(bp_dir.glob("enrich_*.parquet"))
    if not enrichs:
        return None

    def span(p: Path) -> int:
        m = re.search(r"(?:round_)?(\d+).*?(?:to_)?(?:round_)?(\d+)", p.name)
        if m:
            return abs(int(m.group(2)) - int(m.group(1)))
        return 0

    return max(enrichs, key=lambda p: (span(p), p.stat().st_size))


# ---------------------------------------------------------------------------
# Sequence-level analysis primitives (pure, importable)
# ---------------------------------------------------------------------------


def length_distribution(seqs: list[str]) -> tuple[int, float, list[tuple[int, float]]]:
    """Return (mode_length, mode_fraction, top5_lengths_with_fractions)."""
    if not seqs:
        return 0, 0.0, []
    Lc = Counter(len(s) for s in seqs)
    total = sum(Lc.values())
    sorted_L = Lc.most_common()
    mode_L, mode_n = sorted_L[0]
    return mode_L, mode_n / total, [(L, n / total) for L, n in sorted_L[:5]]


def alphabet_violations(seqs: list[str]) -> int:
    return sum(1 for s in seqs if any(c not in ALPHABET for c in s))


def gc_content(seqs: list[str]) -> tuple[float, float]:
    """Return (mean_gc, stddev_gc) across sequences."""
    if not seqs:
        return 0.0, 0.0
    gcs = [(s.count("G") + s.count("C")) / len(s) for s in seqs if s]
    if not gcs:
        return 0.0, 0.0
    mean = sum(gcs) / len(gcs)
    var = sum((g - mean) ** 2 for g in gcs) / len(gcs)
    return mean, math.sqrt(var)


def positional_top_base(seqs: list[str], from_end: bool, n_pos: int) -> list[tuple[str, float]]:
    """Per-position top base + fraction across `seqs` over `n_pos` positions.

    When ``from_end=True``, indexes from the 3' end (variable-length-tolerant).
    """
    counts = [Counter() for _ in range(n_pos)]
    for s in seqs:
        L = len(s)
        for i in range(min(n_pos, L)):
            idx = (n_pos - 1 - i) if from_end else i
            char = s[L - 1 - i] if from_end else s[i]
            counts[idx][char] += 1
    out: list[tuple[str, float]] = []
    for c in counts:
        if not c:
            out.append(("-", 0.0))
            continue
        b = max(c, key=c.get)
        out.append((b, c[b] / sum(c.values())))
    return out


def positional_entropy(seqs: list[str], from_end: bool, n_pos: int) -> list[float]:
    """Per-position Shannon entropy (bits) across `seqs` over `n_pos` positions."""
    counts = [Counter() for _ in range(n_pos)]
    for s in seqs:
        L = len(s)
        for i in range(min(n_pos, L)):
            idx = (n_pos - 1 - i) if from_end else i
            char = s[L - 1 - i] if from_end else s[i]
            counts[idx][char] += 1
    out: list[float] = []
    for c in counts:
        total = sum(c.values())
        if total == 0:
            out.append(0.0)
            continue
        H = -sum((v / total) * math.log2(v / total) for v in c.values() if v > 0)
        out.append(H)
    return out


def top_kmers(
    seqs: list[str], counts: list[int], k: int, top_n: int = 10
) -> list[tuple[str, float]]:
    """Top-N k-mers across the pool, weighted by read count."""
    total: Counter = Counter()
    for s, c in zip(seqs, counts, strict=False):
        if len(s) < k:
            continue
        for i in range(len(s) - k + 1):
            total[s[i : i + k]] += c
    n_total = sum(total.values())
    if n_total == 0:
        return []
    return [(km, n / n_total) for km, n in total.most_common(top_n)]


def truseq_contamination(seqs: list[str], counts: list[int]) -> tuple[float, float, int]:
    """Return (unique_frac, reads_frac, n_reads_with_truseq)."""
    if not seqs:
        return 0.0, 0.0, 0
    unique_hits = sum(1 for s in seqs if TRUSEQ_R1_PREFIX in s)
    reads_hits = sum(c for s, c in zip(seqs, counts, strict=False) if TRUSEQ_R1_PREFIX in s)
    total_reads = sum(counts)
    return (
        (unique_hits / len(seqs)) if seqs else 0.0,
        (reads_hits / total_reads) if total_reads else 0.0,
        reads_hits,
    )


def singleton_fraction(counts: list[int]) -> float:
    """Fraction of unique sequences observed exactly once."""
    if not counts:
        return 0.0
    return sum(1 for c in counts if c == 1) / len(counts)


# ---------------------------------------------------------------------------
# Section checks (each mutates `report` in place)
# ---------------------------------------------------------------------------


def _section_pre(
    bp_dir: Path,
    report: BPReport,
    rounds_data: dict[Path, tuple[list[str], list[int]]],
) -> bool:
    counts = round_parquets(bp_dir, "counts")
    clusters = round_parquets(bp_dir, "clusters")
    enrich = list(bp_dir.glob("enrich_*.parquet"))
    summary = bp_dir / "summary.json"
    cstats = bp_dir / "cluster_stats.json"

    missing = []
    if not counts:
        missing.append("counts parquets")
    if not clusters:
        missing.append("clusters parquets")
    if not enrich:
        missing.append("enrich parquets")
    if not summary.exists():
        missing.append("summary.json")
    if not cstats.exists():
        missing.append("cluster_stats.json")

    if missing:
        report.add("pre", FAIL, f"missing: {', '.join(missing)}")
        return False

    for p in counts:
        try:
            rounds_data[p] = load_parquet_full(p)
        except Exception as exc:
            report.add("pre", FAIL, f"{p.name} read error: {exc}")
            return False

    n_rounds = len(rounds_data)
    n_unique_total = sum(len(s) for s, _ in rounds_data.values())
    n_reads_total = sum(sum(c) for _, c in rounds_data.values())
    report.stats["n_rounds"] = n_rounds
    report.stats["n_unique_total"] = n_unique_total
    report.stats["n_reads_total"] = n_reads_total
    report.stats["n_enrich"] = len(enrich)
    report.add(
        "pre",
        PASS,
        f"{n_rounds} rounds; {n_unique_total:,} unique seqs; {n_reads_total:,} reads",
    )
    return True


def _section_alphabet(
    report: BPReport,
    rounds_data: dict[Path, tuple[list[str], list[int]]],
) -> None:
    total_bad = 0
    total = 0
    bad_rounds = []
    for p, (seqs, _) in sorted(rounds_data.items()):
        bad = alphabet_violations(seqs)
        total_bad += bad
        total += len(seqs)
        if bad > 0:
            bad_rounds.append((p.name, bad, len(seqs)))

    report.stats["alphabet_total_seqs"] = total
    report.stats["alphabet_bad_seqs"] = total_bad

    if total == 0:
        report.add("alphabet", FAIL, "no sequences to check")
        return
    frac = total_bad / total
    if frac > ALPHABET_BAD_FRAC_FAIL:
        details = "; ".join(f"{n}: {b}/{t}" for n, b, t in bad_rounds[:3])
        report.add(
            "alphabet",
            FAIL,
            f"{total_bad:,}/{total:,} ({frac:.3%}) non-ACGT across "
            f"{len(bad_rounds)} rounds: {details}" + (" …" if len(bad_rounds) > 3 else ""),
        )
    elif frac > ALPHABET_BAD_FRAC_WARN:
        report.add(
            "alphabet",
            WARN,
            f"{total_bad:,}/{total:,} ({frac:.3%}) non-ACGT in {len(bad_rounds)} round(s)",
        )
    elif total_bad > 0:
        report.add(
            "alphabet",
            INFO,
            f"{total_bad}/{total:,} ({frac:.4%}) non-ACGT — sequencing-error N's, "
            "filterable downstream",
        )
    else:
        report.add(
            "alphabet",
            PASS,
            f"{total:,} unique seqs across {len(rounds_data)} rounds all pure ACGT",
        )


def _section_lengths(
    random_region_len: int,
    tag: str,
    report: BPReport,
    rounds_data: dict[Path, tuple[list[str], list[int]]],
) -> None:
    rrl = random_region_len or 0

    sorted_rounds = sorted(rounds_data.items())
    per_round = []
    issues = []

    relaxed = tag in (
        "variable_length_undoped",
        "no_r0_post_selection",
        "doped_3pct",
        "doped_8pct",
        "multiplexed_origin",
    )
    min_mode_frac_earliest = 0.30 if relaxed else 0.60

    for i, (p, (seqs, _)) in enumerate(sorted_rounds):
        mode_L, mode_frac, top5 = length_distribution(seqs)
        per_round.append(
            {
                "round": p.name,
                "n_unique": len(seqs),
                "mode_L": mode_L,
                "mode_frac": round(mode_frac, 3),
                "top5_L": [(L, round(f, 3)) for L, f in top5],
            }
        )

        if abs(mode_L - rrl) > MODE_L_TOLERANCE:
            issues.append(f"{p.name} mode_L={mode_L} vs rrl={rrl}")

        if i == 0 and mode_frac < min_mode_frac_earliest:
            issues.append(
                f"{p.name} (earliest) mode covers {mode_frac:.0%} (<{min_mode_frac_earliest:.0%})"
            )

        if rrl:
            for L, f in top5:
                if rrl * 0.5 > L and f > 0.005:
                    issues.append(f"{p.name} short outlier L={L} at {f:.1%}")
                    break

    report.stats["per_round_lengths"] = per_round
    report.stats["rrl_seed"] = rrl

    mode_values = [r["mode_L"] for r in per_round]
    if mode_values and max(mode_values) - min(mode_values) > MODE_L_TOLERANCE:
        issues.append(f"mode_L inconsistent across rounds: {min(mode_values)}-{max(mode_values)}")

    if issues:
        report.add(
            "lengths",
            WARN,
            "; ".join(issues[:3]) + (" …" if len(issues) > 3 else ""),
        )
    else:
        modes = ",".join(str(r["mode_L"]) for r in per_round)
        report.add(
            "lengths",
            PASS,
            f"all rounds within rrl={rrl}±{MODE_L_TOLERANCE}; modes=[{modes}]",
        )


def _section_trim_seq(
    primer_5p: str,
    primer_3p: str,
    report: BPReport,
    rounds_data: dict[Path, tuple[list[str], list[int]]],
) -> None:
    worst_5p = 0.0
    worst_3p = 0.0
    worst_truseq_unique = 0.0
    worst_truseq_reads = 0.0
    total_truseq_reads = 0
    per_round = []

    for p, (seqs, counts) in sorted(rounds_data.items()):
        if not seqs:
            continue
        f5 = f3 = 0.0
        if primer_5p:
            pfx = primer_5p[: min(PRIMER_MATCH_PREFIX_LEN, len(primer_5p))]
            f5 = sum(1 for s in seqs if s.startswith(pfx)) / len(seqs)
        if primer_3p:
            sfx = primer_3p[: min(PRIMER_MATCH_PREFIX_LEN, len(primer_3p))]
            f3 = sum(1 for s in seqs if s.endswith(sfx)) / len(seqs)
        tru_uniq, tru_reads, tru_n = truseq_contamination(seqs, counts)
        per_round.append(
            {
                "round": p.name,
                "frac_starts_5p": round(f5, 5),
                "frac_ends_3p": round(f3, 5),
                "truseq_unique_frac": round(tru_uniq, 5),
                "truseq_reads_frac": round(tru_reads, 5),
                "truseq_n_reads": tru_n,
            }
        )
        worst_5p = max(worst_5p, f5)
        worst_3p = max(worst_3p, f3)
        worst_truseq_unique = max(worst_truseq_unique, tru_uniq)
        worst_truseq_reads = max(worst_truseq_reads, tru_reads)
        total_truseq_reads += tru_n

    report.stats["per_round_trim"] = per_round
    report.stats["worst_starts_5p_frac"] = round(worst_5p, 5)
    report.stats["worst_ends_3p_frac"] = round(worst_3p, 5)
    report.stats["worst_truseq_reads_frac"] = round(worst_truseq_reads, 5)
    report.stats["total_truseq_reads"] = total_truseq_reads

    issues = []
    if worst_5p > PRIMER_RESIDUE_MAX_FRAC:
        issues.append(f"max {worst_5p:.1%} reads start with primer_5p")
    if worst_3p > PRIMER_RESIDUE_MAX_FRAC:
        issues.append(f"max {worst_3p:.1%} reads end with primer_3p")
    if worst_truseq_reads > TRUSEQ_MAX_FRAC_READS:
        issues.append(f"max TruSeq R1 reads={worst_truseq_reads:.1%}")

    if issues:
        report.add("trim_seq", FAIL, "; ".join(issues))
    else:
        report.add(
            "trim_seq",
            PASS,
            f"all {len(per_round)} rounds clean "
            f"(max TruSeq reads={worst_truseq_reads:.4%}, "
            f"max starts_5p={worst_5p:.4%}, max ends_3p={worst_3p:.4%})",
        )


def _section_composition(
    tag: str,
    report: BPReport,
    rounds_data: dict[Path, tuple[list[str], list[int]]],
) -> None:
    earliest = sorted(rounds_data.keys())[0]
    seqs, counts = rounds_data[earliest]
    mode_L, _, _ = length_distribution(seqs)
    mode_seqs = [s for s in seqs if len(s) == mode_L]
    if not mode_seqs:
        report.add("composition", FAIL, "no mode-L sequences")
        return

    head = positional_top_base(mode_seqs, from_end=False, n_pos=10)
    tail = positional_top_base(mode_seqs, from_end=True, n_pos=10)
    gc_mean, gc_std = gc_content(mode_seqs)
    report.stats["composition_earliest_round"] = earliest.name
    report.stats["composition_mode_L_seqs"] = len(mode_seqs)
    report.stats["gc_mean"] = round(gc_mean, 3)
    report.stats["gc_std"] = round(gc_std, 3)
    report.stats["head_topbase"] = [f"{b}:{p:.0%}" for b, p in head]
    report.stats["tail_topbase"] = [f"{b}:{p:.0%}" for b, p in tail]

    if tag == "no_r0_post_selection":
        top5 = sorted(zip(counts, seqs, strict=False), reverse=True)[:5]
        prefixes_6mer = {s[:6] for _, s in top5 if len(s) >= 6}
        report.stats["top5_prefix6_distinct"] = len(prefixes_6mer)
        if len(prefixes_6mer) >= 3:
            report.add(
                "composition",
                PASS,
                f"no_r0; top-5 6-mer prefix diversity = {len(prefixes_6mer)} distinct; "
                f"GC={gc_mean:.0%}±{gc_std:.0%}",
            )
        else:
            report.add(
                "composition",
                WARN,
                f"no_r0; only {len(prefixes_6mer)} distinct top-5 6-mer prefixes — residue suspect",
            )
        return

    is_doped = tag in ("doped_3pct", "doped_8pct")
    threshold = DOPED_POS_TOPBASE_MAX if is_doped else UNDOPED_POS_TOPBASE_MAX

    over = []
    for i, (b, p) in enumerate(head):
        if p > threshold:
            over.append(f"head[{i}]={b}:{p:.0%}")
    for i, (b, p) in enumerate(tail):
        if p > threshold:
            over.append(f"tail[-{10 - i}]={b}:{p:.0%}")

    head_max = max((p for _, p in head), default=0)
    tail_max = max((p for _, p in tail), default=0)

    H_head = positional_entropy(mode_seqs, from_end=False, n_pos=10)
    H_tail = positional_entropy(mode_seqs, from_end=True, n_pos=10)
    H_mean = (
        (sum(H_head) + sum(H_tail)) / (len(H_head) + len(H_tail)) if (H_head or H_tail) else 0.0
    )
    report.stats["composition_entropy_mean_bits"] = round(H_mean, 2)

    detail = (
        f"head_max={head_max:.0%}, tail_max={tail_max:.0%}, "
        f"GC={gc_mean:.0%}±{gc_std:.0%}, H_mean={H_mean:.2f} "
        f"(n_mode_L_seqs={len(mode_seqs):,})"
    )

    if over:
        if not is_doped and H_mean >= UNDOPED_MIN_ENTROPY_BITS:
            report.add(
                "composition",
                WARN,
                f"top-base over {threshold:.0%} at: {'; '.join(over[:4])}"
                f"{' …' if len(over) > 4 else ''} — but H_mean={H_mean:.2f} bits ≥ "
                f"{UNDOPED_MIN_ENTROPY_BITS} = library bias, not residue ({detail})",
            )
        else:
            report.add(
                "composition",
                FAIL,
                f"top-base over {threshold:.0%} at: {'; '.join(over[:4])}"
                f"{' …' if len(over) > 4 else ''} ({detail})",
            )
    elif gc_mean < 0.25 or gc_mean > 0.75:
        report.add("composition", WARN, f"GC={gc_mean:.0%} unusual ({detail})")
    else:
        report.add("composition", PASS, detail)


def _section_diversity(
    tag: str,
    report: BPReport,
    rounds_data: dict[Path, tuple[list[str], list[int]]],
) -> None:
    earliest = sorted(rounds_data.keys())[0]
    seqs, counts = rounds_data[earliest]
    mode_L, _, _ = length_distribution(seqs)
    mode_seqs = [s for s in seqs if len(s) == mode_L]
    if not mode_seqs:
        report.add("diversity", FAIL, "no mode-L sequences")
        return

    H_head = positional_entropy(mode_seqs, from_end=False, n_pos=10)
    H_tail = positional_entropy(mode_seqs, from_end=True, n_pos=10)
    H_head_mean = sum(H_head) / len(H_head) if H_head else 0
    H_tail_mean = sum(H_tail) / len(H_tail) if H_tail else 0
    report.stats["entropy_head_mean_bits"] = round(H_head_mean, 2)
    report.stats["entropy_tail_mean_bits"] = round(H_tail_mean, 2)

    singletons = singleton_fraction(counts)
    report.stats["singleton_frac"] = round(singletons, 3)

    kmers = top_kmers(seqs, counts, TOP_KMER_K, top_n=5)
    report.stats["top_kmers"] = [(k, round(f, 4)) for k, f in kmers]
    top_kmer_frac = kmers[0][1] if kmers else 0.0

    issues = []
    if tag in (
        "undoped",
        "standard",
        "variable_length_undoped",
        "multiplexed_origin",
        "aptasim_excluded",
    ):
        if H_head_mean < 1.5:
            issues.append(f"head entropy={H_head_mean:.2f} bits (low; expected ~2.0)")
        if H_tail_mean < 1.5:
            issues.append(f"tail entropy={H_tail_mean:.2f} bits (low)")
        if top_kmer_frac > TOP_KMER_MAX_FRAC_UNDOPED:
            issues.append(
                f"top 6-mer={kmers[0][0]} at {top_kmer_frac:.1%} "
                f"(>{TOP_KMER_MAX_FRAC_UNDOPED:.0%}; unusual)"
            )

    detail = (
        f"H_head={H_head_mean:.2f}, H_tail={H_tail_mean:.2f}, "
        f"singletons={singletons:.0%}, "
        f"top6mer={kmers[0][0] if kmers else '-'}@{top_kmer_frac:.1%}"
    )
    if issues:
        report.add("diversity", WARN, "; ".join(issues) + f" ({detail})")
    else:
        report.add("diversity", PASS, detail)


def _section_selection(
    bp_dir: Path,
    report: BPReport,
    rounds_data: dict[Path, tuple[list[str], list[int]]],
) -> None:
    enrich = final_enrich_parquet(bp_dir)
    if enrich is None:
        report.add("selection", FAIL, "no enrich parquet")
        return

    try:
        tbl = pq.read_table(enrich)
    except Exception as exc:
        report.add("selection", FAIL, f"enrich read error: {exc}")
        return

    cols = tbl.column_names
    log2fc_col = next((c for c in ("log2FC", "log2_fold_change", "log2fc") if c in cols), None)
    if log2fc_col is None:
        report.add("selection", WARN, f"no log2FC column in {enrich.name}; cols={cols}")
        return

    vals = tbl.column(log2fc_col).to_pylist()
    n_seqs = len(vals)
    if n_seqs == 0:
        report.add("selection", FAIL, "enrich parquet empty")
        return

    max_log2fc = max(vals)
    n_positive = sum(1 for v in vals if v > 0)
    frac_positive = n_positive / n_seqs
    report.stats["enrich_file"] = enrich.name
    report.stats["max_log2fc_final"] = round(max_log2fc, 2)
    report.stats["n_enrich_seqs"] = n_seqs
    report.stats["enrich_frac_positive"] = round(frac_positive, 3)

    rounds_sorted = sorted(rounds_data.keys())
    if len(rounds_sorted) >= 2:
        f_seqs, _ = rounds_data[rounds_sorted[0]]
        l_seqs, l_counts = rounds_data[rounds_sorted[-1]]
        top5_last = sorted(zip(l_counts, l_seqs, strict=False), reverse=True)[:5]
        top5_last_set = {s for _, s in top5_last}
        f_seq_set = set(f_seqs)
        trace = sum(1 for s in top5_last_set if s in f_seq_set)
        report.stats["top5_last_present_in_first"] = trace
    else:
        trace = -1

    issues = []
    if max_log2fc < MIN_MAX_LOG2FC:
        issues.append(f"max_log2FC={max_log2fc:.2f} < {MIN_MAX_LOG2FC}")
    if n_seqs < MIN_POST_NOISE_SEQS:
        issues.append(f"post-noise n_seqs={n_seqs} < {MIN_POST_NOISE_SEQS}")

    detail = (
        f"max_log2FC={max_log2fc:.2f}, n_seqs={n_seqs:,}, "
        f"frac_positive={frac_positive:.0%}, top5_traceback={trace}"
    )
    if issues:
        report.add("selection", WARN, "; ".join(issues) + f" ({detail})")
    else:
        report.add("selection", PASS, detail)


def _section_consistency(
    bp_dir: Path,
    primer_5p: str,
    primer_3p: str,
    report: BPReport,
    rounds_data: dict[Path, tuple[list[str], list[int]]],
) -> None:
    ratios = []
    for p, (seqs, cnts) in sorted(rounds_data.items()):
        n_unique = len(seqs)
        n_reads = sum(cnts)
        ratios.append(
            {
                "round": p.name,
                "n_unique": n_unique,
                "n_reads": n_reads,
                "uniq_per_read": round(n_unique / n_reads, 5) if n_reads else 0.0,
            }
        )
    report.stats["per_round_ratios"] = ratios

    if len(ratios) >= 3:
        first_r = ratios[0]["uniq_per_read"]
        last_r = ratios[-1]["uniq_per_read"]
        report.stats["uniq_per_read_first"] = first_r
        report.stats["uniq_per_read_last"] = last_r
        if last_r > first_r * 1.10:
            report.add(
                "consistency",
                WARN,
                f"unique/reads ratio INCREASES across rounds "
                f"({first_r:.4f}→{last_r:.4f}; selection signature absent)",
            )
            return

    sj_path = bp_dir / "summary.json"
    try:
        sj = json.loads(sj_path.read_text())
    except Exception as exc:
        report.add("consistency", WARN, f"summary.json parse error: {exc}")
        return

    j_5p = sj.get("primer_5p")
    j_3p = sj.get("primer_3p")
    if j_5p is None and j_3p is None:
        report.add(
            "consistency",
            INFO,
            "summary.json primers=None (known stale from cluster/enrich recompute)",
        )
        return

    drift = []
    if j_5p and j_5p != primer_5p:
        drift.append("5p drift: seed≠summary")
    if j_3p and j_3p != primer_3p:
        drift.append("3p drift: seed≠summary")
    if drift:
        report.add("consistency", WARN, "; ".join(drift))
    else:
        report.add("consistency", PASS, "round trend monotone & seed↔summary consistent")


# ---------------------------------------------------------------------------
# Top-level per-BP entry point
# ---------------------------------------------------------------------------


def review_bioproject(
    bp_dir: Path,
    bp_id: str | None = None,
    primer_5p: str = "",
    primer_3p: str = "",
    random_region_len: int | None = None,
    tag: str = "untagged",
    rrl_source: str = "",
) -> BPReport:
    """Run all eight readiness sections on a single processed BP directory.

    `bp_dir` must contain ``round_*.counts.parquet``,
    ``round_*.clusters.parquet``, at least one ``enrich_*.parquet``, plus
    ``summary.json`` and ``cluster_stats.json`` for full coverage.

    `primer_5p`, `primer_3p`, `random_region_len`, `tag` and `rrl_source` are
    seed-derived configuration; the function never reads YAML files itself.
    `tag` drives the composition/diversity thresholds (see module docstring
    for accepted values).
    """
    if bp_id is None:
        bp_id = bp_dir.name
    report = BPReport(bp_id=bp_id, tag=tag)
    rounds_data: dict[Path, tuple[list[str], list[int]]] = {}

    if not _section_pre(bp_dir, report, rounds_data):
        return report

    # Seed authority quick-check (tag-aware)
    primer_5p_ok = bool(primer_5p) and len(primer_5p) >= PRIMER_MIN_LEN
    primer_3p_required = tag not in NO_PRIMER_3P_TAGS
    primer_3p_ok = (bool(primer_3p) and len(primer_3p) >= PRIMER_MIN_LEN) or not primer_3p_required
    rrl_ok = bool(random_region_len)

    src_ok = (
        rrl_source in VALID_RRL_SOURCES_LITERAL
        or rrl_source.startswith("post_rebuild_verified")
        or rrl_source.startswith("empirical_audit_")
    )

    if not (primer_5p_ok and primer_3p_ok and rrl_ok):
        details = []
        if not primer_5p_ok:
            details.append(f"5p missing or <{PRIMER_MIN_LEN} nt")
        if not primer_3p_ok:
            details.append(f"3p missing or <{PRIMER_MIN_LEN} nt (and tag={tag!r} requires it)")
        if not rrl_ok:
            details.append("rrl missing")
        report.add("pre", FAIL, "incomplete seed: " + "; ".join(details))
    elif rrl_source and not src_ok:
        report.add(
            "pre",
            WARN,
            f"seed rrl source = {rrl_source!r} (expected paper, "
            "post_rebuild_verified, or empirical_audit_*)",
        )
    elif not primer_3p and tag in NO_PRIMER_3P_TAGS:
        report.add("pre", INFO, f"5p-only mode (tag={tag}; primer_3p intentionally absent)")

    _section_alphabet(report, rounds_data)
    _section_lengths(random_region_len or 0, tag, report, rounds_data)
    _section_trim_seq(primer_5p, primer_3p, report, rounds_data)
    _section_composition(tag, report, rounds_data)
    _section_diversity(tag, report, rounds_data)
    _section_selection(bp_dir, report, rounds_data)
    _section_consistency(bp_dir, primer_5p, primer_3p, report, rounds_data)

    rounds_data.clear()
    return report


def reports_to_json(reports: list[BPReport]) -> list[dict]:
    """Serialise a list of BPReports into JSON-friendly dicts."""
    return [r.to_dict() for r in reports]
