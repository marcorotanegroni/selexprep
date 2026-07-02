"""Inter-extractor agreement between two dual-extraction JSON outputs.

Usage:
    python compute_agreement.py claude_extractions.json codex_extractions.json \
        [--norm normalization.yaml] [--outdir out]

Each input is a JSON list of per-deposit records (see README.md). Values are
normalized via ``normalize.py`` and compared per field. Writes:
  - ``agreement_summary.tsv`` — per-field counts + agreement rates
  - ``disagreements.tsv``     — every mismatch with both values + both quotes,
                                the worklist for human adjudication

Pre-registered metric definitions (computed per field, over deposits present in
both inputs):
  - ``both_absent``        — both extractors said not_stated (concordant absence)
  - ``agree``              — both stated, normalized values equal
  - ``disagree``           — both stated, normalized values differ
  - ``presence_mismatch``  — one stated, the other not_stated
  - ``raw_agreement``      — (agree + both_absent) / total
  - ``substantive_agreement`` — agree / (agree + disagree + presence_mismatch),
                                i.e. agreement where at least one extractor stated
                                a value (None if that denominator is 0)
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import normalize as N

FIELDS = [
    "study_type",
    "target",
    "target_class",
    "chemistry",
    "n_random",
    "n_rounds",
    "selection_format",
    "counter_selection",
]


def _index(records: list[dict]) -> dict[str, dict]:
    return {r["accession"]: r for r in records}


def _field(rec: dict, field: str) -> dict:
    return (rec.get("fields") or {}).get(field) or {}


def _raw_value(f: dict) -> object:
    return None if f.get("status") == "not_stated" else f.get("value")


def _cite(f: dict) -> str:
    v = _raw_value(f)
    if v is None:
        return "not_stated"
    q = f.get("evidence_quote") or ""
    loc = f.get("location") or ""
    return f'{v}  |  "{q}" @ {loc}'.strip()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("claude", help="JSON list of Claude extraction records")
    ap.add_argument("codex", help="JSON list of Codex extraction records")
    ap.add_argument("--norm", default=None, help="path to normalization.yaml")
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()

    a = _index(json.loads(Path(args.claude).read_text(encoding="utf-8")))
    b = _index(json.loads(Path(args.codex).read_text(encoding="utf-8")))
    accs = sorted(set(a) & set(b))
    missing = set(a) ^ set(b)
    if missing:
        print(f"WARNING: {len(missing)} accession(s) not in both inputs: {sorted(missing)}")
    if not accs:
        print("ERROR: no overlapping accessions between the two inputs.")
        return 1

    stats = {
        f: dict(agree=0, disagree=0, both_absent=0, presence_mismatch=0, unmapped=0) for f in FIELDS
    }
    disagreements: list[list[str]] = []

    for acc in accs:
        for field in FIELDS:
            fa, fb = _field(a[acc], field), _field(b[acc], field)
            va, flag_a = N.normalize(field, _raw_value(fa), table_path=args.norm)
            vb, flag_b = N.normalize(field, _raw_value(fb), table_path=args.norm)
            if "UNMAPPED" in (flag_a, flag_b):
                stats[field]["unmapped"] += 1
            present_a, present_b = va is not None, vb is not None
            if not present_a and not present_b:
                stats[field]["both_absent"] += 1
            elif present_a != present_b:
                stats[field]["presence_mismatch"] += 1
                disagreements.append([acc, field, "presence_mismatch", _cite(fa), _cite(fb)])
            elif va == vb:
                stats[field]["agree"] += 1
            else:
                stats[field]["disagree"] += 1
                disagreements.append([acc, field, "value_mismatch", _cite(fa), _cite(fb)])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summary_path = outdir / "agreement_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(
            [
                "field",
                "n",
                "agree",
                "disagree",
                "presence_mismatch",
                "both_absent",
                "unmapped",
                "raw_agreement",
                "substantive_agreement",
            ]
        )
        tot = dict(agree=0, disagree=0, both_absent=0, presence_mismatch=0, unmapped=0)
        for field in FIELDS:
            s = stats[field]
            for k in tot:
                tot[k] += s[k]
            n = len(accs)
            raw = (s["agree"] + s["both_absent"]) / n if n else 0.0
            denom = s["agree"] + s["disagree"] + s["presence_mismatch"]
            subst = s["agree"] / denom if denom else None
            w.writerow(
                [
                    field,
                    n,
                    s["agree"],
                    s["disagree"],
                    s["presence_mismatch"],
                    s["both_absent"],
                    s["unmapped"],
                    f"{raw:.3f}",
                    "NA" if subst is None else f"{subst:.3f}",
                ]
            )
        n_all = len(accs) * len(FIELDS)
        raw_all = (tot["agree"] + tot["both_absent"]) / n_all if n_all else 0.0
        denom_all = tot["agree"] + tot["disagree"] + tot["presence_mismatch"]
        subst_all = tot["agree"] / denom_all if denom_all else None
        w.writerow(
            [
                "ALL",
                n_all,
                tot["agree"],
                tot["disagree"],
                tot["presence_mismatch"],
                tot["both_absent"],
                tot["unmapped"],
                f"{raw_all:.3f}",
                "NA" if subst_all is None else f"{subst_all:.3f}",
            ]
        )

    dis_path = outdir / "disagreements.tsv"
    with dis_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "field", "kind", "claude", "codex"])
        w.writerows(disagreements)

    print(f"deposits compared: {len(accs)}  |  fields: {len(FIELDS)}")
    print(
        f"overall raw agreement: {raw_all:.3f}  |  substantive: "
        f"{'NA' if subst_all is None else f'{subst_all:.3f}'}"
    )
    print(f"disagreements to adjudicate: {len(disagreements)}")
    print(f"wrote {summary_path} and {dis_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
