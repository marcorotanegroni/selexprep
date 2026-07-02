"""Three-way comparison: existing curation vs the two fresh extractors.

For the deposits present in both extraction files, compares each field across
three sources — the EXISTING curation (``project_metadata.json``, built during
package construction), CLAUDE, and CODEX — normalized via ``normalize.py``. It
shows how deep the original curation was relative to the fresh dual pass, which
tells you whether to tune the prompt or run a supplementary-file pass.

Buckets (per deposit x field), checked in this order:
  - all_absent      none of the three has a value
  - conflict        >=2 distinct present values (real disagreement, incl. existing
                    vs fresh) -> adjudicate; reveals original-curation errors too
  - existing_only   existing present, BOTH fresh extractors absent -> the original
                    dug deeper / the fresh prompt missed it (prompt-tune / supplement signal)
  - fresh_adds      both fresh extractors agree on a value, existing blank -> the
                    fresh pass adds coverage the original curation lacked
  - concordant      the present values agree (all three, or existing + one extractor)

Usage:
  python compare_to_existing.py project_metadata.json claude.json codex.json [--norm ...] [--outdir out]
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
CATS = ["concordant", "fresh_adds", "existing_only", "conflict", "all_absent"]


def _ext_val(fields: dict, f: str) -> object:
    fd = fields.get(f) or {}
    return None if fd.get("status") == "not_stated" else fd.get("value")


def _bucket(e, c, x) -> str:
    present = {k: v for k, v in (("e", e), ("c", c), ("x", x)) if v is not None}
    distinct = set(present.values())
    if not present:
        return "all_absent"
    if len(distinct) >= 2:
        return "conflict"
    if "e" in present and "c" not in present and "x" not in present:
        return "existing_only"
    if "e" not in present and "c" in present and "x" in present:
        return "fresh_adds"
    return "concordant"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("existing", help="project_metadata.json (existing curation)")
    ap.add_argument("claude")
    ap.add_argument("codex")
    ap.add_argument("--norm", default=None)
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()

    E = {r["accession"]: r for r in json.loads(Path(args.existing).read_text(encoding="utf-8"))}
    C = {
        r["accession"]: (r.get("fields") or {})
        for r in json.loads(Path(args.claude).read_text(encoding="utf-8"))
    }
    X = {
        r["accession"]: (r.get("fields") or {})
        for r in json.loads(Path(args.codex).read_text(encoding="utf-8"))
    }
    accs = sorted(set(C) & set(X))

    stats = {f: dict.fromkeys(CATS, 0) for f in FIELDS}
    rows: list[list[str]] = []
    for acc in accs:
        for f in FIELDS:
            e, _ = N.normalize(f, (E.get(acc) or {}).get(f), table_path=args.norm)
            c, _ = N.normalize(f, _ext_val(C[acc], f), table_path=args.norm)
            x, _ = N.normalize(f, _ext_val(X[acc], f), table_path=args.norm)
            cat = _bucket(e, c, x)
            stats[f][cat] += 1
            if cat in ("existing_only", "conflict", "fresh_adds"):
                rows.append(
                    [
                        acc,
                        f,
                        cat,
                        str((E.get(acc) or {}).get(f) or ""),
                        str(_ext_val(C[acc], f) or ""),
                        str(_ext_val(X[acc], f) or ""),
                    ]
                )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "threeway_summary.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["field", *CATS])
        tot = dict.fromkeys(CATS, 0)
        for f in FIELDS:
            w.writerow([f, *[stats[f][c] for c in CATS]])
            for c in CATS:
                tot[c] += stats[f][c]
        w.writerow(["ALL", *[tot[c] for c in CATS]])
    with (outdir / "threeway_cells.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "field", "bucket", "existing", "claude", "codex"])
        w.writerows(rows)

    print(f"deposits: {len(accs)}  fields: {len(FIELDS)}  cells: {len(accs) * len(FIELDS)}")
    for c in CATS:
        print(f"  {c}: {tot[c]}")
    print(f"wrote {outdir}/threeway_summary.tsv and {outdir}/threeway_cells.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
