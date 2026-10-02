# Supplementary Table — public-corpus audit

_selexprep public-corpus audit  ·  30 of 127 INSDC rows audit-eligible · 113 non-INSDC passthrough · 240 catalog total  ·  N=23 sampled / 21 fetchable / 21 with LibraryReport  ·  catalog v0.2.1-snapshot-2026-07-03, seed 42_

**Fetch outcomes (denominator = sampled accessions)**

| Category | Rows |
|---|---|
| OK | 10 |
| EXTRACT_REFUSED | 11 |
| FETCH_FAILED | 2 |

**Inference confidence — LibraryReport.status (denominator = rows with a LibraryReport)**

| Category | Rows |
|---|---|
| HIGH | 4 |
| MEDIUM | 6 |
| UNABLE_TO_INFER | 11 |

**Extraction mode (rows with a LibraryReport)**

| Category | Rows |
|---|---|
| BOTH_PRIMERS_SINGLE_READ | 4 |
| FIVE_PRIME_ONLY | 3 |
| THREE_PRIME_ONLY | 3 |
| UNABLE_TO_EXTRACT | 11 |

**Required action (rows with a LibraryReport)**

| Category | Rows |
|---|---|
| NONE | 10 |
| MANUAL_PRIMERS_REQUIRED | 11 |

**Inference safe-failure rate:** 52% (11/21)
