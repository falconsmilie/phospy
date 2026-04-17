# ADR 0005: Keep Per-Column Candidate Selection Instead of Full-Matrix Vectorization

- Status: Accepted
- Date: 2026-04-14

## Context

Candidate substrate selection in `phospy.prediction.candidates` picks top-scoring
sites per kinase, then applies score-threshold and inclusion rules.

A review recommendation suggested replacing the per-column selection path with a
single full-matrix vectorized pass using `np.argpartition(..., axis=0)`.

The goal was to reduce Python-loop overhead for large kinase panels.

## Decision

PhosPy keeps the per-column implementation:

- one `np.argpartition` per kinase column
- deterministic tie ordering by score, then row position
- existing threshold/inclusion semantics unchanged

The full-matrix vectorized variant is not adopted.

## Evidence Reviewed

The full-matrix variant was implemented and benchmarked against the current
per-column implementation on representative random score matrices. Outputs were
verified equal for each case.

Measured means (local benchmark run on 2026-04-14):

| Rows | Cols | Top | Per-column (ms) | Full-matrix (ms) | Relative |
| --- | --- | --- | ---: | ---: | ---: |
| 2,000 | 100 | 50 | 1.14 | 1.61 | 0.71x |
| 5,000 | 200 | 50 | 3.75 | 6.68 | 0.56x |
| 10,000 | 300 | 50 | 8.00 | 17.08 | 0.47x |
| 5,000 | 500 | 50 | 8.76 | 15.89 | 0.55x |
| 5,000 | 1,000 | 50 | 17.40 | 42.67 | 0.41x |
| 20,000 | 500 | 100 | 25.67 | 59.94 | 0.43x |
| 20,000 | 1,000 | 100 | 51.90 | 176.26 | 0.29x |

Interpretation: the attempted full-matrix approach increased runtime due to
additional whole-matrix ordering and mask/materialization costs.

## Consequences

Positive:

- preserves the faster observed path for current workloads
- avoids shipping a refactor that is "more vectorized" but slower
- keeps candidate semantics stable

Trade-off:

- retains a small explicit per-column loop in Python

## Follow-On Guidance

- Future optimizations must be benchmark-backed before replacing this path.
- If matrix characteristics change materially, rerun targeted candidate-selection
  benchmarks and revisit this ADR.
