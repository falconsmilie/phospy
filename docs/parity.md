# Parity to PhosR

PhosPy parity is intentionally narrow and fixture-backed. Passing a parity test
for one lane does not mean the whole PhosR package is implemented.

PhosPy does not claim global PhosR parity. Use
`docs/scientific-coverage.md` as the single maintained scope matrix for support
status labels (`parity-gated`, `validated PhosPy implementation`,
`experimental`, `open gap`, `deliberate scope difference`, `not planned`).

Scope ownership split:

- [Scientific Coverage](scientific-coverage.md) owns user-facing coverage status,
  intended parity scope, and interpretation limits.
- [Parity](parity.md) owns fixture-level comparison evidence, test locations,
  and comparison mechanics.

## What Parity Means Here

A parity claim must say:

- which input fixture was used
- which output table or metric was compared
- which tolerance or acceptance rule was used
- which PhosPy workflow or stage produced the output

Parity evidence here should be interpreted only for the exact fixture + output +
comparison rule documented by each test lane.

## Active Parity Areas

Current active parity coverage includes:

- differential phosphorylation (`tests/parity/test_differential_analysis_parity.py`)
- differential parity envelope contracts (`tests/parity/test_differential_limma_parity.py`)
- kinase scoring and prediction surfaces (L6/public fixture lanes)
- selected preprocessing and activity-stage behaviours with explicit fixtures
- signalome workflow and clustering backend fixture lanes

## Fixture Scope By Lane

| Lane | Main fixture/evidence scope |
| --- | --- |
| Differential | Two-condition unpaired simple contrasts and related limma-envelope checks (`tests/fixtures/rewrite_parity/differential_r_reference/`, `tests/fixtures/rewrite_parity/differential_limma_envelope/`) |
| Kinase scoring/prediction | L6 and public workflow reference lanes (`tests/fixtures/rewrite_parity/r_reference_l6/`, `tests/fixtures/public_workflow_reference/`) |
| Signalome | Public workflow reference and backend parity lanes (`tests/fixtures/public_workflow_reference/`) |
| Activity parity | Activity-stage parity fixtures and threshold-bearing checks in `tests/parity/test_activity_stage_parity.py` |

Run the parity suite with:

```bash
pytest tests/parity -m parity -s
```

Release decisions should run the full release gate (`make test-release-gate`).
Parity failures in that gate are release-blocking.

Some diagnostic parity tests are informational. Release decisions should use the
threshold-bearing gates and the documented fixture expectations, not visual
inspection alone.

Release-gated parity command in the Makefile is:

```bash
pytest tests/parity -m "parity and not parity_diagnostic" -s
```

The broader CI parity smoke command is:

```bash
pytest tests/parity -m parity -s
```

## Fixture Locations

| Purpose | Location |
| --- | --- |
| Parity tests | `tests/parity/` |
| Shared parity helpers | `tests/support/` |
| Public workflow reference fixtures | `tests/fixtures/public_workflow_reference/` |
| Differential limma parity fixtures | `tests/fixtures/rewrite_parity/differential_r_reference/`, `tests/fixtures/rewrite_parity/differential_limma_envelope/` |
| Regeneration scripts | `scripts/active/` |

## Differential Parity Envelope Notes

- Differential parity claims are feature-scoped. Current limma-backed fixtures
  protect:
  - two-condition unpaired simple contrasts (`B_vs_A`, `A_vs_B`)
  - small-`n` moderated-statistics behavior
  - zero-variance and unequal-variance feature handling
  - Benjamini-Hochberg adjusted p-values and contrast ordering/sign conventions
- Differential parity comparisons use explicit floating-point tolerances in
  parity tests (`rtol=1e-6`, `atol=1e-8`).
- Missing-value handling is an intentional contract difference:
  `AnalysisReadyPhosphoDataset` requires complete matrices, so missing values
  are rejected before differential execution.

## Open Gaps

Open gaps should be described as open gaps, not as partial equivalence. Common
examples include broader organism-specific bundled references, additional PhosR
workflow surfaces, and any method not protected by fixture-backed comparison.
