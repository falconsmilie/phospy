# Archived Script Utilities

These scripts are retained for historical provenance and parity forensics.
They are not part of the supported maintainer workflow.

## Small R Fixture Regeneration (historical)

- `generate_r_fixtures.R`

Why this existed:

- During earlier migration phases, maintainers regenerated a small synthetic
  R-backed fixture family for parity forensics under `tests_legacy/fixtures/`.

Current status:

- This script is archival-only and intentionally excluded from
  `make fixtures-all`.
- Use `make fixtures-r-small-archive` only when explicit historical
  parity-forensics work is needed.

## Prediction Trace Forensics (historical)

- `export_python_prediction_traces.py`
- `diff_prediction_trace_probabilities.py`

Why these existed:

- During early donor-parity debugging, maintainers exported Python-side adaptive
  prediction seam traces and compared them against R trace tables to isolate
  candidate-selection, sampling, and learner deltas.

Current status:

- The active parity-maintenance lane is fixture-family regeneration plus replay
  and parity tests.
- These trace scripts are archived and should not be treated as required for
  routine development or release maintenance.
- Non-gated trace-debug tables promoted during earlier seam work are now
  classified as archival provenance under `tests/fixtures/archive/`.

## Synthetic Adaptive Edge Fixtures (historical)

- `generate_synthetic_adaptive_sampling_edge_fixtures.py`

Why this existed:

- During earlier adaptive prediction seam debugging, maintainers used a small
  synthetic fixture family to pin deterministic edge-case sampling behavior.

Current status:

- This generator is archival parity-debug tooling and not part of the active
  maintainer fixture-regeneration workflow.

## Fragile-Support Seam Fixture Generation (historical)

- `generate_fragile_support_reference.py`

Why this existed:

- During an earlier seam-oriented parity-debug phase, maintainers used this
  script to regenerate `fragile_support_reference` from low-level
  prediction/profile internals.

Current status:

- The script is retained for historical parity forensics only.
- It is not part of the supported fixture bootstrap path or routine maintainer
  workflow.

## L6 Seam-Stress Fixture Generation (historical)

- `generate_l6_seam_stress_reference.py`

Why this existed:

- During an earlier L6 seam-level parity-debug phase, maintainers used this
  script to slice L6 donor-backed references into a seam-stress fixture family
  and probe candidate-selection and score-combination seam behavior.

Current status:

- The script depends on older seam-oriented helper surfaces and historical
  fixture-layout assumptions.
- It is retained as archival parity-maintenance tooling only and is not part
  of the supported current maintainer lane.
