# Archived Script Utilities

These scripts are retained for historical provenance and parity forensics.
They are not part of the supported maintainer workflow.

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

## Synthetic Adaptive Edge Fixtures (historical)

- `generate_synthetic_adaptive_sampling_edge_fixtures.py`

Why this existed:

- During earlier adaptive prediction seam debugging, maintainers used a small
  synthetic fixture family to pin deterministic edge-case sampling behavior.

Current status:

- This generator is archival parity-debug tooling and not part of the active
  maintainer fixture-regeneration workflow.
