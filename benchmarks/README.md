# Benchmarks

This directory contains rewrite-native benchmark scripts that track active code paths.

## Active scripts

- `measure_kinase_scoring_runtime_alignment.py`
  - Measures kinase workflow runtime and memory for filtered motif-lane execution versus unfiltered and diagnostic-table variants.
  - Targets `phospy.workflows.kinase.executor` via `KinaseWorkflow.run`.
- `measure_signalome_prediction_hot_paths.py`
  - Measures signalome module-table, expanded-signalome, and prediction-output science hot paths against legacy-style baseline implementations.
  - Targets `phospy.workflows.signalome.science.build_signalome_module_table`, `phospy.workflows.signalome.science.build_expanded_signalome_table`, and `phospy.workflows.kinase.science.build_prediction_outputs`.

## Drift guard

CI runs `pytest tests/unit/test_benchmark_scripts_smoke.py` to ensure benchmark scripts:

- import cleanly,
- keep rewrite-native local import targets resolvable,
- include a short module header describing what they measure.
