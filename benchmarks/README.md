# Benchmarks

This directory contains rewrite-native benchmark scripts that track active code paths.

## Active scripts

- `measure_kinase_scoring_runtime_alignment.py`
  - Measures kinase workflow runtime and memory for filtered motif-lane execution versus unfiltered and diagnostic-table variants.
  - Targets `phospy.workflows.kinase.executor` via `KinaseWorkflow.run`.
- `measure_signalome_prediction_hot_paths.py`
  - Measures signalome module-table, expanded-signalome, and prediction-output science hot paths against historical-baseline implementations.
  - Targets `phospy.signalomes.science.build_signalome_module_table`, `phospy.signalomes.science.build_expanded_signalome_table`, and `phospy.workflows.kinase.science.build_prediction_outputs`.
- `measure_preprocessing_performance_contracts.py`
  - Measures preprocessing contract hot paths: row-median imputation, median centering, and quantile normalisation.
  - Targets `phospy.datasets.preprocessing.stages.missing_data` and `phospy.datasets.preprocessing.stages.normalisation`.
- `measure_signalome_clustering_contracts.py`
  - Measures signalome module-selection scoring below and above the full-correlation threshold, including peak-memory snapshots and approximation-diagnostics reporting.
  - Targets `phospy.signalomes.clustering.select_module_count_with_diagnostics` and threshold constants in `phospy.signalomes.clustering`.

All benchmark scripts print plain `key=value` metrics for easy CI/log parsing.

## Drift guard

CI runs `pytest tests/unit/test_benchmark_scripts_smoke.py` to ensure benchmark scripts:

- import cleanly,
- keep rewrite-native local import targets resolvable,
- include a short module header describing what they measure,
- print at least one `key=value` metric from `main()`,
- and only write reports under `benchmarks/reports/` if they write files.
