# Benchmarks

This directory contains rewrite-native benchmark scripts that track active code paths.

## Active scripts

- `measure_kinase_scoring_runtime_alignment.py`
  - Measures kinase workflow runtime and memory for filtered motif-lane execution versus unfiltered and diagnostic-table variants.
  - Targets `phospy.workflows.kinase.executor` via `KinaseWorkflow.run`.
- `measure_signalome_prediction_hot_paths.py`
  - Measures signalome module-table, expanded-signalome, prediction-output, and adaptive-prediction hot paths.
  - Targets `phospy.science.signalomes.science.build_signalome_module_table`, `phospy.science.signalomes.science.build_expanded_signalome_table`, `phospy.workflows.kinase.science.build_prediction_outputs`, and `phospy.science.prediction.execution.run_adaptive_ensemble_prediction`.
- `measure_preprocessing_performance_contracts.py`
  - Measures preprocessing contract hot paths: row-median imputation, median centering, and quantile normalisation.
  - Targets `phospy.science.datasets.preprocessing.stages.missing_data` and `phospy.science.datasets.preprocessing.stages.normalisation`.
- `measure_signalome_clustering_contracts.py`
  - Measures exact clustering below the exact-tree guardrail, guard failure above `max_exact_tree_sites`, and candidate-scoring policy behavior (`full` vs `sampled`).
  - Targets `phospy.science.signalomes.clustering.run_signalome_clustering_engine` and scale-guard policy thresholds in `phospy.science.signalomes.clustering`.

All benchmark scripts print plain `key=value` metrics for easy CI/log parsing.

## Drift guard

CI runs `pytest tests/unit/test_benchmark_scripts_smoke.py` to ensure benchmark scripts:

- import cleanly,
- keep rewrite-native local import targets resolvable,
- include a short module header describing what they measure,
- print at least one `key=value` metric from `main()`,
- and only write reports under `benchmarks/reports/` if they write files.
