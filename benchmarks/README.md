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
  - Measures preprocessing contract hot paths: row-median imputation, deterministic KNN imputation, median centering, and quantile normalisation.
  - Targets `phospy.science.datasets.preprocessing.stages.missing_data` and `phospy.science.datasets.preprocessing.stages.normalisation`.
- `measure_dataframe_ownership_copy_policy.py`
  - Measures internal DataFrame borrow copy counts, mutation leaks, and runtime for NumPy-backed and extension-array-backed frames.
  - Targets `phospy.frames.ownership._borrow_dataframe`.
- `measure_repeated_workflow_dataset_snapshot_reuse.py`
  - Measures repeated differential and kinase workflow runtime, tracemalloc peak memory, full-frame deep-copy counts, and dataset internal snapshot construction counts for the same unchanged dataset.
  - Targets `DifferentialAnalysisWorkflow.run`, `KinaseWorkflow.run`, and dataset-owned internal snapshot reuse.
- `measure_release_scale_builder_differential.py`
  - Measures the explicitly invoked 50,000-site x 48-sample public dataset-builder, preprocessing/provenance, and one-contrast differential workload.
  - Targets `AnalysisReadyDatasetBuilder.run` and `DifferentialAnalysisWorkflow.run`.
- `measure_signalome_clustering_contracts.py`
  - Measures exact clustering below the exact-tree guardrail, guard failure above `max_exact_tree_sites`, and candidate-scoring policy behavior (`full` vs `sampled`).
  - Targets `phospy.science.signalomes.clustering.run_signalome_clustering_engine` and scale-guard policy thresholds in `phospy.science.signalomes.clustering`.

All benchmark scripts print plain `key=value` metrics for easy CI/log parsing.
The optional release-scale builder+differential benchmark and repeated workflow
snapshot reuse benchmark can also write JSON reports under
`benchmarks/reports/`; those scratch reports remain ignored by git.

Retained release-scale benchmark evidence lives separately under
`benchmarks/evidence/`. The retained 50,000 x 48 observation is
[release-scale-builder-differential-2026-08-12.json](evidence/release-scale-builder-differential-2026-08-12.json).
It records Python, selected dependency versions, machine metadata, runtime
timings, peak RSS or an explicit unavailable state, metrics, scientific-summary
digest, output fingerprints, and source provenance for the executed checkout.
This is dated machine-specific benchmark evidence only; it is not a portable
performance guarantee, a release gate, scientific validation, or external
parity evidence.

To retain a new evidence report, run the workload explicitly and provide a
dated evidence path:

```bash
python benchmarks/measure_release_scale_builder_differential.py \
  --report-path benchmarks/evidence/release-scale-builder-differential-YYYY-MM-DD.json
```

Do not commit ordinary local scratch reports from `benchmarks/reports/`.

## Drift guard

CI runs `pytest tests/unit/test_benchmark_scripts_smoke.py` to ensure benchmark scripts:

- import cleanly,
- keep rewrite-native local import targets resolvable,
- include a short module header describing what they measure,
- print at least one `key=value` metric from `main()`,
- keep default scratch reports under `benchmarks/reports/`,
- and keep the heavy release-scale workload explicitly invoked rather than
  adding it to pytest, CI, or release gates.
