# Performance Contracts

These contracts are practical guardrails, not promises that every machine will
finish at the same speed. CI thresholds are intentionally loose enough to avoid
noise but strict enough to catch major regressions.

## Dataset Preprocessing

The builder keeps the final `AnalysisReadyPhosphoDataset` complete and numeric.
Most default preprocessing is light because the default lane performs no
transform, no normalisation, no imputation, and no total-protein correction.

Cost increases when you enable:

- quantile normalisation
- row-median imputation over large matrices
- total-protein correction with large mapping tables
- site-matrix construction with duplicate-site aggregation

## Quantile Normalisation

Quantile normalisation sorts dense sample columns and creates additional float64
matrix copies. Expect memory use to scale with rows times samples. For very wide
or very tall matrices, run it deliberately and keep a raw input copy outside the
workflow.

CI includes a quantile-normalisation performance contract with a loose runtime
and memory ceiling.

## Kinase and Motif Scoring

Runtime depends mainly on:

- number of dataset sites
- number of kinases retained after reference overlap filtering
- number and quality of site sequences
- motif-library size

Large raw kinase-substrate maps are acceptable, but the expensive work happens
after filtering to overlapping dataset sites and valid sequences.

## Signalome Clustering

Signalome is the main scale-sensitive lane.

Default guards in `SignalomeConfig` are:

- `max_exact_tree_sites=2000`
- `max_full_candidate_scoring_sites=2000`
- `candidate_scoring_policy="full"`
- `clustering_engine="scipy_hierarchical"`

`clustering_engine="exact_python"` remains available for reference/debug checks.
Both engines still honour the exact tree guard.

`candidate_scoring_policy="sampled"` avoids materialising a full site-by-site
correlation matrix for candidate module-count scoring. It does not avoid exact
tree construction, so `max_exact_tree_sites` still matters.

If `module_count` is provided explicitly, candidate module-count scoring is
skipped, but final module assignment can still require exact tree construction.

Inspect what happened in a completed run:

```python
scale_guard = result.provenance.workflow_parameters["scale_guard"]
print(scale_guard["tree_generation_mode"])           # full_exact_tree_construction
print(scale_guard["tree_generation_is_approximate"]) # False
print(scale_guard["candidate_scoring_mode"])         # full / sampled / not_evaluated
print(scale_guard["candidate_scoring_is_approximate"])
print(scale_guard["candidate_scoring_sampled_site_total"])
print(scale_guard["candidate_scoring_sampled_pair_count"])
print(scale_guard["max_exact_tree_sites"])
print(scale_guard["max_full_candidate_scoring_sites"])
```

## Recommended Ranges

| Lane | Practical guidance |
| --- | --- |
| First run | 2 to 100 sites; rat bundled references; activity optional |
| Routine kinase scoring | usually limited by reference overlap and valid sequences |
| Quantile normalisation | check memory before running on very large dense matrices |
| Signalome full candidate scoring | stay at or below the configured full-scoring guard unless you have profiled locally |
| Signalome sampled candidate scoring | useful when candidate scoring is the bottleneck, but not a bypass for exact-tree limits |
| Very large signalome runs | reduce interpreted sites, set explicit `module_count`, or profile before raising guards |

## Benchmark Scripts

Benchmark helpers live in `benchmarks/`. They are for local review and release
work, not for ordinary first use.

Useful scripts include:

```bash
python benchmarks/measure_preprocessing_performance_contracts.py
python benchmarks/measure_kinase_scoring_runtime_alignment.py
python benchmarks/measure_signalome_prediction_hot_paths.py
python benchmarks/measure_signalome_clustering_contracts.py
```

Generated reports should stay local unless a release process asks for them.
`benchmarks/reports/` is ignored by git.

## CI Policy

Performance tests live under `tests/performance` and are marked with
`@pytest.mark.performance`. They use deterministic inputs and broad thresholds
for preprocessing, motif scoring, reference filtering, signalome backend runs,
and signalome workflow execution.

Run them locally with:

```bash
pytest tests/performance -m performance
```
