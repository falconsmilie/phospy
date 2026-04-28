# Performance Contracts

PhosPy is designed for common phosphoproteomics workflow sizes, not unlimited-scale distributed execution. This page documents explicit performance contracts for the current release and the lightweight CI checks that guard against accidental regressions.

## Dataset preprocessing

Most preprocessing stages are dense `pandas`/`numpy` transforms over the phospho matrix and scale roughly with:

```text
O(n_sites x n_samples)
```

This includes:

- missing-data checks,
- row-median imputation,
- log2 transform,
- median centering,
- total-protein correction when matrices are already aligned.

Site-matrix construction and duplicate-site resolution also depend on metadata size and duplicate grouping.

You can inspect executed preprocessing operations from the dataset boundary:

```python
report = dataset.preprocessing_report
print(report.operations)
```

## Quantile normalisation

Quantile normalisation is dense and sort-heavy.

```text
Time: O(n_sites x n_samples x log(n_sites))
Memory: additional dense float64 copies of approximately n_sites x n_samples
```

Practical guidance:

Quantile normalisation is reasonable for ordinary phosphoproteomics matrices with thousands to low tens of thousands of sites and tens of samples. For much larger matrices, expect memory pressure and prefer running on a machine with enough RAM for several dense float64 matrix copies.

PhosPy does not guarantee an exact runtime for this stage.

## Kinase scoring and motif scoring

Motif scoring scales roughly with:

```text
O(n_scored_sites x n_eligible_kinases x motif_window_width)
```

Workflow contract details:

- the kinase workflow filters the motif lane to kinases that are eligible from profile-overlap scoring,
- this prevents reference-only off-lane kinases from dominating motif work when they have no dataset-overlapping substrates,
- enabling diagnostic scoring tables (`motif_scores`, `score_fusion_weights`) retains additional intermediate outputs and can increase memory/runtime overhead.

Supported downstream method names in this lane:

- `rank_weighted_fusion_scores`
- `score_fusion_weights`
- `thresholded_substrate_mean_activity`

## Large kinase-substrate references

Large reference maps are supported. Runtime depends more on post-filter overlap than on raw reference-map size.

Expected behavior:

- large maps are allowed,
- filtering by dataset/reference overlap defines eligible kinases for scoring,
- off-lane reference-only kinases should not dominate scoring outputs,
- enabling diagnostic scoring tables retains more intermediate tables and can increase memory usage.

## Signalome clustering and module selection

Current module-selection scoring contract constants:

```python
MAX_FULL_CORRELATION_SITE_COUNT = 2000
MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER = 256
SIGNALOME_MAX_EXACT_CLUSTER_TREE_SITES_DEFAULT = 2000
SIGNALOME_MAX_FULL_CORRELATION_SITES_DEFAULT = 2000
```

Behavior:

- `SignalomeConfig.cluster_tree_backend="exact"` builds the exact cluster tree,
- `SignalomeConfig.clustering_backend` selects the numerical implementation
  (`"exact_python"` or `"scipy_hierarchical"`) for that exact tree,
- exact cluster-tree construction is hard-guarded by
  `SignalomeConfig.max_exact_cluster_tree_sites` (default `2000`): runs above
  this limit fail with `SignalomeScaleError`,
- `SignalomeConfig.candidate_scoring_backend="full"` uses full site-by-site
  candidate correlations and is hard-guarded by
  `SignalomeConfig.max_full_correlation_sites`,
- `SignalomeConfig.candidate_scoring_backend="sampled"` uses sampled
  within-cluster correlation estimates for module-count candidate scoring,
- sampled candidate scoring only changes candidate module-count evaluation; it
  does not change exact cluster-tree construction or final module assignment,
- sampled scoring is deterministic and order-invariant in the current
  implementation,
- both clustering backends use the same module-selection semantics, scale guards,
  and provenance fields,
- candidate-scoring mode and exact-tree construction are recorded in
  provenance (`workflow_parameters.scale_guard`),
- provenance records requested backend separately as
  `candidate_scoring_requested_backend`,
- provenance also records whether candidate scoring was evaluated
  (`candidate_scoring_evaluated`) and why it was skipped when not evaluated
  (`candidate_scoring_skip_reason`, for example `explicit_module_count`),
- sampled runs also record `candidate_scoring_sampling` provenance
  (`sampling_cap`, `sampling_method`, `deterministic_seed_policy`,
  `actual_sampled_pair_count`, and per-cluster sample-count summary),
- explicit `module_count` runs skip candidate scoring entirely and record
  `candidate_scoring_sampling=None`,
- `tests/performance/test_performance_contracts.py` keeps a lightweight contract
  test that intentionally stubs the internal exact cluster-tree builder to isolate
  module-selection correlation-path behavior,
- `tests/performance/test_signalome_clustering_benchmark.py` adds real
  agglomerative-tree coverage for the internal exact tree-construction path on a
  deterministic medium matrix (`500` sites x `8` kinase columns) and asserts
  loose runtime/memory bounds plus tree-shape invariants.

Example:

```python
diagnostics = result.module_selection_diagnostics
print(diagnostics.reason)
```

Exact full-correlation memory floor:

```text
n_sites x n_sites x 8 bytes
```

Example: a `2000 x 2000` float64 correlation matrix is about `32 MiB` before pandas/object overhead and additional intermediate arrays.

Important limit notes:

Exact cluster-tree construction remains pairwise/agglomerative and can become expensive as site count grows. Approximate candidate scoring avoids full correlation-matrix materialisation, but does not make tree construction free.

Exact-mode guard failures are intentional scientific/runtime boundaries. Do not
silently reinterpret them as automatic approximation; choose approximation or
site filtering explicitly and record that choice in provenance.

Practical memory guidance for larger runs:

- candidate-scoring memory can be kept bounded with
  `candidate_scoring_backend="sampled"` when site count is above
  `max_full_correlation_sites`,
- exact-tree construction itself is still the dominant cost at higher site
  counts; plan for substantial transient memory and runtime growth as site count
  approaches the configured exact-tree guard,
- for multi-thousand-site runs, use a machine with headroom for dense score
  matrices, temporary clustering intermediates, and workflow tables in addition
  to raw matrix size.

### Choosing a clustering backend

`SignalomeConfig.clustering_backend` controls implementation details, not
scientific output schema:

| Backend | Practical guidance | Notes |
| --- | --- | --- |
| `exact_python` | Good default when SciPy is unavailable or when keeping dependencies minimal matters more than runtime. | Pure Python/Numpy implementation of the exact Ward tree path. |
| `scipy_hierarchical` | Prefer for medium-to-larger exact-tree runs when SciPy is available. | Uses SciPy hierarchical routines for the same exact-tree semantics; often faster in practice on larger fixtures. |

Both backends respect the same `max_exact_cluster_tree_sites` and
`max_full_correlation_sites` guards. Exceeding those guards raises
`SignalomeScaleError`; there is no silent backend fallback.

### Benchmark reporting

Signalome benchmark scripts remain manually runnable and emit parse-friendly
metrics:

- `benchmarks/measure_signalome_clustering_contracts.py` prints stable JSONL
  records (plus `key=value` header lines) for fixed fixture names:
  - `signalome_small_deterministic_v1`
  - `signalome_medium_realistic_v1`
  - `signalome_near_exact_tree_limit_v1`
  - `signalome_full_correlation_guard_v1`
  - `signalome_sampled_candidate_scoring_v1`
- each clustering record includes fixture name, backend name, site/kinase counts,
  selected module count (when available), candidate-scoring mode, sampled/skipped
  flags, exact-tree construction flag, runtime, and peak memory,
- `benchmarks/measure_signalome_prediction_hot_paths.py` uses the same JSONL
  style for hot-path timing/memory comparisons.

## Recommended ranges

| Area | Current recommended range | Behaviour above range |
| --- | ---: | --- |
| Signalome exact cluster-tree construction | up to configured `max_exact_cluster_tree_sites` (default `2,000`) | workflow fails with `SignalomeScaleError` |
| Signalome full candidate scoring | up to configured `max_full_correlation_sites` (default `2,000`) | workflow fails with `SignalomeScaleError` |
| Signalome sampled candidate scoring | low thousands of sites | runtime grows quickly because clustering is still pairwise |
| Quantile normalisation | thousands to low tens of thousands of sites, tens of samples | dense copies and sorting may become memory-heavy |
| Motif scoring | thousands of sites x hundreds of kinases | cost grows with scored sites, eligible kinases, and window width |
| Kinase-substrate reference maps | large maps are acceptable | overlap filtering should prevent off-lane kinases dominating |

These are performance contracts and operational guidance, not scientific-validity cutoffs.

## CI benchmark policy

Performance contracts are guarded by lightweight synthetic tests in `tests/performance/` and a dedicated CI job (`performance-contracts`).

Run the same suite locally with:

```bash
pytest tests/performance -m performance -q
```

Policy:

- normal unit tests remain non-timing-sensitive,
- performance tests use fixed synthetic inputs and loose thresholds/ratios,
- tests assert both behaviour and bounded runtime/memory for hot paths,
- benchmark scripts under `benchmarks/` remain manually runnable and emit `key=value` metrics for log parsing.

Threshold guidance:

- treat runtime ceilings as regression guards, not micro-benchmarks,
- keep them loose (roughly `5x` to `10x` observed baseline on CI hardware),
- if CI runner class changes, re-measure and adjust thresholds in
  `tests/performance/` while preserving enough strictness to catch accidental
  scaling blowups.
