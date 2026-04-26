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
```

Behavior:

- for `n_sites <= 2000`, module-count scoring uses a full site-by-site correlation matrix,
- for `n_sites > 2000`, module-count scoring uses sampled within-cluster correlation estimates,
- sampled scoring is deterministic and order-invariant in the current implementation,
- approximation use is reported in `SignalomeModuleSelectionDiagnostics.reason`.

Example:

```python
diagnostics = result.module_selection_diagnostics
print(diagnostics.reason)
```

Full-correlation memory floor:

```text
n_sites x n_sites x 8 bytes
```

Example: a `2000 x 2000` float64 correlation matrix is about `32 MiB` before pandas/object overhead and additional intermediate arrays.

Important limit note:

`build_cluster_tree()` remains pairwise/agglomerative and can become expensive as site count grows. Approximate correlation scoring avoids full correlation-matrix materialisation, but does not make clustering free.

## Recommended ranges

| Area | Current recommended range | Behaviour above range |
| --- | ---: | --- |
| Signalome full correlation | up to 2,000 sites | sampled correlation estimates activate |
| Signalome clustering | low thousands of sites | runtime grows quickly because clustering is pairwise |
| Quantile normalisation | thousands to low tens of thousands of sites, tens of samples | dense copies and sorting may become memory-heavy |
| Motif scoring | thousands of sites x hundreds of kinases | cost grows with scored sites, eligible kinases, and window width |
| Kinase-substrate reference maps | large maps are acceptable | overlap filtering should prevent off-lane kinases dominating |

These are performance contracts and operational guidance, not scientific-validity cutoffs.

## CI benchmark policy

Performance contracts are guarded by lightweight synthetic tests in `tests/performance/` and a dedicated CI job (`performance-contracts`).

Policy:

- normal unit tests remain non-timing-sensitive,
- performance tests use fixed synthetic inputs and loose thresholds/ratios,
- tests assert both behavior and bounded runtime/memory for hot paths,
- benchmark scripts under `benchmarks/` remain manually runnable and emit `key=value` metrics for log parsing.
