# Performance Contracts and Local Benchmarks

PhosPy uses two different performance mechanisms:

- Bounded pytest performance contracts under `tests/performance/`, which are
  suitable for routine automation and remain part of `make test-performance`
  and `make release-check`.
- Optional local benchmark scripts under `benchmarks/`, which report
  machine-dependent runtime and capacity observations but do not block releases.

Performance results are operational guardrails and profiling evidence. They are
not guarantees of identical runtime on every machine.

Validated dataset construction installs one private dataset-owned immutable
frame store. Ordinary kinase, differential, and Signalome workflow runs still
receive independent `DatasetInternalView` instances, but those views return
workflow-local wrappers over the dataset-owned snapshots. Shareable NumPy-backed
columns use immutable backing buffers whose writeability cannot be restored;
object, extension, and otherwise unshareable columns are copied per wrapper.
Repeated runs against the same unchanged dataset therefore avoid rebuilding full
internal numeric-table snapshots, while public exports remain detached copies.

## Target Dataset Scales

PhosPy currently targets these practical execution scales:

| Scale | Phosphosites | Samples | Conditions | Missingness in raw phospho input | Reference bundle workload | Kinase scoring workload | Signalome graph/network workload |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Smoke performance contracts | ~800 | ~8 | 2 | ~8% | up to ~60 eligible kinases plus off-lane background map | score/predict from ~800 x 60 site-kinase support matrix | module/network outputs from ~150-300 interpreted sites and up to ~40 retained kinases |
| Medium performance contracts | ~3,000 | ~12 | 4 | ~18% | up to ~100 eligible kinases plus large off-lane background map | score/predict from ~2,000 to ~3,000 sites x ~100 candidate kinases | module/network outputs from a few hundred interpreted sites and up to ~100 retained kinases |
| Optional local release-scale benchmark | 50,000 | 48 | 2 | ~3% | no bundled reference lookup; public dataset builder plus differential workflow | not exercised | not exercised |

The first two scales are designed to catch major regressions in CI without
changing scientific semantics, validation requirements, or provenance capture.
The 50,000 x 48 scale is intentionally opt-in because it is too sensitive to
machine capacity, operating system, Python build, dependency versions, and
runner contention to serve as a release gate.

## Contract Table

| Operation | Primary input dimensions | Complexity (known/expected) | Recommended input size | Guardrail threshold | Memory-heavy behavior | Fallback / approximation behavior | Failure behavior when limits are exceeded |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Dataset preprocessing (pipeline orchestration) | sites x samples; optional metadata/total rows | Sum of stage costs | 2,000 x 8 to 5,000 x 12 | None (pipeline-level) | Stage-dependent; quantile and correction dominate | No hidden fallback; runs configured stages only | Stage-specific validation or scale errors propagate |
| Site matrix building | input rows; duplicate-site groups; samples | O(rows x samples) plus grouping | <= 5,000 rows, <= 12 samples | None | Duplicate resolution and grouping materialize intermediate tables | Policy-driven duplicate handling only after explicit non-error `duplicate_site_policy` selection | Raises `PhosPyInputError` for unsupported policy, missing required metadata, or zero retained rows |
| Duplicate-site resolution | duplicate groups, rows per group, samples | Groupby/sort path: roughly O(rows log rows) | Moderate duplicate groups; keep group size small | None | Aggregate policies allocate grouped matrices | `aggregate_mean` / `aggregate_median` deliberately collapse duplicate `site_key` evidence; `first` and `max_mean_signal` deliberately select one source row | Default `duplicate_site_policy='error'` fails fast on duplicate `site_key` rows |
| Missing data handling | sites x samples; missing fraction | Row-median/MinProb: O(rows x columns). KNN: O(rows_with_missing x retained_rows x columns) distance work | Row-median/MinProb: 2,000 x 8 to 5,000 x 12. KNN: sparse 12-sample and moderate 24-sample retained missing-target workloads up to 50,000 retained sites under the KNN guardrails | KNN retained rows <= 50,000, samples <= 64, distance-feature operations <= 2,000,000,000 | Imputation copies matrix and row-audit state. KNN additionally chunks pairwise target-row distance matrices | Policy-driven (`forbid`, `impute_row_median`, `impute_minprob`, `impute_knn`); no hidden fallback between policies | Validation errors for policy/parameter mismatches; strict policies reject missing values; KNN rejects impractical retained shapes with actionable `PhosPyInputError` |
| Quantile normalisation | sites x samples (dense numeric) | ~O(samples x sites log sites) due to per-column sort | 5,000 x 12 (CI contract fixture) | None | Sorting and rank-averaging create additional dense float arrays | None | Numeric/shape validation failures propagate |
| Total protein correction | phospho rows, total rows, samples, identity mapping size | O(matched rows x samples) plus mapping resolution | <= 5,000 rows, <= 12 samples | None | Produces corrected phospho copy and diagnostics hashes | Unmatched-row policy can retain uncorrected rows (`allow_uncorrected`) | Raises `PhosPyInputError` for identity mismatches, missing total rows, unresolved mapping, or invalid scale |
| Differential workflow | sites x samples; design samples; conditions; contrasts | Core fit is roughly O(sites x design columns^2) with per-site moderation/testing | 800 x 8 (2 conditions) to 3,000 x 12 (4 conditions) | Validation contract enforces balanced/estimable design and minimum replicates | Stores per-contrast full output tables (`logFC`, `t`, `P.Value`, `adj.P.Val`) across all sites | No hidden approximations in moderated-statistics path | Raises `WorkflowValidationError` for unsupported/misaligned design, insufficient replicates, missing values, or invalid contrasts |
| Repeated workflow internal dataset reads | unchanged validated dataset frames across repeated workflow runs | First internal access per frame pays one owner-detached snapshot copy; later views/runs use immutable-buffer wrappers for shareable NumPy-backed columns and per-wrapper copies for unshareable columns | Same dataset reused across differential/kinase/Signalome runs | Unit instrumentation asserts one phospho/site-metadata snapshot across repeated differential and kinase runs | Retains one private immutable copy per accessed dataset frame; optional frames are snapshotted only when internally read; unshareable columns allocate per returned wrapper | No global cache and no workflow-instance mutable frame cache | Public exports still copy; workflow computation may still make explicit per-run derived/result tables when scientifically required |
| Optional release-scale builder plus differential benchmark | 50,000 sites x 48 samples with realistic site/sample metadata, required `site_sequence`, log2 transform, median centering, row-median imputation, preprocessing provenance, table fingerprints, and one two-condition differential contrast | Sum of request preparation, public builder execution, preprocessing/provenance fingerprinting, differential fitting, result-table export, and benchmark summary generation | 50,000 x 48 | None; runtime and RSS are informational local observations | Dense input/output matrices, metadata tables, preprocessing reports, provenance fingerprints, and one full 50,000-row differential result table | No hidden approximation; this benchmark uses public builder/workflow entrypoints | Fails only when the scientific workflow errors, required invariants are violated, output dimensions are wrong, expected preprocessing/differential outputs are absent, or a requested report cannot be produced |
| ssGSEA substrate enrichment activity | finite ranked sites x kinases x profiles x seeded permutations | Observed score pass is O(sites x kinases x profiles); permutation work is O(kinases x profiles x permutations x selected substrates), with reusable null-score constants cached by equivalent background size, substrate count, and tie-block structure | 720 sites x 32 kinases x 6 profiles x 48 permutations (CI contract fixture) | None | Does not materialize a sites x kinases x profiles x permutations cube; p/q matrices scale with kinases x profiles | No approximation; seeded per-kinase/profile permutation streams and tie-block midrank semantics are preserved | Validation/status diagnostics report insufficient substrates, empty finite background, or all-substrate backgrounds |
| Motif scoring | dataset sites; eligible kinases; sequence window width | Approximately O(sites x eligible kinases) after reference filtering | 2,000 sites x 100 kinases | None | Motif-library and score matrices scale with kinase count | Kinases without valid motif support are naturally excluded | Validation errors for malformed sequence/reference inputs |
| Profile scoring | sites x samples; kinase substrate supports | Dominated by correlation computations; typically O(sites x kinases x samples) | 2,000 to 5,000 sites; 8 to 12 samples | No hard scale guard in scoring stage | Dense downstream score matrices can be large | Profile-only fallback remains available when motif evidence is absent | Boundary errors when no eligible scoring/prediction candidates remain |
| Adaptive prediction | prediction score matrix (sites x kinases); candidate substrates per kinase; ensemble runs | Roughly O(candidate kinases x ensemble runs x sites x kinases) | 2,000 x 100 with fixed seed | No explicit size guard; bounded by config (`adaptive_ensemble_runs`, `n_iterations`) | Repeated train/test allocations per ensemble run | Deterministic seeded sampling (`prediction_config.random_state`) | Raises workflow-stage/boundary errors for empty candidates, missing random seed, or dependency issues |
| Signalome clustering (exact tree path) | interpreted sites x kinases | Exact-tree cost grows superlinearly; practical behavior near O(sites^2) memory/time | <= 2,000 sites by default | `performance.max_exact_tree_sites` (default `2000`) | Exact tree and correlation paths can allocate O(sites^2) structures | None for tree construction; still exact when candidate scoring is sampled | Raises `SignalomeScaleError` when `n_sites > max_exact_tree_sites` |
| Signalome candidate scoring | candidate cluster range; sites x kinases | `full`: O(sites^2); `sampled`: reduced by deterministic per-cluster subsampling | `full` at/below 2,000 sites; sampled above that | `performance.max_full_candidate_scoring_sites` (default `2000`) for `full` policy | Full mode can materialize full site-by-site correlation matrix | `candidate_scoring_policy='sampled'` uses seeded, order-invariant subsampling for candidate-count evaluation only | `SignalomeScaleError` for `full` policy above `max_full_candidate_scoring_sites` |
| Protein module derivation | clustered sites; site->protein mapping; proteins | Crosstab/grouping cost roughly O(sites x proteins_nonzero) | Same scale as clustering input | None | Crosstab expands to cluster-by-protein membership matrix | No approximation path | Raises `ValueError` when site->protein mappings are missing for clustered sites |
| Bundle writing | number/size of workflow output tables | O(total cells written) + serialization overhead | Representative signalome outputs from ~200-site workflows | None | File IO across many tables; manifest/config JSON serialization | None | IO/validation errors propagate |
| Provenance hashing | table rows x columns; dtype/structure complexity | O(total cells) for value normalization + hashing | 5,000 x 12 numeric plus large metadata tables | None | Converts table cells to object normalization payloads before digest updates | None | Hashing/serialization errors propagate for unsupported payload values |

## Signalome Guardrails

Default signalome guardrails are configured through `SignalomeConfig.performance`:

- `max_exact_tree_sites=2000`
- `max_full_candidate_scoring_sites=2000`

Guardrail enforcement is owned by
`src/phospy/science/signalomes/clustering/scale_guards.py`.

## KNN Missing-Data Guardrails

`missing_data.policy="impute_knn"` is a custom deterministic preprocessing
implementation in
`src/phospy/science/datasets/preprocessing/stages/missing_data/knn.py`.

Current KNN execution budgets:

- retained site rows: `<= 50,000`
- sample columns: `<= 64`
- estimated distance-feature operations:
  `rows_with_missing_values x retained_rows x sample_columns <= 2,000,000,000`
- pairwise distance chunk target: `48 MiB` per target-by-donor matrix
- release-check peak memory budget: `< 384 MiB`

Release-check KNN benchmark fixtures measure both sparse and moderate retained
missing-target workloads. They intentionally do not claim broad random
missingness across all retained rows.

| Tier | Sites | Samples | Rows With Missing Values | Missing Cells Per Target Row | Runtime Budget | Peak-Memory Budget |
| --- | --- | --- | --- | --- | --- | --- |
| Sparse | 10,000 | 12 | 96 | 1 | < 5 seconds | < 384 MiB |
| Sparse | 25,000 | 12 | 96 | 1 | < 8 seconds | < 384 MiB |
| Sparse | 50,000 | 12 | 96 | 1 | < 12 seconds | < 384 MiB |
| Moderate | 10,000 | 24 | 256 | 2 | < 8 seconds | < 384 MiB |
| Moderate | 25,000 | 24 | 512 | 2 | < 15 seconds | < 384 MiB |
| Moderate | 50,000 | 24 | 768 | 2 | < 30 seconds | < 384 MiB |

These benchmarks validate deterministic donor semantics, chunk-equivalent
output, and peak-memory behavior for supported retained target densities.
Broader random missingness at 25,000-50,000 retained sites can exceed the
distance-work budget and is intentionally rejected. In that case, reduce
retained missing rows, lower `missing_data.max_missing_fraction_per_row`,
pre-filter low-value features, or use `missing_data.policy="impute_row_median"`
when its scientific semantics are acceptable.

## CI Performance Contract Ownership

- Performance thresholds and representative fixture sizes are centralized in
  `tests/support/performance_contracts.py`.
- CI performance tests live in `tests/performance/`.
- The bounded ssGSEA activity contract covers ranked sites, kinases, profiles,
  seeded permutations, and allocation sanity without using the optional
  release-scale workload.
- `make test-performance` selects
  `pytest tests/performance -m "performance or release_gate"`.
- `make release-check` includes `make test-performance`.
- Performance CI jobs retain JUnit/duration reports from `build/reports/`, so
  current runtimes can be reviewed with the budget constants in
  `tests/support/performance_contracts.py`.
- Failing bounded performance contracts block release until fixed, waived, or
  intentionally updated with matching test and documentation changes.

The 50,000 x 48 end-to-end workload is not owned by pytest, CI, or release
checks. It is owned by the benchmark script described below.

## Optional Repeated Workflow Snapshot Reuse Benchmark

The explicit local command is:

```bash
python benchmarks/measure_repeated_workflow_dataset_snapshot_reuse.py
```

The benchmark builds one unchanged analysis-ready dataset, runs two
differential workflow executions and two kinase workflow executions against
that same dataset, and reports:

- dataset dimensions
- frame dtypes
- first and repeated workflow runtime
- tracemalloc peak memory
- full-frame deep-copy counts for dataset-owned frames
- dataset internal snapshot construction counts
- Python/platform metadata
- selected dependency versions

It can also write a JSON report under `benchmarks/reports/`:

```bash
python benchmarks/measure_repeated_workflow_dataset_snapshot_reuse.py --write-report
```

Observed local run on 2026-08-08, Windows 11
(`Windows-11-10.0.26200-SP0`), Python 3.13.12, NumPy 2.4.4, pandas 2.3.3,
SciPy 1.17.1, and PhosPy 1.5.2:

| Metric | Observed value |
| --- | ---: |
| Dataset phospho dimensions | 240 x 8 |
| Dataset site-metadata columns | 10 |
| Differential first run | 0.438641 s |
| Differential repeated run | 0.327899 s |
| Differential first-run tracemalloc peak | 0.670 MiB |
| Differential repeated-run tracemalloc peak | 0.556 MiB |
| Kinase first run | 1.486004 s |
| Kinase repeated run | 1.479663 s |
| Kinase first-run tracemalloc peak | 1.263 MiB |
| Kinase repeated-run tracemalloc peak | 0.942 MiB |
| Full-frame deep-copy counts | `{"dataset.phospho": 1, "dataset.site_metadata": 5}` |
| Snapshot construction counts | `{"dataset.phospho internal snapshot": 1, "dataset.site_metadata internal snapshot": 1}` |

The copy-count observation is the relevant ownership signal: the unchanged
numeric phospho table is deep-copied once for the dataset-owned snapshot across
all four workflow executions, and the dataset phospho/site-metadata snapshots
are each constructed once. The runtime and memory values are local observations
for this machine and dependency set, not portable performance guarantees.

## Optional Release-Scale Local Benchmark

The explicit local command is:

```bash
make benchmark-release-scale
```

It invokes:

```bash
python benchmarks/measure_release_scale_builder_differential.py
```

The benchmark preserves the representative release-scale workload:

- 50,000 phosphosites
- 48 samples
- deterministic approximately 3% missingness
- realistic phosphosite metadata with required `site_sequence`
- realistic sample metadata
- public `AnalysisReadyDatasetBuilder`
- log2 transformation
- median centering
- row-median imputation
- preprocessing provenance
- table fingerprint generation
- one two-condition differential contrast
- all 50,000 differential result rows

The normal invocation runs the workload once. It does not launch a second
tracemalloc subprocess. It reports `key=value` metrics including total runtime,
request preparation, builder execution, preprocessing, preprocessing report
assembly, provenance/fingerprinting, differential analysis, result-table
assembly, result-table fingerprinting, output dimensions, tested feature count,
original/final missing-cell counts, scientific-summary digest, and process RSS
peak when the platform exposes a truthful current-process peak. If process RSS
is unavailable, it reports `process_rss_peak_mib=unavailable`.

To also write a JSON report:

```bash
python benchmarks/measure_release_scale_builder_differential.py --write-report
```

Reports are written only under `benchmarks/reports/`. The report records the
Python executable/version, selected dependency versions (`phospy`, `numpy`,
`pandas`, and `scipy` when installed metadata is available), machine/platform
metadata, runtime timings, peak RSS when the platform exposes it, metrics,
scientific-summary payload, and output fingerprints/digests for the benchmarked
tables and provenance trace.

Benchmark observations are informational. A slow result is evidence for local
profiling or same-machine comparison; it is not a release failure and is not
evidence that all supported machines meet a fixed runtime or memory envelope.
The benchmark fails only if the scientific workflow raises an error, required
output invariants are violated, expected preprocessing or differential outputs
are absent, output dimensions are wrong, or a requested report cannot be
produced.

## Execution and Release Policy

- `tests/performance/` are release-check confidence checks for bounded
  performance contracts.
- They are excluded from default local unit/integration pytest runs.
- They are selected by `make test-performance` and included in
  `make release-check`.
- The optional 50,000 x 48 benchmark is excluded from `make test-performance`,
  `make release-check`, GitHub Actions workflows, release workflows, and
  publication targets.
- Do not move the optional benchmark into nightly, scheduled, manual,
  tag-only, release, or non-blocking GitHub Actions jobs.
