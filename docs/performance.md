# Performance Contracts

This page defines explicit, CI-tested performance contracts for realistic
phosphoproteomics workloads. These are operational guardrails, not guarantees of
identical runtime on every machine.

## Target Dataset Scale Contract

PhosPy currently targets three practical execution scales for release-checked
scientific workflows:

| Scale | Phosphosites | Samples | Conditions | Missingness in raw phospho input | Reference bundle workload | Kinase scoring workload | Signalome graph/network workload |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Smoke (CI sanity) | ~800 | ~8 | 2 | ~8% | up to ~60 eligible kinases plus off-lane background map | score/predict from ~800 x 60 site-kinase support matrix | module/network outputs from ~150-300 interpreted sites and up to ~40 retained kinases |
| Medium (release contract) | ~3,000 | ~12 | 4 | ~18% | up to ~100 eligible kinases plus large off-lane background map (hundreds of extra kinases) | score/predict from ~2,000 to ~3,000 sites x ~100 candidate kinases | module/network outputs from a few hundred interpreted sites and up to ~100 retained kinases; dense edge tables can reach thousands of rows |
| End-to-end release scale | 50,000 | 48 | 2 | ~3% | no bundled reference lookup; public dataset builder plus differential workflow | not exercised in this contract | not exercised in this contract |

These targets are intentionally conservative for CI stability and are designed
to catch major regressions without changing scientific semantics, validation
requirements, or provenance capture.

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
| End-to-end release-scale builder plus differential | 50,000 sites x 48 samples with realistic site/sample metadata, log2 transform, median centering, row-median missing-data imputation, provenance/fingerprinting, and one two-condition differential contrast | Sum of request preparation, public builder execution, preprocessing/provenance fingerprinting, differential fitting, and result-table export over all retained rows | 50,000 x 48 | Ordinary production runtime < 1,200 seconds; tracemalloc peak < 4 GiB; instrumented memory probe completion timeout 1,800 seconds | Dense input/output matrices, metadata tables, preprocessing reports, provenance fingerprints, and one full differential result table | No hidden approximation; this contract intentionally uses public builder/workflow entrypoints | Fails if shape, provenance completeness, row status, ordinary runtime, memory-probe completion, or tracemalloc peak-memory contract is violated |
| Motif scoring | dataset sites; eligible kinases; sequence window width | Approximately O(sites x eligible kinases) after reference filtering | 2,000 sites x 100 kinases | None | Motif-library and score matrices scale with kinase count | Kinases without valid motif support are naturally excluded | Validation errors for malformed sequence/reference inputs |
| Profile scoring | sites x samples; kinase substrate supports | Dominated by correlation computations; typically O(sites x kinases x samples) | 2,000 to 5,000 sites; 8 to 12 samples | No hard scale guard in scoring stage | Dense downstream score matrices can be large | Profile-only fallback remains available when motif evidence is absent | Boundary errors when no eligible scoring/prediction candidates remain |
| Adaptive prediction | prediction score matrix (sites x kinases); candidate substrates per kinase; ensemble runs | Roughly O(candidate kinases x ensemble runs x sites x kinases) | 2,000 x 100 with fixed seed | No explicit size guard; bounded by config (`adaptive_ensemble_runs`, `n_iterations`) | Repeated train/test allocations per ensemble run | Deterministic seeded sampling (`prediction_config.random_state`) | Raises workflow-stage/boundary errors for empty candidates, missing random seed, or dependency issues |
| Signalome clustering (exact tree path) | interpreted sites x kinases | Exact-tree cost grows superlinearly; practical behavior near O(sites^2) memory/time | <= 2,000 sites by default | `performance.max_exact_tree_sites` (default `2000`) | Exact tree and correlation paths can allocate O(sites^2) structures | None for tree construction; still exact when candidate scoring is sampled | Raises `SignalomeScaleError` when `n_sites > max_exact_tree_sites` |
| Signalome candidate scoring | candidate cluster range; sites x kinases | `full`: O(sites^2); `sampled`: reduced by deterministic per-cluster subsampling | `full` at/below 2,000 sites; sampled above that | `performance.max_full_candidate_scoring_sites` (default `2000`) for `full` policy | Full mode can materialize full site-by-site correlation matrix | `candidate_scoring_policy='sampled'` uses seeded, order-invariant subsampling for candidate-count evaluation only | `SignalomeScaleError` for `full` policy above `max_full_candidate_scoring_sites` (when exact-tree guard allows entry) |
| Protein module derivation | clustered sites; site->protein mapping; proteins | Crosstab/grouping cost roughly O(sites x proteins_nonzero) | Same scale as clustering input | None | Crosstab expands to cluster-by-protein membership matrix | No approximation path | Raises `ValueError` when site->protein mappings are missing for clustered sites |
| Bundle writing | number/size of workflow output tables | O(total cells written) + serialization overhead | Representative signalome outputs from ~200-site workflows | None | File IO across many tables; manifest/config JSON serialization | None | IO/validation errors propagate (invalid path, write failure) |
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

## CI Benchmark Ownership

- Performance thresholds and representative fixture sizes are centralized in
  `tests/support/performance_contracts.py`.
- CI performance tests live in `tests/performance/`.
- The end-to-end release-scale contract lives in
  `tests/performance/test_end_to_end_release_scale_contract.py` and emits
  ordinary production wall-clock, segmented phase timings, separately measured
  tracemalloc peak-memory, final-shape, tested-feature-count, and RSS diagnostics
  through pytest reporting and
  `build/reports/release-scale-performance-contract.json`.
- The release-scale runtime budget is applied only to the uninstrumented
  production execution. The separately identified tracemalloc execution reports
  its runtime for diagnosis and must complete within its explicit CI timeout,
  but that tracing-instrumented runtime is not compared with the production
  runtime threshold.
- The release-scale memory probe runs in a subprocess so the parent test can
  sample process RSS where the platform exposes it (`/proc` on Linux CI and
  Windows process counters locally). RSS is reported for capacity review; the
  required memory gate remains Python-tracked tracemalloc peak memory.
- Local benchmark scripts live in `benchmarks/` and report plain `key=value` or
  JSONL metrics without affecting production logic.
- DataFrame ownership copy behavior is tracked by
  `benchmarks/measure_dataframe_ownership_copy_policy.py`, which reports
  shallow/deep copy counts and owner-mutation leak counts for representative
  borrowed-frame operations.

## End-to-End Release-Scale Baseline Policy

The release-scale gate preserves 50,000 sites, 48 samples, approximately 3%
missingness, public builder construction, log2 transformation, median centering,
row-median missing-data handling, provenance/fingerprinting, and one differential
contrast.

Current enforced budgets:

| Metric | Budget | Applies to |
| --- | --- | --- |
| Production runtime | < 1,200 seconds | Uninstrumented request preparation + builder + differential + result-table export |
| Python-tracked peak memory | < 4,096 MiB | Separate subprocess run with `tracemalloc` enabled |
| Instrumented completion timeout | 1,800 seconds | Tracemalloc subprocess completion only; not a production-runtime budget |
| Process RSS peak | Reported when available | Parent-side subprocess RSS sampling; no release threshold until a stable cross-platform baseline is established |

The release-scale report records these phase timings:

- `dataset_request_preparation_seconds`
- `builder_execution_seconds`
- `preprocessing_execution_seconds`
- `preprocessing_report_assembly_seconds`
- `provenance_fingerprinting_seconds`
- `differential_execution_seconds`
- `serialization_report_assembly_seconds`

Latest local supported-version sanity observation, not a replacement for the
required supported-CI artifact set: on 2026-07-27, Windows Python 3.12.10 passed
the split contract with ordinary production runtime 299.484 seconds,
tracemalloc peak 427.406 MiB, and sampled process RSS peak 1,288.086 MiB.

Budget changes require retained reports from Python 3.10, 3.11, and 3.12 on the
supported CI job, with two consecutive successful CI executions after the change.
The budget must include a documented margin over the slowest supported-version
ordinary production runtime and the highest supported-version tracemalloc peak.
The per-version observations are the retained
`performance-contracts-py3.10`, `performance-contracts-py3.11`, and
`performance-contracts-py3.12` CI artifacts; unsupported-interpreter local
measurements may explain investigation work but must not set release budgets.

## Execution and Release Policy

- `tests/performance/` are release-check confidence checks.
- They are excluded from default local unit/integration pytest runs.
- They are not manual-only checks.
- They should run in dedicated CI/release validation jobs or the explicit
  release command (`make release-check`).
- The release-check selector is
  `pytest tests/performance -m "performance or release_gate"`.
- Performance CI jobs publish pytest duration summaries and retain JUnit reports
  from `build/reports/`, so current runtimes can be reviewed with the budget
  constants in `tests/support/performance_contracts.py`.
- Failing performance contracts block release until fixed, formally waived, or
  intentionally updated with matching test and documentation changes.
