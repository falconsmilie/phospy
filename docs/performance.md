# Performance Contracts

This page defines explicit, CI-tested performance contracts for realistic
phosphoproteomics workloads. These are operational guardrails, not guarantees of
identical runtime on every machine.

## Target Dataset Scale Contract

PhosPy currently targets two practical execution scales for release-gated
scientific workflows:

| Scale | Phosphosites | Samples | Conditions | Missingness in raw phospho input | Reference bundle workload | Kinase scoring workload | Signalome graph/network workload |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Smoke (CI sanity) | ~800 | ~8 | 2 | ~8% | up to ~60 eligible kinases plus off-lane background map | score/predict from ~800 x 60 site-kinase support matrix | module/network outputs from ~150-300 interpreted sites and up to ~40 retained kinases |
| Medium (release contract) | ~3,000 | ~12 | 4 | ~18% | up to ~100 eligible kinases plus large off-lane background map (hundreds of extra kinases) | score/predict from ~2,000 to ~3,000 sites x ~100 candidate kinases | module/network outputs from a few hundred interpreted sites and up to ~100 retained kinases; dense edge tables can reach thousands of rows |

These targets are intentionally conservative for CI stability and are designed
to catch major regressions without changing scientific semantics, validation
requirements, or provenance capture.

## Contract Table

| Operation | Primary input dimensions | Complexity (known/expected) | Recommended input size | Guardrail threshold | Memory-heavy behavior | Fallback / approximation behavior | Failure behavior when limits are exceeded |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Dataset preprocessing (pipeline orchestration) | sites x samples; optional metadata/total rows | Sum of stage costs | 2,000 x 8 to 5,000 x 12 | None (pipeline-level) | Stage-dependent; quantile and correction dominate | No hidden fallback; runs configured stages only | Stage-specific validation or scale errors propagate |
| Site matrix building | input rows; duplicate-site groups; samples | O(rows x samples) plus grouping | <= 5,000 rows, <= 12 samples | None | Duplicate resolution and grouping materialize intermediate tables | Policy-driven duplicate handling (`max_mean_signal`, `first`, `aggregate_*`) | Raises `PhosPyInputError` for unsupported policy, missing required metadata, or zero retained rows |
| Duplicate-site resolution | duplicate groups, rows per group, samples | Groupby/sort path: roughly O(rows log rows) | Moderate duplicate groups; keep group size small | None | Aggregate policies allocate grouped matrices | `aggregate_mean` / `aggregate_median` collapse duplicates; `first` and `max_mean_signal` select one row | `duplicate_site_policy='error'` fails fast on duplicate constructed site IDs |
| Missing data handling | sites x samples; missing fraction | O(rows x columns) | 2,000 x 8 to 5,000 x 12 | None | Imputation copies matrix and row-audit state | Policy-driven (`forbid` vs `impute_row_median`) | Validation errors for policy/parameter mismatches; strict policies reject missing values |
| Quantile normalisation | sites x samples (dense numeric) | ~O(samples x sites log sites) due to per-column sort | 5,000 x 12 (CI contract fixture) | None | Sorting and rank-averaging create additional dense float arrays | None | Numeric/shape validation failures propagate |
| Total protein correction | phospho rows, total rows, samples, identity mapping size | O(matched rows x samples) plus mapping resolution | <= 5,000 rows, <= 12 samples | None | Produces corrected phospho copy and diagnostics hashes | Unmatched-row policy can retain uncorrected rows (`allow_uncorrected`) | Raises `PhosPyInputError` for identity mismatches, missing total rows, unresolved mapping, or invalid scale |
| Differential workflow | sites x samples; design samples; conditions; contrasts | Core fit is roughly O(sites x design columns^2) with per-site moderation/testing | 800 x 8 (2 conditions) to 3,000 x 12 (4 conditions) | Validation contract enforces balanced/estimable design and minimum replicates | Stores per-contrast full output tables (`logFC`, `t`, `P.Value`, `adj.P.Val`) across all sites | No hidden approximations in moderated-statistics path | Raises `WorkflowValidationError` for unsupported/misaligned design, insufficient replicates, missing values, or invalid contrasts |
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

## CI Benchmark Ownership

- Performance thresholds and representative fixture sizes are centralized in
  `tests/support/performance_contracts.py`.
- CI performance tests live in `tests/performance/`.
- Local benchmark scripts live in `benchmarks/` and report plain `key=value` or
  JSONL metrics without affecting production logic.

## Execution and Release Policy

- `tests/performance/` are release-gate checks.
- They are excluded from default local unit/integration pytest runs.
- They are not manual-only checks.
- They should run in dedicated CI/release validation jobs or explicit
  release-validation commands (`make test-release-gate`).
- Failing performance contracts block release until fixed, formally waived, or
  intentionally updated with matching test and documentation changes.

