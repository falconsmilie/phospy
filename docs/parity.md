# Parity to PhosR

PhosPy parity is intentionally narrow and fixture-backed. Passing a parity test
for one lane does not mean the whole PhosR package is implemented.

PhosPy does not claim global PhosR parity. Use
`docs/scientific-coverage.md` as the single maintained scope matrix for support
status labels (`parity-gated`, `validated PhosPy implementation`,
`experimental`, `open gap`, `deliberate scope difference`, `not planned`).

Scope ownership split:

- [Scientific Coverage](scientific-coverage.md) owns user-facing coverage status,
  intended parity scope, and interpretation limits.
- [Parity](parity.md) owns fixture-level comparison evidence, test locations,
  and comparison mechanics.

## What Parity Means Here

A parity claim must say:

- which input fixture was used
- which output table or metric was compared
- which tolerance or acceptance rule was used
- which PhosPy workflow or stage produced the output

Parity evidence here should be interpreted only for the exact fixture + output +
comparison rule documented by each test lane.

Evidence categories are intentionally separated:

- External parity requires an independently generated external output whose
  semantics match the compared PhosPy output.
- PhosPy-generated expected values are regression fixtures, not external parity.
- Closed-form planted fixtures are synthetic validation evidence, not empirical
  validation and not external parity.
- Real-world empirical validation must name the independent empirical dataset or
  study and the comparison result; no current PhosPy-owned fixture supplies that
  evidence category.

Parity evidence is also release-scoped. A public scientific claim should be made
only after the maintainer release checks pass for the tagged source and freshly
built artifacts. A local parity-only pass, default pytest pass, or copied result
from another commit or distribution is not sufficient release confidence.

## Active Parity Areas

Current active parity coverage includes:

- differential phosphorylation (`tests/parity/test_differential_analysis_parity.py`)
- differential parity envelope contracts (`tests/parity/test_differential_limma_parity.py`)
- large-feature differential trend parity against R/limma
  (`tests/parity/test_differential_limma_trend_large.py`)
- kinase scoring and prediction surfaces:
  `tests/parity/test_kinase_workflow_parity.py`,
  `tests/parity/test_prediction_science_parity.py`,
  `tests/parity/test_l6_prediction_parity.py`, and
  `tests/parity/test_public_predmat_parity.py`
- selected preprocessing behaviours with explicit fixtures
- activity-stage behaviours in `tests/parity/test_activity_stage_parity.py`
- signalome workflow and clustering backend fixture lanes:
  `tests/parity/test_signalome_workflow_parity.py` and
  `tests/parity/test_signalome_clustering_backend_parity.py`

## Non-Parity Support and Scope Differences

Some implemented PhosPy features are not PhosR parity claims. They are
validated PhosPy implementations:

- `EnrichmentWorkflow` is implemented as offline ORA over caller-supplied
  gene-set or PTM-set collections. It is not a PhosR enrichment parity lane and
  does not imply GSEA, ssGSEA, PTM-SEA, Enrichr, gseapy, or clusterProfiler
  support.
- `MappedPhosphositeTableImporter`, `MaxQuantPhosphositeImporter`, and
  `FragPipePTMProphetImporter` are input-preparation adapters that feed
  dataset-builder validation. They are not PhosR parity lanes, broad vendor
  parity, Spectronaut/DIA-NN support, or upstream statistical result import.
- `linear_residualize_batch` is limited fixed-effect residualisation under the
  dataset preprocessing `batch_correction` config group. It is not native
  SPS/RUV-style correction, not PhosR-equivalent batch correction, not ComBat,
  not limma `removeBatchEffect` parity, not `duplicateCorrelation`, and not
  mixed-effects modelling.
- `SpsRuvBatchCorrectionConfig` exposes native SPS/RUV-style preprocessing
  correction with explicit controls, protected design metadata, missingness
  policy, factor count, diagnostics, and provenance. It estimates unwanted
  factors from eligible control-site residuals after protected-design handling.
  Batch terms are resolved for validation and diagnostics, not directly
  residualized as fixed effects by the native correction. It is not
  PhosR-equivalent SPS/RUV-III parity.
- Differential fixed-effect batch covariates are ordinary model terms. They are
  not a data-cleaning batch-correction step. Differential
  `duplicate_correlation` fixtures are a separate limma-envelope lane covering
  one REML-estimated consensus compound-symmetry correlation and GLS, not a
  general mixed-effects claim.
- Release-validation fixtures under
  `tests/fixtures/release_validation_regression/` are PhosPy regression
  contracts. They cover adverse differential designs, peptide evidence
  resolution, sparse kinase support, and signalome safety behavior without
  claiming external parity.
- Synthetic validation fixtures under
  `tests/fixtures/release_validation_regression/sps_ruv_planted_unwanted_factor/`,
  `tests/fixtures/release_validation_regression/peptide_site_bias_regimes/`,
  `tests/fixtures/release_validation_regression/kinase_activity_known_membership/`,
  and
  `tests/fixtures/release_validation_regression/signalome_planted_modules/`
  are closed-form known-truth checks. They are not PhosR parity and not
  empirical validation.
- `ruv_readiness` is diagnostic/report-only metadata readiness reporting. It is
  not RUV/SPS/RUV-III correction support.

## Fixture Scope by Lane

| Lane | Main fixture/evidence scope |
| --- | --- |
| Differential | Two-condition unpaired simple contrasts and related limma-envelope checks (`tests/fixtures/rewrite_parity/differential_r_reference/`, `tests/fixtures/rewrite_parity/differential_limma_envelope/`), duplicate-correlation paired-design fixtures (`tests/fixtures/rewrite_parity/differential_duplicate_correlation/`), plus the large-feature trend fixture (`tests/fixtures/rewrite_parity/differential_limma_trend_large/`) |
| Release-validation regression | PhosPy-owned regression fixtures for evidence resolution, sparse kinase support, and signalome safety (`tests/fixtures/release_validation_regression/`) |
| Kinase scoring/prediction | L6 and public workflow reference lanes (`tests/fixtures/rewrite_parity/r_reference_l6/`, `tests/fixtures/public_workflow_reference/`) plus sparse-support regression fixtures under `tests/fixtures/release_validation_regression/kinase_sparse_support/` |
| Signalome | Public workflow reference and backend parity lanes (`tests/fixtures/public_workflow_reference/`) plus safety regression fixtures under `tests/fixtures/release_validation_regression/signalome_safety/` |
| Activity parity | Activity-stage parity fixtures and threshold-bearing checks in `tests/parity/test_activity_stage_parity.py` |
| Synthetic scientific validation | Planted unwanted-factor/protected-signal, peptide-to-site bias regimes, known-membership kinase activity, and planted signalome modules under `tests/fixtures/release_validation_regression/` |
| Importer edge-case regression | Targeted MaxQuant and FragPipe/PTMProphet fixture index under `tests/fixtures/release_validation_regression/importer_edge_cases/`, referencing checked-in importer fixture bytes |

Run the blocking parity suite with:

```bash
pytest tests/parity -m "parity and not parity_diagnostic" -s
```

Release decisions should run the maintainer release checks (`make release-check`).
Parity failures in that check are release-blocking, and
performance, release/golden/reproducibility, checked-in reference,
packaged-reference, installed wheel/sdist verification, metadata, and
threshold-bearing parity checks must all pass before a public scientific
release. Packaged reference validation blocks release when a bundled reference
is missing file hashes, fails hash verification, lacks license/organism/namespace
metadata, or declares a non-release-eligible `redistribution_status`. Installed
verification separately installs both built artifacts outside the checkout,
checks installed import origins, verifies bundled resources, and exercises
representative public workflow contracts.
`unresolved` bundled references block release, and `external_only` references
must not be shipped as bundled data. `approved` requires verified evidence in
the manifest for the exact packaged files; developers and Codex agents must not
use optimistic wording or source-lineage notes as a substitute for that
evidence.

This parity guidance is part of normal CI/build confidence. Final release
verification, including Git-backed staged-byte checks and built wheel/sdist
verification, is described in [Maintenance](maintenance.md).

Some diagnostic parity tests are informational. Release decisions should use the
threshold-bearing gates and the documented fixture expectations, not visual
inspection alone. Tests marked `parity_diagnostic` remain non-blocking unless
they are intentionally promoted into the release selector.

Release-gated parity command in the Makefile is:

```bash
pytest tests/parity -m "parity and not parity_diagnostic" -s
```

The broader CI parity smoke command is:

```bash
pytest tests/parity -m "parity and not parity_diagnostic" -s
```

Informational parity diagnostics run separately in CI with:

```bash
pytest tests/parity -m "parity_diagnostic" -s
```

## Fixture Locations

| Purpose | Location |
| --- | --- |
| Parity tests | `tests/parity/` |
| Shared parity helpers | `tests/support/` |
| Public workflow reference fixtures | `tests/fixtures/public_workflow_reference/` |
| Differential limma parity fixtures | `tests/fixtures/rewrite_parity/differential_r_reference/`, `tests/fixtures/rewrite_parity/differential_limma_envelope/`, `tests/fixtures/rewrite_parity/differential_limma_trend_large/` |
| PhosPy release-validation regression fixtures | `tests/fixtures/release_validation_regression/` |
| Regeneration scripts | `scripts/active/` |

## Differential Parity Envelope Notes

- Current differential analysis is not full PhosR or limma parity unless a
  specific behavior is explicitly documented by fixture-backed parity tests.
- Supported designs are limited to tested design and contrast envelopes; outside
  those envelopes, results should be treated as unsupported unless public docs
  and tests say otherwise.
- Differential parity claims are feature-scoped. Current limma-backed fixtures
  protect:
  - two-condition unpaired simple contrasts (`B_vs_A`, `A_vs_B`)
  - small-`n` moderated-statistics behavior
  - unequal-variance feature handling
  - the >1,024-feature empirical-Bayes trend branch through a 1,600-feature
    R/limma fixture with unbalanced 5/7 condition groups
  - Benjamini-Hochberg adjusted p-values for fixtures whose rows are all tested
    and contrast ordering/sign conventions
- PhosPy withholds all-constant feature rows before differential model fitting.
  Those rows are reported with `result_status="withheld_all_constant"` and
  missing `logFC`, `t`, `P.Value`, and `adj.P.Val` values.
- Withheld rows are not included in multiple-testing correction. Adjustment is
  computed per contrast over tested rows only, then withheld rows are expanded
  back into the result table with missing statistics.
- Differential parity comparisons that use fixtures containing withheld rows are
  performed only on rows with `result_status="tested"`. Withheld all-constant
  rows document a PhosPy safety policy, not a limma parity surface.
- This tested-row comparison scope does not claim that PhosPy reproduces all
  limma edge-case behavior for constant or otherwise untestable features.
- Small differential parity fixtures use explicit floating-point tolerances in
  parity tests (`rtol=1e-6`, `atol=1e-8`).
- The large trend fixture compares exact condition coefficients at
  `rtol=1e-10`, `atol=1e-10` because the fitted condition coefficients are
  scientifically equivalent ordinary least-squares quantities. Moderated
  quantities use drift-envelope checks because limma and PhosPy use different
  trend smoothers: moderated-t correlation > 0.995, negative-log10-p
  correlation > 0.98, standard-error correlation > 0.94,
  log-prior-variance correlation > 0.90, median p-value absolute difference
  < 0.015, p99 p-value absolute difference < 0.08, p99 adjusted-p-value
  absolute difference < 0.15, and p99 standard-error absolute difference
  < 0.035.
- Missing-value handling is an intentional contract difference:
  `AnalysisReadyPhosphoDataset` requires complete matrices, so missing values
  are rejected before differential execution.
- Upstream-imputed datasets are rejected by default. The explicit
  `withhold_imputed_features` policy is a PhosPy safety contract, not a limma
  parity surface.

### Large Differential Trend Fixture Provenance

The checked-in large trend fixture was generated by:

```bash
Rscript scripts/active/generate_large_differential_limma_trend_fixture.R --outdir tests/fixtures/rewrite_parity/differential_limma_trend_large --seed 20260724 --timestamp 2026-07-24T00:00:00Z --n_features 1600
```

Its manifest records R version `4.5.2`, limma version `3.66.0`, seed
`20260724`, output file SHA-256 hashes, generator SHA-256, the explicit
`~0 + condition` design, and contrast `B_vs_A = B - A`. Only the exported
limma scientific result columns are external-reference comparison targets;
simulation diagnostics are sanity metadata.

The PhosPy-owned regression fixtures were generated by:

```bash
python scripts/active/generate_release_validation_regression_fixtures.py --outdir tests/fixtures/release_validation_regression --timestamp 2026-07-24T00:00:00Z --seed 20260724
```

Those fixtures are classified as `regression` in their manifests and are not
external parity evidence.

The additional known-truth fixture families generated by the same command are
classified as `synthetic_validation` in their manifests. They validate planted
factor recovery, protected-signal retention, peptide-to-site adverse-regime
bias, kinase activity direction/coverage sensitivity, and signalome planted
module recovery. They do not broaden PhosR equivalence claims and must not be
used as parameter-tuning targets.

The manifest-governed fixture families use project-standard text bytes: UTF-8, LF
line endings, and a final newline. Manifest hashes are raw byte hashes of the
checked-in files; parity and regression tests must not repair or normalize
line endings before hashing. Use
`tests/release/test_manifest_fixture_byte_reproducibility.py` to verify that
regeneration into a temporary directory reproduces the checked-in bytes.

## Open Gaps

Open gaps should be described as open gaps, not as partial equivalence. Common
examples include broader organism-specific bundled references, additional PhosR
workflow surfaces, and any method not protected by fixture-backed comparison.

PhosPy is not currently parity-equivalent with PhosR for SPS/RUV-III
correction. There are no SPS control-selection fixtures, no native RUV-III
correction-kernel parity fixtures, and no PhosR `RUVphospho` corrected-output
parity fixtures. Native SPS/RUV-style preprocessing correction is a validated
PhosPy implementation, not current PhosR parity.
