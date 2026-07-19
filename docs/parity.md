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

Parity evidence is also artifact-scoped. A public scientific claim is supported
only for the exact release artifact that passed the release gate with that
fixture evidence. A local parity-only pass, default pytest pass, or copied result
from another commit or distribution is not sufficient release evidence.

## Active Parity Areas

Current active parity coverage includes:

- differential phosphorylation (`tests/parity/test_differential_analysis_parity.py`)
- differential parity envelope contracts (`tests/parity/test_differential_limma_parity.py`)
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
  not a data-cleaning batch-correction step.
- `ruv_readiness` is diagnostic/report-only metadata readiness reporting. It is
  not RUV/SPS/RUV-III correction support.

## Fixture Scope by Lane

| Lane | Main fixture/evidence scope |
| --- | --- |
| Differential | Two-condition unpaired simple contrasts and related limma-envelope checks (`tests/fixtures/rewrite_parity/differential_r_reference/`, `tests/fixtures/rewrite_parity/differential_limma_envelope/`) |
| Kinase scoring/prediction | L6 and public workflow reference lanes (`tests/fixtures/rewrite_parity/r_reference_l6/`, `tests/fixtures/public_workflow_reference/`) |
| Signalome | Public workflow reference and backend parity lanes (`tests/fixtures/public_workflow_reference/`) |
| Activity parity | Activity-stage parity fixtures and threshold-bearing checks in `tests/parity/test_activity_stage_parity.py` |

Run the parity suite with:

```bash
pytest tests/parity -m parity -s
```

Release decisions should run the full release gate (`make test-release-gate`).
Parity failures in that gate are release-blocking, and performance,
reproducibility, golden, reference-manifest, and threshold-bearing parity gates
must all pass before a public scientific release. The same gate also scans
packaged reference manifests and blocks release when a bundled reference is
missing file hashes, fails hash verification, lacks license/organism/namespace
metadata, or declares a non-release-eligible `redistribution_status`.
`unresolved` bundled references block release, and `external_only` references
must not be shipped as bundled data. `approved` requires verified evidence in
the manifest for the exact packaged files; developers and Codex agents must not
use optimistic wording or source-lineage notes as a substitute for that
evidence.

Release-gate source reports are written under:

```text
build/reports/
```

Treat those files as source-suite evidence, not as publication authorization.
The publication audit record is `release-attestation.json`; it is written only
after the policy-required source reports, build manifest, wheel, sdist, and
installed-artifact verification reports agree.

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
| Differential limma parity fixtures | `tests/fixtures/rewrite_parity/differential_r_reference/`, `tests/fixtures/rewrite_parity/differential_limma_envelope/` |
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
- Differential parity comparisons use explicit floating-point tolerances in
  parity tests (`rtol=1e-6`, `atol=1e-8`).
- Missing-value handling is an intentional contract difference:
  `AnalysisReadyPhosphoDataset` requires complete matrices, so missing values
  are rejected before differential execution.

## Open Gaps

Open gaps should be described as open gaps, not as partial equivalence. Common
examples include broader organism-specific bundled references, additional PhosR
workflow surfaces, and any method not protected by fixture-backed comparison.

PhosPy is not currently parity-equivalent with PhosR for SPS/RUV-III
correction. There are no SPS control-selection fixtures, no native RUV-III
correction-kernel parity fixtures, and no PhosR `RUVphospho` corrected-output
parity fixtures. Native SPS/RUV-style preprocessing correction is a validated
PhosPy implementation, not current PhosR parity.
