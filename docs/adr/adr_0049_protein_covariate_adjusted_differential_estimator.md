# ADR-0049: Protein-Covariate-Adjusted Differential Estimator

## Status

- **ADR ID:** ADR-0049
- **Title:** Protein-Covariate-Adjusted Differential Estimator
- **Status:** Accepted
- **Date:** 2026-08-31
- **Decision Type:** Scientific and Software Contract

## Context

ADR-0025 records protein-aware modelling preparation as a current
preprocessing capability and identifies downstream differential modelling as a
gap to address before broader kinase and activity-method expansion.

The current dataset-owned preparation contract already supports
`DatasetProteinAwarePreparationConfig(policy="prepare_model_inputs")`.
That stage produces typed, non-mutating phosphosite/protein model inputs:
matched phosphosite/protein-row pairs, a sample-aligned protein covariate
matrix, site eligibility, mapping diagnostics, sample-alignment diagnostics,
transformation-state diagnostics, limitations, and preparation provenance when
available. It does not alter phosphosite intensities and does not run
differential modelling.

`DifferentialAnalysisWorkflow` does not currently consume
`AnalysisReadyPhosphoDataset.protein_aware_preparation`. That boundary remains
correct until a later implementation ticket explicitly activates an opt-in
lane. `DifferentialAnalysisRequest` must remain free of caller-supplied protein
matrices, mapping tables, preparation results, or other injected
protein-aware scientific inputs.

The ordinary differential model currently uses one shared design decomposition
for all phosphosite rows in a request. A protein covariate is feature-specific:
different phosphosite rows can map to different total-protein rows. A valid
protein-aware estimator therefore cannot use one global augmented design for
all sites.

This ADR freezes the version-1 estimator contract before runtime implementation
begins. It does not add public exports, runtime support, result fields,
validation behavior, feature documentation, or a support claim beyond this
accepted design record.

## Decision

### Public Selection and Claim

The public opt-in field name is frozen as
`DifferentialAnalysisConfig.protein_aware_model`.

The public advanced configuration type is frozen as
`DifferentialProteinAwareModelConfig`.

The only version-1 method identifier is
`protein_covariate_adjusted_moderated_linear_model_v1`.

The initial scientific claim is `experimental`. Later implementation and
documentation must not describe this method as production validated,
parity-gated, MSstatsPTM-compatible, or a general protein-normalized
differential workflow.

The public workflow request shape remains:

```python
DifferentialAnalysisRequest(
    dataset=dataset,
    design=design,
    contrasts=contrasts,
    config=config,
)
```

The selected lane consumes only the
`ProteinAwarePreparationResult` owned by `request.dataset`. The workflow
request must not accept caller-supplied protein matrices, matched-pair tables,
or preparation sidecars.

### Model Definition

For phosphosite row `s`, mapped to total-protein row `p(s)`, version 1 fits:

```text
y_s = X beta_s + z_p(s) gamma_s + error
```

Definitions:

- `y_s` is the phosphosite abundance vector on an established log2 scale.
- `X` is PhosPy's existing explicit design matrix for condition terms,
  declared fixed covariates, and, when selected, complete `fixed_block` terms.
- `z_p(s)` is the matched total-protein abundance vector on an established
  log2 scale.
- `z_p(s)` is mean-centred over the exact execution samples, after any
  supported deterministic sample subsetting/reordering, and is not
  standardized to unit variance.
- `gamma_s` is a site-specific nuisance protein coefficient.
- Every requested condition contrast is extended with a zero coefficient for
  the nuisance protein coefficient.

The reported `logFC` is the requested condition contrast conditional on the
measured matched total-protein abundance under this model. It is not
phosphosite/protein subtraction, not a stoichiometry or occupancy estimate, and
not a causal decomposition of protein abundance and phosphorylation
regulation.

### Execution Population and Fitting

Version 1 has no within-lane fallback. A site selected into the
protein-aware lane must either be tested by the protein-aware estimator or
retained as a typed withheld row. It must not be silently routed through the
ordinary phosphosite-only estimator.

Fitting is grouped by distinct `total_protein_row_key`. Sites sharing the same
valid protein row share one augmented design decomposition. Sites mapped to
different protein rows require separate augmented design decompositions.
The implementation must preserve authoritative input `site_key` order in
returned tables regardless of grouping order.

Residual variances from all successfully fitted protein-aware sites are pooled
into one existing empirical-Bayes moderation step. The moderation population is
the tested adjusted-site set only.

Multiple testing is applied independently per contrast over successfully tested
protein-aware sites only. Withheld rows and ordinary-lane rows are excluded
from the protein-aware adjusted-p-value population.

The nuisance protein coefficient is reported in protein-aware diagnostics. It
is not a separate inferential result and must not receive its own public
per-site *p*-value or adjusted *p*-value claim in version 1.

### Global Errors and Per-Site Withholding

The selected protein-aware request fails globally before fitting when any of
these conditions holds:

- `dataset.protein_aware_preparation` is absent;
- the sidecar type, schema version, or preparation policy is unsupported;
- the sidecar cannot be proven to belong to the current dataset;
- the current dataset has no total-protein matrix;
- sidecar and dataset sample identities are inconsistent;
- phosphosite values are not on an established log2 scale;
- total-protein values are not on an established log2 scale;
- the dataset has already applied `subtract_log_total` or another incompatible
  protein subtraction;
- actual technical-replicate aggregation is required;
- `paired_design_policy="duplicate_correlation"` is selected;
- no site remains testable after all eligibility and augmented-design checks.

A malformed or cross-attached sidecar is a contract error, not ordinary row
attrition.

For requests that pass global validation, every authoritative `site_key` remains
in the full result index. Per-site status follows this deterministic primary
precedence:

1. existing phosphosite numeric and imputation eligibility;
2. protein-aware preparation and mapping eligibility;
3. protein-covariate numeric and variance eligibility;
4. augmented-design rank, conditioning, residual-degrees-of-freedom, and
   contrast eligibility;
5. defensive post-fit numeric eligibility.

Version 1 must include typed statuses sufficient to distinguish at least:

- `withheld_protein_preparation_ineligible`;
- `withheld_protein_covariate_invalid`;
- `withheld_protein_augmented_design_invalid`;
- `withheld_protein_contrast_non_estimable`.

Specific reason codes must be stable beneath those statuses. Protein-aware
failures must not collapse into `withheld_other`.

### Version-1 Compatibility Matrix

| Capability | Ordinary Lane | Protein-Aware Version 1 |
| --- | --- | --- |
| Fixed condition design | Existing behavior | Supported |
| Declared fixed covariates | Existing behavior | Supported when each augmented design is admissible |
| `fixed_block` | Existing behavior | Supported when each augmented design is admissible |
| `duplicate_correlation` | Existing behavior | Rejected |
| No-op technical-replicate policy | Existing behavior | Supported |
| Actual technical-replicate aggregation | Existing behavior | Rejected |
| `allow_design_subset=True` | Existing behavior | Supported through exact deterministic protein-covariate subsetting/reordering |
| Existing phosphosite imputation policy | Existing behavior | Preserved |
| Protein-covariate imputation | Not applicable | Not performed; non-finite covariates are withheld |
| Prior `subtract_log_total` | Existing ordinary behavior | Rejected |
| Protein-aware fallback to ordinary fitting | Not applicable | Not allowed |

### Diagnostics, Provenance, Fingerprints, and Caveats

The version-1 result must expose typed diagnostics and provenance for:

- selected method and claim status;
- preparation and mapping policies;
- protein reference/source context when genuinely available;
- exact execution sample order;
- mean-centring and no-standardization policy;
- no-protein-imputation and no-fallback policies;
- total site count, ordinary-eligible count, preparation-eligible count, tested
  count, and reason counts;
- base-design rank and base residual degrees of freedom;
- expected augmented rank and common tested augmented residual degrees of
  freedom;
- augmented-design condition-number summaries across protein-row groups;
- per-site mapping, protein-row, covariate, design, fit, and nuisance
  coefficient diagnostics;
- exact, order-sensitive fingerprints of the phosphosite input matrix,
  matched-pair table, resolved protein-covariate matrix, eligibility table,
  design matrix, and contrast matrix;
- all scientific limitations and unsupported claims.

Scalar design diagnostics must be labelled as base-design values or documented
aggregates. One protein group's singular values or condition number must never
be presented as if it represented every augmented design.

### Required Documentation Language

User-facing documentation for this method must use this claim boundary:

```text
protein_covariate_adjusted_moderated_linear_model_v1 is an experimental,
opt-in differential estimator. It reports the requested phosphosite condition
contrast conditional on measured matched total-protein abundance under
y_s = X beta_s + z_p(s) gamma_s + error. It requires established log2
phosphosite and total-protein inputs; centres, but does not standardize,
normalize, subtract, or impute, the total-protein covariate; withholds
inadmissible sites rather than falling back to the ordinary estimator; and is
not stoichiometry, occupancy, causal decomposition, MSstatsPTM parity, or
MSstatsPTM-style joint PTM/protein inference.
```

Feature documentation must not be added until executable support exists.

## External-Method Distinction

MSstatsPTM is a materially different workflow. It models post-translational
modification (PTM) and protein abundance information through separate PTM-site
and global-protein model outputs, then performs adjusted PTM inference using
those outputs (Kohler et al., 2023). That workflow targets PTM changes in the
context of corresponding protein changes and has its own data processing,
summarization, modelling, and adjusted-inference contracts.

The PhosPy estimator defined here is a single phosphosite-row model with a
matched total-protein abundance vector included as a row-specific nuisance
covariate. It must not be described as MSstatsPTM-like, MSstatsPTM parity,
MSstatsPTM-compatible, or a substitute for MSstatsPTM's separate PTM/protein
model-output workflow.

## Consequences

- The reported condition effect is a conditional model coefficient, not
  phosphosite/protein subtraction.
- Treatment-induced protein abundance can be part of the biological pathway.
  Adjustment changes the estimand and is not a causal decomposition of protein
  abundance and phosphorylation regulation.
- Sites without an admissible matched covariate are retained as withheld rows
  rather than silently analyzed by the ordinary lane.
- Multiple protein rows require multiple augmented design decompositions.
- Technical-replicate aggregation and `duplicate_correlation` remain future
  design work for any protein-aware lane.
- Runtime implementation can proceed in later tickets against one frozen field
  name, type name, method ID, eligibility contract, diagnostics contract, and
  documentation claim boundary.

## Rejected Alternatives

1. Implement version 1 as `subtract_log_total`. Rejected because subtraction
   changes phosphosite values before modelling, fixes the protein relationship
   implicitly, and is already a distinct preprocessing policy. The
   protein-aware lane must preserve phosphosite and protein inputs separately.
2. Fix the protein coefficient at 1. Rejected because version 1 estimates a
   site-specific nuisance coefficient `gamma_s`; a fixed coefficient would
   impose a different scientific model.
3. Silently fall back to phosphosite-only fitting for unmatched or inadmissible
   sites. Rejected because fallback mixes estimands within one selected lane and
   hides protein-aware attrition.
4. Add caller-supplied matrices or matched-pair tables to
   `DifferentialAnalysisRequest`. Rejected because the dataset owns the typed
   preparation sidecar and the workflow request must remain an explicit
   dataset/design/contrast/config contract.
5. Fit separate PTM and protein models and combine their contrast estimates.
   Rejected for version 1 because that is a different method family requiring a
   separate design, result, validation, and parity/evidence contract.
6. Claim immediate MSstatsPTM parity. Rejected because PhosPy is defining a
   different estimator and has no external MSstatsPTM parity evidence for this
   lane.
7. Add general mixed-effects, random-intercept, random-slope, or arbitrary
   covariance support. Rejected because version 1 uses the existing fixed-design
   linear-model family plus one row-specific protein covariate.
8. Automatically preprocess total-protein values inside differential analysis.
   Rejected because normalization, transformation, imputation, batch
   correction, and protein mapping are dataset/preprocessing responsibilities,
   not differential workflow responsibilities.

## Implementation Notes

- Later tickets may add configuration contracts, validation, fitting kernels,
  result diagnostics, provenance, tests, and user-facing documentation.
- This ADR itself adds no runtime behavior and must not be used as evidence
  that the method is currently executable.
- Related ADRs: [ADR-0019](adr_0019_experimental_design_and_contrast_contract.md),
  [ADR-0025](adr_0025_competitive_phosphoproteomics_workflow_coverage.md),
  [ADR-0032](adr_0032_differential_declared_scale_eligibility_override.md),
  [ADR-0033](adr_0033_result_caveats.md), and
  [ADR-0048](adr_0048_duplicate_correlation_scientific_contract.md).

## References

- Kohler, D., Tsai, T.-H., Verschueren, E., Huang, T., Hinkle, T., Phu, L.,
  Choi, M., & Vitek, O. (2023). MSstatsPTM: Statistical relative quantification
  of posttranslational modifications in bottom-up mass spectrometry-based
  proteomics. *Molecular & Cellular Proteomics, 22*(1), 100477.
  https://doi.org/10.1016/j.mcpro.2022.100477
