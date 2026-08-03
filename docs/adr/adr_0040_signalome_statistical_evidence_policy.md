# ADR-0040: Signalome Statistical Evidence Policy

## Status

- **ADR ID:** ADR-0040
- **Title:** Signalome Statistical Evidence Policy
- **Status:** Accepted
- **Date:** 2026-07-23
- **Decision Type:** Scientific Safety and Workflow Contract

## Context

Signalome network output is built from pairwise Pearson correlations between
kinase downstream score profiles. A two-point Pearson correlation is always
`+1` or `-1` when both profiles are nonconstant, so accepting two-observation
pairs as ordinary network edges overstates the evidential strength of small
sample profiles.

Signalome clustering also previously filled fully unobserved kinase/dimension
columns with zero. That made absence of observation look like an observed
numerical value and could affect clustering distance calculations without a
clear audit trail.

## Decision

New signalome network execution requires at least three paired finite
observations for any candidate edge, and the built-in default is five paired
finite observations. Public config rejects explicit values below three. The
production mode and `SignalomeConfig.production()` require the effective
minimum to be at least five. Network science code enforces the effective
minimum and classifies pairs below it with the existing
insufficient-observations status instead of accepting an edge.

Production mode is the default and recommended Signalome entry point. It also
requires present site-level localisation evidence with the production
probability threshold. Historical exploratory behavior is available only
through the explicitly named `exploratory_compatibility` mode via
`SignalomeConfig.compatibility()`, where network paired-observation settings
may use the public floor of three.

Signalome grouping identity is named `protein_group_id`. It belongs to the
Signalome workflow/domain because it controls module and protein-site context
grouping. The legacy `protein_id` column is accepted only as a migration alias;
payloads or datasets containing both names with conflicting values fail
validation. Core protein identity remains under the site-identity domain as
`organism`, `protein_namespace`, and `protein_identifier`; Signalome must not
reinterpret `protein_identifier` as grouping identity.

Every accepted public network edge includes `valid_observations`. Network
provenance records the requested threshold, effective threshold, and stable
policy identifier:
`signalome_network_min_paired_finite_observations_floor3_default5_v1`.
Workflow provenance also records `signalome_mode` and Signalome grouping
identity metadata, including whether the legacy alias was used.

Signalome network correlations are descriptive score-profile associations.
They are not inferential evidence, causal evidence, directional evidence, or
experimental validation. This ADR does not add p-values, confidence intervals,
or multiple-testing correction. Any future inferential network statistics must
be governed by a separate statistical ADR.

Signalome clustering preparation drops fully missing kinase/dimension columns
before clustering and median-imputes only partially missing values in retained
columns. The stable preparation policy identifier is
`drop_fully_missing_then_column_median_impute`. Preparation diagnostics and
provenance record retained labels, dropped fully missing labels, dropped-cell
counts, per-column imputation counts, total imputed retained cells, and the
exact prepared-matrix fingerprint. If no dimension remains after dropping fully
missing columns, execution fails before tree construction.

## Historical Payloads

Historical signalome result bundles remain readable. If a stored result records
an effective network threshold of two, bundle reconstruction preserves that
historical threshold instead of rewriting its meaning.

Historical configuration is not permission for new execution. Replay or
re-execution using a legacy threshold of two must fail validation with an
actionable migration message requiring
`config.output.network_min_paired_finite_observations >= 3`.

New snapshots must preserve requested and effective threshold values when the
schema distinguishes them. Production defaults record mode `production`,
requested `5`, and effective `5`. Legacy snapshots that omit mode or
localisation policy are loaded as explicit exploratory compatibility payloads
rather than being silently upgraded to production semantics.

## Consequences

Positive consequences:

- Two-observation Pearson correlations cannot become ordinary signalome
  network edges in new executions.
- Public edge tables expose the observation count needed to interpret an edge.
- Production and exploratory compatibility behavior are visible in typed config
  and provenance instead of being hidden behind permissive defaults.
- Signalome grouping identity cannot be confused with the core
  `protein_identifier` identity field.
- Missing clustering dimensions are not silently converted into numerical zero.
- Clustering preparation is reproducible through an exact prepared-matrix
  fingerprint and deterministic diagnostics.
- Historical results remain interpretable under the policy recorded when they
  were produced.

Negative consequences:

- Small signalome executions with fewer than the effective minimum paired
  finite observations may now produce no accepted network edges and fail at the
  existing network boundary.
- Golden/parity fixtures that include network edges or signalome provenance
  require fixture updates.

## Related Records

- [ADR-0025: Competitive Phosphoproteomics Workflow Coverage Roadmap](adr_0025_competitive_phosphoproteomics_workflow_coverage.md)
- [ADR-0033: Result Caveats and Scientific Warning Ownership](adr_0033_result_caveats.md)
- [ADR-0034: Quantitative State, Motif Scoring Semantics, and Reference Context](adr_0034_quantitative_state_motif_semantics_and_reference_context.md)
