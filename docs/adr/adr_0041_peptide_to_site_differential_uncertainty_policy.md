# ADR: Peptide-to-Site Differential Uncertainty Policy

## Document Control

- **ADR ID:** ADR-0041
- **Title:** Peptide-to-Site Differential Uncertainty Policy
- **Status:** Accepted
- **Date:** 2026-07-29
- **Decision Type:** Architecture Decision Record

## Context

Peptide-level differential outputs are tempting to combine into site-level
results after model fitting, but same-experiment peptide estimates usually come
from the same samples and are correlated. Treating finite-degree-of-freedom
peptide t-statistics as independent normal observations can overstate site-level
evidence and changes the meaning of p-values.

PhosPy already has a safer peptide-evidence lane: resolve peptide evidence into
site-level sample intensities during dataset building, then fit the core
site-level differential model.

## Decision

The preferred PhosPy-origin lane aggregates or resolves peptide evidence before
differential model fitting:

1. `phospy.science.evidence` and dataset-building preprocessing resolve
   peptide evidence into site-level sample-intensity rows.
2. `DifferentialAnalysisWorkflow` fits the existing moderated fixed-effect
   model on those resolved site rows.
3. Provenance records peptide mapping, multi-site handling, attrition, and the
   downstream differential model policy separately.

PhosPy also supports a narrow advanced post-hoc estimate-combination lane through
`PeptideDifferentialEstimateTable` and `PeptideToSiteAggregator`. This lane is
not used for PhosPy-origin peptide evidence by default. It requires a typed input
model with:

- effect estimate;
- standard error;
- original statistic;
- original p-value;
- residual degrees of freedom;
- moderated degrees of freedom;
- source experiment/run identifier;
- dependence policy;
- peptide-to-site mapping policy.

The supported multi-estimate post-hoc methods require independent source
experiments or runs. Same-experiment peptide estimates are rejected because
same-sample peptide dependence is not modelled by the current lane.

## Supported Statistical Methods

Supported methods are:

- `single_estimate_passthrough`: used for sites with one estimate. This is not a
  meta-analysis. It preserves the original effect, standard error, statistic,
  finite-df p-value, residual degrees of freedom, and moderated degrees of
  freedom.
- `stouffer_signed_p_independent`: combines independent source estimates by
  converting each original two-sided p-value to a signed z value and combining
  weighted z values. The sign comes from the original statistic, falling back to
  the effect sign only when the statistic is exactly zero.
- `fixed_effect_inverse_variance_independent`: combines independent source
  estimates using the typed standard error as inverse-variance weight.

Finite-degree-of-freedom t evidence must not be treated as `z = t`. When a
z-scale input is needed, conversion is through signed two-sided p-values.

## Supported Dependence Assumptions

The current supported post-hoc dependence policy is:

- `independent_sources`: estimates for the same site must come from distinct
  source experiment/run identifiers and must be declared independent.

The current lane does not implement generalized least squares, known covariance
matrices, robust clustered covariance, mixed effects, duplicate-correlation
style methods, or any same-sample peptide dependence model.

## Minimum Evidence

`min_estimates_per_site` defaults to `1`. Sites below the configured minimum are
emitted with missing statistics and explicit attrition/provenance. The default
therefore lets a single estimate pass through exactly, but callers can require
more independent source estimates when appropriate.

## Multiple Testing

Post-hoc site p-values are adjusted with the shared configurable multiple-testing
correction domain. The post-hoc lane must not hardcode Benjamini-Hochberg outside
the configured correction policy.

## Multi-Site Mapping

Multi-site peptide allocation or retention remains explicit input provenance.
The typed estimate table records `peptide_to_site_mapping_policy` per estimate
and optional mapping uncertainty. Output tables and run provenance report the
observed mapping policies and multi-site estimate counts.

## Output Interpretation

Output tables record:

- aggregation level;
- dependence assumption;
- uncertainty method;
- p-value method;
- statistic distribution;
- correction method;
- source experiment/run identifiers;
- peptide-to-site mapping policy.

A post-hoc combined site row is a site-level summary under the recorded typed
uncertainty and dependence assumptions. It is not evidence that same-experiment
peptides were independent, and it is not a replacement for sample-level
peptide-to-site resolution followed by the core differential model.

## Unsupported Claims

This ADR does not claim:

- no post-hoc same-experiment peptide meta-analysis;
- independence of peptide estimates from the same samples;
- no limma `duplicateCorrelation`, mixed-effects, or clustered-covariance modelling;
- no broad upstream statistical result import lane;
- equivalence between single-estimate pass-through and meta-analysis;
- that finite-df t statistics can be used directly as z statistics.

## Responsibility Audit

Ownership boundaries are:

- validators validate eligibility and request shape only;
- interpreters resolve user policy and routing only;
- result assemblers attach already-computed outputs and provenance only;
- numerical uncertainty handling lives in
  `phospy.science.differential.aggregation` or the core
  `phospy.science.differential` implementation;
- peptide evidence sample-intensity resolution remains owned by
  `phospy.science.evidence` and dataset preprocessing/builder logic.
