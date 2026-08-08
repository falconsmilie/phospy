# ADR: Peptide-to-Site Differential Uncertainty Policy

## Document Control

- **ADR ID:** ADR-0041
- **Title:** Peptide-to-Site Differential Uncertainty Policy
- **Status:** Accepted
- **Date:** 2026-07-29
- **Amended:** 2026-08-05
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

After the original ADR text, implementation review found that the post-hoc
peptide differential estimate-combination lane accepted mapping policies that it
did not execute coherently:

- `exclude_from_statistical_model` could be accepted while included in the
  calculation;
- `split_equal_weight` could be accepted without splitting or weighting;
- a typed input could produce a concrete combined estimate and p-value despite
  unresolved combined-effect, inference, and mapping semantics.

Withdrawal is the smallest scientifically safe correction. Completing a new
combination model requires a separate scientific design decision.

## Decision

The preferred and supported PhosPy-origin lane remains:

1. `phospy.science.evidence` and dataset-building preprocessing resolve
   peptide evidence into site-level sample-intensity rows using the
   ADR-0020 production policy
   `peptide_to_site_linear_abundance_fractional_allocation_arithmetic_mean_v1`.
2. `DifferentialAnalysisWorkflow` fits the existing moderated fixed-effect
   model on those resolved site rows.
3. Provenance records peptide mapping, multi-site handling, attrition, and the
   downstream differential model policy separately.

ADR-0041 does not redefine sample-level site construction. Its supported input
boundary is the final analysis-ready site matrix produced by ADR-0020 or an
equivalent externally resolved site-level dataset with explicit provenance.

The post-hoc peptide differential estimate-combination lane is withdrawn from
public support. It must not be exported as production functionality from
supported public facades, documented as supported, or silently executed through a
compatibility path. Its retained compatibility shell is internal/experimental and
fails closed with an error explaining that coherent combined effect/inference
semantics and executable peptide-to-site mapping semantics are not implemented.

The withdrawal status is
`unsupported_withdrawn_posthoc_estimate_combination_v1`.

## Amendment: Validation Evidence Does Not Reintroduce Post-Hoc Support

**Date:** 2026-08-08

The release validation fixture
`tests/fixtures/release_validation_regression/peptide_site_bias_regimes/`
quantifies duplicate, ambiguous, missing, and localisation-error bias for the
sample-intensity peptide-to-site resolution route. It is
`synthetic_validation`, not post-hoc differential estimate-combination support,
not external parity, and not empirical validation.

This evidence documents adverse-regime sensitivity for the supported
pre-differential route. It does not satisfy the required future conditions for
reintroducing post-hoc peptide-level differential estimate combination.

## Required Future Conditions Before Reintroduction

Future public support requires a new ADR-backed scientific model and executable
implementation that defines:

- peptide-to-site mapping semantics, including ambiguous evidence, equal
  splitting, statistical-model exclusion, and weighted allocation;
- the site-level combined estimand and effect interpretation;
- the inferential result, uncertainty statistic, p-value meaning, and degrees of
  freedom or asymptotic approximation policy;
- dependence handling for same-experiment peptide estimates and any
  independent-source assumptions;
- multiple-testing semantics and correction domain;
- provenance semantics that make the mapping, dependence, estimand, inference,
  and attrition policies auditable;
- public API, documentation, and tests proving unsupported mapping or dependence
  states fail closed.

Until those conditions are met, mapping policies must not be labelled
experimental while executing silently as ordinary evidence.

## Retained Internal Source

Science-owned source under `phospy.science.differential.aggregation` may be
retained only as internal future-work material. Retention is not public support.
If retained internal source is exercised in experiments, its provenance must
record the unsupported/withdrawn status and must not be presented as production
PhosPy functionality.

## Responsibility Audit

Ownership boundaries are:

- validators validate eligibility and request shape only;
- public request DTOs remain passive and do not execute combination semantics;
- interpreters resolve user policy and routing only;
- result assemblers attach already-computed outputs and provenance only;
- withdrawal enforcement belongs at the public compatibility/export boundary;
- numerical combination source, if retained for future work, remains
  science-owned and must not move into workflow orchestration, validation, or
  request DTOs;
- peptide evidence sample-intensity resolution remains owned by
  `phospy.science.evidence` and dataset preprocessing/builder logic;
- the core `DifferentialAnalysisWorkflow` remains independent of the withdrawn
  feature.

## Non-Claims

This ADR does not authorize or provide:

- a public post-hoc peptide-to-site differential estimate-combination route;
- coherent post-hoc combined site-level effects or inferential results;
- executable mapping semantics for split, exclude, keep-joint, or weighted
  post-hoc policies;
- independence of peptide estimates from the same samples;
- no limma `duplicateCorrelation`, no mixed-effects modelling,
  no clustered-covariance modelling, and no covariance-aware peptide
  combination modelling;
- equivalence between single-estimate pass-through and meta-analysis;
- that finite-degree-of-freedom t statistics can be used directly as z
  statistics.
