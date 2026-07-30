# ADR-0042: Signalome Module-Selection Stability Diagnostics

## Status

- **ADR ID:** ADR-0042
- **Title:** Signalome Module-Selection Stability Diagnostics
- **Status:** Accepted
- **Date:** 2026-07-30
- **Decision Type:** Scientific Safety and Workflow Contract

## Context

Signalome module-count auto-selection is deterministic for a fixed input matrix,
configuration, and implementation. Determinism is necessary for reproducible
execution, but it does not show that the selected module count is scientifically
stable under small input or policy perturbations.

Returning only the selected count can therefore hide boundary cases where nearby
thresholds or tiny score changes select a different count. Selecting the modal
count alone is also insufficient because disagreement is the diagnostic signal.

## Decision

Automatic signalome module-count selection returns a typed stability report
alongside the selected count. The report is descriptive sensitivity analysis, not
an inferential statistical result.

The supported method is versioned as
`seeded_score_perturbation_and_threshold_grid` version `1`. It records:

- evaluation method and version;
- seed policy and concrete seed when available;
- number of seeded score perturbations;
- selected-count frequency across perturbations;
- assignment similarity metric and summary;
- threshold-grid sensitivity;
- `stable`, `unstable`, or `not_computable` status;
- limitations and a `not_computable` reason when applicable.

Seeded perturbations must be reproducible. Callers may provide a fixed seed, or
the science layer derives a deterministic seed from the scoring matrix and method
version. Reports where no perturbation was applicable record a not-applicable
seed policy. Unseeded resampling is not allowed.

Insufficient sample structure returns `not_computable`. The implementation must
not fabricate a stability score when too few phosphosite profiles, too few score
dimensions, or too few non-degenerate profiles are available.

The report must not contain p-values, confidence intervals, confidence
probabilities, or multiple-testing claims. Selected-count frequencies are counts
across seeded perturbations, not calibrated probabilities. Assignment similarity
is a descriptive partition-agreement score, not a confidence probability.

## Responsibility Boundary

Stability computation belongs in `phospy.science.signalomes.clustering`.
Workflow code may coordinate execution and pass typed diagnostics forward.
Result assembly may attach the already-computed report to result contracts and
provenance. Validators may check request shape and eligibility, but they must not
perform clustering, perturbation, threshold-grid, or module-stability
calculations.

## Output Interpretation

Signalome modules and kinase-network correlations remain descriptive
score-derived summaries. A stable report means the automatic module-count choice
was unchanged under the implemented seeded perturbation and threshold-grid checks.
It does not establish biological validity, inferential evidence, causality,
directionality, or experimental validation.

An unstable report means at least one implemented perturbation or threshold-grid
point changed the selected count or assignment partition. Downstream module
summaries remain usable as descriptive output, but the reported instability must
travel with the result and provenance.

A not-computable report means the implemented stability method was not applicable
to the observed input or execution path. This is distinct from stability and must
not be silently converted into a stable or modal selected count.

## Consequences

Positive consequences:

- Automatic module-count selection exposes deterministic reproducibility and
  scientific sensitivity separately.
- Boundary cases surface selected-count disagreement instead of hiding it behind
  a modal count.
- Result bundles and provenance can replay the stability interpretation from
  typed serialized diagnostics.
- Small or degenerate inputs report `not_computable` rather than fabricated
  stability evidence.

Negative consequences:

- Automatic signalome clustering performs additional deterministic perturbation
  work for eligible inputs.
- Large inputs may receive a guarded `not_computable` report until a bounded
  stability design is implemented for that scale.
- Golden fixtures and bundle schema tests must include the new typed report.

## Related Records

- [ADR-0017: Stochastic Reproducibility and Adaptive Prediction Seed Governance](adr_0017_stochastic_reproducibility_and_adaptive_prediction_seed_governance.md)
- [ADR-0033: Result Caveats and Scientific Warning Ownership](adr_0033_result_caveats.md)
- [ADR-0035: Provenance Immutability and Stable Serialization](adr_0035_provenance_immutability_and_stable_serialization.md)
- [ADR-0040: Signalome Statistical Evidence Policy](adr_0040_signalome_statistical_evidence_policy.md)
