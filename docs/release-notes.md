# PhosPy Release Notes

## Version 1.7.2

Release date: 2026-09-08.

These notes describe the changes since Version 1.7.1.

## Release Overview

PhosPy 1.7.2 is an experimental protein-aware differential-analysis release. It
adds an explicit opt-in lane that can fit phosphosite condition contrasts with a
matched total-protein covariate prepared during dataset building.

The ordinary phosphosite differential workflow remains the default. A dataset
may carry protein-aware preparation sidecars, but `DifferentialAnalysisWorkflow`
uses them only when the caller supplies `DifferentialAnalysisConfig` with a
`DifferentialProteinAwareModelConfig`. No breaking public API change is
introduced relative to 1.7.1.

## Kinase Scientific-Policy Versions

The current implementation owns these policy and schema versions:

| Policy | Implemented version |
| --- | ---: |
| KSEA activity policy | 5 |
| Membership-selection policy | 4 |
| Inferential policy | 4 |
| Membership payload schema | 2 |
| Membership-independence policy | 2 |

These versions govern the KSEA scientific contract recorded in provenance and
bundles: membership evidence, whether substrate membership was selected
independently of the tested matrix, whether ordinary KSEA p/q output is
eligible, and compatibility for persisted membership and provenance payloads.
They are compatibility and interpretation contract identifiers, not empirical
proof of scientific validity.

## Compatibility and Migration

- The new protein-aware estimator is selected with
  `DifferentialProteinAwareModelConfig(method="protein_covariate_adjusted_moderated_linear_model_v1")`
  inside `DifferentialAnalysisConfig`.
- Dataset inputs for this lane must be prepared with
  `DatasetProteinAwarePreparationConfig(policy="prepare_model_inputs")` during
  dataset building. This preparation records aligned phosphosite/protein model
  inputs and diagnostics; it does not transform the phosphosite matrix.
- The estimator requires established log2 phosphosite and total-protein values.
  Requests with prior `subtract_log_total` correction are rejected because
  subtraction and covariate adjustment are different analyses.
- `fixed_block` remains supported when each augmented protein-aware design is
  admissible. `duplicate_correlation` and actual technical-replicate
  aggregation are rejected for the protein-aware lane in this release.
- Unsupported protein-aware requests fail closed before fitting. Withheld rows
  remain visible in the full result index with typed status and reason fields.

## Major Additions

- Added the experimental
  `protein_covariate_adjusted_moderated_linear_model_v1` differential method.
  It reports the requested phosphosite condition contrast conditional on the
  matched measured total-protein covariate.
- Added dataset-owned protein-aware preparation inputs, sidecar binding, sample
  alignment checks, transformation-state checks, mapping eligibility checks,
  and no-fallback validation for protein-aware execution.
- Added grouped augmented linear modelling, per-site attrition, mapping/model
  diagnostics, post-fit diagnostics, protein-aware provenance, and quantitative
  input fingerprints.
- Added public documentation and an executable example for the preparation and
  differential workflow.
- Added synthetic scientific-validation fixtures and an independent-oracle
  estimator check for the implemented method. These are PhosPy-owned validation
  contracts, not external empirical validation or parity evidence.

## Fixes and Hardening

- Protein-aware eligibility now uses the resolved reliability-profile replicate
  threshold, so production and exploratory replicate policies apply
  consistently.
- Post-fit numerical failures are withheld for the affected site instead of
  aborting valid sites. Invalid post-fit sites are excluded from
  empirical-Bayes moderation and multiple-testing families.
- Protein-row grouping and duplicate validation avoid unnecessary quadratic
  scans while preserving the same supported semantics.
- Ordinary differential behavior is protected by regression coverage. The
  default non-protein-aware path remains unchanged when the protein-aware model
  is not requested.

## Scientific Scope

The protein-aware differential estimator is experimental. It is a matched
protein-covariate adjustment lane for a single phosphosite differential
workflow, not MSstatsPTM-style joint PTM/protein inference, not MSstatsPTM
parity, not production validation, not phosphorylation stoichiometry, not
occupancy estimation, and not causal separation of protein abundance and
phosphorylation regulation.

`subtract_log_total` remains a direct preprocessing transformation:
`log2_phospho - log2_total`. The protein-aware lane is different: it keeps the
phosphosite response intact and includes the matched measured total-protein
value as a covariate in the differential model.

The broader PhosPy scientific boundaries from 1.7.1 remain in force:
differential analysis is limited to tested design and contrast envelopes;
bundled runtime references are rat-only; enrichment is offline ORA over
caller-supplied collections; native SPS/RUV-style correction is not
PhosR-equivalent RUV/SPS/RUV-III parity; `duplicate_correlation` is a narrow
paired-design GLS route rather than a mixed-effects framework; and
kinase/signalome outputs should be interpreted through the documented workflow
assumptions and caveats.

Next: [Quickstart](quickstart.md) or [API Guide](api/guide.md).
