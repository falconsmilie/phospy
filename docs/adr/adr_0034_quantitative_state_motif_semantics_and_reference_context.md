# ADR-0034: Quantitative State, Motif Scoring Semantics, and Reference Context

## Status

- **ADR ID:** ADR-0034
- **Title:** Quantitative State, Motif Scoring Semantics, and Reference Context
- **Status:** Accepted
- **Date:** 2026-07-09
- **Decision Type:** Scientific Architecture and Workflow Contract

## Context

Recent workflow work introduced explicit scientific state for quantitative
matrices, kinase motif scoring modes, profile self-inclusion policy, and
reference-context compatibility. These changes are implemented in typed domain
models and workflow validators, but future contributors need one decision
record that explains how these pieces fit together.

The governing principle is that scientific state must be explicit, typed, and
validated at the correct boundary. ADR text must not be used to excuse weak
validation or to make invalid inputs appear acceptable.

## Decision

PhosPy records quantitative state and reference context as first-class domain
models, not unstructured dictionaries or inferred prose.

The primary models are:

- `IntensityScaleState`, `IntensityScaleEstablishmentProvenance`,
  `IntensityScaleEvidenceLevel`, and `QuantitativeMeaning`
- `InputIntensityScaleEvidence`
- `ReferenceContext`
- `KinaseScoringModeInputContract`
- `ProfileSelfInclusionPolicy`

Workflow validators must continue to compose shared validation with
workflow-specific checks. Private dataset validators remain internal validation
support and must not be promoted through `phospy.api` or the root package.
Request DTOs may enforce narrow local type checks, but they must not become the
owner of dataset validation.

## Kinase Scoring Mode Semantics

Kinase scoring modes are explicit workflow contracts.
`kinase_scoring_mode_input_contract(...)` defines whether a scoring mode
requires site sequences, centered sequence context, substrate/reference overlap,
a `KinaseLibraryResource`, and profile construction. Every current kinase
scoring mode requires `site_sequence` evidence. This matches the
`AnalysisReadyPhosphoDataset` boundary, which requires complete
`site_sequence` metadata.

`kinase_library_contextual_motif` is contextual motif mode. It requires a
caller-supplied `KinaseLibraryResource`, normal kinase workflow reference
resolution, display-ID projection, site-sequence support, eligible
kinase-substrate-map context, substrate/reference overlap, and profile
construction. Its authoritative downstream score matrix is
`kinase_library_motif_scores`; the profile context remains part of workflow
eligibility and diagnostics.

`kinase_library_motif_only` is motif-only mode. It requires a
`KinaseLibraryResource`, validated site identity, and centered sequence context,
but it does not construct substrate-derived kinase profiles and does not require
quantified known-substrate profile overlap. Its authoritative downstream score
matrix is `kinase_library_motif_scores`. Result caveats must preserve the fact
that motif-only scores use sequence motif evidence only and should not be
interpreted as kinase activity or causal kinase assignment.

`combined_profile_motif` uses the workflow profile context and caller-supplied
motif scores together, with `combined_profile_motif_scores` as the
authoritative matrix.

The profile self-inclusion policy is explicit. `profile_self_inclusion_policy`
defaults to `allow`, preserving historical profile-scoring behavior: a known
substrate site may contribute to the kinase profile used to score that same
site. This condition is caveated because it can inflate exploratory support
scores for known substrates. The `leave_one_out` policy recomputes applicable
profile scores after excluding the scored site from its own kinase profile and
records diagnostics for cells that cannot be scored after exclusion.

## Intensity-Scale Evidence

Intensity-scale evidence is recorded separately from the numeric scale itself.
`IntensityScaleEvidenceLevel` has these implemented values:

- `observed_transformation`: PhosPy observed and owned the transformation event.
- `declared_by_user`: the user declaration established or preserved the scale.
- `inferred_from_metadata`: scale was inferred from metadata rather than a
  directly observed numeric transformation.
- `unknown`: evidence exists only as an unknown or legacy-compatible state.

Downstream workflows consume the evidence payload through
`input_intensity_scale_evidence_payload(...)` and preserve it in workflow
provenance. When user-declared scale evidence materially affects result
interpretation, workflows attach the shared
`input_intensity_scale_declared_by_user` result caveat. Caveats do not replace
validation: workflows that require established intensity scale must still reject
unestablished or incompatible state.

Unknown evidence is not a permission to proceed silently. It is an explicit
state that validators and result assembly must treat according to workflow
policy. Scientific caveats must not be removed to make tests pass.

## Quantitative Matrix Meaning

`QuantitativeMeaning` records the scientific interpretation of matrix values.
It is distinct from `IntensityScaleKind`, though the two must be coherent.

The implemented meanings are:

| Category | Implemented value | Scale rule | Meaning |
|---|---|---|---|
| Abundance | `phosphosite_abundance` | `linear` | Phosphosite abundance matrix on linear scale. |
| Log abundance | `phosphosite_log_abundance` | `log2` | Phosphosite abundance matrix after log2-scale establishment. |
| Total-corrected log ratio | `phospho_total_log_ratio` | `log2` | Log-scale phosphosite values after total-protein correction. |
| Contrast/effect | `contrast_log2_fold_change` | `log2` | Contrast matrix, not abundance input. |
| Contrast/effect | `differential_effect_size` | `log2` | Effect-size matrix, not abundance input. |
| Activity score | `activity_score` | `linear` or `log2` | Kinase activity or activity-like score matrix. |
| Mixed | `mixed_phospho_total_log_ratio_and_phosphosite_log_abundance` | `log2` | Mixed corrected and uncorrected phosphosite rows. |
| Unknown | `unknown` | `linear` or `log2` | Scale may be known, but scientific meaning is not established. |

Workflow quantitative validators decide which meanings are executable for a
workflow. Phosphosite-abundance workflows accept phosphosite abundance/log
abundance and total-corrected log ratio; mixed total-protein meaning requires an
explicit workflow opt-in. Differential analysis requires
`phosphosite_log_abundance`. Unknown quantitative meaning is rejected by default
unless a workflow adds an explicit policy for it.

## Reference Context

`site_key` does not include reference version.

`site_key` encodes protein-scoped biological row identity: organism, protein
namespace, protein identifier, residue, position, and optional isoform. It must
stay stable enough for analysis-ready row identity, joins, saved outputs, and
workflow tables. Adding reference version to `site_key` would make row identity
change whenever a reference bundle is rebuilt, even when the biological protein
coordinate is the same.

Reference version, source name, proteome version, and table hash are provenance
and compatibility context, not row-key fields. That is why `ReferenceContext`
exists. It records comparable reference identity fields and derives a
`reference_context_id` from the identity payload.

Compatibility checks are required whenever a workflow combines artifacts whose
site identity or reference-derived semantics may come from different reference
contexts, including:

- an input dataset and an explicit or resolved `ReferenceBundle`;
- a kinase result and downstream signalome execution;
- reference-derived dataset construction or bundle restoration paths where both
  sides have reference context; and
- result caveat assembly when unknown context was explicitly allowed earlier.

`validate_reference_context_compatibility(...)` is the shared validator.
Mismatched known contexts fail under every policy. Unknown contexts fail by
default through `ReferenceContextCompatibilityPolicy.REQUIRE_KNOWN_MATCH`.
Kinase and signalome workflows may proceed with unknown reference context only
when their public config explicitly selects
`ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT`. In that case,
the workflow must emit a warning-severity result caveat that records the policy,
workflow scope, operation, missing side or sides, and available left/right
`reference_context_id` values.

`ReferenceContext` is not license metadata and does not approve redistribution.
Reference manifests and reference validation remain responsible for license and
redistribution status. Contributors must not mark references approved without
verified evidence for the exact packaged files.

## Consequences

Positive consequences:

- Quantitative state is explicit and auditable at dataset and workflow
  boundaries.
- Kinase scoring modes cannot silently change their required inputs.
- Motif-only scoring remains clearly separated from contextual/profile-aware
  scoring.
- Reference compatibility can be checked without destabilizing `site_key`.
- Caveats carry interpretation risk without weakening validators.

Negative consequences:

- Adding a new quantitative meaning or kinase scoring mode requires typed model
  updates, validator policy, result/provenance review, and tests.
- Reference-context checks may reject workflows that previously relied on
  ambiguous provenance.

Neutral consequences:

- `site_key` remains a biological row identity, while `ReferenceContext`
  carries reference-version compatibility.
- Declared scale remains acceptable only through established, auditable lanes
  and may still produce result caveats.

## Implementation Notes

- Quantitative state and evidence enums:
  `src/phospy/science/transformations/models.py`
- Workflow-visible intensity-scale evidence and declared-scale caveats:
  `src/phospy/workflows/intensity_scale_evidence.py`
- Provenance payload model:
  `src/phospy/provenance/models.py`
- Kinase scoring config and public scoring-mode strings:
  `src/phospy/contracts/configs/kinase.py`
- Internal kinase scoring-mode input contracts:
  `src/phospy/workflows/kinase/scoring_mode_contracts.py`
- Kinase scoring execution branches:
  `src/phospy/workflows/kinase/scoring_runner.py`
- Kinase result caveat assembly:
  `src/phospy/workflows/kinase/caveats.py`
- Scoring policy enums:
  `src/phospy/science/scoring/policy_models.py`
- Protein-scoped `site_key` encoding:
  `src/phospy/science/sites/site_keys.py`
- Reference context and reference bundle models:
  `src/phospy/science/references/models.py`
- Shared identity and reference-context validation:
  `src/phospy/validation/identity_contracts.py`
- Kinase and signalome workflow validators:
  `src/phospy/workflows/kinase/validator.py` and
  `src/phospy/workflows/signalome/validator.py`
- Analysis-ready dataset boundary:
  `src/phospy/science/datasets/models.py`

## Related Records

- [ADR-0004: Reference Resolution Strategy and `ReferenceBundle` Contract for PhosPy](adr_0004_reference_resolution_strategy_and_referencebundle_contract.md)
- [ADR-0006: Intensity-Scale and Processing-State Contract for PhosPy Datasets](adr_0006_transformation_state_and_transformer_contract.md)
- [ADR-0007: Validation Domain Architecture for PhosPy](adr_0007_validation_domain_architecture.md)
- [ADR-0015: Reference and Fixture Data Policy for PhosPy](adr_0015_reference_and_fixture_data_policy.md)
- [ADR-0024: Protein-Scoped Phosphosite Row Identity](adr_0024_protein_scoped_phosphosite_row_identity.md)
- [ADR-0033: Result Caveats and Scientific Warning Ownership](adr_0033_result_caveats.md)
