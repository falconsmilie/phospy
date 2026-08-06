# ADR-0034: Quantitative State, Motif Scoring Semantics, and Reference Context

## Status

- **ADR ID:** ADR-0034
- **Title:** Quantitative State, Motif Scoring Semantics, and Reference Context
- **Status:** Amended
- **Date:** 2026-07-09
- **Decision Type:** Scientific Architecture and Workflow Contract

Amended 2026-07-29 to require explicit kinase reliability intent and make
kinase activity execution opt-in.

Amended 2026-08-02 to require method-owned quantitative input contracts for
kinase scoring and kinase activity methods.

Amended 2026-08-05 to make ssGSEA-style substrate enrichment equal-value ties
method-owned and row-order invariant.

Amended 2026-08-06 to require distinct kinase reference-projection attrition
before source substrate identifiers are projected into dataset `site_key`
identity.

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
  `IntensityScaleEvidenceLevel`, `QuantitativeMeaning`,
  `QuantitativeMeaningTransitionProvenance`, and
  `QuantitativeMeaningEvidenceMode`
- `InputIntensityScaleEvidence`
- `ReferenceContext`
- `KinaseScoringModeInputContract`
- `MethodQuantitativeInputContract` and
  `ResolvedMethodQuantitativeInputContract`
- `KinaseReliabilityProfile`
- `ProfileSelfInclusionPolicy`

Workflow validators must continue to compose shared validation with
workflow-specific checks. Private dataset validators remain internal validation
support and must not be promoted through `phospy.api` or the root package.
Request DTOs may enforce narrow local type checks, but they must not become the
owner of dataset validation.

## Sequence Readiness Contract

This ADR distinguishes three related states:

- Analysis-ready sequence evidence: `AnalysisReadyPhosphoDataset` requires
  `site_sequence` for every row. The dataset boundary checks that the value is a
  non-empty, plausible amino-acid context aligned to phosphosite identity.
- Workflow-specific sequence-context readiness: sequence-aware workflows resolve
  the selected dataset/reference sequence and validate it against the active
  workflow and scoring-mode contract before execution.
- Kinase/motif-specific centered window readiness: Kinase Library resource-backed
  motif modes require the selected sequence to match the
  `KinaseLibraryResource.sequence_window` exactly, including
  upstream-plus-site-plus-downstream length, center index, `S/T/Y` center
  residue matched to the site, known sequence source, supported alphabet, and
  accepted padding/lowercase/modified-residue policy.

Base dataset validation is intentionally plausibility-level. Dataset
construction happens before a caller selects a workflow, scoring mode,
reference bundle, sequence conflict policy, or motif resource. It must therefore
require `site_sequence` evidence without implying that every value is biologically
correct, centered on the phosphosite, sourced for motif scoring, or compatible
with a future resource-specific window. Stricter centered-sequence validation
belongs to the sequence-aware workflow validators.

For peptide-evidence input, ADR-0020 owns the earlier collapse from peptide rows
to one site-level `site_sequence` value. Workflow dataset/reference conflict
policies must not be used to choose among contradictory supplied peptide
contexts for the same resolved site; those conflicts are rejected before
`AnalysisReadyPhosphoDataset` construction.

## Kinase Scoring Mode Semantics

Kinase scoring modes are explicit workflow contracts.
`kinase_scoring_mode_input_contract(...)` defines whether a scoring mode
requires site sequences, centered sequence context, substrate/reference overlap,
a `KinaseLibraryResource`, and profile construction. Every current kinase
scoring mode requires `site_sequence` evidence and centered sequence context.
The `site_sequence` evidence requirement matches the
`AnalysisReadyPhosphoDataset` boundary, but the centered-context requirement is
stricter and remains workflow-owned. Passing dataset construction does not make
a sequence motif-ready.

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

Kinase scoring reliability profiles are public policy labels and caller intent.
Public kinase workflow requests must not manufacture an implicit scoring
profile. Callers choose one of:

- `KinaseScoringConfig.exploratory()` for the historical permissive scoring
  preset.
- `KinaseScoringConfig.production(...)` for production reliability; callers
  must supply study-specific non-zero reference-overlap, sequence-support, and
  scored-site attrition floors.
- direct `KinaseScoringConfig(..., reliability_profile="custom")` construction
  for caller-defined settings.

Direct `KinaseScoringConfig()` construction without a reliability profile is
invalid. Modified direct construction is not inferred as custom unless the
caller explicitly selects `custom`. `production` requires at least five
substrates, leave-one-out profile scoring, production site-level localisation,
and caller-selected non-zero attrition thresholds with error-on-violation
behavior. Production is never inferred from strict localisation alone.

`KinaseWorkflowRequest.activity_config` defaults to `None`. Activity-like
execution is opt-in by providing a `KinaseActivityConfig`; omitting the field
must not run the historical simplified weighted activity stage. Provenance must
record the selected scoring reliability profile and `activity_config=None`
when the caller did not request activity.

## Kinase Method Quantitative Input Contracts

Each kinase scoring and activity method owns its own quantitative input
contract. A method contract declares accepted `IntensityScaleKind` values,
accepted `QuantitativeMeaning` values, required centring or standardisation,
missing-value treatment, profile-axis requirements, statistical interpretation,
and p-value interpretation when applicable. Shared validators may enforce only
the contract they are given; they must not contain a single global kinase scale
or meaning policy.

Profile-dependent kinase scoring modes may consume declared abundance or
total-corrected log-ratio profiles when the dataset state explicitly records
those meanings. They must not accept contrast or effect matrices as abundance
profiles. Those inputs are scale-sensitive: a linear matrix and its
log2-transformed matrix are valid as separate declared inputs, but they are not
equivalent and scoring does not transform one into the other. Motif-only
scoring records the resolved dataset scale and meaning for audit but does not
consume phospho values for sequence-motif score calculation. Contrast/effect
activity analyses in the single-dataset kinase workflow must use a scoring mode
whose scoring step does not consume quantitative profiles, such as motif-only
scoring, or a future explicit separate activity input.

Kinase activity contracts are method-specific:

- Simplified weighted substrate activity accepts sample-level or
  condition-summary abundance semantics on declared linear or log2 abundance
  scales. It computes heuristic substrate-supported means and does not produce
  p-values.
- KSEA-style z-score activity accepts log2 sample abundance,
  total-corrected log-ratio, log2 contrast fold-change, or pre-standardised
  effect semantics. It rejects linear abundance. It computes z-scores for
  eligible substrate sets, but ordinary normal-approximation p-values and
  optional Benjamini-Hochberg q-values are emitted only when typed
  substrate-membership provenance shows that membership was selected
  independently of the tested quantitative matrix. Adaptive profile-derived,
  fused profile/motif, and leave-one-out profile-derived membership reports
  descriptive z-scores with p/q values unavailable.
- ssGSEA-style substrate enrichment accepts only log2 contrast fold-change or
  pre-standardised effect semantics. It produces no p-values unless seeded
  permutations are explicitly requested, in which case p-values are empirical
  substrate-label permutation p-values with optional Benjamini-Hochberg q-values.
  Equal finite effect values are handled inside the activity method as tie
  blocks, never by dataset order, workflow pre-sorting, stable-sort accident, or
  lexical site-label meaning. Version 2 uses a midrank block-expectation policy:
  a block containing `h` substrates and `m` non-substrates contributes the
  expected rank-walk area over all within-block orders,
  `b * running_before + ((b + 1) / 2) * (h / n_substrates - m / n_non_substrates)`,
  where `b = h + m`; the walk then advances by that same block delta. Blocks
  containing only substrates or only non-substrates therefore contribute the
  ordinary uninterrupted hit or miss segment, while mixed blocks give all
  equal-valued sites equivalent treatment. Seeded permutation p-values use the
  same tie-block score definition. Their active deterministic RNG policy is
  profile-based: method ID, method version, `profile_id`, kinase, stream name,
  and caller-supplied seed define the scientific stream identity.

Invalid scale/meaning/typed activity-semantics combinations must fail before
workflow execution. Activity methods also guard their direct method boundary so
internal callers cannot silently relabel abundance as effect input, or effect
input as abundance. Provenance must record the resolved method contract used by
the run, including the observed scale, meaning, activity profile axis, and
activity quantitative semantics when available.

## Kinase Method-Specific Site Universes

Kinase workflow interpretation must resolve named site universes instead of
using one shared scoring index for all downstream methods. The typed resolved
contract separates:

- measured quantitative sites, which match `dataset.phospho.index` and define
  the quantitative activity/background rows;
- sequence-supported scoring sites, which are the only rows eligible for
  centered-sequence motif/profile scoring;
- reference-supported membership sites, which are projected reference
  kinase-substrate rows in the measured site-key namespace;
- predicted membership sites, which are the prediction matrix rows emitted by
  the scoring/prediction stage or supplied by a typed fixed-membership source,
  and are bounded by measured site identity rather than by sequence support;
- KSEA background sites, which default to the measured quantitative universe
  and are independent of sequence support unless a future explicit KSEA
  background policy says otherwise;
- ssGSEA effect-ranking sites, which are the declared contrast/effect input
  rows and are not inherited from motif-scoring sequence support.

Scoring science owns sequence eligibility and centered sequence-context
validation. Motif scoring must continue to reject invalid centered sequence
context before scoring. Activity science owns background, effect-ranking,
membership-intersection, finite-value, and inferential-eligibility rules. KSEA
and ssGSEA must not require sequence context or infer sequence eligibility.

Workflow interpretation composes these contracts without conflating them.
`activity_phospho_matrix.index` must match the measured quantitative universe
and is not required to equal `scoring_site_index`. The executor must consume
named resolved matrices/maps, such as scoring phospho rows for scoring, KSEA
background rows for KSEA, and ssGSEA effect-ranking rows for ssGSEA; it must not
silently choose a background by reusing the motif-scoring universe. Membership
maps are explicitly intersected with the selected method universe before the
activity method is called.

Reference projection is a separate provenance transition that happens before
`reference_supported_membership_sites` exists. The source universe is
`references.kinase_substrate_map.substrate_site`, whose values may be reference
display IDs, dataset `site_key` values, or unmatched source identifiers. The
projected membership universe remains dataset `site_key` rows only. Provenance
must therefore serialize source-reference projection attrition separately from
later site-key universe filters; it must not infer projection loss from the
already-projected membership universe.

Kinase universe provenance must distinguish:

- source-reference projection attrition;
- sequence-support attrition;
- membership intersection;
- finite-value attrition; and
- activity-background or effect-universe attrition.

One-to-many display projection is projection ambiguity/diagnostic evidence, not
attrition. A source display identifier that maps to multiple dataset `site_key`
rows remains matched; it must not be counted as removed.

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

`IntensityScaleEstablishmentProvenance` keeps scale, establishment mode, source,
evidence level, transformer, trace ID, and warnings as typed top-level fields.
Its `parameters` mapping is JSON-like evidence and follows ADR-0035 recursive
immutability: constructor input is frozen without key stringification, and
`to_payload()` returns fresh ordinary JSON containers.

Unknown evidence is not a permission to proceed silently. It is an explicit
state that validators and result assembly must treat according to workflow
policy. Scientific caveats must not be removed to make tests pass.

For trusted `AnalysisReadyPhosphoDataset.from_trusted_tables(...)`
construction, intensity-scale evidence is also a separate
`TrustedDatasetConstructionAssertions.intensity_scale` dimension. It is not
satisfied by the quantitative-meaning assertion; callers must explicitly record
how the already-established `IntensityScaleState` was known.
Direct `AnalysisReadyPhosphoDataset(...)` construction is sealed; advanced
trusted reconstruction must use `from_trusted_tables(...)` with complete
assertions.

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

Quantitative meaning is a provenance-bearing scientific fact, not a mutable
label attached to an intensity scale. A scale-compatible target value does not
prove that the corresponding transformation occurred. Any meaning change must
use the internal transformation authority and must carry
`QuantitativeMeaningTransitionProvenance` with source meaning, target meaning,
stable operation identifier, producer/component identifier, evidence mode,
immutable parameters, trace ID when available, deterministic caveat codes, and
data table fingerprints for derived transitions.

`DatasetBuildRequest.quantitative_meaning` is retained only as an explicit
caller declaration for supplied input matrices. The exact caller-declarable set
is:

- `unknown`
- `phosphosite_abundance`
- `phosphosite_log_abundance`

Those declarations record `declared_by_caller` evidence and may add a
`quantitative_meaning_user_declared` caveat. They are assertions about the input
matrix, not proof that PhosPy or an upstream tool performed a transformation.
The builder rejects caller declarations for meanings that require an operation
not performed by the builder, including `phospho_total_log_ratio`,
`mixed_phospho_total_log_ratio_and_phosphosite_log_abundance`,
`contrast_log2_fold_change`, `differential_effect_size`, and
`activity_score`.

Default base meanings inferred solely from an established linear or log2 scale
record `inferred_from_scale_contract` evidence. Successful subtract-log-total
protein correction records a derived transition to
`phospho_total_log_ratio` or the mixed total-protein meaning, with stage/event
identifier, parameters, consumed phospho/total input fingerprints, and produced
phospho output fingerprint.

## Reference Context

`site_key` does not include reference version.

`site_key` encodes protein-scoped biological row identity: organism, protein
namespace, protein identifier, residue, position, and optional isoform. It must
stay stable enough for analysis-ready row identity, joins, saved outputs, and
workflow tables. Adding reference version to `site_key` would make row identity
change whenever a reference bundle is rebuilt, even when the biological protein
coordinate is the same.

The organism encoded in `site_key` is a shared `Organism` value, not an
arbitrary string. Analysis-ready construction normalizes supported aliases
and rejects contradictions between dataset-level organism, row-level organism,
and decoded site-key organism. This keeps reference-context compatibility from
having to compensate for mixed or ambiguous dataset organism state.

Reference version, source name, proteome version, and table hash are provenance
and compatibility context, not row-key fields. That is why `ReferenceContext`
exists. It records comparable reference identity fields and derives a
`reference_context_id` from the identity payload.

Update note (2026-07-17, resolved organism identity): `ReferenceContext`
stores organism as the shared resolved `Organism` value, not as a free string.
Serialized payloads continue to emit the standard string value. Construction
and deserialization must reject unsupported organisms and contradictions among
dataset organism, row/site-key organism, run-provenance reference context, and
selected reference bundle provenance.

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

Reference organism compatibility remains separate from reference-context
version compatibility, but it follows the same single-organism policy:
`ReferencePreset.AUTO`, explicit presets, and explicit bundles compare against
the same resolved dataset `Organism`. Resolver-time checks are defense in depth;
dataset, provenance, and bundle construction boundaries must already prevent
internally contradictory organism/reference identity.

## Row-Attrition Provenance

Typed row-attrition records are causal execution facts, not explanations guessed
from initial-versus-final absence. A workflow may emit a `row_attrition` record
only when the stage that performed a filter has captured its immediate input
site-row index and output site-row index. The record stores counts and a capped,
deterministic sample of removed identifiers; public provenance must not
serialize complete removed-ID lists.

Compatibility metrics remain available through `row_attrition_metrics`, but
metrics do not create machine-readable stage causality. When a metric and a
typed causal record describe the same stage count, they must agree.

Kinase attrition threshold evidence is stored separately as
`KinaseWorkflowAttritionProvenance`. Kinase science and workflow components own
the calculation of attrition counts and policy outcomes; the public result
contract only validates and freezes the already-calculated `metrics`, `policy`,
and `policy_violations` payloads under ADR-0035.

Current classification:

- `AnalysisReadyPhosphoDataset` still requires `site_sequence`. Missing or
  structurally invalid sequence context is a validation failure, not a
  synthetic kinase or signalome attrition record.
- Explicit localisation requirements are validation preconditions for kinase
  and signalome. Missing or below-threshold localisation metadata fails before
  execution when the workflow contract requires it.
- Signalome `protein_group_id` grouping metadata is required by the Signalome
  validator. Legacy `protein_id` is accepted only as a migration alias. Missing
  grouping metadata is not inferred as a later execution filter, and core
  `protein_identifier` identity is not reinterpreted as grouping identity.
- Kinase reference/resource overlap is required for profile-dependent modes,
  but motif-only kinase scoring does not require substrate-reference overlap
  and must not emit a reference-overlap site-row record.
- Site-row attrition and site/kinase-pair attrition are distinct. Pair-only
  scoring loss remains pair attrition in metrics and must not be encoded as a
  site-row removal.
- Signalome score preconditioning and any actual downstream scoring/clustering
  site-retention filter are execution-stage filters and must emit typed records
  when they remove site rows.

## Signalome Matrix and Network Provenance

Signalome network provenance must preserve both the requested and effective
paired finite observation threshold when the bundle schema distinguishes them.
`None` in public config is caller intent, not an effective execution value; the
interpreter resolves it and provenance records the resolved minimum plus the
stable policy identifier that governed edge eligibility.

Signalome clustering preparation is a workflow-specific matrix policy, not a
generic dataset preprocessing operation. Provenance must record retained
dimension labels, fully missing dimensions dropped before clustering, partial
median-imputation counts, and the exact prepared-matrix fingerprint consumed by
tree construction or module-count selection.

## Consequences

Positive consequences:

- Quantitative state is explicit and auditable at dataset and workflow
  boundaries.
- Kinase scoring modes cannot silently change their required inputs.
- Kinase scoring and activity methods cannot silently reinterpret linear,
  log2, abundance, contrast, or effect inputs as interchangeable.
- Motif-only scoring remains clearly separated from contextual/profile-aware
  scoring.
- Reference compatibility can be checked without destabilizing `site_key`.
- Caveats carry interpretation risk without weakening validators.

Negative consequences:

- Adding a new quantitative meaning, kinase scoring mode, or activity method
  requires typed model updates, method contract declaration, provenance review,
  documentation generated or checked from that declaration, and tests.
- Reference-context checks may reject workflows that previously relied on
  ambiguous provenance.

Neutral consequences:

- `site_key` remains a biological row identity, while `ReferenceContext`
  carries reference-version compatibility.
- Declared scale remains acceptable only through established, auditable lanes
  and may still produce result caveats.

## Implementation Notes

- Quantitative state and evidence enums:
  `src/phospy/science/transformations/policy.py`
- Transformation value models and state transitions:
  `src/phospy/science/transformations/scale_values.py` and
  `src/phospy/science/transformations/scale_state.py`; the historical
  `src/phospy/science/transformations/models.py` route is identity-preserving
  only.
- Workflow-visible intensity-scale evidence and declared-scale caveats:
  `src/phospy/workflows/intensity_scale_evidence.py`
- Provenance payload model:
  `src/phospy/provenance/models/workflows.py`
- Kinase scoring config and public scoring-mode strings:
  `src/phospy/contracts/configs/kinase.py`
- Internal kinase scoring-mode input contracts:
  `src/phospy/workflows/kinase/scoring_mode_contracts.py`
- Method quantitative input contract models:
  `src/phospy/science/quantitative_method_contracts.py`
- Kinase activity method quantitative contracts:
  `src/phospy/science/activities/method_contracts.py`
- Shared method-contract validation primitive:
  `src/phospy/validation/workflows/method_quantitative.py`
- Kinase scoring execution branches:
  `src/phospy/workflows/kinase/scoring_runner.py`
- Kinase result caveat assembly:
  `src/phospy/workflows/kinase/caveats.py`
- Scoring policy enums:
  `src/phospy/science/scoring/policy_models.py`
- Protein-scoped `site_key` encoding:
  `src/phospy/science/sites/site_keys.py`
- Reference context and reference bundle models:
  `src/phospy/science/references/models.py`; reference manifest value models
  are split across `src/phospy/science/references/manifest_policy.py`,
  `src/phospy/science/references/manifest_files.py`,
  `src/phospy/science/references/redistribution.py`, and
  `src/phospy/science/references/manifest_model.py`, with
  `src/phospy/science/references/manifest.py` retained as an identity-preserving
  route.
- Shared phosphosite identity and reference-context contract implementation:
  `src/phospy/science/sites/identity_contracts.py`; the validation-package
  route `src/phospy/validation/identity_contracts.py` is a compatibility
  re-export.
- Workflow phosphosite sequence-context contract implementation:
  `src/phospy/science/sites/sequence_context.py`
- Kinase and signalome workflow validators:
  `src/phospy/workflows/kinase/validator.py` and
  `src/phospy/workflows/signalome/validator.py`
- Analysis-ready dataset boundary:
  `src/phospy/science/datasets/construction/analysis_ready.py`

## Related Records

- [ADR-0004: Reference Resolution Strategy and `ReferenceBundle` Contract for PhosPy](adr_0004_reference_resolution_strategy_and_referencebundle_contract.md)
- [ADR-0006: Intensity-Scale and Processing-State Contract for PhosPy Datasets](adr_0006_transformation_state_and_transformer_contract.md)
- [ADR-0007: Validation Domain Architecture for PhosPy](adr_0007_validation_domain_architecture.md)
- [ADR-0015: Reference and Fixture Data Policy for PhosPy](adr_0015_reference_and_fixture_data_policy.md)
- [ADR-0024: Protein-Scoped Phosphosite Row Identity](adr_0024_protein_scoped_phosphosite_row_identity.md)
- [ADR-0033: Result Caveats and Scientific Warning Ownership](adr_0033_result_caveats.md)
- [ADR-0047: ssGSEA Tie-Block Policy](adr_0047_ssgsea_tie_block_policy.md)
