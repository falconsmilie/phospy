# ADR-0006: Intensity-Scale and Processing-State Contract for PhosPy Datasets

## Status

- **ADR ID:** ADR-0006
- **Title:** Intensity-Scale and Processing-State Contract for PhosPy Datasets
- **Status:** Accepted
- **Date:** 2026-05-13

## Context

`AnalysisReadyPhosphoDataset` must carry explicit scientific state at the
dataset boundary. Earlier wording around `TransformationState` blended two
different concerns:

- quantitative matrix scale (`linear` vs `log2`)
- preprocessing-policy state (missing-data handling, normalisation,
  total-protein correction, site-matrix policy, and comparison-building policy)

That blending made ownership and validation ambiguous.

Recent implementation changes also established explicit intensity-scale
establishment modes, so the ADR must record the supported provenance model.

Update note (2026-07-21, quantitative-meaning provenance):
`IntensityScaleState` now separates numerical scale establishment from
quantitative-meaning establishment. Changing scientific meaning, such as moving
from abundance to phospho/total ratio, fold change, differential effect, or
activity score, requires an authority-gated semantic transition with its own
`QuantitativeMeaningTransitionProvenance`. Scale-compatible enum values are not
evidence that the scientific operation happened.

Update note (2026-07-13, enrichment identifier-set provenance):
`EnrichmentWorkflowRequest` may carry typed selected/background identifier-set
provenance. For PhosPy-derived quantitative identifier sets, the provenance must
include the shared `InputIntensityScaleEvidence` model and the enrichment-owned
typed derived-set provenance model described in
[ADR-0045](adr_0045_enrichment_derived_identifier_set_provenance.md). This
keeps enrichment provenance aligned with dataset and workflow intensity-scale
evidence without inferring scale from values, column names, diagnostics, or
labels. Manual and raw identifier lists remain valid without intensity-scale
evidence or derived-set provenance. Declared quantitative scale evidence
produces a role-specific enrichment caveat; observed transformation evidence is
recorded without that declared-only caveat.

## Decision

PhosPy uses two required boundary models on `AnalysisReadyPhosphoDataset`:

1. `intensity_scale_state: IntensityScaleState`
2. `processing_state: DatasetProcessingState`

`IntensityScaleState` is narrow and authoritative for quantitative scale and
quantitative meaning, but it keeps those two facts as separate provenance
contracts. `DatasetProcessingState` is authoritative for preprocessing-policy
summary.

Intensity-scale establishment is evidence-backed and explicit. Establishment is
tracked with `IntensityScaleEstablishmentMode`:

- `declared`: incoming matrix scale was explicitly declared
- `transformed`: scale changed through a scale-changing transformation path
- `identity`: pass-through preservation of an already established declared
  state without changing numeric values
- `derived`: scale established through a supported derived path

In addition to mode, provenance records an explicit establishment source:

- `transformed_by_phospy`: PhosPy performed/owned the establishing transformation lane
- `declared_by_user`: user declaration was preserved as declared input scale
- `restored_from_trusted_provenance`: state was reconstructed from trusted bundle/provenance payload

Configured intent alone must not assign `log2`. Established state must come
from a supported builder/transformer/bundle path with structured establishment
provenance.
Identity/pass-through paths must not establish `linear`/raw scale from
undeclared inputs; they may only preserve an explicit trusted declaration.

Declared-scale diagnostics are safeguards and audit aids only. They can flag
suspicious declared values (for example ranges that resemble raw linear values
when declared as `log2`), but diagnostics do not prove scientific truth of the
declaration and do not silently change scale.

Quantitative meaning is established or transitioned separately from intensity
scale. Initial base meanings may be declared by the caller or inferred from an
already established scale contract. Operation-derived meanings must be produced
by the operation that performed the scientific transformation and must record
source meaning, target meaning, operation and producer identifiers, immutable
parameters, relevant input/output table fingerprints, trace ID when available,
and deterministic caveat codes. Bundle reconstruction may restore trusted
serialized semantic provenance or explicitly migrate legacy payloads as
unverified; it must not reinterpret a missing record as derived evidence.

Corrective clarification (2026-08-01, additive preprocessing scale policy):
Until PhosPy has a general transition model for additive versus multiplicative
abundance operations, additive preprocessing may not produce datasets labelled
as established linear abundance. Median centring, fixed-effect batch
residualisation (`linear_residualize_batch`), and native SPS/RUV-style
residualisation require an established log2 phosphosite-abundance scale at the
point where the additive operation runs. Linear abundance users must first apply
the supported log2 transform or provide already-log2 abundance data with an
explicit `input_intensity_scale='log2'` declaration. A separate future
multiplicative median-scaling operation may support linear abundance by
division rather than subtraction.

Update note (2026-08-01, operation-level quantitative contracts):
Preprocessing stages that can affect quantitative values now carry an explicit
`QuantitativeOperationContract` resolved from stage-owned metadata. The contract
is the semantic authority for accepted input scale kinds, accepted quantitative
meanings, output scale transition, output meaning transition,
abundance-preservation status, negative-domain behavior, required evidence,
reversibility, and information-loss category. Registry construction rejects any
quantitative preprocessing stage that omits this contract.

The contract layer is intentionally split by responsibility:

- stage execution owns numerical logic and typed operation evidence;
- `science.transformations` owns the typed transition vocabulary and semantic
  state model;
- `science.datasets.preprocessing.state_builder` folds executed contracts into
  `IntensityScaleState` and mints quantitative-meaning provenance;
- workflows and public configuration DTOs request operations but do not own or
  mint semantic transitions.

State builders must consume typed contract metadata and typed stage evidence.
They must not infer scientific meaning from diagnostic text. Diagnostics may
mirror the resolved contract output for reporting and audit, but a diagnostic
value that disagrees with the typed transition is an error. This makes
total-protein correction a normal operation-derived meaning transition rather
than a special state-builder branch.
Contracts that affect pre-execution scale/meaning folding without minting a
separate operation-derived quantitative-meaning provenance event must declare
that explicitly.

Update note (2026-08-01, final numeric-semantic coherence):
Final analysis-ready dataset construction now validates the observed finite
numeric sign domain against the established intensity scale and quantitative
meaning. This is a meaning-aware dataset-boundary rule, not a generic
table-schema positivity rule. Linear phosphosite abundance and linear total
protein abundance must be non-negative. Signed scientific quantities remain
valid when their established meaning permits signed values, including centred
log abundance, phospho/total log ratios, fold changes, differential effect
sizes, activity scores, and mixed log-ratio/log-abundance matrices. Unknown
quantitative meaning has no numeric-domain contract and therefore cannot be
promoted to analysis-ready merely because the values look plausible.

The rule consumes only typed transformation state and matrix values. It must
not infer scientific meaning from column names, diagnostics text, or informal
labels. Table wrappers continue to own structural numeric validity
(DataFrame-ness, numeric dtype, finiteness, shape, and alignment); dataset
validation owns the scientific coherence between scale, quantitative meaning,
and observed numeric domain.

This ADR supersedes the old transformation-state wording but does not remove
compatibility behavior where historical names appear in non-contract internals.

## Consequences

- **Positive**
  - Dataset-boundary semantics are clearer and auditable.
  - Validation can check scale establishment and policy-state coherence
    independently.
  - Bundle/publisher payloads can report both scale and processing policy with
    less ambiguity.
- **Negative**
  - Existing docs/tests/fixtures required contract wording and key updates.
  - Some internal historical symbols may remain until touched, even though
    preferred terminology is now intensity-scale/processing-state.
- **Neutral**
  - Existing lanes remain supported; this change primarily clarifies and
    hardens contract semantics.

## Alternatives Considered

1. Keep `TransformationState` as the primary public model.
   Rejected because it conflates scale metadata with preprocessing policy.
2. Keep only `IntensityScaleState` and leave policy state implicit.
   Rejected because policy-level dataset meaning would stay fragmented.
3. Merge `DatasetPreprocessingReport` into boundary state.
   Rejected because operational diagnostics are larger and serve a different
   purpose than compact boundary contract state.

## Implementation Notes

- Dataset boundary enforcement is in
  `src/phospy/science/datasets/models.py` and
  `src/phospy/validation/transformations/state.py`.
- Intensity-scale establishment modes and structured establishment provenance
  are defined in `src/phospy/science/transformations/models.py`.
- Quantitative-meaning transition authority and provenance are defined in
  `src/phospy/science/transformations/_authority.py` and
  `src/phospy/science/transformations/models.py`. Public DTOs may request a
  declaration, but they do not mint transition authority.
- Operation-level quantitative contracts are defined in
  `src/phospy/science/transformations/quantitative_contracts.py`, declared by
  stage-owned contracts under
  `src/phospy/science/datasets/preprocessing/stages/`, validated by the
  preprocessing stage registry, pre-folded by the preprocessing pipeline and
  dataset-build executor before numerical execution, and folded into processing
  state by
  `src/phospy/science/datasets/preprocessing/state_builder.py`.
- Numeric intensity transformations (including default preprocessing log2) must
  execute through transformer implementations in
  `src/phospy/science/transformations/transformers/`; preprocessing stages
  orchestrate and report but do not own transformation science.
- Final numeric-semantic coherence is implemented in
  `src/phospy/science/transformations/state_coherence.py` and composed by the
  private dataset construction boundary. This validator checks observed
  numeric sign domain against established `IntensityScaleState.quantity`.
- Additive preprocessing scale eligibility is privately enforced by
  `src/phospy/science/datasets/preprocessing/quantitative_scale_policy.py` and
  composed by the dataset builder before preprocessing execution. Workflow
  validators do not reimplement this rule.
- Intensity-scale preservation is expressed as transformer capability metadata
  (`preserves_input_scale_state`, `changes_numeric_values`,
  `requires_established_input_state`); resolvers consume these capabilities and
  must not branch on concrete transformer classes.
- Bundle metadata and reconstruction use explicit `intensity_scale_state` keys
  under `src/phospy/io/bundles/_shared/intensity_scale_state.py`.
- Workflow validators consume established boundary state; they do not establish
  intensity scale or repair missing quantitative-meaning lineage.

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R.,
& Yang, J. Y. H. (2019). Kinase activity inference from quantitative
phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14),
i349-i356. https://doi.org/10.1093/bioinformatics/btz306

YangLab. (n.d.). *PhosR* (Version release) [Computer software]. GitHub.
https://github.com/PYangLab/PhosR
