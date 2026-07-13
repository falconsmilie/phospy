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

Update note (2026-07-13, enrichment identifier-set provenance):
`EnrichmentWorkflowRequest` may carry typed selected/background identifier-set
provenance. For PhosPy-derived quantitative identifier sets, the provenance
must include the shared `InputIntensityScaleEvidence` model. This keeps
enrichment provenance aligned with dataset and workflow intensity-scale
evidence without inferring scale from values, column names, diagnostics, or
labels. Manual and raw identifier lists remain valid without intensity-scale
evidence. Declared quantitative scale evidence produces a role-specific
enrichment caveat; observed transformation evidence is recorded without that
declared-only caveat.

## Decision

PhosPy uses two required boundary models on `AnalysisReadyPhosphoDataset`:

1. `intensity_scale_state: IntensityScaleState`
2. `processing_state: DatasetProcessingState`

`IntensityScaleState` is narrow and authoritative for quantitative scale and
quantitative meaning; `DatasetProcessingState` is authoritative for
preprocessing-policy summary.

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
- Numeric intensity transformations (including default preprocessing log2) must
  execute through transformer implementations in
  `src/phospy/science/transformations/transformers/`; preprocessing stages
  orchestrate and report but do not own transformation science.
- Intensity-scale preservation is expressed as transformer capability metadata
  (`preserves_input_scale_state`, `changes_numeric_values`,
  `requires_established_input_state`); resolvers consume these capabilities and
  must not branch on concrete transformer classes.
- Bundle metadata and reconstruction use explicit `intensity_scale_state` keys
  under `src/phospy/io/bundles/_shared/intensity_scale_state.py`.
- Workflow validators consume established boundary state; they do not establish
  intensity scale.

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R.,
& Yang, J. Y. H. (2019). Kinase activity inference from quantitative
phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14),
i349-i356. https://doi.org/10.1093/bioinformatics/btz306

YangLab. (n.d.). *PhosR* (Version release) [Computer software]. GitHub.
https://github.com/PYangLab/PhosR
