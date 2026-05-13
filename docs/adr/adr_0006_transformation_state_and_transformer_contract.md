# ADR: Intensity-Scale and Processing-State Contract for PhosPy Datasets

## Document Control

- **ADR ID:** ADR-0006
- **Title:** Intensity-Scale and Processing-State Contract for PhosPy Datasets
- **Status:** Accepted
- **Date:** 2026-04-26
- **Decision Type:** Architecture Decision Record

## Abstract

This ADR defines the scientific state contract carried by
`AnalysisReadyPhosphoDataset`.

The previous label `TransformationState` was too broad for what the model
actually represented. The old object only encoded quantitative intensity scale
(`linear` or `log2`) for phospho and optional total matrices. It did not
describe broader preprocessing policy state.

PhosPy now separates these concerns explicitly:

- `IntensityScaleState`: narrow quantitative scale state only.
- `DatasetProcessingState`: compact summary of preprocessing policy state at
  the analysis-ready boundary.

## Status

Accepted.

This ADR remains aligned with ADR-0003 (dataset boundary), ADR-0007
(validation architecture), and ADR-0011 (builder public contract). It
supersedes earlier transformation-state wording.

## Context and Problem Statement

The name `TransformationState` made the dataset boundary easy to misunderstand.
Readers could incorrectly assume it represented full preprocessing history.

In reality, a built dataset may have crossed multiple preprocessing policies:

- missing-data policy and optional imputation
- normalisation policy
- total-protein correction policy
- site-matrix construction and duplicate-site policy
- comparison-building policy

Those policy choices should be explicit at the dataset boundary, but should not
be conflated with matrix scale.

## Decision

PhosPy uses two distinct state objects on `AnalysisReadyPhosphoDataset`:

1. `intensity_scale_state: IntensityScaleState` (required)
2. `processing_state: DatasetProcessingState` (required in supported builder
   lane)

`IntensityScaleState` is authoritative for quantitative scale interpretation.
`DatasetProcessingState` is authoritative for preprocessing-policy summary.

Intensity-scale establishment is evidence-backed, not expectation-backed.
Configured target policy alone must never mint a scale label. A `log2` state is
valid only when backed by either:

- an executed scale-changing transformation path, or
- an explicit trusted declaration that the incoming matrix is already `log2`.

Identity pass-through may preserve declared state, but it must not upgrade
unknown/raw input to `log2` solely because a caller expects `log2`.

Intensity-scale establishment also records explicit establishment mode:

- `declared`: user declared the incoming scale
- `transformed`: PhosPy executed a scale-changing transformation
- `identity`: PhosPy identity pass-through established the state without changing values
- `derived`: a supported transformer-derived establishment path

Declaration remains auditable metadata, not proof of scientific correctness.
PhosPy records declaration diagnostics/warnings for suspicious declared scales,
but does not silently override the user declaration.

The old public name `TransformationState` is no longer the preferred dataset
contract model.

## State Responsibilities

`IntensityScaleState` answers:
- Are quantitative values `linear` or `log2`?
- Was that scale established through a supported path?
- Which establishment mode produced that state (`declared`/`transformed`/`identity`/`derived`)?

`DatasetProcessingState` answers:
- Which preprocessing-policy state crossed the analysis-ready boundary?
- Was imputation applied?
- Which imputation method/assumptions were used and on which cells/rows?
- Was total-protein correction applied?
- Was a site matrix constructed?
- How were comparisons configured?

`DatasetPreprocessingReport` remains separate and answers:
- What happened operationally during preprocessing, with diagnostics/tables?

## Model Direction

### Intensity Scale Model

`IntensityScaleState` is intentionally narrow and includes:

- phospho matrix scale (`linear` or `log2`)
- optional total matrix scale (`linear` or `log2`)
- establishment metadata

### Processing-State Model

`DatasetProcessingState` is a compact dataclass summary containing:

- `intensity_scale`
- `missing_data`
- `normalisation`
- `total_protein_correction`
- `site_matrix`
- `comparisons`

This model is not a replacement for the full preprocessing report.

## Builder and Workflow Responsibilities

Builders and preprocessing paths must:

- establish `IntensityScaleState` via supported establishment paths
- ensure establishment evidence matches the resulting declared scale
- construct `DatasetProcessingState` from the active preprocessing plan
- ensure dataset boundary coherence between both states

Workflows may:

- consume both states as established boundary metadata

Workflows may not:

- establish intensity scale
- infer missing processing state
- reinterpret absent policy metadata heuristically

## Validation Direction

Validation domain responsibilities include:

- intensity-scale type and establishment checks
- coherence between `dataset.intensity_scale_state` and
  `dataset.processing_state.intensity_scale`
- strict analysis-ready guardrails (for example complete-matrix expectation at
  the public dataset boundary)
- missing-data/imputation compatibility checks against intensity-scale
  assumptions (for example MinProb requiring log2-scale state)

## Bundle and Publisher Contract

Bundle manifests and workflow publisher metadata use explicit keys:

- `intensity_scale_state`
- `processing_state`

Publisher workflow metadata uses:

- `intensity_scale`
- `processing_state` (structured payload)

`transformation_state` is no longer the preferred outward key for this
boundary.

## Consequences

### Positive

- Scientific contract clarity improves.
- Dataset boundary meaning is less ambiguous.
- Scale state and preprocessing policy state can evolve independently.
- Validation and documentation become clearer.

### Negative

- Public/internal naming changed in multiple modules.
- Existing tests/docs/fixtures required updates.

### Neutral

- Error class names may retain historical identifiers where not contractually
  harmful, but messages should reflect intensity-scale wording.

## Rejected Alternatives

### Alternative 1: Keep `TransformationState` as the Public Dataset State Name

Rejected because it conflates narrow scale metadata with broader preprocessing
state.

### Alternative 2: Keep Only `IntensityScaleState`, Leave Preprocessing State Implicit

Rejected because policy-level analysis-ready state would remain fragmented and
harder to reason about.

### Alternative 3: Merge Preprocessing Report Into Dataset State Object

Rejected because report payloads are operational/diagnostic and should remain
separate from compact contract state.

## Resolved Decisions

1. `IntensityScaleState` is the explicit narrow replacement for the previous
   transformation-state model.
2. `DatasetProcessingState` is required as the broader analysis-ready policy
   state summary in the supported builder lane.
3. Builder output datasets must expose both `intensity_scale_state` and
   `processing_state`.
4. Bundle and publisher manifests use `intensity_scale_state` and
   `processing_state` keys.
5. Docs must not imply that intensity scale alone captures full preprocessing
   history.

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi, R.,
& Yang, J. Y. H. (2019). Kinase activity inference from quantitative
phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14),
i349-i356.

YangLab. (n.d.). *PhosR* [Computer software]. GitHub.
https://github.com/PYangLab/PhosR
