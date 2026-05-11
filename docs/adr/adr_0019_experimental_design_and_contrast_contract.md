# ADR: Experimental Design and Contrast Contract

## Document Control

- **ADR ID:** ADR-0019
- **Title:** Experimental Design and Contrast Contract
- **Status:** Accepted
- **Date:** 2026-05-11
- **Decision Type:** Architecture Decision Record

## Context

Differential workflows previously accepted numeric design and contrast matrices
directly. That shape allowed valid models, but it also left workflow intent
implicit:

- condition labels were only indirectly encoded in matrix columns
- replicate semantics were not explicit
- batch/block metadata had no typed home
- request auditability depended on external conventions

Future workflow lanes (replicate-aware filtering, batch-aware modeling,
condition-aware imputation, paired designs) require a typed design contract
that is validated before statistical execution.

## Decision

PhosPy adopts a typed experimental-design contract for differential workflows:

1. `ExperimentalDesign` with `SampleDesignRecord` entries.
2. Typed `Contrast` definitions (`numerator_condition`,
   `denominator_condition`).
3. Standalone design-contract validation in
   `phospy.validation.workflows.differential`.

`DifferentialAnalysisRequest` now requires:

- `design: ExperimentalDesign`
- `contrasts: tuple[Contrast, ...]`

The validator resolves validated matrix-ready representations before the
interpreter/executor stages.

## Validation Responsibilities

The design validator owns:

- dataset/design sample alignment
- duplicate sample ID rejection
- empty condition-label rejection
- unknown contrast-condition rejection
- minimum replicate-count checks
- optional field alignment checks (`batch`, `block`)
- explicit unsupported-feature errors for currently non-executable design
  features (batch-aware and block/paired modeling in this release)

Executor logic does not parse or infer design semantics.

## Separation From Raw Metadata

`AnalysisReadyPhosphoDataset.sample_metadata` remains raw tabular metadata.
`ExperimentalDesign` is interpreted statistical intent for one analysis run.
Keeping these separate avoids overloading generic datasets with all statistical
semantics.

## Consequences

### Positive

- Differential runs require an auditable, typed design input.
- No condition inference from sample names.
- Request validation fails early with explicit diagnostics.
- Future extensions can add batch/block/paired/time-course support by extending
  design validation and design-to-model translation without breaking the request
  contract.

### Tradeoffs

- Existing matrix-only differential request call sites must migrate to typed
  design/contrast objects.
- Some previously tolerated ambiguous inputs now fail at validation boundaries.

## Scope Boundaries

This ADR does not:

- implement batch-adjusted or paired statistical fitting in the current release
- auto-derive design conditions from sample labels
- move design parsing into statistical execution stages

## References

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B (Methodological), 57*(1), 289-300.

Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research, 43*(7), e47.

Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. *Statistical Applications in Genetics and Molecular Biology, 3*(1), Article 3.
