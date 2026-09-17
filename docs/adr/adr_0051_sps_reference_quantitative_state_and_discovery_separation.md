# ADR-0051: SPS Reference Quantitative State and Discovery Separation

## Status

- **ADR ID:** ADR-0051
- **Title:** SPS Reference Quantitative State and Discovery Separation
- **Status:** Accepted
- **Date:** 2026-09-17
- **Decision Type:** Scientific Architecture and Input Governance
- **Refines:** ADR-0027, ADR-0029, and ADR-0034

## Context

ADR-0027 preserved the distinction between the planned native SPS/RUV-style
lane and true RUV-III. ADR-0029 required governed controls, missingness,
diagnostics, and provenance for correction. ADR-0034 established explicit
quantitative-state authority. SPS discovery and replicate-aware
`ruv_iii_style` now exist, with narrow pinned external comparison evidence, so
the quantitative boundary between reference preparation, discovery, and target
correction must be recorded explicitly.

PhosR `getSPS` relies on caller convention to provide biologically meaningful
condition-relative measurements. Numeric matrices do not carry enough evidence
for PhosPy to infer scale, baseline, or quantitative meaning safely.

## Decision

PhosPy treats reference preparation, SPS discovery, and RUV correction as
separate scientific operations.

An SPS reference accepted by `SpsDiscoveryWorkflow` must have:

1. established `IntensityScaleKind.LOG2` scale;
2. established `QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE` meaning;
3. zero representing a caller-established biological reference/control
   baseline; and
4. scale and meaning-establishment evidence bound to the submitted matrix
   fingerprint.

`SpsReferenceDataset.from_condition_relative_log2(...)` is the supported
external assertion boundary. It copies and records the supplied values; it does
not choose a baseline, centre abundance, transform values, or infer meaning.
Absolute abundance, linear scale, unknown state, and absent, incompatible, or
unbound establishment evidence fail before ranking.

This requirement is scoped to SPS reference datasets. Target experiments
corrected by `sps_ruv_style` or `ruv_iii_style` continue to use normal governed
PhosPy input and preprocessing contracts; they do not have to be represented as
SPS-style condition-relative measurements.

SPS discovery returns a governed `ControlSiteSet`. Correction consumes that
object and does not require source reference matrices. Discovery provenance
records reference identity, quantitative state and establishment evidence,
matrix fingerprints, algorithm/version, contribution and attrition counts, and
selected controls. Correction provenance separately records method, control
provenance, replicate definition where relevant, `k`, missingness strategy,
diagnostics, and fingerprints. Neither provenance record must embed complete
SPS reference matrices.

The method identities remain distinct:

- `sps_ruv_style` is the existing PhosPy-native protected-design control-based
  correction and is not retrospectively renamed RUV-III.
- `ruv_iii_style` is the replicate-aware estimator. Its replicate-set mapping
  participates directly in estimation, so it is numerically distinct from
  `sps_ruv_style` and may legitimately produce different output.

PhosPy's explicit state rejection is a deliberate input-governance difference
from PhosR's caller convention, not evidence of a numerical algorithm
difference.

## Consequences

Positive consequences:

- invalid or ambiguous SPS evidence fails before it influences ranking;
- correction targets are not burdened with an unrelated SPS-reference state;
- discovered controls can be persisted and reused with bounded provenance; and
- documentation and parity claims can distinguish native behavior from narrow
  external numerical evidence.

Costs and limitations:

- callers must establish and document a defensible biological baseline before
  constructing references;
- PhosPy cannot rescue an arbitrary abundance matrix by guessing its meaning;
  and
- technical reuse of a `ControlSiteSet` does not prove biological transferability
  across contexts.

## External Evidence Boundary

Pinned fixtures compare consensus SPS intermediates and selection with PhosR
1.13.1 commit `1be74902b775833c64f5833e70538eaf843cf6a5`, and complete-input
finite-`k` RUV-III corrected matrices with `ruv` 0.9.7.1. Fixture preparation
starts from synthetic processed absolute log2 abundance, subtracts the
pre-established baseline-condition mean per site, and passes the identical
condition-relative matrices to both SPS implementations.

The evidence excludes missing-value RUV-III execution, partial-reference SPS
contributions, the historical `sps_ruv_style` estimator, the full PhosR
`RUVphospho` workflow, and empirical biological validation. It therefore does
not establish full PhosR equivalence.

### Compatibility-boundary clarification

The accepted ADR originally recorded the evidence exclusions above but omitted
two existing compatibility differences from this boundary. This clarification
does not change the decision or the implementation:

- PhosPy requires every RUV-III replicate set to contain at least two samples
  and rejects singleton replicate groups. The pinned `ruv` reference is more
  permissive.
- PhosPy rejects a requested `k` that exceeds the estimable latent rank. It does
  not silently cap or reduce `k`, so requested and effective `k` remain equal
  for a successful fit; the pinned `ruv` reference may reduce the effective
  factor count.

Both are conservative PhosPy input-validation contracts. They narrow external
comparison to inputs supported by both implementations and do not contradict
the demonstrated complete-input corrected-matrix parity within that domain.

## Related Records

- [ADR-0027: Target Future Native PhosR-Style SPS/RUV-III Correction](adr_0027_target_future_native_phosr_style_sps_ruv_iii_correction.md)
- [ADR-0029: Native SPS/RUV-Style Batch Correction Prerequisites](adr_0029_native_sps_ruv_style_batch_correction_prerequisites.md)
- [ADR-0034: Quantitative State, Motif Scoring Semantics, and Reference Context](adr_0034_quantitative_state_motif_semantics_and_reference_context.md)
- [SPS Discovery and RUV-Style Correction](../api/sps-ruv.md)
- [Scientific Coverage](../scientific-coverage.md)
- [Comparison With PhosR](../parity.md)
