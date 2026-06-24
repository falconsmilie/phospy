# ADR-0027: Target Future Native PhosR-Style SPS/RUV-III Correction

## Status

- **ADR ID:** ADR-0027
- **Title:** Target Future Native PhosR-Style SPS/RUV-III Correction
- **Status:** Accepted
- **Date:** 2026-06-14
- **Decision Type:** Architecture and Scientific Roadmap
- **Refined By:** ADR-0029
- **Implemented By:** `SpsRuvBatchCorrectionConfig` native SPS/RUV-style
  preprocessing correction. The implementation is not a PhosR-equivalent
  SPS/RUV-III parity claim.

## Decision

PhosPy targets native SPS/RUV-style correction as a preprocessing/normalisation
capability only after the required representation, provenance, validation, and
test contracts exist.

The supported public implementation is `SpsRuvBatchCorrectionConfig`, which is
a validated PhosPy preprocessing lane with explicit controls, protected design
metadata, missingness policy, unwanted-factor count, diagnostics, and
provenance. It is not PhosR-equivalent SPS/RUV-III parity.

ADR-0029 refines this record with prerequisite details:
[ADR-0029: Native SPS/RUV-Style Batch Correction Prerequisites](adr_0029_native_sps_ruv_style_batch_correction_prerequisites.md).

## Scope

The capability must live before downstream analysis workflows consume the
quantitative matrix. It belongs to preprocessing/normalisation, and must not be
owned by differential analysis, kinase analysis, enrichment, or signalome
execution.

`linear_residualize_batch` remains a limited fixed-effect residualisation
method. It is not native SPS/RUV-style correction, not PhosR-equivalent batch
correction, and not equivalent to PhosR-style RUV/SPS correction.

`ruv_readiness` and similarly named diagnostics are report-only readiness
signals. They must not select controls, estimate unwanted factors, apply
correction, or alter the analysis-ready matrix.

## Required Future Constraints

Future implementation must define explicit contracts for:

- SPS/control phosphosite selection;
- unwanted-factor count;
- temporary imputation rules;
- whether to restore missingness after correction;
- whether to carry an observation/imputation mask forward;
- whether to flag imputed positions in downstream outputs;
- whether to withhold features from downstream statistical testing; and
- provenance for every correction decision.

Temporary imputation for correction mechanics is not equivalent to observed
biological evidence. Any future design must make that distinction auditable
before corrected values reach differential, kinase, enrichment, or signalome
workflows.

## Non-Goals

This ADR does not implement SPS/RUV/RUV-III correction, add public workflow
flags, bundle control phosphosite resources, or claim PhosR batch-correction
parity.
