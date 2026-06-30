# ADR-0019: Experimental Design, Contrast, and Replicate Contract

## Status

- **ADR ID:** ADR-0019
- **Title:** Experimental Design, Contrast, and Replicate Contract
- **Status:** Accepted
- **Date:** 2026-05-13

## Context

Differential workflows previously allowed matrix-first design/contrast surfaces
that left intent implicit. That made condition semantics, replicate handling,
and auditability fragile.

Implementation now enforces typed design and explicit technical-replicate policy
before design-matrix assembly and model fitting.

Dataset `sample_metadata` remains a passive, aligned metadata table. It is not
the scientific design contract for differential analysis, and adding metadata
columns does not implicitly define conditions, replicate structure, batches, or
blocks.

Current parity-protected differential execution remains deliberately narrow:
two-condition unpaired simple contrasts, with unsupported design features
rejected at validation boundaries.

Differential outputs include a `logFC` column. That quantitative interpretation
is only valid when phospho intensities are established as log2-scale before
differential interpretation/execution.

## Decision

Differential workflows use a typed contract:

1. `ExperimentalDesign` with `SampleDesignRecord` entries.
2. Typed `Contrast` definitions
   (`numerator_condition`, `denominator_condition`).
3. Contract validation in
   `phospy.validation.workflows.differential.ExperimentalDesignContractValidator`.
4. Explicit technical-replicate handling via
   `DifferentialAnalysisConfig.technical_replicate_policy`:
   - default `TechnicalReplicatePolicy.REJECT`
   - explicit aggregation modes `TechnicalReplicatePolicy.MEAN` and
     `TechnicalReplicatePolicy.MEDIAN`
5. Structured differential policy provenance on result objects through
   `DifferentialAnalysisResult.policy_provenance`, recording:
   - design formula and matrix summary
   - typed contrast definitions
   - replicate requirements and technical-replicate lineage
   - empirical-Bayes settings
   - p-value and adjusted p-value methods
   - missing-value handling policy
   - intentionally rejected unsupported design features (`batch`, `block`)

Technical-replicate handling is owned by the public
`DifferentialAnalysisWorkflow.run(...)` pipeline. Normal users should call that
workflow rather than internal replicate helpers. The workflow validator uses
`phospy.workflows.differential.replicates.TechnicalReplicateAggregationPlanner`
to enforce the explicit `technical_replicate_policy` contract and construct a
`TechnicalReplicateAggregationPlan` before design-matrix assembly. When that
plan requires aggregation, the workflow interpreter uses
`phospy.workflows.differential.replicates.TechnicalReplicateAggregator` to apply
the plan, aggregate the supported matrices, record technical-replicate lineage,
and revalidate the resolved design before execution.

`phospy.workflows.differential.replicates.TechnicalReplicateResolver` remains a
backward-compatible wrapper around the planner and aggregator. It is not the
current ownership point for differential workflow enforcement.

Differential workflow eligibility requires dataset phospho intensity scale to be
both:

- established/trusted through supported `IntensityScaleState` establishment, and
- `IntensityScaleKind.LOG2`.

Raw/linear user input is still supported at dataset-building/preprocessing
boundaries, where PhosPy may apply log2 transformation and establish scale
state before differential execution.

Repeated `biological_replicate_id` values within condition groups are treated as
technical replicates and require explicit aggregation policy. Aggregation also
requires `biological_replicate_id` on every design sample and consistent
optional group fields (`batch`, `block`) within each
`condition + biological_replicate_id` group.

No condition inference from sample names is allowed.
No condition/replicate/batch/block inference from dataset sample metadata is
allowed.

## Consequences

- **Positive**
  - Differential inputs are auditable and explicit.
  - Replicate policy is deterministic and validated before modeling.
  - Ambiguous replicate structure fails early with clear errors.
  - Future extensions can add richer modeling without breaking request shape.
- **Negative**
  - Matrix-only call sites must migrate to typed design/contrast objects.
  - Previously tolerated ambiguous replicate inputs now fail unless explicitly
    resolved.
- **Neutral**
  - Statistical parity claims remain scoped to documented executable lanes.

## Alternatives Considered

1. Continue allowing ad hoc design-matrix-only requests.
   Rejected because intent and replicate semantics stayed implicit.
2. Auto-detect and collapse technical replicates with hidden defaults.
   Rejected because silent aggregation can change scientific interpretation.
3. Allow executors to infer design semantics downstream.
   Rejected because boundary validation must enforce correctness before
   execution.

## Implementation Notes

- Request/config contracts:
  `src/phospy/api/requests.py` and
  `src/phospy/api/configs/differential.py`.
- Technical replicate policy enum:
  `src/phospy/science/differential/policy_models.py`.
- Differential structured policy provenance model:
  `src/phospy/science/differential/models/provenance.py`.
- Planner, aggregator, and compatibility-resolver behavior:
  `src/phospy/workflows/differential/replicates.py`.
- Differential request validation pipeline:
  `src/phospy/workflows/differential/validator.py`.
- Differential provenance construction:
  `src/phospy/workflows/differential/provenance.py`.
- Ownership registry alignment:
  `docs/validation-ownership.md` (Design matrix validity, Contrast validity,
  Replicate policy).
- Scope guardrails:
  `docs/scientific-coverage.md` and `docs/workflow_contracts.md`.

## References

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A
practical and powerful approach to multiple testing. *Journal of the Royal
Statistical Society: Series B (Methodological), 57*(1), 289-300.
https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth,
G. K. (2015). limma powers differential expression analyses for RNA-sequencing
and microarray studies. *Nucleic Acids Research, 43*(7), e47.
https://doi.org/10.1093/nar/gkv007

Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing
differential expression in microarray experiments. *Statistical Applications in
Genetics and Molecular Biology, 3*(1), Article 3.
https://doi.org/10.2202/1544-6115.1027
