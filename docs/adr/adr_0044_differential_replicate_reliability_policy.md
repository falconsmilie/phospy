# ADR-0044: Differential Replicate Reliability Policy

## Status

- **ADR ID:** ADR-0044
- **Title:** Differential Replicate Reliability Policy
- **Status:** Accepted
- **Date:** 2026-08-02

## Context

Differential linear-model execution can sometimes compute coefficients and
moderated statistics for a contrasted condition represented by only one
biological replicate. Computability is not the same as production-supported
inference.

The previous public surface allowed callers to lower
`minimum_condition_replicates` directly. That made the exploratory lane look like
an ordinary integer tuning parameter and risked hiding a scientific support
change.

Technical replicate metadata creates a related risk: repeated measurements of
the same biological unit must not be counted as independent biological
replicates.

## Decision

Differential analysis has an explicit reliability profile:

- `reliability_profile="production"` is the default supported inference lane.
  It requires at least two biological replicates for every contrasted condition.
- `reliability_profile="exploratory_single_replicate"` is the only supported
  single-biological-replicate opt-in. It resolves the replicate threshold to one
  and records the run as exploratory.

Lowering `minimum_condition_replicates` below the production floor while staying
in `production` mode is rejected. The named reliability profile is the override;
the generic integer is not.

Biological replicate counts are evaluated by experimental-design validation
before numerical execution:

- declared biological replicate IDs are counted by unique
  `biological_replicate_id` within each condition;
- declared technical replicate IDs require biological replicate IDs on all
  affected samples;
- technical replicates never inflate the biological replicate count;
- fixed-block and paired designs keep their existing rank, residual degrees of
  freedom, block coverage, and contrast-estimability checks.

Exploratory single-replicate results carry a structured top-level
`ResultCaveat` with code
`differential_exploratory_single_replicate`. The caveat records the reliability
profile, production support flag, computability flag, resolved replicate
threshold, observed condition replicate counts, and contrasted conditions below
the production minimum. Result payload and bundle-like persistence surfaces must
preserve this caveat.

## Consequences

- Production differential inference fails early for one-versus-many or
  one-versus-one contrasted conditions.
- Exploratory execution remains possible for users who intentionally opt in and
  accept that the output is not production-supported inferential evidence.
- Numerical model code remains responsible for linear algebra, fitting, and
  defensive estimability checks only. It does not own replicate semantics.
- Existing tests that intentionally exercise low-residual-DoF computability must
  identify themselves through the exploratory profile.

## Alternatives Considered

1. Change the default minimum only.
   Rejected because it would not make exploratory execution explicit and would
   not protect callers who lower the integer threshold.
2. Keep `minimum_condition_replicates` as the single override.
   Rejected because a generic integer hides the reliability change.
3. Count technical replicate rows as independent samples.
   Rejected because repeated technical measurements do not create independent
   biological evidence.
4. Disable residual degrees-of-freedom or estimability checks for exploratory
   execution.
   Rejected because exploratory reliability does not make numerically
   non-estimable designs executable.

## Implementation Notes

- Public config contract:
  `src/phospy/contracts/configs/differential.py`.
- Reliability-profile resolution:
  `src/phospy/workflows/differential/reliability.py`.
- Biological-replicate counting:
  `src/phospy/validation/workflows/differential_design_rules.py`.
- Workflow validation composition:
  `src/phospy/workflows/differential/validator.py`.
- Policy provenance:
  `src/phospy/workflows/differential/provenance.py` and
  `src/phospy/science/differential/models/provenance.py`.
- Structured caveat construction:
  `src/phospy/workflows/differential/caveats.py`.
- Ownership audit:
  `docs/validation-ownership.md`.

## References

- [ADR-0019: Experimental Design, Contrast, and Replicate Contract](adr_0019_experimental_design_and_contrast_contract.md)
- [ADR-0033: Result Caveats and Scientific Warning Ownership](adr_0033_result_caveats.md)
