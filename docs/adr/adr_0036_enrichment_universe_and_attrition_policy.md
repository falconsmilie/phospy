# ADR-0036: Enrichment Universe and Attrition Policy

## Status

- **ADR ID:** ADR-0036
- **Title:** Enrichment Universe and Attrition Policy
- **Status:** Accepted
- **Date:** 2026-07-16
- **Decision Type:** Scientific Workflow Contract

## Context

Offline enrichment uses a selected foreground, an explicit background universe,
and caller-supplied enrichment sets. If identifiers are outside the background,
PhosPy can either reject the request or drop those identifiers before ORA.

Silent dropping is risky for foreground identifiers. A namespace mismatch, an
incorrect background universe, or accidental gene/site mixing can otherwise
become a smaller statistical analysis without the caller noticing. Set members
have a different practical shape: caller-supplied reference collections are
often broader than a measured experiment-specific universe, so intersecting set
members with the explicit background can be a legitimate configured behavior.

## Decision

`EnrichmentConfig` owns enrichment universe policy.

Selected foreground identifiers outside `background_universe` default to
`"error"`. Callers who intentionally want intersection behavior must set
`selected_outside_background_policy="drop"`.

Reference/set-member outside-background behavior is separately configurable via
`set_member_outside_background_policy`. Its default is `"drop"` because broad
local set collections commonly exceed an experiment-specific background, but
callers can set `"error"` to reject any set member outside the background.

`minimum_retained_foreground_fraction` optionally enforces the fraction of
workflow-normalized selected identifiers that must remain after background
intersection. When selected dropping is enabled, a foreground that becomes empty
after filtering is always rejected.

Workflow validation owns these policy checks before interpretation or
execution. The enrichment interpreter may translate the resolved request policy
into method-engine configuration, but it must not invent, override, or hardcode
an attrition policy absent from the resolved request.

## Diagnostics and Provenance

Successful enrichment results record universe policy in diagnostics and run
provenance. Recorded facts include:

- selected and set-member outside-background policies;
- configured minimum retained foreground fraction;
- selected identifier count before background intersection;
- selected identifiers retained in the background;
- selected identifiers outside the background;
- retained foreground fraction;
- set-member identifiers outside the background, with deterministic bounded
  preview metadata; and
- set-size filtering counts when configured.

Validation errors include counts and offending identifiers where practical so
namespace mistakes fail loudly before a smaller analysis is executed.

## Consequences

Positive consequences:

- Namespace and universe mistakes are rejected by default for foreground
  identifiers.
- Intentional foreground intersection remains available and auditable.
- Broad reference collections can still be used against experiment-specific
  backgrounds without requiring users to pre-trim every set.
- Retained-foreground thresholds make excessive attrition a validation failure
  rather than a post-hoc warning.

Negative consequences:

- Some legacy requests with selected identifiers outside the background now
  require `selected_outside_background_policy="drop"` or a corrected universe.

Neutral consequences:

- ORA statistics still run on identifiers inside the explicit background.
- This decision does not add online enrichment resources, ranked enrichment,
  GSEA, ssGSEA, PTM-SEA, or pathway activity inference.

## Implementation Notes

The public policy fields live in `src/phospy/contracts/configs/enrichment.py`.
Workflow validation is implemented in
`src/phospy/validation/workflows/enrichment.py`. Interpreter propagation lives
in `src/phospy/workflows/enrichment/interpreter.py`, and result diagnostics plus
run provenance are assembled in `src/phospy/workflows/enrichment/executor.py`.

## Related Records

- [ADR-0030: Ranked-List Enrichment Prerequisites](adr_0030_ranked_list_enrichment_prerequisites.md)
- [ADR-0035: Provenance Immutability and Stable Serialization](adr_0035_provenance_immutability_and_stable_serialization.md)
- [Enrichment Workflow](../api/enrichment.md)
- [Workflow Contracts](../workflow_contracts.md)
