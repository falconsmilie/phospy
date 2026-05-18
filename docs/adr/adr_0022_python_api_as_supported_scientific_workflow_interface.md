# ADR: Python API as the Supported Scientific Workflow Interface

## Document Control

- **ADR ID:** ADR-0022
- **Title:** Python API as the Supported Scientific Workflow Interface
- **Status:** Accepted
- **Date:** 2026-05-18
- **Decision Type:** Architecture Decision Record

## Context

PhosPy scientific workflows require explicit, typed configuration and strict
validation across dataset building, workflow execution, and scientific policy
boundaries. Workflow outputs also require provenance that captures material
configuration and execution context for reproducibility and auditability.

The prior command-line workflow path could not safely and consistently express
the required scientific configuration surface and provenance obligations for the
supported public workflow contracts.

## Decision

PhosPy no longer supports CLI workflow execution as a public scientific
interface.

The supported scientific workflow execution interface is the Python API.

Reintroducing CLI workflow execution in a future release requires a separate
architecture decision.

Any future CLI workflow interface must satisfy all of the following before it
can be treated as supported:

- complete coverage of required scientific workflow configuration
- parity with Python API validation and boundary error behavior
- provenance capture sufficient for reproducibility and contract auditing
- documentation and test coverage that reflect the supported public surface

## Consequences

### Positive

- Public workflow execution is now aligned with the explicit Python API contract.
- Scientific configuration and provenance expectations are clearer and more
  enforceable.
- Documentation and tests can focus on one supported scientific execution path.

### Tradeoffs

- Users who prefer command-line workflow execution must use Python API paths at
  this time.
- Any future CLI reintroduction requires explicit design and governance work.

## Scope Boundaries

This ADR does not define a future CLI design.

This ADR does not prohibit future CLI support. It requires a separate decision
before CLI workflow execution can be considered supported again.
