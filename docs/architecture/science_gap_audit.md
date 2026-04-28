# Scientific Coverage Audit

This page captures governance boundaries for scientific claims.

## Scope

PhosPy does not claim whole-package equivalence with PhosR.
Coverage is lane-specific and evidence-backed through committed fixtures and tests.

Use [Parity to PhosR](../parity.md) as the detailed inventory and threshold source.

## How To Read Coverage Claims

- Coverage tiers describe confidence and regression protection level.
- Status labels describe governance state.
- A closed status does not automatically mean broad parity.

The active vocabulary and current lane-by-lane inventory live in:

- [Parity to PhosR](../parity.md)

## Evidence Anchors

Primary evidence sources in this repository:

- `tests/parity/`
- `tests/integration/`
- `tests/unit/`
- `tests/fixtures/rewrite_parity/`
- `tests/fixtures/public_workflow_reference/`

## Update Rule

When scientific behaviour, parity gates, or thresholds change, update this page
and [Parity to PhosR](../parity.md) in the same change.
