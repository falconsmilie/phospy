# ADR: Validation Domain Architecture for PhosPy

## Document Control

- **ADR ID:** ADR-0007
- **Title:** Validation Domain Architecture for PhosPy
- **Status:** Accepted
- **Date:** 2026-05-02
- **Decision Type:** Architecture Decision Record

## Abstract

This ADR defines ownership boundaries for shared validation primitives and
workflow/domain validators. Validation remains an internal architecture domain,
but ownership is now explicit enough to prevent drift and generic-helper sprawl.
The executable ownership registry is maintained in
`docs/validation-ownership.md`.

## Status

Accepted.

This ADR supersedes earlier generic statements about reusable validation.

Update note (2026-05-13): invariant-level ownership is now documented in
`docs/validation-ownership.md` (owner, enforcement point, exclusions, and test
coverage expectations).

## Context and Problem Statement

Current code has shared DataFrame validation primitives, workflow request
validators, and domain validators. Without explicit ownership rules:

- shared modules become dumping grounds
- scientific/domain rules get misplaced into generic validators
- workflow validators start owning execution behavior

The boundary must be explicit because it affects public-config behavior,
workflow correctness, and review governance.

## Decision Drivers

1. Keep reusable validation primitives reusable and narrow.
2. Keep domain/scientific rules near domain ownership.
3. Keep workflow validators as boundary enforcers, not executors.
4. Ensure public config presets still pass the same validation boundaries.

## Decision

### Core Ownership Rules

1. Workflow validators compose validation primitives.
2. Validators do not execute scientific workflows.
3. Generic DataFrame/index checks belong in shared validation primitives.
4. Domain-specific rules belong in domain validation modules.
5. Shared validation modules must not become dumping grounds.
6. Site identity/coherence validation is owned by phosphosite/dataset validation
   boundaries, not by generic DataFrame policy decisions.
7. Config presets and public config objects must still pass through validation.
8. Validators may return the original request/config when enforcing boundaries
   rather than constructing new execution objects.

### `validation/common/` Boundary

validation/common/ is for generic structural primitives only.

This includes rules such as:

- DataFrame type/shape checks
- required columns
- uniqueness checks
- finite/missing constraints
- strict generic alignment primitives

The following must live in narrower domain modules, not in
`validation/common/` governance:

- phosphosite identity rules
- reference identity/organism compatibility rules
- workflow-specific alignment boundaries
- scientific policy checks

### Workflow Validator Boundary

Workflow validators must:

- enforce request boundary correctness
- compose shared/domain validators
- reject invalid requests with explicit validation errors
- return unchanged request/config when appropriate

Workflow validators must not:

- run scoring, prediction, clustering, or provenance execution logic
- perform scientific transformations
- become catch-all utility modules

## Consequences

### Positive Consequences

- Validation ownership is reviewable and enforceable.
- Shared primitives stay generic and reusable.
- Domain logic stays close to scientific context.
- Workflow validators remain small and predictable.

### Negative Consequences

- Some existing helper placement may require future cleanup when touched.
- Reviewers must enforce module ownership instead of accepting convenience moves.

## Affected Modules

- `src/phospy/validation/common/dataframes.py`
- `src/phospy/validation/workflows/configs.py`
- `src/phospy/validation/references/bundle.py`
- `src/phospy/validation/references/compatibility.py`
- `src/phospy/workflows/kinase/validator.py`
- `src/phospy/workflows/signalome/validator.py`
- `src/phospy/tables/datasets.py`
- `docs/validation-ownership.md`

## Scope Boundaries

This ADR defines validation ownership and composition governance. It does not
define module-splitting rules (ADR-0010) or stochastic reproducibility policy
(ADR-0017).

## Validation and Review Criteria

Future changes must satisfy all of the following:

1. Is the rule generic structure or domain/scientific policy?
2. If generic, is it in `validation/common/` and still structurally generic?
3. If domain-specific, is it in a domain module with clear ownership?
4. Does the workflow validator compose validation without executing science?
5. Do presets/config objects still pass normal validation boundaries?
6. Does the change align with `docs/validation-ownership.md` for primary owner
   and enforcement point?

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi,
R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative
phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14),
i349-i356.

YangLab. (n.d.). *PhosR* [Computer software]. GitHub.
https://github.com/PYangLab/PhosR
