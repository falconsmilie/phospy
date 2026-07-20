# ADR-0038: Release Verification and Artifact Reproducibility

## Document Control

- **ADR ID:** ADR-0038
- **Title:** Release Verification and Artifact Reproducibility
- **Status:** Accepted
- **Date:** 2026-07-16
- **Decision Type:** Release and Reproducibility Governance

## Context

Release confidence for a scientific package must be tied to a clean supported
environment, the exact Git checkout, and the exact built artifacts. A partial
test pass, a single unsupported interpreter, or a wheel-only check can make the
release evidence look stronger than it is.

PhosPy supports Python 3.10, 3.11, and 3.12. Release tooling therefore needs to
prove installation and verification across those supported versions, while also
showing that source and binary distributions reproduce the committed
reference-bundle manifests and file hashes.

## Decision

The authoritative release-verification command is `make test-release-gate`.
Release CI and publish CI must invoke that same documented command rather than
a private variant.

Release verification requires:

- a clean constrained `[dev,test]` installation on Python 3.10, 3.11, and 3.12;
- the full default pytest suite on Python 3.10, 3.11, and 3.12;
- the full release gate on Python 3.10, 3.11, and 3.12;
- explicit Git-index validation through
  `python scripts/validate_reference_bundle_index.py`;
- wheel and sdist builds through `make build`;
- wheel and sdist reference-bundle validation against the committed Git index;
- wheel and sdist installation checks on every supported Python version; and
- retained duration, JUnit, release metadata, and performance-budget reports.
- installed-artifact `public-boundary-integrity` verification with passing
  detail outcomes for public signatures, dataset provenance binding, public
  DataFrame ownership, and public JSON immutability.

The Pyright development requirement in `pyproject.toml` and the CI constraint
in `constraints/ci.txt` must move together. CI type checking must install the
constrained development environment instead of bypassing constraints.

## Consequences

Release evidence is heavier, but it is reproducible and auditable. A release can
no longer be justified by a single interpreter, a default-only pytest run, or a
wheel that validates only against its own embedded manifest.

Maintainers must update `Makefile`, `.github/workflows/ci.yml`,
`.github/workflows/publish.yml`, `constraints/ci.txt`, and contributor docs
together when release verification changes. Performance budget changes must also
update `docs/performance.md` and the constants in
`tests/support/performance_contracts.py`.

## Validation

The release-policy tests under `tests/release/` enforce the documented command,
the supported Python matrix, the Pyright constraint alignment, the Git-index
reference-bundle check, the wheel/sdist artifact checks, and the
public-boundary integrity attestation details.

This ADR amends ADR-0014's release-gate policy and ADR-0015's reference-data
release rules by making clean supported-version verification and artifact
reproduction mandatory release evidence.
