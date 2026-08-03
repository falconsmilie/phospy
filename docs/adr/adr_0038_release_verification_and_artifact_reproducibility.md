# ADR-0038: Release Verification and Artifact Reproducibility

## Document Control

- **ADR ID:** ADR-0038
- **Title:** Release Verification and Artifact Reproducibility
- **Status:** Superseded
- **Date:** 2026-07-16
- **Decision Type:** Release and Reproducibility Governance

## Supersession Note

Superseded on 2026-07-20 by
[ADR-0039: Lightweight Solo-Maintainer Release Process](adr_0039_lightweight_solo_maintainer_release_process.md).
The retained evidence aggregation, exact selected-source identity, artifact
verification matrix, final attestation, and attested publication directory
described below are historical and are no longer active release requirements.
PhosPy deliberately keeps normal tests, constrained builds, metadata checks,
packaged-reference validation, trusted publishing, and installed wheel/sdist
execution verification, while no longer providing formal
exact-source/exact-artifact attestation.
The active release command is `make release-check`; release/golden tests are
selected through `make test-release-gates`, not the superseded
`make test-release-gate` command named in this historical decision.
Current CI also keeps a supported-version release-science matrix and a
dedicated Python 3.11 minimum-dependency lane under ADR-0039. Those current
checks are active release confidence controls even though this ADR's former
exact-source/exact-artifact attestation system remains superseded.
ADR-0039 now treats the 50,000 x 48 end-to-end workload as an opt-in local
benchmark, not an active release gate or CI responsibility.

## Context

Release confidence for a scientific package must be tied to a clean supported
environment, the exact selected source form, and the exact built artifacts. A
partial test pass, a single unsupported interpreter, or a wheel-only check can
make the release evidence look stronger than it is.

PhosPy supports Python 3.10, 3.11, and 3.12. Release tooling therefore needs to
prove installation and verification across those supported versions, while also
showing that source and binary distributions reproduce the reference-bundle
manifests and file hashes selected by the source identity.

## Decision

The superseded release-verification command was `make test-release-gate`.
Release CI and publish CI were expected to invoke that same documented command
rather than a private variant. ADR-0039 replaces this with `make release-check`.

This superseded ADR historically required:

- a clean constrained `[dev,test]` installation on Python 3.10, 3.11, and 3.12;
- the full default pytest suite on Python 3.10, 3.11, and 3.12;
- the full release gate on Python 3.10, 3.11, and 3.12;
- source-form-aware reference-bundle validation tied to a source-identity
  report;
- wheel and sdist builds through `make build`;
- wheel and sdist reference-bundle validation against the exact source identity
  evidence;
- wheel and sdist installation checks on every supported Python version; and
- retained duration, JUnit, release metadata, and performance-budget reports.

ADR-0039 intentionally does not retain the source-identity report, legacy
public-boundary attestation, or attested publication directory as active
release requirements.

The Pyright development requirement in `pyproject.toml` and the CI constraint
in `constraints/ci.txt` must move together. CI type checking must install the
constrained development environment instead of bypassing constraints.
The current minimum-dependency lane is intentionally separate from
`constraints/ci.txt`; it uses `constraints/minimum.txt` to test declared
lower-bound runtime/test dependencies rather than the current reproducible CI
stack.

## Consequences

Release evidence is heavier, but it is reproducible and auditable. A release can
no longer be justified by a single interpreter, a default-only pytest run, or a
wheel that validates only against its own embedded manifest.

Git worktrees still validate reference resources against staged Git index blobs.
Source archive identities additionally record a digest of the extracted source
tree used for the build, and source-tree identities recompute that normalized
tree digest during validation. Missing Git metadata is not a release failure
when equivalent deterministic source identity evidence is present.

Maintainers must update `Makefile`, `.github/workflows/ci.yml`,
`.github/workflows/publish.yml`, `constraints/ci.txt`, and contributor docs
together when release verification changes. Bounded performance-contract budget
changes must also update `docs/performance.md` and the constants in
`tests/support/performance_contracts.py`; optional local benchmark observations
do not define release budgets.

## Validation

The active release-policy tests under `tests/release/` enforce the documented
`make release-check` command, supported Python matrix, Pyright constraint
alignment, archive-level packaged-reference checks, installed wheel/sdist
verification, and the absence of the obsolete release-scale workload from
required gates.

This ADR amends ADR-0014's release-gate policy and ADR-0015's reference-data
release rules by making clean supported-version verification and artifact
reproduction mandatory release evidence.

## Amendment: Release-Scale Benchmark No Longer Release Evidence (2026-07-29)

The 50,000-site x 48-sample end-to-end workload is no longer a release gate or
CI responsibility. It is retained as an opt-in local benchmark because its
runtime and memory cost are disproportionate for routine hosted CI.

Historical references in this superseded ADR to retained performance-budget
reports do not require Python 3.10-3.12 CI execution, two consecutive CI
matrices, tracemalloc subprocess measurements, or release-blocking thresholds
for that workload.
