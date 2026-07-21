# ADR-0039: Lightweight Solo-Maintainer Release Process

## Document Control

- **ADR ID:** ADR-0039
- **Title:** Lightweight Solo-Maintainer Release Process
- **Status:** Accepted
- **Date:** 2026-07-20
- **Decision Type:** Release and Packaging Governance

## Context

PhosPy previously used a bespoke retained-evidence release process. That system
improved auditability, but it required more machine policy, report aggregation,
artifact cross-checking, and publication staging than a solo maintainer can
operate reliably.

The project still needs release confidence for scientific behavior, reference
metadata, packaged resources, and Python packaging quality.

## Decision

The maintainer release command is:

```bash
make release-check
```

That command runs lint, type checking, the default non-parity test suite,
threshold-bearing parity tests excluding `parity_diagnostic`, performance
contracts, release/golden/reproducibility tests, checked-in reference-bundle
validation, and a fresh distribution build.

Default `pytest` remains a local development command. Its configured
`testpaths` omit `tests/release`, `tests/golden`, and `tests/performance`.
`make test-release-gates` selects `tests/release` and `tests/golden` explicitly
with the marker expression `release_gate or golden or reproducibility` and
clears global pytest `addopts`, so global marker defaults cannot change the
release-gate selector.

`make build` starts from an empty `dist/`, builds one wheel and one sdist using
the constrained no-isolation build policy, runs metadata checks, and validates
the packaged reference manifests and declared file hashes inside both archives.
It must work from a copied source tree without Git metadata.

Publishing uses GitHub trusted publishing. The publish workflow builds once on
the checked-out tag, uploads the freshly built `dist/` directory, and publishes
that same workflow artifact to TestPyPI for `tv*` tags or PyPI for `v*` tags.

## Consequences

The release process is easier to understand and operate. It preserves ordinary
CI, scientific tests, reference metadata checks, clean builds, package metadata
checks, packaged-reference validation, and a lightweight installed-wheel smoke
test.

The deliberate trade-off is lower auditability. PhosPy no longer provides
formal exact-source/exact-artifact attestation, retained report aggregation, or
machine-enforced binding between human approval records and publication.

License, attribution, source, redistribution status, reference identity,
namespace, organism, and declared reference-file hashes remain release-relevant
metadata and must not be weakened to simplify publishing.

## Validation

Release policy tests check the Makefile command flow, CI and publish workflow
shape, dependency constraints, selector coverage, and packaged-reference build
checks. The selector coverage audit uses collection-only pytest subprocesses to
compare actual node IDs and effective markers against the authoritative release
targets. Scientific runtime invariants remain protected by focused unit,
integration, parity, golden, release, validation, workflow, architecture, and
performance tests.
