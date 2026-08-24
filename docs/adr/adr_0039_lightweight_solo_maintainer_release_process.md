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
external-consumer public API contract tests, threshold-bearing parity tests
excluding `parity_diagnostic`, performance contracts,
release/golden/reproducibility tests, checked-in reference-bundle validation,
and a fresh distribution build.

CI expands the release-science evidence beyond the single local aggregate
command by running the non-parity suite, external-consumer public API contract
suite, threshold-bearing parity suite, bounded performance contracts, and
release/golden gates on Python 3.11 and 3.12. Documentation validation remains
a separate maintenance command through `make docs-build`; it is not part of the
package release gate. Distribution building and archive-level packaged-reference
validation remain a dedicated single-build job so wheel publication is not
duplicated. The uploaded wheel and sdist from that job are then installed and
executed by a separate installed-distribution verifier matrix on Python 3.11 and
3.12.

The performance-contract CI job is a release blocker for bounded tests under
`tests/performance/`. The 50,000-site x 48-sample end-to-end workload is no
longer part of that job, `make test-performance`, or `make release-check`.
It is retained as an opt-in local benchmark under `benchmarks/` and is invoked
explicitly with `make benchmark-release-scale`. Its runtime and process-memory
observations are informational and do not establish release budgets.

CI also includes a Python 3.11 minimum-dependency lane. That lane uses
`constraints/minimum.txt`, not `constraints/ci.txt`, installs the project with
test dependencies under declared lower-bound pins, runs `pip check`, and then
runs the non-parity suite plus release/golden/reproducibility selectors that do
not require external scientific tools.

Default `pytest` remains a local development command. Its configured
`testpaths` omit `tests/contract`, `tests/release`, `tests/golden`, and
`tests/performance`. `make test-contract` selects `tests/contract` explicitly
and clears global pytest `addopts`, so the external-consumer public API contract
does not depend on local default collection. `make test-release-gates` selects
`tests/release` and `tests/golden` explicitly with the marker expression
`release_gate or golden or reproducibility` and clears global pytest `addopts`,
so global marker defaults cannot change the release-gate selector.

`make build` starts from an empty `dist/`, builds one wheel and one sdist using
the constrained no-isolation build policy, runs metadata checks, and validates
the packaged reference manifests and declared file hashes inside both archives.
It must work from a copied source tree without Git metadata. This archive-level
validation remains separate from installed execution.

`make verify-installed-distributions` installs exactly one wheel and exactly
one sdist from `dist/` into separate temporary environments outside the
checkout. Its installed probe runs Python with isolation enabled, asserts
`phospy.__file__` resolves inside the installed environment rather than the
source tree, imports the supported public package surface, loads the bundled
rat `l6_native` reference manifest, verifies every manifest-declared bundled
resource and SHA-256 digest, and exercises representative dataset,
differential, kinase, and resolved public-boundary contracts. The verifier must
not import repository tests, fixtures, or `conftest.py`.

Publishing uses GitHub trusted publishing. The publish workflow builds once on
the checked-out tag, uploads the freshly built `dist/` directory, waits for the
installed-distribution verifier matrix to pass against that uploaded artifact
set, and publishes that same workflow artifact to TestPyPI for `tv*` tags or
PyPI for `v*` tags.

## Consequences

The release process is easier to understand and operate. It preserves ordinary
CI, scientific tests, reference metadata checks, clean builds, package metadata
checks, archive-level packaged-reference validation, and substantive installed
wheel/sdist execution verification.

The deliberate trade-off is lower auditability. PhosPy no longer provides
formal exact-source/exact-artifact attestation, retained report aggregation, or
machine-enforced binding between human approval records and publication.

License, attribution, source, redistribution status, reference identity,
namespace, organism, and declared reference-file hashes remain release-relevant
metadata and must not be weakened to simplify publishing.

## Validation

Release policy tests check the Makefile command flow, CI and publish workflow
shape, dependency constraints, minimum-dependency lane, supported Python
release-science and installed-distribution matrices, public-consumer contract
reachability, standalone `docs-build` target shape, selector coverage,
archive-level packaged-reference build checks, and installed-distribution
verifier source constraints. They also audit release-reachable Make/workflow
command and local-helper import closure for effective 50,000 x 48 workload
dimensions, including renamed scripts/targets, simple arithmetic, positional
and keyword aliases, dictionaries, and configuration objects.
The selector coverage audit uses collection-only pytest subprocesses to compare
actual node IDs and effective markers against the authoritative release targets.
Scientific runtime invariants remain protected by focused unit, integration,
contract, parity, golden, release, validation, workflow, architecture, and
bounded performance tests. Installed artifact behavior is checked by standalone
release tooling rather than pytest source-tree behavior, with regression tests
that damage manifest-declared resources in both wheel and sdist artifacts and
assert clear installed-probe failures. Machine-dependent scale observations are
collected through explicit local benchmark scripts.

## Amendment: Release-Scale Scientific Summary Equality (2026-07-27)

Superseded by the 2026-07-29 amendment below. Historically, the 50,000-site x
48-sample release-scale performance contract treated
ordinary and tracemalloc-instrumented executions as equivalent only when their
compact scientific summaries match exactly. The summary records exact
dimensions, exact original/final missing-cell counts, the accepted
preprocessing stage sequence, compact stable input/output table fingerprint
records, preprocessing stage hash traces, processing-state completeness flags,
contrast names, tested-feature count, a stable differential result-table digest,
and relevant policy/workflow provenance digests. It uses PhosPy stable
hashing/serialization helpers and never Python's process-randomized `hash()`.

After the ordinary production run is asserted and summarized, the parent keeps
only timings and the summary before launching the traced child subprocess. The
full ordinary dataset, differential result, and result table are not retained
across the child workload. No new attestation or approval subsystem is added.
Runtime, tracemalloc memory, subprocess timeout, and RSS-reporting policies
remain separate; traced runtime remains diagnostic and is not compared with the
ordinary runtime budget.

The subprocess resource sampler spools child stdout/stderr to temporary files
while sampling RSS so the child can return JSON metrics, including the summary,
without blocking on platform pipe-buffer limits.

## Amendment: Release-Scale Workload Moved to Local Benchmark (2026-07-29)

The 50,000-site x 48-sample end-to-end workload is no longer a release gate or
CI responsibility. It is retained as an opt-in local benchmark because its
runtime and memory cost are disproportionate for routine hosted CI.

The benchmark lives at
`benchmarks/measure_release_scale_builder_differential.py` and is invoked
through `make benchmark-release-scale`. It preserves the full 50,000 x 48
builder, preprocessing, provenance/fingerprinting, and one-contrast
differential workload, but it runs the workload once and reports runtime,
capacity, dimension, missingness, and tested-feature metrics as informational
observations. There is no duplicate tracemalloc subprocess run and no
release-blocking runtime, memory, or child-timeout threshold.

`make release-check`, `make test-performance`, GitHub Actions workflows,
release workflows, tag workflows, scheduled workflows, and publication targets
must not invoke this benchmark.
The release-policy audit detects equivalent 50,000 x 48 dimensions through
required helper modules and transitive Make/workflow reachability; current
bounded 50,000 x 12 and 50,000 x 24 contracts remain permitted.

## Amendment: Python 3.10 Support Removed (2026-08-03)

Python 3.10 is no longer supported. The active supported interpreter matrix is
Python 3.11 and Python 3.12, with package metadata declaring
`Requires-Python: >=3.11,<3.13`.

The minimum-dependency lane remains part of the release process and now runs on
Python 3.11 as the lowest supported interpreter. Python 3.10-specific
compatibility dependencies and branches, including the `tomli` dependency and
the dynamic `tomllib`/`tomli` parser fallback, have been removed. Environment
provenance uses the Python 3.11 standard-library `tomllib` module directly.
