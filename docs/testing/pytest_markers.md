# Pytest Marker Model

This project uses pytest markers to describe suite intent and release-blocking
policy.

## Marker Intent

- `unit`: Focused tests for a single module/function contract with isolated collaborators.
- `integration`: Multi-component workflow tests across public boundaries.
- `parity`: Python vs R/PhosR parity checks against approved reference fixtures.
- `performance`: Performance contract checks for key scientific paths.
- `slow`: Long-running tests that are useful for deeper local validation.
- `reproducibility`: Deterministic replay/provenance checks required for release confidence.
- `golden`: Golden fixture contract checks for stable scientific/provenance outputs.
- `release_gate`: Extended release-confidence tests for critical checks.
- `activity_parity`: Hard activity-stage parity fixture and kernel gate.
- `parity_diagnostic`: Informational parity reporting that is not release-blocking.

## Release-Check Suites

Public releases must run the maintainer release command, `make release-check`.
Default `pytest` is a local development check and is not sufficient for
publishing. The configured default `testpaths` deliberately omit
`tests/release`, `tests/golden`, and `tests/performance`; release-only suites are
selected by explicit Makefile targets. The release check blocks release on:

| Marker or suite | Release status |
| --- | --- |
| Lint | Blocking through `ruff check .`. |
| Type checking | Blocking through `python scripts/run_pyright.py`. |
| Default non-parity suite | Blocking through `pytest -m "not parity"` over configured `testpaths`. |
| Checked-in reference bundles | Blocking through `python scripts/validate_reference_bundle_index.py --repo-root .`. |
| Built distributions | Blocking through `make build`, metadata checks, and packaged-reference validation. |
| Installed distributions | Blocking through `make verify-installed-distributions`, which installs and executes the built wheel and sdist outside the checkout. |
| `release_gate`, `golden`, `reproducibility` in `tests/release` or `tests/golden` | Blocking through `make test-release-gates`, which runs `pytest -o addopts= tests/release tests/golden -m "release_gate or golden or reproducibility"`. |
| `parity` | Blocking through `pytest tests/parity -m "parity and not parity_diagnostic" -s`. |
| `activity_parity` | Blocking because the activity parity file is also marked `parity` and is not marked `parity_diagnostic`; CI also has a dedicated activity parity gate. |
| `performance` | Blocking through `pytest tests/performance -m "performance or release_gate"`. |
| `parity_diagnostic` | Explicitly excluded from the blocking parity target unless intentionally promoted out of the exclusion. |
| `slow` | Not selected solely by marker for release; it runs only when also collected by a blocking selector. |

The optional 50,000 x 48 release-scale benchmark lives under `benchmarks/`, not
`tests/`, and is invoked only with `make benchmark-release-scale`. It is not a
pytest marker category, is not collected by pytest, and is excluded from
`make test-performance`, `make release-check`, and CI.

CI runs the default non-parity suite, threshold-bearing parity suite, release
and golden gates, and performance contracts on each supported Python version:
3.11 and 3.12. Build and archive-level packaged-reference validation remain a
dedicated single-build artifact job; the uploaded wheel and sdist are then
installed and executed by a separate Python 3.11 and 3.12 verifier matrix.
Manifest-governed fixture byte integrity also runs on both Ubuntu and Windows
to catch checkout newline conversion regressions.

CI also has a Python 3.11 minimum-dependency lane using
`constraints/minimum.txt`. That lane installs the project with test
dependencies under declared lower-bound pins, runs `pip check`, then runs the
non-parity suite plus release/golden/reproducibility selectors that do not
require external tools.

## Local Command Conventions

The current pytest defaults are configured in `pyproject.toml` with:

- `testpaths = ["tests/unit", "tests/integration", "tests/parity", "tests/workflows", "tests/validation", "tests/science", "tests/architecture"]`
- `addopts = -m "not parity"`

Based on that configuration:

- Default local run: `pytest`
- Blocking parity validation: `pytest tests/parity -m "parity and not parity_diagnostic" -s`
- Exclude slow tests in local loops: `pytest -m "not parity and not slow"`
- Performance-only validation: `pytest tests/performance -m "performance or release_gate"`
- Optional local release-scale benchmark: `make benchmark-release-scale`
- Release/golden validation: `make test-release-gates`
- Installed wheel/sdist validation: `make verify-installed-distributions`
- Full release-check command: `make release-check`

The default local run deliberately omits release tests, golden tests,
threshold-bearing parity tests, and performance contracts unless they are
selected separately through the release check. `make release-check` is the
authoritative aggregate command. This process provides normal CI/build
confidence, not formal exact-source/exact-artifact attestation.

The release policy test suite includes a collection-only selector audit using
`tools/testing/release_selector_coverage.py`. It compares actual collected node
IDs and effective markers against the authoritative Makefile selectors, so a new
release-blocking node cannot be added without either being selected by a release
target or failing the release-policy test.

This page documents marker usage and the release-check command path used by
release CI.
