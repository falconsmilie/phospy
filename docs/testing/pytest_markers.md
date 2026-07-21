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
publishing. The release check blocks release on:

| Marker or suite | Release status |
| --- | --- |
| Lint | Blocking through `ruff check .`. |
| Type checking | Blocking through `python scripts/run_pyright.py`. |
| Default non-parity suite | Blocking through `pytest -m "not parity"` over configured `testpaths`. |
| Checked-in reference bundles | Blocking through `python scripts/validate_reference_bundle_index.py --repo-root .`. |
| Built distributions | Blocking through `make build`, metadata checks, and packaged-reference validation. |
| `parity` | Blocking through `pytest tests/parity -m parity -s`. |
| `activity_parity` | Blocking because the activity parity file is also marked `parity`; CI also has a dedicated activity parity gate. |
| `performance` | Blocking through `pytest tests/performance -m "performance or release_gate"`. |
| `parity_diagnostic` | Non-blocking diagnostic unless intentionally promoted out of the exclusion. |
| `slow` | Not selected solely by marker for release; it runs only when also collected by a blocking selector. |

## Local Command Conventions

The current pytest defaults are configured in `pyproject.toml` with:

- `testpaths = ["tests/unit", "tests/integration", "tests/parity", "tests/workflows", "tests/validation", "tests/science", "tests/architecture"]`
- `addopts = -m "not parity"`

Based on that configuration:

- Default local run: `pytest`
- Include parity locally: `pytest tests/parity -m parity -s`
- Exclude slow tests in local loops: `pytest -m "not parity and not slow"`
- Performance-only validation: `pytest tests/performance -m "performance or release_gate"`
- Full release-check command: `make release-check`

The default local run deliberately omits release tests, threshold-bearing parity
tests, and performance contracts unless they are selected separately through
the release check. This process provides normal CI/build confidence, not formal
exact-source/exact-artifact attestation.

This page documents marker usage and the release-check command path used by
release CI.
