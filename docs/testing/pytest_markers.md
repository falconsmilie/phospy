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
- `release_gate`: Release sign-off tests for critical confidence checks.
- `activity_parity`: Hard activity-stage parity fixture and kernel gate.
- `parity_diagnostic`: Informational parity reporting that is not release-blocking.

## Release-Gated Suites

Public releases must run the authoritative release-gate command,
`make test-release-gate`. Default `pytest` is a local development check and is
not sufficient for release. The release gate blocks release on:

| Marker or suite | Release status |
| --- | --- |
| Git-index reference bundle check | Blocking through `python scripts/validate_reference_bundle_index.py` before release manifest tests. |
| Default non-parity suite | Blocking through `pytest -m "not parity and not performance and not release_gate"` over configured `testpaths`. |
| `release_gate` + `reproducibility` / `golden` | Blocking through explicit provenance and golden paths, including `tests/golden`. |
| `tests/release` reference manifest gates | Blocking through `pytest tests/release -m "release_gate"`. |
| `parity` | Blocking through `pytest tests/parity -m "parity and not parity_diagnostic" -s`. |
| `activity_parity` | Blocking because the activity parity file is also marked `parity`; CI also has a dedicated activity parity gate. |
| `performance` | Blocking through `pytest tests/performance -m "performance or release_gate" -q`. |
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
- Full release-gate validation command: `make test-release-gate`

The default local run deliberately omits release tests, threshold-bearing parity
tests, and performance contracts unless they are selected separately through the
release gate. Release-gate pytest commands include duration reporting and write
JUnit reports under `build/reports/`; CI retains those reports for release
review.

This page documents marker usage and the release-gate command path used by
release CI.
