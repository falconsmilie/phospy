# Pytest Marker Model

This project uses pytest markers to describe suite intent. This document defines the
intended marker model without changing current CI behavior.

## Marker Intent

- `unit`: Focused tests for a single module/function contract with isolated collaborators.
- `integration`: Multi-component workflow tests across public boundaries.
- `parity`: Python vs R/PhosR parity checks against approved reference fixtures.
- `performance`: Performance contract checks for key scientific paths.
- `slow`: Long-running tests that are useful for deeper local validation.
- `reproducibility`: Deterministic replay/provenance checks required for release confidence.
- `golden`: Golden fixture contract checks for stable scientific/provenance outputs.
- `release_gate`: Release sign-off tests for critical confidence checks.

## Release-Gated Suites

`parity`, `performance`, and explicit `release_gate` checks are release-gated suites.
In practice this means they should be part of explicit release validation runs, even
if excluded from default fast local runs.

## Local Command Conventions

The current pytest defaults are configured in `pyproject.toml` with:

- `testpaths = ["tests/unit", "tests/integration", "tests/parity"]`
- `addopts = -m "not parity"`

Based on that configuration:

- Default local run: `pytest`
- Include parity locally: `pytest -m "parity"`
- Exclude slow tests in local loops: `pytest -m "not parity and not slow"`
- Performance-only validation (explicit path, because performance is outside `testpaths`):
  `pytest tests/performance -m "performance or release_gate"`
- Reproducibility/golden release checks:
  `pytest tests/unit/test_provenance_regressions.py tests/integration/test_kinase_workflow_integration.py::test_kinase_public_predmat_provenance_matches_golden_contract tests/integration/test_signalome_workflow_integration.py::test_signalome_l6_provenance_matches_golden_contract -m "release_gate and (reproducibility or golden)"`
- Full release-gate validation command:
  `make test-release-gate`

This page documents marker usage and the release-gate command path used by release CI.
