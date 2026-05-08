# Pytest Marker Model

This project uses pytest markers to describe suite intent. This document defines the
intended marker model without changing current CI behavior.

## Marker Intent

- `unit`: Focused tests for a single module/function contract with isolated collaborators.
- `integration`: Multi-component workflow tests across public boundaries.
- `parity`: Python vs R/PhosR parity checks against approved reference fixtures.
- `performance`: Performance contract checks for key scientific paths.
- `slow`: Long-running tests that are useful for deeper local validation.
- `release_gate`: Release sign-off tests for critical confidence checks.

## Release-Gated Suites

`parity` and `performance` are intended to be release-gated suites. In practice this
means they should be part of explicit release validation runs, even if excluded from
default fast local runs.

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
- Release-gated validation (explicit paths):
  `pytest tests/unit tests/integration tests/parity tests/performance -m "parity or performance or release_gate"`

This page documents intended usage only; it does not change CI behavior in this ticket.
