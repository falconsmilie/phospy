# Type Checking (Pyright)

Pyright is the required static type-checking gate for PhosPy’s core scientific and
public API modules.

The intent is scientific and API reliability: type boundaries prevent invalid
states from reaching numerical kernels, workflow orchestration, and result
assembly/export paths.

## Why Pyright

- Fast, deterministic static checks in CI.
- Strong support for typed dataclasses, protocols, enums, and structured models.
- Clear diagnostics for boundary mismatches before runtime.

## Current Checked Scope

Pyright currently checks:

- `src/phospy/api`
- `src/phospy/datasets`
- `src/phospy/prediction`
- `src/phospy/signalomes/clustering`
- `src/phospy/workflows`

This scope is configured in `[tool.pyright]` in `pyproject.toml`.

Strict-mode enforcement is enabled where currently practical (a targeted subset of
non-pandas-heavy modules) while the full scope is enforced in standard mode.

## Exclusions (Initial Pass)

The initial gate excludes generated/cache/build artifacts (`__pycache__`, `.venv`,
`build`, `dist`, egg-info, etc.) and does not attempt whole-repository strict
typing in one pass.

Compatibility/migration boundaries may remain looser as long as they normalize
external payloads into typed internal models before entering core scientific paths.

## Run Locally

```bash
pip install -e ".[dev]"
pyright
```

If your local shell cannot resolve a Python interpreter for Pyright automatically
(for example on some Windows setups), run:

```bash
pyright --pythonpath .venv\Scripts\python.exe
```

## Ignore Policy

Use ignores only when a runtime invariant is already validated elsewhere and cannot
be represented cleanly to the type checker.

Rules:

- Prefer fixing types over suppressing diagnostics.
- Avoid file-level blanket ignores in checked core modules.
- Every ignore must be specific and justified inline.

Example:

```python
value = payload["field"]  # pyright: ignore[reportUnknownVariableType] - validated by manifest validator upstream
```

Avoid broad suppression like:

```python
value = payload["field"]  # type: ignore
```

## CI Gate

CI runs `pyright` as a required job in `.github/workflows/ci.yml`. Any Pyright
error fails the workflow.
