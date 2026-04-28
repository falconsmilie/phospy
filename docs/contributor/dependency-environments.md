# Dependency Environments and Reproducibility

PhosPy uses three dependency layers on purpose:

1. **Minimum supported versions**: defined in `pyproject.toml` with lower bounds.
2. **Tested versions**: the Python and scientific stack versions validated in CI.
3. **CI/release constrained versions**: pinned in `constraints/ci.txt` for deterministic CI and release verification.

## Current Tested Matrix

Last updated: `2026-04-28`

| Component | Versions used in CI |
| --- | --- |
| Python | 3.10, 3.11, 3.12 |
| NumPy | 2.2.6 |
| pandas | 2.2.3 |
| SciPy | 1.15.2 |
| scikit-learn | 1.6.1 |

Tooling versions used in CI/release jobs are also pinned in `constraints/ci.txt`.

## Why This Split Exists

- Lower bounds in package metadata keep installation flexible for users.
- Pinned CI/release constraints keep scientific checks reproducible and audit-friendly.
- Scientific parity/regression outputs are interpreted against the constrained set, not against whichever newest packages happen to resolve on a given day.

## Reproduce the CI Environment Locally

Create an environment and install with the same constraints CI uses:

```bash
python -m venv .venv-ci
source .venv-ci/bin/activate
python -m pip install --upgrade pip
python -m pip install -c constraints/ci.txt -e ".[dev,test]"
```

Windows PowerShell:

```powershell
python -m venv .venv-ci
.venv-ci\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -c constraints/ci.txt -e ".[dev,test]"
```

Run the same scientific regression gates:

```bash
pytest tests/parity/test_activity_stage_parity.py -m "parity and activity_parity" -s
pytest tests/parity -m parity -s
```

## Update Policy

When bumping the tested scientific stack:

1. Update `constraints/ci.txt`.
2. Update the version table on this page.
3. Run CI (including parity gates) under the new constraints.
4. Include the dependency-set update in release notes/changelog entries for traceability.
