# Validation and Parity

This is the quickest guide to what PhosPy validates and what its current parity claim actually means.

## The Three Checks That Matter

PhosPy uses three practical layers of evidence.

### 1. Non-Parity Tests

These tests do not depend on R. They cover core package behaviour, validation rules, preprocessing, native workflow
components, and the documented example smoke path.

```bash
pytest -m "not parity"
```

### 2. Parity Tests

These tests compare Python outputs against committed R/PhosR-generated reference tables.

The claim is deliberately narrow: parity is asserted only for the tested seam.

```bash
pytest -m parity
```

### 3. Lint and Formatting

These checks keep the release surface tidy and consistent.

```bash
pre-commit run --all-files
```

## Release Gate for 1.0.0

From a clean checkout, the 1.0.0 release gate is simply:

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest -m parity
```

That is the whole gate.

## What Parity Means Here

In PhosPy, parity means:

- committed fixtures exist for the seam being discussed
- automated tests compare Python results against those fixtures
- the claim stays limited to that seam

Parity does **not** mean:

- the whole package is numerically identical to PhosR
- every PhosR option or workflow branch is implemented
- every native Python path should match the R implementation exactly

## `KinaseWorkflow` and Parity

`KinaseWorkflow` is part of the supported 1.0.0 API, but it is still a **native Python workflow**.

PhosPy includes fixture-backed validation at selected seams within that workflow. `svm_mode="r_parity"` exists for a
narrower learner-seam comparison against committed references. The default `svm_mode="default"` remains the preferred
Python-native path and is **not** a claim of package-wide numerical equivalence to PhosR.

## Typical Commands

From the repo root:

```bash
python -m pip install --upgrade pip
pip install -e ".[test]"

pytest -m "not parity"
pytest -m parity
pytest
```

Use the split commands when you want the release gate to stay easy to reason about.

## Regenerating R Fixtures

You only need R when you want to regenerate or extend the committed fixtures.

```bash
Rscript scripts/generate_r_fixtures.R
Rscript scripts/generate_r_l6_fixtures.R
```

For the fixture and trace directory layout, see [`docs/fixtures.md`](fixtures.md).
