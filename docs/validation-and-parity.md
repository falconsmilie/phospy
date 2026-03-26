# Validation and Parity Guide

This document is the practical entry point for understanding what PhosPy validates, what its current parity claims
actually mean, and which checks form the v1 release gate.

## Validation Layers

PhosPy uses three complementary layers of evidence.

### 1. Core Python Tests

These tests do not depend on R.

They cover:

- schema and request validation
- preprocessing rules
- matrix construction
- public API behaviour
- native workflow components
- the documented example smoke path

Run them with:

```bash
pytest -m "not parity"
```

### 2. Fixture-Backed Parity Tests

These tests compare Python outputs against committed reference tables generated from R/PhosR.

They cover selected seams rather than the package as a whole.

Run them with:

```bash
pytest -m parity
```

### 3. Lint and Formatting Checks

These checks keep the release surface tidy and consistent.

Run them with:

```bash
pre-commit run --all-files
```

## v1 Release Gate

PhosPy v1 is ready to cut when all three checks are green from a clean checkout:

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest -m parity
```

That is the whole gate. Anything beyond that is useful, but not required for the basic v1 release contract.

## What Parity Means Here

Parity in PhosPy means:

- committed fixtures exist for the seam being discussed
- automated tests compare the Python result against those fixtures
- the corresponding documentation names the seam and keeps the claim narrow

Parity does **not** mean:

- the whole package is numerically identical to PhosR
- every PhosR option, corner case, or workflow branch is implemented
- every native Python path should match the R implementation exactly

## `KinaseWorkflow` Parity Wording

`KinaseWorkflow` is part of the supported v1 API, but its parity claim is deliberately narrower than the older
preprocessing path.

For v1, the correct wording is:

- `KinaseWorkflow` is a **native Python workflow** for profile construction, motif scoring, score combination,
  candidate selection, and adaptive SVM prediction.
- PhosPy includes **fixture-backed validation at selected seams** of that workflow.
- `svm_mode="r_parity"` exists for narrower learner-seam comparison against committed references.
- The default `svm_mode="default"` is the preferred Python-native path and is **not** a claim of package-wide
  numerical equivalence to PhosR.

That wording is intentionally modest. It is still strong enough for a real v1.

## Typical Test Commands

From the repo root:

```bash
python -m pip install --upgrade pip
pip install -e ".[test]"

pytest -m "not parity"
pytest -m parity
pytest
```

`pytest` runs the whole collected suite. The split commands are more useful when you want to keep the release gate easy to
reason about.

## Regenerating R Fixtures

You only need R when you want to regenerate or extend the committed fixtures.

```bash
Rscript scripts/generate_r_fixtures.R
Rscript scripts/generate_r_l6_fixtures.R
```

For fixture and trace directory details, see [`docs/fixtures.md`](fixtures.md).

## When to Update the Docs

When you change parity-backed behaviour, update at least one of these in the same line of work:

- the fixtures
- the parity tests
- this guide or [`docs/parity.md`](parity.md)

Do not silently broaden parity claims in the README or release notes without adding matching evidence.
