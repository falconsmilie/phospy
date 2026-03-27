# Parity Notes

This document explains what **parity** means in PhosPy and where the native kinase workflow fits into that claim.

For the bigger picture and the release gate, start with
[`docs/validation-and-parity.md`](validation-and-parity.md).

## What Parity Means

In this repository, parity means:

- Python outputs are compared against committed reference tables generated from R/PhosR
- those comparisons are automated in the parity-marked test suite
- the claim stays limited to the tested seam

So parity here is:

- fixture-backed
- seam-level
- narrower than full package equivalence

It does **not** mean that PhosPy is a complete behavioural, numerical, or feature-level replacement for PhosR.

## What the Parity Suite Covers

The current parity layer covers:

- deterministic preprocessing and matrix-building seams backed by small synthetic fixtures
- downstream kinase-analysis summaries backed by committed R-generated fixtures
- selected native workflow seams backed by committed L6 reference tables
- prediction-stage debugging through committed R and Python trace exports
- a curated fragile-support dataset used to widen evidence beyond the main L6 path

For fixture and trace directory details, see [`docs/fixtures.md`](fixtures.md).

## `KinaseWorkflow` and Parity

`KinaseWorkflow` is part of the supported 1.0.0 public API, but that does **not** turn the whole workflow into a broad
PhosR-equivalence claim.

The practical wording is:

- PhosPy provides a native Python workflow for profile construction, motif scoring, score combination, candidate
  selection, and adaptive SVM prediction.
- The repository includes fixture-backed validation for selected seams within that workflow.
- `svm_mode="r_parity"` is available when you want a closer comparison to the PhosR learner seam.
- The default `svm_mode="default"` is still the preferred Python-native mode.

## Running the Parity Suite

```bash
pytest -m parity
```

A few useful variations:

```bash
pytest -m parity -rs
pytest -m parity -vv
pytest -m parity --maxfail=1
pytest -m parity -k l6
```

## Optional Trace Regeneration

Most users do not need this. It is mainly useful when you are debugging a divergence between the committed R traces and
Python prediction behaviour.

Regenerate the committed R fixture sets:

```bash
Rscript scripts/generate_r_fixtures.R
Rscript scripts/generate_r_l6_fixtures.R
```

Regenerate the committed Python trace exports:

```bash
python scripts/export_python_prediction_traces.py --trace-kinases PRKAA1,MAPK1 --svm-mode r_parity --debug-top-n 10 --outdir tests/fixtures/python_reference_l6/prediction_trace
```

Replay the R sampling rows in Python:

```bash
python scripts/export_python_prediction_traces.py --trace-kinases PRKAA1,MAPK1 --svm-mode r_parity --sampling-trace-dir tests/fixtures/r_reference_l6/prediction_trace --outdir tests/fixtures/python_reference_l6/prediction_trace
```
