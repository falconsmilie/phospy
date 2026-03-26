# Parity Notes

This document explains what the repository currently means by **parity** and where the native kinase workflow fits into
that claim.

For the broader validation model and the v1 release gate, start with
[`docs/validation-and-parity.md`](validation-and-parity.md).

## What Parity Means in This Repository

In PhosPy, parity means:

- Python outputs are compared against committed reference tables generated from R/PhosR
- those comparisons are automated in the parity-marked test suite
- the claim stays limited to the tested seam

Parity is therefore:

- fixture-backed
- seam-level
- narrower than full package equivalence

It does **not** mean that the repository as a whole is a complete behavioural, numerical, or feature-level replacement
for PhosR.

## What the Parity Suite Covers

The current parity layer covers:

- deterministic preprocessing and matrix-building seams backed by small synthetic fixtures
- downstream kinase-analysis summaries backed by R-generated fixtures
- selected native kinase workflow seams backed by committed L6 reference tables
- seam-level prediction debugging through committed R and Python trace exports
- a committed curated fragile-support dataset used to widen evidence beyond the main L6 path

For fixture and trace directory details, see [`docs/fixtures.md`](fixtures.md).

## `KinaseWorkflow` and Parity

`KinaseWorkflow` is part of the supported v1 public API.

That does **not** turn the whole workflow into a blanket PhosR-equivalence claim.

The current parity wording for `KinaseWorkflow` is:

- PhosPy provides a native Python workflow for profile construction, motif scoring, score combination, candidate
  selection, and adaptive SVM prediction.
- The repository includes fixture-backed validation for selected seams within that workflow.
- `svm_mode="r_parity"` is available when you want a closer comparison to the PhosR learner seam.
- The default `svm_mode="default"` remains the preferred Python-native mode and should not be described as numerically
  identical to PhosR across all datasets or settings.

That distinction matters. The package is strongest when it is explicit about which parts are validated against committed
references and which parts are deliberately native.

## Running the Parity Suite

```bash
pytest -m parity
```

Useful variations:

```bash
pytest -m parity -rs
pytest -m parity -vv
pytest -m parity --maxfail=1
pytest -m parity -k l6
pytest -m parity --collect-only -q
```

## Optional Diagnostic Output

The metrics-oriented parity tests can print extra diagnostic output through environment variables.

### Base Switch for Parity Metrics

```bash
PHOSPY_SHOW_PARITY=1 pytest -m parity -s
```

<details>
<summary>Windows (PowerShell)</summary>

```powershell
$env:PHOSPY_SHOW_PARITY=1
pytest -m parity -s
```

</details>

### Profile-Construction Metrics

```bash
PHOSPY_SHOW_PARITY=1 PHOSPY_SHOW_PROFILE_CONSTRUCTION=1 pytest -m parity -s
```

<details>
<summary>Windows (PowerShell)</summary>

```powershell
$env:PHOSPY_SHOW_PARITY=1
$env:PHOSPY_SHOW_PROFILE_CONSTRUCTION=1
pytest -m parity -s
```

</details>

### Prediction-Mode Comparison Metrics

```bash
PHOSPY_SHOW_PARITY=1 PHOSPY_SHOW_PREDICTION_MODE_COMPARISON=1 pytest -m parity -s
```

<details>
<summary>Windows (PowerShell)</summary>

```powershell
$env:PHOSPY_SHOW_PARITY=1
$env:PHOSPY_SHOW_PREDICTION_MODE_COMPARISON=1
pytest -m parity -s
```

</details>

### Replayed Prediction-Mode Comparison Metrics

```bash
PHOSPY_SHOW_PARITY=1 PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON=1 pytest -m parity -s
```

<details>
<summary>Windows (PowerShell)</summary>

```powershell
$env:PHOSPY_SHOW_PARITY=1
$env:PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON=1
pytest -m parity -s
```

</details>

## Regenerating the Python-Side Prediction Trace

This is useful for seam debugging and trace comparison, but it is **not required** for the parity suite itself.

```bash
python scripts/export_python_prediction_traces.py --trace-kinases PRKAA1,MAPK1 --svm-mode r_parity --debug-top-n 10 --outdir tests/fixtures/python_reference_l6/prediction_trace
```

If you want Python to replay the committed R sampling rows so the remaining delta is model-side only:

```bash
python scripts/export_python_prediction_traces.py --trace-kinases PRKAA1,MAPK1 --svm-mode r_parity --sampling-trace-dir tests/fixtures/r_reference_l6/prediction_trace --outdir tests/fixtures/python_reference_l6/prediction_trace
```

## Clean Regeneration Flow

```bash
rm -rf tests/fixtures/r_reference
rm -rf tests/fixtures/r_reference_l6
rm -rf tests/fixtures/python_reference_l6/prediction_trace

Rscript scripts/generate_r_fixtures.R
Rscript scripts/generate_r_l6_fixtures.R

python scripts/export_python_prediction_traces.py --trace-kinases PRKAA1,MAPK1 --svm-mode r_parity --debug-top-n 10 --outdir tests/fixtures/python_reference_l6/prediction_trace

pytest
```

<details>
<summary>Windows (PowerShell)</summary>

```powershell
Remove-Item -Recurse -Force tests\fixtures\r_reference -ErrorAction SilentlyContinue
Remove-Item -Recurse -Force tests\fixtures\r_reference_l6 -ErrorAction SilentlyContinue
Remove-Item -Recurse -Force tests\fixtures\python_reference_l6\prediction_trace -ErrorAction SilentlyContinue

Rscript scripts/generate_r_fixtures.R
Rscript scripts/generate_r_l6_fixtures.R

python scripts/export_python_prediction_traces.py --trace-kinases PRKAA1,MAPK1 --svm-mode r_parity --debug-top-n 10 --outdir tests/fixtures/python_reference_l6/prediction_trace

pytest
```

</details>
