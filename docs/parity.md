# Parity Model

This document covers the fixture-backed test suite: what parity means in this repository, how to run it, and which
options change the pytest output you see.

Unless otherwise noted, commands below assume:

- **Linux**
- **repo root**
- a shell that understands standard `bash` syntax

macOS uses the same commands unless a section says otherwise. Windows is only shown where the syntax changes.

## What Parity Means Here

In this repository, parity means that a specific Python path has been compared against outputs generated from the R
package and shown to agree within the limits of the committed fixture-backed tests for that path.

Parity here is:

- explicit
- test-backed
- limited to named workflow seams
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

## Running the Fixture-Backed Test Suite

### 1) Install the Python Test Dependencies

```bash
python -m pip install --upgrade pip
pip install -e ".[test]"
```

### 2) Generate the Small Synthetic R Fixtures

```bash
Rscript scripts/generate_r_fixtures.R
```

This writes into:

```text
tests/fixtures/r_reference
```

### 3) Generate the L6 R Fixtures, Including the R Prediction Trace

```bash
Rscript scripts/generate_r_l6_fixtures.R
```

This writes into:

```text
tests/fixtures/r_reference_l6
tests/fixtures/r_reference_l6/prediction_trace
```

### 4) Run the Python Tests That Do Not Need Parity Fixtures

```bash
pytest -m "not parity"
```

### 5) Run the Parity Tests

```bash
pytest -m parity
```

### 6) Run the Whole Python Test Suite

```bash
pytest
```

## One Clean Command Sequence

```bash
python -m pip install --upgrade pip
pip install -e ".[test]"

Rscript scripts/generate_r_fixtures.R
Rscript scripts/generate_r_l6_fixtures.R

pytest -m "not parity"
pytest -m parity
pytest
```

## Useful Pytest Options

The parity suite does **not** currently define custom pytest command-line flags. The output is controlled by:

- standard pytest options
- a small set of environment variables used by `tests/test_parity-with_metrics.py`

### Show Skip Reasons

```bash
pytest -m parity -rs
```

Effect on output:

- prints skip reasons in the summary
- useful when fixture files are missing
- helpful because the parity tests call `pytest.skip(...)` with the missing filenames

### Show Verbose Test Names

```bash
pytest -m parity -vv
```

Effect on output:

- prints each selected parity test name
- makes it easier to see which seam failed
- useful when you want a more detailed test log

### Quiet Output

```bash
pytest -m parity -q
```

Effect on output:

- reduces pytest noise
- useful when you only care whether the parity suite passed

### Stop on the First Failure

```bash
pytest -m parity --maxfail=1
```

Effect on output:

- stops the run after the first failing parity test
- useful when debugging one broken seam at a time

### Restrict Test Selection

```bash
pytest -m parity -k l6
```

Effect on output:

- only runs parity tests whose names match `l6`
- useful when you want to focus on the L6 fixture path

### Show Collection Without Running Tests

```bash
pytest -m parity --collect-only -q
```

Effect on output:

- prints which parity tests would run
- does not execute the tests

## Optional Diagnostic Output

The metrics-oriented parity tests can print extra diagnostic output, but that output is controlled by environment
variables rather than pytest flags.

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

Effect on output:

- enables printed metric summaries from `tests/test_parity-with_metrics.py`
- `-s` is important here because it disables stdout capture and lets the printed metrics appear during a passing run

Typical output may include:

- mean per-kinase Pearson or Spearman correlation
- mean and maximum absolute differences
- ranked-overlap summaries for prediction outputs

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

Effect on output:

- adds optional profile-construction metrics
- only works when `PHOSPY_SHOW_PARITY=1` is also set

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

Effect on output:

- adds optional comparison output between prediction modes such as `default` and `r_parity`
- only works when `PHOSPY_SHOW_PARITY=1` is also set

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

Effect on output:

- adds optional comparison output for the replayed prediction-trace path
- only works when `PHOSPY_SHOW_PARITY=1` is also set

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

If you want to regenerate everything from scratch first:

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

## R-Side Packages

The R scripts in this repository expect at least:

- `PhosR`
- `SummarizedExperiment`
- `e1071`

The small synthetic fixture script also uses:

- `readr`
- `dplyr`
- `tidyr`
- `tibble`
- `janitor`

## Maintenance Rule

Treat the parity tests as the executable definition of the repository’s current parity contract.

When a parity-backed workflow changes, update at least one of these in the same line of work:

- the fixtures
- the tests
- the documented scope in this file or [`docs/fixtures.md`](fixtures.md)

Do not silently broaden parity claims in the README or other project documentation beyond the fixture-backed seams
described here without adding corresponding evidence.