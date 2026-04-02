# Parity

This is the detailed parity guide for PhosPy.

Read [`docs/validation-and-parity.md`](validation-and-parity.md) first if you want the short version of validation
rules and parity scope. Read this page when you want to run parity tests, turn on extra metrics, or interpret the
results.

## What Parity Means Here

In this repository, parity means all three of these are true:

- a committed fixture exists for the seam being discussed
- an automated parity test compares Python output with that fixture
- the claim stays limited to that seam

So parity here is:

- **seam-level**, not package-wide
- **fixture-backed**, not anecdotal
- **specific**, not a blanket claim that all PhosR behaviour is reproduced

Parity does **not** mean:

- the whole package is numerically identical to PhosR
- every PhosR branch, option, or edge case is implemented
- every native Python path should match the R implementation exactly

## What Is Covered Today

The current parity layer covers:

- core preprocessing outputs backed by committed R-generated fixtures
- downstream kinase-analysis summaries backed by committed R-generated fixtures
- selected native workflow seams backed by committed L6 reference tables
- prediction-stage debugging through committed R and Python trace exports
- a curated fragile-support fixture that widens evidence beyond the main L6 path
- smaller seam-stress and synthetic adaptive-sampling fixtures used for narrower replay-style checks

For the fixture and trace directory layout, see [`docs/fixtures.md`](fixtures.md).

## `KinaseWorkflow` and `svm_mode`

`KinaseWorkflow` is part of the supported public API, but it remains a native Python workflow.

The practical rule is simple:

- use `svm_mode="default"` for the normal Python-native path
- use `svm_mode="r_parity"` when you want a closer learner-seam comparison against the committed parity fixtures

Configuration example:

```python
from phospy import KinaseWorkflow

native = KinaseWorkflow(svm_mode="default")
comparison = KinaseWorkflow(svm_mode="r_parity")
```

What `svm_mode="r_parity"` does **not** mean is “the whole workflow now matches PhosR”. It narrows one comparison
seam. It does not widen the package claim.

## Running the Parity Suite

From the repository root:

```bash
pytest -m parity
```

Useful variations:

```bash
pytest -m parity -rs
pytest -m parity -vv
pytest -m parity --maxfail=1
pytest -m parity -k l6
```

Repository shortcuts:

```bash
make test-parity
make test-seams
```

- `make test-parity` runs `pytest -m parity -s` with `PHOSPY_SHOW_PARITY=1` and
  `PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON=1`
- `make test-seams` runs the seam-focused parity files: `tests/test_reference_workflow_seams.py` and
  `tests/test_fragile_support_reference.py`

## Optional Parity Metrics Output

`tests/test_parity-with_metrics.py` can print extra comparison summaries while you investigate a seam. To see those
summaries in the terminal, run pytest with `-s` or `--capture=no`.

Available environment variables:

- `PHOSPY_SHOW_PARITY`: master switch for parity metrics output
- `PHOSPY_SHOW_PROFILE_CONSTRUCTION`: adds the optional profile-construction summary
- `PHOSPY_SHOW_PREDICTION_MODE_COMPARISON`: adds default-versus-`r_parity` prediction comparison metrics
- `PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON`: adds replayed prediction comparison metrics

The more specific flags do nothing unless `PHOSPY_SHOW_PARITY` is also enabled. Truthy values are case-insensitive and
include `1`, `true`, `yes`, and `on`.

If you enable all four flags and run the full parity suite, PhosPy prints every available metrics block reached by
those tests. If you narrow the run with `-k`, you only see the summaries for the matching tests.

### Linux or macOS examples

```bash
PHOSPY_SHOW_PARITY=1 pytest -m parity -s
PHOSPY_SHOW_PARITY=1 PHOSPY_SHOW_PROFILE_CONSTRUCTION=1 pytest -m parity -k l6 -s
PHOSPY_SHOW_PARITY=1 PHOSPY_SHOW_PREDICTION_MODE_COMPARISON=1 pytest -m parity -k comparison -s
PHOSPY_SHOW_PARITY=1 PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON=1 pytest -m parity -k replayed -s
PHOSPY_SHOW_PARITY=1 PHOSPY_SHOW_PROFILE_CONSTRUCTION=1 PHOSPY_SHOW_PREDICTION_MODE_COMPARISON=1 PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON=1 pytest -m parity -s
```

### Windows PowerShell examples

```powershell
$env:PHOSPY_SHOW_PARITY = "1"; pytest -m parity -s
$env:PHOSPY_SHOW_PARITY = "1"; $env:PHOSPY_SHOW_PROFILE_CONSTRUCTION = "1"; pytest -m parity -k l6 -s
$env:PHOSPY_SHOW_PARITY = "1"; $env:PHOSPY_SHOW_PREDICTION_MODE_COMPARISON = "1"; pytest -m parity -k comparison -s
$env:PHOSPY_SHOW_PARITY = "1"; $env:PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON = "1"; pytest -m parity -k replayed -s
$env:PHOSPY_SHOW_PARITY = "1"; $env:PHOSPY_SHOW_PROFILE_CONSTRUCTION = "1"; $env:PHOSPY_SHOW_PREDICTION_MODE_COMPARISON = "1"; $env:PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON = "1"; pytest -m parity -s
```

### Windows Command Prompt examples

```bat
set PHOSPY_SHOW_PARITY=1 && pytest -m parity -s
set PHOSPY_SHOW_PARITY=1 && set PHOSPY_SHOW_PROFILE_CONSTRUCTION=1 && pytest -m parity -k l6 -s
set PHOSPY_SHOW_PARITY=1 && set PHOSPY_SHOW_PREDICTION_MODE_COMPARISON=1 && pytest -m parity -k comparison -s
set PHOSPY_SHOW_PARITY=1 && set PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON=1 && pytest -m parity -k replayed -s
set PHOSPY_SHOW_PARITY=1 && set PHOSPY_SHOW_PROFILE_CONSTRUCTION=1 && set PHOSPY_SHOW_PREDICTION_MODE_COMPARISON=1 && set PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON=1 && pytest -m parity -s
```

## How to Read the Metrics

These summaries help you interpret the parity claim. They do not widen it.

- **Profile-construction, profile-scoring, and combined-score metrics** show how closely numeric tables match the
  committed reference outputs.
- **Prediction parity metrics** focus on ranking agreement and overlap because the prediction seam is better judged by
  ranked output than by strict element-for-element equality.
- **Prediction mode comparison metrics** compare the default learner path with `svm_mode="r_parity"` on the same
  bundled fixtures.
- **Replayed prediction trace metrics** go one level deeper by checking how closely Python follows committed trace data
  for selected kinases.

In the bundled fixtures, the deterministic profile-building seams are effectively numerically identical to the reference
outputs. Prediction-stage parity is intentionally described with rank agreement, overlap, and trace-replay statistics.

## Example Output From the Bundled Fixtures

The block below is reference output from the bundled repository fixtures. It is useful context when you are reading the
parity docs or debugging a seam, but it is not a promise that every dataset or every platform will produce identical
numbers.

```text
tests/test_parity-with_metrics.py ...
Optional profile-construction parity metrics:
  kinases compared: 44
  profile matrix shape: (44, 12) vs (44, 12)
  mean per-kinase Pearson correlation: 100.00%
  mean absolute difference: 4.25993e-16
  max absolute difference: 6.43929e-15
.
Profile-scoring parity metrics:
  sites compared: 589
  kinases compared: 44
  mean per-kinase Pearson correlation: 100.00%
  mean per-kinase Spearman correlation: 100.00%
  mean absolute difference: 2.85064e-16
  max absolute difference: 1.44329e-15
.
Combined-score parity metrics:
  sites compared: 589
  kinases compared: 28
  mean per-kinase Pearson correlation: 100.00%
  mean per-kinase Spearman correlation: 100.00%
  mean absolute difference: 2.81253e-16
  max absolute difference: 9.99201e-16
  mean weight absolute difference: 1.39695e-15
..
Prediction parity metrics:
  svm_mode: default
  kinases compared: 28
  mean Spearman rank agreement: 96.41%
  mean top-10 overlap: 83.21%
  mean top-20 overlap: 88.75%
  mean top-30 overlap: 89.64%
  kinases with top-10 overlap >= 70%: 26/28
.
Prediction parity mode comparison:
  default mean Spearman rank agreement: 96.41%
  r_parity mean Spearman rank agreement: 96.46%
  default mean top-10 overlap: 83.21%
  r_parity mean top-10 overlap: 80.71%
  default mean top-20 overlap: 88.75%
  r_parity mean top-20 overlap: 87.32%
  default mean top-30 overlap: 89.64%
  r_parity mean top-30 overlap: 89.05%
.
Replayed prediction trace parity metrics:
  svm_mode: r_parity
  kinases compared: 2
  trace kinases used: MAPK1, PRKAA1
  trace kinases skipped: none
  initial negative exact matches: 600/600
  iteration sample exact matches: 6000/6000
  iteration prob class-1 Pearson correlation: 99.952%
  iteration prob class-2 Pearson correlation: 99.952%
  iteration decision class-1 Pearson correlation: 100.000%
  iteration decision class-1 mean absolute difference: 2.14513e-15
  iteration prob class-1 mean absolute difference: 0.00675895
  final top-site matches: 198/200
  final top class-1 mean absolute difference: 0.00325916
.
Replayed prediction parity mode comparison:
  default trace kinases used: MAPK1, PRKAA1
  r_parity trace kinases used: MAPK1, PRKAA1
  trace kinases skipped: none
  default iteration prob class-1 Pearson correlation: 99.953%
  r_parity iteration prob class-1 Pearson correlation: 99.952%
  default iteration prob mean absolute difference: 0.00670807
  r_parity iteration prob mean absolute difference: 0.00675895
  default final top-site matches: 183/200
  r_parity final top-site matches: 198/200
  default final top class-1 mean absolute difference: 0.00326616
  r_parity final top class-1 mean absolute difference: 0.00325916
```

## Trace Replay and Debugging

Most users will not need this. It is useful when you are trying to answer a narrower question: is the delta coming
from candidate sampling, or from the later learner path?

A common pattern is:

1. compare against the committed R trace tables
2. export Python traces for the same kinases
3. optionally replay the R sampling rows in Python to isolate the remaining difference

## Optional Trace Regeneration

You only need R when you want to regenerate or extend the committed fixture sets.

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
