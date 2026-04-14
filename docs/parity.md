# Parity to PhosR

PhosPy is inspired by `PhosR`, but the parity claim is intentionally narrow.

## What Parity Means Here

Parity in this repository means:

- fixture-backed
- seam-level
- selective

It does not mean:

- full package equivalence with `PhosR`
- every `PhosR` feature is implemented
- every native Python path should match `PhosR` numerically

## Fixture Families at a Glance

Trace directories such as `prediction_trace/` belong to the fixture family listed below. They are not separate fixture families.

- `tests/fixtures/r_reference`: small R-generated fixtures protecting preprocessing, site-matrix construction, and the downstream wrapper flow
- `tests/fixtures/r_reference_l6`: main L6 reference set protecting downstream kinase-analysis outputs, native prediction seams, ranking agreement, and replay against committed R sampling traces
- `tests/fixtures/fragile_support_reference`: curated support-screening fixtures protecting boundary conditions around support and inclusion rules
- `tests/fixtures/r_reference_l6_seam_stress`: smaller seam-stress fixtures protecting combined-score and replay boundary behaviour
- `tests/fixtures/synthetic_adaptive_sampling_edge`: synthetic fixtures protecting deterministic adaptive-sampling edge cases
- `tests/fixtures/public_workflow_reference`: committed benchmark outputs for the public `SimpleKinaseWorkflow` and `SignalomeWorkflow` demos

For fixture rebuild commands, see [`fixtures.md`](fixtures.md).

## `svm_mode`

PhosPy exposes two public prediction presets:

- `svm_mode="default"` for the recommended stable native path
- `svm_mode="r_parity"` for the supported parity-oriented learner, sampling, and final-scoring preset

Using `svm_mode="r_parity"` does not make the full workflow equivalent to `PhosR`. It is the narrower preset used when parity-sensitive prediction checks need the closest supported reference-oriented path.

The explicit support decision is recorded in [ADR 0002](adr/0002-r-parity-public-preset.md).

## What Is Covered

The current parity layer covers selected seams, including:

- core preprocessing outputs
- downstream kinase-analysis outputs
- selected native workflow seams
- selected prediction trace and replay checks
- end-to-end benchmark fixtures for the documented `SimpleKinaseWorkflow` and `SignalomeWorkflow` demos

## Release Thresholds by Mode

The release bar is deliberately different for `default` and `r_parity` because the presets have different jobs.

### `default`

`default` is acceptable for release when all of the following still hold:

- the public workflow benchmarks still match `tests/fixtures/public_workflow_reference`
- the L6 ranking floor in `tests/test_parity-with_metrics.py` still holds:
  - mean Spearman rank agreement `>= 0.96`
  - mean top-20 overlap `>= 0.85`
  - mean top-30 overlap `>= 0.88`
  - kinases with top-10 overlap of at least 70%: `>= 20`
- non-parity tests and parity tests both pass

`default` is not required to meet the replay-trace bar used to justify `r_parity`.

### `r_parity`

`r_parity` is acceptable for release only when all of the following still hold:

- the public workflow benchmarks still match `tests/fixtures/public_workflow_reference`
- on the L6 ranking benchmark, it remains at least as strong as `default` for:
  - mean Spearman rank agreement
  - mean top-10 overlap
  - mean top-20 overlap
  - mean top-30 overlap
- on the same ranking benchmark, mean top-10 overlap remains `>= 0.82`
- on the replayed L6 sampling-trace benchmark, it continues to meet the protected replay floor:
  - initial negative rows: exact match
  - iteration sample rows: exact match
  - iteration decision class-1 Pearson correlation `>= 0.999999`
  - iteration decision mean absolute difference `<= 1e-12`
  - iteration probability class-1 Pearson correlation `>= 0.998`
  - iteration probability mean absolute difference `<= 0.015`
  - final top-site matches: exact match

## Release Review Checklist

When a change touches prediction policy, sampling, scoring, fixture generation, or the public workflow examples, release review should check:

1. `pytest -m "not parity"`
2. `pytest -m parity`
3. `pytest tests/test_readme_smoke.py tests/test_end_to_end_parity.py`
4. regenerate affected fixtures if the change is intentional and explain the contract change in the pull request
5. regenerate and review the mode-comparison benchmark when prediction policy or sampling changes:

   ```bash
   python benchmarks/compare_prediction_modes.py --repeats 1
   ```

6. confirm the public docs still describe `default` as the recommended native mode and `r_parity` as the supported parity-oriented mode

## Run the Parity Suite

```bash
pytest -m parity
```

Useful variants:

```bash
pytest -m parity -rs
pytest -m parity -vv
pytest -m parity -k l6
```

Make targets:

```bash
make test-parity
make test-seams
```

## Benchmark the Public Prediction Modes

Run the reproducible mode-comparison harness from the repository root:

```bash
python benchmarks/compare_prediction_modes.py
```

By default this writes two review artifacts under `benchmarks/reports/latest/`:

- `compare_prediction_modes.json`
- `compare_prediction_modes.md`

The harness uses:

- `tests/fixtures/r_reference_l6` for ranking parity and replayed sampling-trace fidelity
- `tests/fixtures/public_workflow_reference` for the documented `SimpleKinaseWorkflow` and `SignalomeWorkflow` demo outputs

Useful variants:

```bash
python benchmarks/compare_prediction_modes.py --repeats 1
python benchmarks/compare_prediction_modes.py --stdout-only
```

## Optional Debug Output

Some parity tests can print extra summaries when environment flags are enabled.

Available flags:

- `PHOSPY_SHOW_PARITY`
- `PHOSPY_SHOW_PROFILE_CONSTRUCTION`
- `PHOSPY_SHOW_PREDICTION_MODE_COMPARISON`
- `PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON`
