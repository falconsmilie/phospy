# Fixtures and Traces

This page is the quick map for the committed fixture and trace directories.

Most commands assume the repository root.

## Rebuild Everything

```bash
make fixtures-all
```

## Main Fixture Families

### `tests/fixtures/r_reference`

Small R-backed fixtures for deterministic preprocessing, site-matrix construction, and the downstream wrapper flow.

Generate with:

```bash
Rscript scripts/generate_r_fixtures.R
make fixtures-r-small
```

### `tests/fixtures/r_reference_l6`

The main committed L6 reference set for downstream kinase-analysis parity and selected native workflow seams.

Generate with:

```bash
Rscript scripts/generate_r_l6_fixtures.R
make fixtures-r-l6
```

### `tests/fixtures/fragile_support_reference`

A curated reference set used to widen seam coverage beyond the main L6 path.

Generate with:

```bash
python scripts/generate_fragile_support_reference.py
make fixtures-fragile
```

### `tests/fixtures/r_reference_l6_seam_stress`

A smaller L6-derived seam-stress fixture family.

Generate with:

```bash
python scripts/generate_l6_seam_stress_reference.py --outdir tests/fixtures/r_reference_l6_seam_stress
make fixtures-r-l6-seam-stress
```

### `tests/fixtures/synthetic_adaptive_sampling_edge`

A fully synthetic set for deterministic adaptive-sampling edge cases.

Generate with:

```bash
python scripts/generate_synthetic_adaptive_sampling_edge_fixtures.py --outdir tests/fixtures/synthetic_adaptive_sampling_edge
make fixtures-synthetic-edge
```

### `tests/fixtures/public_workflow_reference`

Small committed benchmark outputs for the public `PredMatWorkflow` and `SignalomeWorkflow` demo paths.

These fixtures are intentionally narrow. They protect the documented example-driven workflow contracts rather than every internal numeric seam.

Generate with:

```bash
PYTHONPATH=src python examples/predmat_workflow_demo.py
PYTHONPATH=src python examples/signalome_workflow_demo.py
```

Then promote the resulting workflow tables into `tests/fixtures/public_workflow_reference` when the benchmark contract changes intentionally.

## Trace Directories

### `tests/fixtures/r_reference_l6/prediction_trace`

Committed R prediction traces for seam-level debugging.

Generate with:

```bash
Rscript scripts/generate_r_l6_fixtures.R \
  --outdir tests/fixtures/r_reference_l6 \
  --trace_kinases PRKAA1,MAPK1 \
  --trace_top_n 10
make traces-r
```

### `tests/fixtures/python_reference_l6/prediction_trace`

Committed Python prediction traces used alongside the R traces.

Generate with:

```bash
python scripts/export_python_prediction_traces.py \
  --trace-kinases PRKAA1,MAPK1 \
  --svm-mode r_parity \
  --debug-top-n 10 \
  --outdir tests/fixtures/python_reference_l6/prediction_trace
make traces-python
```

Replay mode:

```bash
python scripts/export_python_prediction_traces.py \
  --trace-kinases PRKAA1,MAPK1 \
  --svm-mode r_parity \
  --sampling-trace-dir tests/fixtures/r_reference_l6/prediction_trace \
  --outdir tests/fixtures/python_reference_l6/prediction_trace
make traces-python-replay
```

## Temporary Debug Output

`tmp_trace_out` is scratch output for local debugging. It is not part of the committed parity contract unless you deliberately promote files into a fixture directory.

## Promotion Rule

If a temporary trace becomes part of the supported seam contract:

1. promote it into a committed fixture directory
2. update the related tests
3. update this page or [`parity.md`](parity.md)
