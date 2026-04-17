# Fixtures and Traces

This page is the quick map for committed fixture and trace directories.

For fixture intent, protected seams, metrics, and mode coverage, see [`parity.md`](parity.md). This page focuses on provenance, rebuild commands, and promotion rules.

Most commands assume the repository root.

## Rebuild Everything

```bash
make fixtures-all
```

## Main Fixture Families

### `tests/fixtures/r_reference`

Small R-backed fixtures for deterministic preprocessing, site-matrix construction, and downstream wrapper flow.

```bash
Rscript scripts/generate_r_fixtures.R
make fixtures-r-small
```

### `tests/fixtures/r_reference_l6`

Main committed L6 reference data for downstream kinase-analysis parity, native prediction seams, and committed R trace replay.

```bash
Rscript scripts/generate_r_l6_fixtures.R
make fixtures-r-l6
```

### `tests/fixtures/fragile_support_reference`

Curated L6-derived support-screening fixtures used to widen coverage around support boundaries and inclusion cut-offs.

```bash
python scripts/generate_fragile_support_reference.py
make fixtures-fragile
```

### `tests/fixtures/r_reference_l6_seam_stress`

Smaller L6-derived seam-stress fixtures.

```bash
python scripts/generate_l6_seam_stress_reference.py --outdir tests/fixtures/r_reference_l6_seam_stress
make fixtures-r-l6-seam-stress
```

### `tests/fixtures/synthetic_adaptive_sampling_edge`

Fully synthetic fixtures for deterministic adaptive-sampling edge cases.

```bash
python scripts/generate_synthetic_adaptive_sampling_edge_fixtures.py --outdir tests/fixtures/synthetic_adaptive_sampling_edge
make fixtures-synthetic-edge
```

### `tests/fixtures/public_workflow_reference`

Small committed benchmark outputs for the supported `SimpleKinaseWorkflow` and `SignalomeWorkflow` demo paths.

```bash
PYTHONPATH=src python examples/simple_workflow_demo.py
PYTHONPATH=src python examples/signalome_workflow_demo.py
```

Promote the resulting workflow tables into `tests/fixtures/public_workflow_reference` only when the public workflow contract changed intentionally.

### `tests/fixtures/public_workflow_reference/output_bundle`

Optional reproducibility bundle fixture directory for
`save_simple_kinase_workflow_bundle(...)` /
`load_simple_kinase_workflow_bundle(...)` round-trip checks.

Use this only when the output-bundle manifest contract or table inventory changes intentionally.

## Trace Directories

### `tests/fixtures/r_reference_l6/prediction_trace`

Committed R prediction traces used by the `r_reference_l6` fixture family.

```bash
Rscript scripts/generate_r_l6_fixtures.R \
  --outdir tests/fixtures/r_reference_l6 \
  --trace_kinases PRKAA1,MAPK1 \
  --trace_top_n 10
make traces-r
```

### `tests/fixtures/python_reference_l6/prediction_trace`

Committed Python prediction traces used alongside the R traces for export and replay debugging.

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
3. update [`parity.md`](parity.md) if the contract changes
4. update this page if the rebuild or promotion workflow changes
