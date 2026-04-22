# Fixtures

This page is the quick map for committed fixture directories used by the
supported public lanes.

> Audience: contributors and maintainers working on tests, parity, and provenance.
> If you are using PhosPy as an end user, this page is optional.

For parity intent and protected seams, see [`parity.md`](parity.md).

## Canonical Active Roots

- Rewrite-owned fixture families: `tests/fixtures/rewrite_parity/*`
- Public workflow regression family: `tests/fixtures/public_workflow_reference`

## Maintainer Regeneration Defaults

The default output/input roots used by `Makefile` and fixture scripts are:

- `scripts/generate_r_fixtures.R` -> `tests/fixtures/rewrite_parity/r_reference`
- `scripts/generate_r_l6_fixtures.R` -> `tests/fixtures/rewrite_parity/r_reference_l6`
- `scripts/export_python_prediction_traces.py` -> `tests/fixtures/rewrite_parity/python_reference_l6/prediction_trace`
- `scripts/generate_fragile_support_reference.py` -> `tests/fixtures/rewrite_parity/fragile_support_reference`
- `scripts/generate_l6_seam_stress_reference.py` -> `tests/fixtures/rewrite_parity/r_reference_l6_seam_stress`
- `scripts/diff_prediction_trace_probabilities.py` compares `tests/fixtures/rewrite_parity/r_reference_l6/prediction_trace` (R trace) against `tests/fixtures/rewrite_parity/python_reference_l6/prediction_trace` (Python trace)
- `scripts/generate_signalome_public_workflow_reference.py` -> `tests/fixtures/public_workflow_reference`

## Active Public Fixture Families

### `tests/fixtures/rewrite_parity/r_reference_l6`

Current L6 parity inputs/expectations used by active tests in
`tests/parity/`.

Included files currently cover:

- activity-stage fixture checks (`predMat`, `ksea_*`, `kinase_*`)
- kinase profile scoring checkpoints (`native_profile_scores`)
- shared L6 phospho input (`l6_phospho_matrix`)

Activity parity tests in `tests/parity/test_activity_stage_parity.py` consume
these committed files via `tests/support/rewrite_fixture_data.py` and do not run
live legacy activity code.

Provenance is documented in:

- `tests/fixtures/rewrite_parity/r_reference_l6/PROVENANCE.md`

### `tests/fixtures/public_workflow_reference`

Committed workflow regression expectations for the public
`SignalomeWorkflow` lane:

- `signalome_rewrite_l6_module_assignments.csv`
- `signalome_rewrite_l6_modules.csv`
- `signalome_rewrite_l6_network_nodes.csv`
- `signalome_rewrite_l6_network_edges.csv`
- `signalome_rewrite_l6_expanded_signalome.csv`
- `signalome_rewrite_l6_contract.json`

Supported regeneration path:

```bash
make fixtures-public-workflow-reference
```

The generator (`scripts/generate_signalome_public_workflow_reference.py`) uses a
maintainer-owned helper in `scripts/support/public_workflow_reference.py`
instead of test-only support modules.

## Historical Archive

Legacy fixture trees remain in `tests_legacy/fixtures/` for provenance and
historical traceability. Active tests should resolve from
`tests/fixtures/` as their normal source.

The following top-level roots are historical/legacy carryovers and are not used
as default maintainer output locations:

- `tests/fixtures/r_reference`
- `tests/fixtures/r_reference_l6`
- `tests/fixtures/fragile_support_reference`
- `tests/fixtures/r_reference_l6_seam_stress`
- `tests/fixtures/python_reference_l6`

## Public Example Regeneration

The public examples remain runnable via:

```bash
PYTHONPATH=src python examples/dataset_builder_demo.py
PYTHONPATH=src python examples/kinase_workflow_demo.py
PYTHONPATH=src python examples/signalome_workflow_demo.py
```

## Where Next

- Governance policy around fixture use: [Parity to PhosR](parity.md)
- Maintainer navigation hub: [Contributor and maintainer docs](contributor/index.md)
