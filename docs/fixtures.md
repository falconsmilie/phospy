# Fixtures

This page is the quick map for committed fixture directories used by the
supported public lanes.

> Audience: contributors and maintainers working on tests, parity, and provenance.
> If you are using PhosPy as an end user, this page is optional.

For parity intent and protected seams, see [`parity.md`](parity.md).

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

## Historical Archive

Legacy fixture trees remain in `tests_legacy/fixtures/` for provenance and
historical traceability. Active tests should resolve from
`tests/fixtures/` as their normal source.

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
