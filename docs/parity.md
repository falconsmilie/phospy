# Parity to PhosR

PhosPy parity is intentionally narrow and fixture-backed. The rewrite does not
claim full package equivalence with PhosR.

## What Parity Means Here

Parity in this repository is:

- seam-level
- selective
- tied to committed fixtures

Parity here does not mean:

- every PhosR feature is implemented
- every Python path must numerically match PhosR

## Active Parity Coverage

The parity suite currently protects three rewrite-era slices:

- activity-stage outputs from fixed `predMat` + phospho inputs
- selected kinase-scoring/prediction points on the supported L6 lane
- selected signalome regression contracts on the supported L6 lane:
  `module_assignments`, `signalome_modules`, `kinase_network.nodes`,
  `kinase_network.edges`

## Fixture Locations

### Rewrite-owned parity inputs and expectations

- `tests/fixtures/rewrite_parity/r_reference_l6/`
- provenance and promotion history:
  `tests/fixtures/rewrite_parity/r_reference_l6/PROVENANCE.md`

These files are the normal source for active parity tests in `tests/parity/`
and helpers in `tests/support/rewrite_fixture_data.py`.

### Rewrite workflow regression expectations

- `tests/fixtures/public_workflow_reference/signalome_rewrite_l6_*.csv`
- `tests/fixtures/public_workflow_reference/signalome_rewrite_l6_contract.json`

### Historical reference archive

- `tests_legacy/fixtures/` is retained for provenance and archival material.
- Active rewrite parity tests should not resolve fixtures from this path as their
  normal source.

## Run the Parity Suite

```bash
pytest -m parity
```

Useful variants:

```bash
pytest -m parity -rs
pytest -m parity -vv
pytest tests/parity/test_signalome_workflow_parity.py -vv
```
