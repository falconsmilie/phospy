# Rewrite Parity Fixture Provenance (`r_reference_l6`)

These files are rewrite-owned parity inputs/expectations used by active tests in
`tests/parity/` via loaders in `tests/support/rewrite_fixture_data.py`.

## Source Donor Material

Promoted from:

- `tests_legacy/fixtures/r_reference_l6/` (archival donor fixtures)

Donor snapshot lineage in this repository:

- `fe3b07f` (legacy fixture tree introduction)
- `ced835e` (`legacy_archive/phospy_legacy` archival snapshot)

Promotion commits for rewrite-owned fixture path:

- `d55b164` (activity science support landing)
- `6e10739` (activity parity hard-gate alignment)

Fixture ownership: rewrite parity maintainers (`tests/parity/` + science-port
review owners) are responsible for this directory as a blocking regression
asset.

## Included Files (Active Slice)

- `l6_phospho_matrix.csv`
- `native_profile_scores.csv`
- `predMat.csv`
- `kinase_activity_matrix.csv`
- `ksea_scores.csv`
- `ksea_counts.csv`
- `kinase_target_counts.csv`
- `kinase_target_table.csv`

## Generation Procedure

Activity expected outputs were materialized from the donor fixture inputs
(`predMat.csv` + `l6_phospho_matrix.csv`) using the archived legacy activity
kernels and request validator through:

- `scripts/archive/generate_activity_donor_snapshot.py`

Active parity tests do not run that donor code path. The script is archival
fixture-regeneration tooling only.

## Environment Baseline (Fixture Regeneration / Verification)

- Python `3.14.0` (`.venv/Scripts/python.exe`)
- `pandas==3.0.2`
- `numpy==2.4.4`

## Normalization Applied

- Index normalization:
  `kinase_activity_matrix.csv` and `ksea_scores.csv` are normalized to index
  name `kinase`.
- Series normalization:
  `ksea_counts.csv` uses series name `n_substrates`;
  `kinase_target_counts.csv` uses series name `n_targets`.
- Target-table normalization:
  `kinase_target_table.csv` is committed as explicit columns
  (`site_id`, `kinase`, `score`) for deterministic equality checks.

## Why Stable for Parity

- Inputs are fixed, committed fixtures under `tests/fixtures/rewrite_parity/`.
- Activity runtime parameters are fixed in active parity tests
  (`threshold=0.6`, `min_substrates=3`, `top_n_substrates=20`).
- Rewrite parity assertions compare full tables/series against these committed
  expectations.
- CI enforces this lane through `activity-parity-gate`
  (`pytest tests/parity/test_activity_stage_parity.py -m "parity and activity_parity"`).
