# Rewrite Parity Fixture Provenance (`r_reference_l6`)

These files are the rewrite-owned parity input/expectation subset used by active
tests in `tests/parity/` and helpers in `tests/support/rewrite_fixture_data.py`.

## Source

Promoted from:

`tests_legacy/fixtures/r_reference_l6/`

on 2026-04-18.

Fixture ownership: rewrite parity maintainers (`tests/parity/` + science-port
review owners) are responsible for this directory as a blocking regression
asset.

## Why This Directory Exists

- Active rewrite parity tests should resolve fixtures from rewrite-owned paths.
- Legacy fixture layout remains archived for provenance and historical reference.

## Included Files (Active Slice)

- `l6_phospho_matrix.csv`
- `native_profile_scores.csv`
- `predMat.csv`
- `kinase_activity_matrix.csv`
- `ksea_scores.csv`
- `ksea_counts.csv`
- `kinase_target_counts.csv`
- `kinase_target_table.csv`

`kinase_target_table.csv` is generated from the promoted `predMat.csv` fixture
using the archived legacy target-table kernel
(`legacy_archive/phospy_legacy/activities/scoring.py::build_kinase_target_table`)
with `threshold=0.6`.

## Promotion Policy (Hard Gate)

- Activity fixture files above are locked as a regression baseline for
  `tests/parity/test_activity_stage_parity.py`.
- Any behavior-changing update must be intentional: regenerate fixtures from the
  approved source, update this provenance note, and land together with parity
  test updates.
- CI keeps an explicit activity parity gate; removing or bypassing this fixture
  family is treated as a regression.
