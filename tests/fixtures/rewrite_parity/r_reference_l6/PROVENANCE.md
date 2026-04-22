# Rewrite Parity Fixture Provenance (`r_reference_l6`)

These files are rewrite-owned parity inputs and expectations used by active
tests in `tests/parity/` via `tests/support/rewrite_fixture_data.py`.

## Source

- Promoted from historical project snapshots.
- Promotion commits for the current rewrite-owned fixture path:
  - `d55b164`
  - `6e10739`

## Fixture ownership

- Rewrite parity maintainers own this fixture family as a blocking regression
  asset for activity-stage parity.

## Included Files

- `l6_phospho_matrix.csv`
- `native_profile_scores.csv`
- `predMat.csv`
- `kinase_activity_matrix.csv`
- `ksea_scores.csv`
- `ksea_counts.csv`
- `kinase_target_counts.csv`
- `kinase_target_table.csv`

## Stability Notes

- Inputs are fixed, committed fixtures under `tests/fixtures/rewrite_parity/`.
- Activity runtime parameters are fixed in active parity tests.
- CI enforces this lane through `activity-parity-gate`.
