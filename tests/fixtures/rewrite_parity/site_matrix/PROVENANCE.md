# Site Matrix Rewrite Parity Fixture Provenance

- Promoted into rewrite-owned fixture lane on 2026-04-21.
- Source material was promoted from historical project history and committed in
  this folder for active parity use.
- Purpose:
  - lock site-matrix construction parity for supported
    `site_matrix.policy="build_from_metadata"` behavior
    (`drop_any_missing` + `max_mean_signal`).
- Active rewrite tests consume this folder directly.
- Active parity-tier gate:
  - `tests/parity/test_preprocessing_science_parity.py::test_site_matrix_build_from_metadata_matches_rewrite_reference_fixture`

## Fixture Files

- `reference_phospho_corrected.csv`
- `reference_expected_matrix.csv`
- `reference_expected_input.csv`
