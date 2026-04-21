# Site Matrix Rewrite Parity Fixture Provenance

- Promoted into rewrite-owned fixture lane on 2026-04-21.
- Source donor material:
  - `tests_legacy/fixtures/r_reference/df_phospho_corrected.csv`
  - `tests_legacy/fixtures/r_reference/mat_phospho_corrected.csv`
  - `tests_legacy/fixtures/r_reference/phosr_input.csv`
- Purpose:
  - lock site-matrix construction parity for supported
    `site_matrix.policy="build_from_metadata"` default legacy-equivalent policy
    surface (`drop_any_missing` + `max_mean_signal`).
- Active rewrite tests should consume this rewrite fixture folder rather than
  reading from `tests_legacy/fixtures/` directly.
