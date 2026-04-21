# Rewrite Parity Fixture Provenance (`protein_correction`)

This directory contains rewrite-owned promoted donor material for legacy
total/protein correction parity checks in the supported builder preprocessing
lane.

## Source Donor Material

- `tests_legacy/fixtures/r_reference/df_phospho_corrected.csv`
- `tests_legacy/fixtures/r_reference/df_phospho_filtered.csv`
- `tests_legacy/fixtures/r_reference/df_total_filtered.csv`

`legacy_r_reference_corrected_matrix.csv` captures the donor-corrected phospho
matrix slice (`phospho_corrected_1..6`) mapped to rewrite sample columns
(`p_group1..6`) and canonical `site_id` formatting (`GENE;SITE;`).

## Rewrite-Owned Parity Inputs

The active preprocessing parity gate uses rewrite-owned input fixtures in this
directory:

- `legacy_r_reference_input_phospho.csv`
- `legacy_r_reference_input_site_metadata.csv`
- `legacy_r_reference_input_total.csv`

These fixtures encode the supported
`total_protein_correction.policy="ratio_to_total"` lane directly, so ordinary
parity execution does not depend on runtime reads from `tests_legacy/`.
