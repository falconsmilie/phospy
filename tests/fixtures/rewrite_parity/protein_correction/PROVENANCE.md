# Rewrite Parity Fixture Provenance (`protein_correction`)

This directory contains promoted reference material for
`total_protein_correction.policy="subtract_log_total"` parity checks in the
supported builder preprocessing lane.

## Source Material

The baseline source material was promoted from historical project history and
is now committed directly in this folder.

`reference_corrected_matrix.csv` captures the corrected phospho matrix slice
(`phospho_corrected_1..6`) mapped to rewrite sample columns (`p_group1..6`)
and canonical `site_id` formatting (`GENE;SITE;`), using:

`log2(phospho + 1.0) - log2(total + 1.0)`

## Rewrite-Owned Parity Inputs

The active preprocessing parity gate uses rewrite-owned input fixtures in this
directory:

- `reference_input_phospho.csv`
- `reference_input_site_metadata.csv`
- `reference_input_total.csv`

These fixtures encode the supported
`intensity_transform.policy="log2"` +
`total_protein_correction.policy="subtract_log_total"` lane directly, so
ordinary parity execution does not depend on archived trees.
