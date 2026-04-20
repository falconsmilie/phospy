# Rewrite Parity Fixture Provenance (`fragile_support_reference`)

These files are rewrite-owned copies used by active parity tests for restored
prediction science seams (motif scoring, combined profile/motif weighting, and
candidate selection).

## Source

Promoted from:

`tests_legacy/fixtures/fragile_support_reference/`

on 2026-04-18.

Additional motif-lock fixtures were promoted/reconstructed on 2026-04-20 from:

- `tests_legacy/fixtures/r_reference_l6/native_motif_scores.csv`
- `tests_legacy/fixtures/r_reference_l6/l6_site_sequences.csv`

for selected kinases (`AKT1`, `IRAK1`, `LCK`, `MAPK1`, `PRKAA1`, `PRKAA2`).
`motif_frequency_matrices/*.csv` are rewrite-owned reconstructed matrices that
exactly reproduce the promoted full motif score table under the current motif
scoring kernel.

## Included Files

- `phospho_matrix.csv`
- `site_sequences.csv`
- `substrate_map.csv`
- `motif_sequences.csv`
- `profile_scores.csv`
- `motif_scores.csv`
- `motif_scores_full.csv`
- `motif_site_sequences_full.csv`
- `motif_frequency_matrices/*.csv`
- `profile_sizes.csv`
- `motif_sizes.csv`
- `combined_scores.csv`
- `combined_weights.csv`
- `candidate_substrates.csv`
