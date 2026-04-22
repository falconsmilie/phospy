# Comparison Building Fixture Provenance

- Fixture purpose:
  lock comparison-building output naming and values for a canonical one-pair
  reference case (`p_sample_a_sample_b = 3.0`).
- Promotion date: 2026-04-21.

## Rewrite-Owned Parity Inputs

This directory contains rewrite-owned comparison-building input fixtures for
active parity-tier regression gates:

- explicit pair lane:
  - `reference_pairwise_input_phospho.csv`
  - `reference_pairwise_input_sample_metadata.csv`
  - `reference_pairwise_expected.csv`
- inferred all-pairs lane:
  - `inferred_pairs_input_phospho.csv`
  - `inferred_pairs_input_sample_metadata.csv`
  - `inferred_pairs_expected.csv`

These fixtures pin supported `comparisons.policy="sample_metadata_pairs"`
behavior (pair membership, identity/order, and output values) without runtime
dependence on archived paths.
