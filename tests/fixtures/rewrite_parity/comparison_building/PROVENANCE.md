# Comparison Building Fixture Provenance

- Source donor expectation:
  `tests_legacy/test_preprocessing.py::test_add_pairwise_comparisons_uses_schema_group_names`
- Fixture purpose:
  lock rewrite comparison-building output naming and values for a canonical
  one-pair legacy-style case (`p_sample_a_sample_b = 3.0`).
- Promotion date: 2026-04-21

## Rewrite-Owned Parity Inputs

This directory also contains rewrite-owned comparison-building input fixtures
for active parity-tier regression gates:

- explicit pair lane:
  - `legacy_pairwise_input_phospho.csv`
  - `legacy_pairwise_input_sample_metadata.csv`
  - `legacy_pairwise_expected.csv`
- inferred all-pairs lane:
  - `inferred_pairs_input_phospho.csv`
  - `inferred_pairs_input_sample_metadata.csv`
  - `inferred_pairs_expected.csv`

These fixtures pin supported
`comparisons.policy="sample_metadata_pairs"` behavior (pair membership,
identity/order, and output values) without depending on live legacy paths.
