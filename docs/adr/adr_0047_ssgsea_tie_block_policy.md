# ADR-0047: ssGSEA Tie-Block Policy

## Status

- **ADR ID:** ADR-0047
- **Title:** ssGSEA Tie-Block Policy
- **Status:** Accepted
- **Date:** 2026-08-05
- **Decision Type:** Scientific Method Policy

## Context

The ssGSEA-style substrate enrichment activity method ranked finite
contrast/effect values and then walked an ordered hit/miss vector. Equal-valued
rows inherited their order from the input table through stable sorting. That made
activity scores and seeded empirical p-values depend on arbitrary row order when
a tie block crossed kinase substrate membership.

This responsibility belongs inside
`phospy.science.activities.methods.ssgsea_substrate_enrichment`. Dataset
construction, workflow interpretation, reference projection, and callers must
not pre-sort rows to make ssGSEA results deterministic.

## Decision

`ssgsea_substrate_enrichment_activity_v1` method policy version 2 uses a
method-owned tie-block rule named `midrank_block_expectation`.

For each profile, finite values are grouped into ordered equal-value blocks
according to the configured ranking direction. Untied blocks are equivalent to
the previous one-row rank walk. A tied block with `h` substrates and `m`
non-substrates contributes the expected rank-walk area over all within-block
orders:

`b * running_before + ((b + 1) / 2) * (h / n_substrates - m / n_non_substrates)`

where `b = h + m`. The walk then advances by the same block delta.

Blocks containing only substrates or only non-substrates therefore contribute
the ordinary uninterrupted hit or miss segment. Blocks containing both
substrates and non-substrates give all equal-valued sites equivalent treatment
without using row position or lexical site identity as a scientific tie break.

Seeded empirical substrate-label permutations use the same tie-block score
definition. Because ssGSEA permutation streams include method version in seed
material, version 2 may change permutation p-values and q-values even for
non-tied inputs, while non-tied deterministic enrichment scores remain
unchanged.

## Consequences

Positive consequences:

- Permuting rows with identical values does not change ssGSEA activity scores,
  p-values, q-values, or method diagnostics.
- Tied rows crossing substrate and non-substrate membership no longer create
  arbitrary label/order dependence.
- The policy is explicit in scientific policy provenance and per-run
  statistics-table diagnostics.

Negative consequences:

- Tied-input scores can differ from version 1 because the score definition is
  now the block expectation rather than one arbitrary within-tie order.
- Seeded permutation p/q values can differ from version 1 because method version
  is part of the deterministic child RNG seed material.

## Implementation Notes

- Method implementation:
  `src/phospy/science/activities/methods/ssgsea_substrate_enrichment.py`
- Policy provenance:
  `src/phospy/science/activities/scientific_policies.py`
- Method quantitative contract:
  `src/phospy/science/activities/method_contracts.py`

## Related Records

- [ADR-0017: Stochastic Reproducibility and Adaptive Prediction Seed Governance](adr_0017_stochastic_reproducibility_and_adaptive_prediction_seed_governance.md)
- [ADR-0034: Quantitative State, Motif Scoring Semantics, and Reference Context](adr_0034_quantitative_state_motif_semantics_and_reference_context.md)
