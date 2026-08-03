# ADR-0045: Enrichment Derived Identifier-Set Provenance

## Status

- **ADR ID:** ADR-0045
- **Title:** Enrichment Derived Identifier-Set Provenance
- **Status:** Accepted
- **Date:** 2026-08-03
- **Decision Type:** Scientific Workflow Contract

## Context

Offline enrichment consumes a foreground identifier set, a background universe,
and caller-supplied enrichment collections. External lists are often curated
outside PhosPy and can be legitimately described only by source type, label, and
identifier count.

Identifier sets derived from PhosPy quantitative results have a different
reproducibility burden. A publishing-grade ORA input needs to identify the exact
source result, the profile or contrast, the threshold rule, missing-value
handling, identifier namespace, quantitative scale and meaning, and producing
software version. Without those facts, the ORA result cannot be traced back to a
specific upstream quantitative result and selection rule.

## Decision

External/manual identifier lists remain accepted without provenance. When
provenance is supplied for manual or raw lists, it remains minimal and typed.

Identifier lists declared as `phospy_derived_quantitative` must carry typed
`EnrichmentDerivedQuantitativeSetProvenance` plus typed
`InputIntensityScaleEvidence`. The derived provenance records:

- source result fingerprint;
- source result kind, currently profile or contrast;
- source profile or contrast name;
- identifier namespace;
- threshold and direction;
- missing-value rule;
- quantitative scale and meaning; and
- producing software version.

The enrichment validator owns consistency checks over these claims. It verifies
that the derived identifier namespace matches the enrichment request
`identifier_kind`. When `input_table` is the selected identifier source, it also
checks the declared source-result fingerprint against the supplied table.

The enrichment workflow does not own differential threshold calculation. It does
not infer selected identifiers from quantitative columns, re-run differential
filters, infer thresholds from current data, or convert identifiers between
namespaces. Selection remains caller/upstream-workflow owned; enrichment records
and validates the resulting typed provenance.

`EnrichmentConfig.publishing()` provides explicit configurable set-size bounds
for publication-oriented ORA runs. The preset uses `min_set_size=5` and
`max_set_size=500` by default, with caller overrides for study-specific
collections.

## Consequences

Positive consequences:

- PhosPy-derived ORA inputs become reproducible from source result fingerprint
  and explicit selection metadata.
- Manual and arbitrary external lists are not overburdened with provenance
  fields they cannot truthfully provide.
- Enrichment continues to reject namespace mismatches and fingerprint mismatches
  without becoming a thresholding owner.
- Publishing-oriented set-size bounds are opt-in and auditable.

Negative consequences:

- Existing callers that mark sets as `phospy_derived_quantitative` must now
  provide a typed derivation record.

Neutral consequences:

- ORA statistics, background-universe policy, and multiple-testing correction
  behavior are unchanged.
- No online enrichment resources, GSEA, ssGSEA, PTM-SEA, or hidden identifier
  conversion are introduced.

## Implementation Notes

The typed provenance contract lives in
`src/phospy/contracts/enrichment_identifier_sets.py`.
Workflow validation lives in `src/phospy/validation/workflows/enrichment.py`.
Result run provenance serializes the nested identifier-set provenance through
`src/phospy/workflows/enrichment/provenance.py`.

The publishing config preset lives in
`src/phospy/contracts/configs/enrichment.py` and reuses the existing set-size
filtering diagnostics and provenance path.

## Related Records

- [ADR-0030: Ranked-List Enrichment Prerequisites](adr_0030_ranked_list_enrichment_prerequisites.md)
- [ADR-0035: Provenance Immutability and Stable Serialization](adr_0035_provenance_immutability_and_stable_serialization.md)
- [ADR-0036: Enrichment Universe and Attrition Policy](adr_0036_enrichment_universe_and_attrition_policy.md)
- [Enrichment Workflow](../api/enrichment.md)
- [Workflow Contracts](../workflow_contracts.md)
