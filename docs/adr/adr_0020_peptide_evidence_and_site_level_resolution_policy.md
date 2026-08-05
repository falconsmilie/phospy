# ADR: Peptide Evidence and Site-Level Resolution Policy

## Document Control

- **ADR ID:** ADR-0020
- **Title:** Peptide Evidence and Site-Level Resolution Policy
- **Status:** Accepted
- **Date:** 2026-05-12
- **Amended:** 2026-08-05
- **Decision Type:** Architecture Decision Record

## Context

PhosPy workflows downstream of dataset building consume site-level rows.
However, some upstream pipelines provide peptide-level evidence where one
observation can map to one or more phosphosites. Treating that ambiguity as
already resolved silently changes scientific meaning and prevents reproducible
audit of row-level decisions.

## Decision

Dataset building now supports two explicit input modes:

1. `site_level_resolved`
2. `peptide_evidence`

When `peptide_evidence` is used, `multi_site_policy` is required and must be
one of:

- `reject`
- `exclude_from_sequence_scoring`
- `keep_joint`
- `split`

Policy mapping reuses the existing `phospy.science.evidence` models and multi-site
resolution logic (`PeptideEvidenceTable`, `MultiSiteHandlingConfig`,
`SiteEvidenceMapping`) instead of introducing parallel ambiguity models.

## Amendment: Differential Uncertainty Scope

**Date:** 2026-07-29

ADR-0020 covers dataset-construction resolution from peptide evidence to
analysis-ready site-level intensity rows. It remains the preferred
PhosPy-origin lane for peptide-to-site differential analysis: resolve peptide
evidence at the sample-intensity level first, then fit the core site-level
differential model.

Post-hoc peptide-level differential estimate combination is governed separately
by [ADR-0041: Peptide-to-Site Differential Uncertainty Policy](adr_0041_peptide_to_site_differential_uncertainty_policy.md).
That route is withdrawn from public support and its compatibility shell fails
closed because coherent combined effect/inference semantics and executable
peptide-to-site mapping semantics are not implemented. ADR-0020's
sample-intensity resolution lane is the only current lane with an explicit
mapping-weight signal allocation model. Future public support for post-hoc
combination requires executable mapping semantics, a coherent combined
estimand, and an inferential result.

Scientific rationale:

- same-experiment peptide-level differential statistics are not independent
  external studies and cannot be treated as independent production
  meta-analysis input;
- finite-degree-of-freedom t statistics must not be used directly as z
  statistics;
- a single-estimate pass-through is not meta-analysis;
- multi-site allocation remains explicit provenance, not hidden statistical
  repair.

### Explicit Sample-Level Peptide-to-Site Estimand

The semantics in this section are dataset-construction signal-allocation and
site-resolution semantics. They are not post-hoc differential statistic
meta-analysis semantics.

Peptide-to-site aggregation is a typed, run-specific quantitative policy owned
by `phospy.science.evidence.dataset_resolution`. The dataset builder may
orchestrate this policy and record its payload, but it must not duplicate or
infer the peptide aggregation science. Downstream workflows consume the
resulting site-level matrix and provenance.

The supported production policy is
`peptide_to_site_linear_abundance_fractional_allocation_arithmetic_mean_v1`.
It defines these contract fields:

- **Supported input scales:** `linear`, `log2`
- **Supported input quantitative meanings:** `peptide_abundance` for `linear`,
  `peptide_log2_abundance` for `log2`
- **Output scales:** the declared input scale is preserved
- **Output quantitative meanings:** `phosphosite_abundance` for `linear`,
  `phosphosite_log_abundance` for unit-mapped `log2`
- **Mapping-weight source policy:**
  `explicit_mapping_weight_when_supplied_else_equal_fraction_per_resolved_site`
- **Mapping-weight normalisation policy:**
  `sum_to_one_per_peptide_evidence_row`
- **Mapping-weight semantics:**
  `unitless_fraction_of_one_peptide_evidence_row_allocated_to_each_resolved_site`
- **Signal allocation policy:** `multiply_peptide_signal_by_mapping_fraction`
- **Allocation domain:** `linear_abundance` when any non-unit fraction is
  present; `declared_scale_unit_mapping_passthrough` for `log2` input only when
  every mapping fraction is `1.0`
- **Site summarisation policy:** `arithmetic_mean_of_allocated_signals`
- **Missing-value policy:**
  `mean_finite_allocated_values_per_site_sample_preserve_missing_if_none_finite`
- **Duplicate evidence policy:**
  `retain_duplicate_peptide_evidence_rows_as_separate_observations`
- **Mixed ambiguity policy:**
  `combine_ambiguous_and_unambiguous_allocated_signals_in_site_mean`
- **Localisation summary policy:**
  `descriptive_mean_of_finite_reported_localisation_confidence_values`
- **Localisation summary semantics:**
  `descriptive_arithmetic_mean_not_calibrated_posterior_probability`
- **Signal conservation policy:**
  `not_signal_conserving_after_per_site_arithmetic_mean`
- **Legacy aggregation-policy alias:**
  `legacy_alias_for_arithmetic_mean_of_allocated_signals`

The mapping fraction column remains `site_mapping.mapping_weight`, exposed after
validation as `mapping_fraction`. Mapping weights must be positive finite
unitless values that sum to `1.0` for each `peptide_row_id`. When explicit
weights are absent, equal unitless fractions are derived across the resolved
sites for that peptide row.

#### Linear-abundance allocation formula

For a run with `input_intensity_scale="linear"` and
`input_quantitative_meaning="peptide_abundance"`, fractional allocation is
defined in linear abundance units:

```text
a[p,s,j] [linear abundance units]
  = w[p,s] [unitless peptide-row allocation fraction]
    * x[p,j] [linear peptide-abundance units]

y[s,j] [linear phosphosite-abundance estimate units]
  = arithmetic_mean(
      a[p,s,j] [linear abundance units]
      over finite retained peptide evidence rows p mapped to site s for sample j
    )
```

where:

- `p` is a retained peptide evidence row.
- `s` is a resolved phosphosite.
- `j` is a sample.
- `x[p,j]` is the peptide-row abundance in linear abundance units for sample
  `j`.
- `w[p,s]` is the unitless allocation fraction for peptide row `p` and site
  `s`.
- `a[p,s,j]` is the allocated peptide-row abundance in linear abundance units.
- `y[s,j]` is the site-level abundance estimate in linear phosphosite-abundance
  estimate units.

#### Log2 unit-mapping formula

For a run with `input_intensity_scale="log2"` and
`input_quantitative_meaning="peptide_log2_abundance"`, only unit mappings are
supported. In that case there is no fractional allocation and no inverse
transformation:

```text
a[p,s,j] [log2 peptide-abundance units]
  = x[p,j] [log2 peptide-abundance units]
    because w[p,s] = 1.0 [unitless peptide-row allocation fraction]

y[s,j] [log2 phosphosite-abundance estimate units]
  = arithmetic_mean(
      a[p,s,j] [log2 peptide-abundance units]
      over finite retained peptide evidence rows p mapped to site s for sample j
    )
```

Peptide evidence declared on a non-linear scale such as `log2` must fail closed
before allocation whenever a non-unit mapping fraction would be applied,
including derived equal split fractions and explicit `site_mapping.mapping_weight`
values. PhosPy does not assume `2**x` is an invertible recovery of linear
abundance when pseudocount, transformation direction, censoring, or
transformation provenance is unknown. A future scale-aware log-domain estimator
requires an ADR update, complete typed transformation provenance, numerical
round-trip tests, and scientific validation.

`w[p,s]` is an allocation fraction, not an inverse-variance weight,
localisation-confidence weight, posterior probability, or statistical evidence
weight. The final arithmetic mean is taken over allocated evidence-row signals.
This is not equivalent in general to `sum(w * x)` or
`sum(w * x) / sum(w)` when `x` is in linear abundance units and `w` is
unitless.

This policy does not generally conserve total signal after site-level
summarisation because the allocated signals are averaged per resolved site.
Duplicate retained rows affect the arithmetic mean because they are treated as
separate evidence observations. Ambiguous and unambiguous evidence rows mapped
to the same site are summarised by the same arithmetic mean over allocated
signals. Non-conservation is intentional for this retained estimator and is
recorded in provenance as
`not_signal_conserving_after_per_site_arithmetic_mean`.

Localisation output from peptide evidence is descriptive, not inferential. The
site metadata includes `localisation_confidence_descriptive_mean` and
`localisation_confidence_summary_semantics`. The legacy-compatible
`localisation_confidence` alias may be retained for existing localisation
threshold configuration, but provenance must state that the value is a
descriptive arithmetic mean of finite reported confidence values and not a
calibrated posterior localisation probability.

This ADR documents the currently supported policy. It does not claim that this
policy is universally optimal.

#### Evaluated estimator and rejected alternatives

The selected estimator is evaluated by:

- synthetic ground-truth fixtures covering unit, equal fractional, unequal
  explicit fractional, duplicate, missing-value, mixed ambiguous/unambiguous,
  contradictory-localisation, row-order-invariance, and log2 unit-mapping
  round-trip cases; and
- a checked-in realistic/reference regression fixture under
  `tests/fixtures/release_validation_regression/evidence_resolution`, with an
  independently generated expected site matrix.

Rejected alternatives for the current production contract:

- **Summation of allocated signals:** conserves allocated linear signal better
  than per-site means, but changes the row statistic with peptide-evidence row
  count and needs a separate abundance-estimand decision.
- **Conventional normalized weighted mean (`sum(w * x) / sum(w)`, with `x` in
  linear abundance units and `w` unitless):** has clearer weighted-estimator
  semantics but would treat mapping fractions as estimator weights rather than
  allocation fractions, conflicting with the current evidence-row contract.
- **Arithmetic mean of unallocated peptide signals:** ignores explicit
  fractional mapping and can overstate split ambiguous evidence.
- **Best peptide / highest localisation confidence:** discards usable
  quantitative evidence and would require a justified peptide-selection rule.
- **Median or robust estimators:** may reduce outlier sensitivity but require a
  defined minimum evidence count, uncertainty semantics, and validation.
- **Hierarchical or mixture models:** scientifically preferable for some
  ambiguous evidence settings, but no production model, priors, likelihood,
  diagnostics, or uncertainty calibration is implemented.
- **Log2 fractional allocation via `2**x`:** rejected because invertibility is
  not guaranteed without complete transformation provenance including
  pseudocount and prior transformations.

### Protein Identity Metadata

Peptide-evidence protein_accession is row-identity metadata. It must be
preserved as protein_accession or explicit protein_namespace/protein_identifier
metadata. It must not be rewritten into protein_id, which remains available for
legacy grouping-alias compatibility. New Signalome grouping metadata uses
protein_group_id.

If peptide-evidence rows collapse to one resolved site, the resolved site must
not aggregate multiple distinct accessions silently. Conflicting accessions for
one resolved site are rejected until a future ADR defines explicit aggregation
semantics.

## Peptide vs Site Ambiguity

- **Peptide-level ambiguity:** one peptide evidence row references multiple site
  tokens (for example `S10,T12`).
- **Site-level data contract:** downstream lanes require explicit site rows with
  explicit provenance of how ambiguity was resolved.

Ambiguity is not hidden by implicit defaults inside normalisation, kinase
scoring, or signalome scoring.

## Where Resolution Happens

Ambiguity resolution is owned by evidence/preprocessing boundary logic during
dataset interpretation, before `AnalysisReadyPhosphoDataset` construction.

`dataset.preprocessing_report` and dataset run provenance record the decision
and counts.

## Sequence Context Policy

Peptide-evidence resolution may normalise sequence strings by trimming
surrounding whitespace and uppercasing letters, but it must not repair
biological sequence context. In particular, the resolver must never rewrite the
centre residue of a provided `site_sequence` to match a resolved site token.

When multiple peptide-evidence rows collapse to one resolved site, supplied
`site_sequence` context is resolved as a set, not by row order. All non-null,
non-blank supplied values for that resolved site are normalized with the same
rules. They must normalize to exactly one unique valid value, and there must be
no invalid supplied non-blank values in the same group.

Conflicting valid contexts are rejected with `PhosPyInputError`. PhosPy must not
guess by choosing the first row, the most common value, the lexicographically
smallest value, or any other implicit precedence. Mixed evidence containing one
valid supplied value and one or more invalid supplied values is also rejected,
because reducing that group to the valid value would hide source-evidence
disagreement.

Split multi-site peptide-context derivation is a fallback for absent or
unusable split-site context when existing split-domain rules can derive the
resolved-site window from `peptide_sequence` plus `site_string`. It is not a
precedence rule over conflicting valid supplied contexts. When all supplied
contexts for a split target are invalid, derivation may succeed only if it
produces exactly one distinct normalized sequence; multiple derived candidates
or no derived candidate fail deterministically.

For unambiguous centred sequence windows, centre-residue mismatch is a hard
input error. Users must remove the peptide-evidence `site_sequence` value to
allow trusted reference derivation, or correct the upstream evidence. A future
repair path would require an explicit ADR and machine-readable provenance for
the repair decision.

## Provenance Requirements

The dataset report/provenance includes:

- peptide-to-site aggregation policy ID
- supported input scales and quantitative meanings
- input intensity scale and input quantitative meaning
- output intensity scale and output quantitative meaning
- allocation domain and whether fractional mappings were present
- peptide observations received
- mapped peptide observations
- site-mapping rows and allocated evidence rows
- unique analysis-ready `site_key` rows produced
- ambiguous observations
- unambiguous observations
- excluded observations
- split observations (when applicable)
- fractional mapping row count and unit mapping row count
- selected multi-site policy
- mapping-weight source policy and actual source observed for the run
- mapping-weight normalisation policy
- mapping-weight semantics
- signal allocation policy
- site summarisation policy
- missing-value policy
- legacy aggregation-policy alias + formula
- duplicate evidence policy + duplicate peptide row count
- mixed-ambiguity handling policy
- localisation aggregation policy
- localisation summary policy, output column, compatibility alias, and
  descriptive-not-posterior semantics
- signal-conservation policy
- uncertainty limitations

## Consequences

### Positive

- Input boundary is explicit about whether site-level resolution is already done.
- Multi-site handling is reproducible and audit-ready.
- Downstream workflows can trust dataset rows without owning peptide collapse.

### Tradeoffs

- Call sites that provide peptide evidence must now provide
  `peptide_evidence_sample_intensity_columns` and `multi_site_policy`.
- `keep_joint` rows preserve ambiguous site tokens and must be interpreted as
  ambiguous by consumers.
- Explicit mapping-weight contracts reject malformed mappings where per-peptide
  weights do not sum to `1.0`.

## Responsibility Audit

This ADR explicitly keeps ownership boundaries as:

- evidence resolution: `phospy.science.evidence` + dataset preprocessing/builder
- dataset model: no peptide collapse
- kinase workflow: no peptide ambiguity resolution
- signalome workflow: no peptide ambiguity resolution
- dataset builder: invokes evidence resolution and records provenance, but does
  not duplicate aggregation policy logic
- downstream workflows: consume resolved site-level data after provenance is
  attached
