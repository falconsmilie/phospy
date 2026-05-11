# ADR: Phosphosite Identity and Localisation Policy

## Document Control

- **ADR ID:** ADR-0018
- **Title:** Phosphosite Identity and Localisation Policy
- **Status:** Accepted
- **Date:** 2026-05-11
- **Decision Type:** Architecture Decision Record

## Context

Phosphoproteomics tables frequently provide `gene_symbol` and a site token (for
example, `S123`). Those fields are necessary, but not always sufficient for
biological interpretation at the site level. Ambiguous protein mapping,
site-token inconsistencies, and uncertain phosphosite localisation can produce
apparently precise rows that carry different confidence levels.

PhosPy already enforces an analysis-ready boundary with required
`site_sequence`. This ADR extends that boundary by making identity assumptions
and localisation confidence explicit and policy-visible.

## Decision

PhosPy adopts a two-layer phosphosite policy:

1. **Analysis-ready dataset identity contract** in dataset/table validation.
2. **Workflow-specific strictness policy** in workflow validators (not executors).

### Analysis-Ready Metadata Classes

At `AnalysisReadyPhosphoDataset` boundary:

- **Required:** `gene_symbol`, `site`, `site_sequence`
- **Strongly recommended (optional):** `protein_id`, `residue`,
  `site_position` (or legacy `position`), `localisation_probability`

### Identity Rules at Analysis-Ready Boundary

For `dataset.site_metadata`:

- `site` must be parseable as `<residue><position>` by default.
- `residue` and `site_position`/`position`, when provided, must agree with
  parsed `site`.
- `site_sequence`, when sequence-format checks are possible, must have a
  central phosphorylatable residue (`S`, `T`, or `Y`) compatible with site
  residue metadata.
- Duplicate rows remain governed by existing duplicate-resolution policy;
  identity validation does not introduce silent row dropping.

### Localisation Representation

`localisation_probability` is modeled as optional numeric metadata:

- valid values: `0.0 <= x <= 1.0`
- unknown values: missing (`NA`/blank)
- invalid values (for example `-0.1`, `1.2`, `"high"`, `"unknown"`) are rejected

Unknown localisation is not interpreted as high confidence.

### Workflow Policy Object

Workflow validation uses explicit `LocalisationRequirement`:

- `allow_unknown` (default)
- `require_present`
- `require_threshold` (via `minimum_probability`)

This allows statements like:

- dataset is analysis-ready, localisation unknown
- workflow requires localisation present
- workflow requires localisation >= threshold

### Workflow Responsibilities

- Workflow validators enforce workflow-required metadata (`protein_id`,
  localisation policy).
- Executors do not repair identity metadata late.
- Signalome/kinase execution consumes already validated assumptions.

## Consequences

### Positive

- Preserves permissive ingestion while strengthening analysis-ready integrity.
- Makes localisation uncertainty explicit instead of implicit.
- Keeps biologically meaningful strictness where it belongs: validator boundary.
- Produces row-context diagnostics for missing/invalid localisation metadata.

### Tradeoffs

- Some legacy placeholder sequences/site tokens that are not biologically
  interpretable now fail earlier at validation boundaries.
- Workflows that need strict localisation/protein identity must opt into policy
  and handle validation failures explicitly.

## Scope Boundaries

This ADR does not:

- force every raw input file to provide all recommended metadata columns
- add organism-specific hard-coded phosphosite rules to generic DataFrame helpers
- introduce implicit filtering based on localisation thresholds

If filtering is later introduced, it must be explicit and provenance-recorded.

## References

Macek, B., Mann, M., & Olsen, J. V. (2009). Global and site-specific quantitative phosphoproteomics: Principles and applications. *Annual Review of Pharmacology and Toxicology, 49*, 199-221.

Olsen, J. V., Blagoev, B., Gnad, F., Macek, B., Kumar, C., Mortensen, P., & Mann, M. (2006). Global, in vivo, and site-specific phosphorylation dynamics in signaling networks. *Cell, 127*(3), 635-648.

Sharma, K., D’Souza, R. C. J., Tyanova, S., Schaab, C., Wiśniewski, J. R., Cox, J., & Mann, M. (2014). Ultradeep human phosphoproteome reveals a distinct regulatory nature of Tyr and Ser/Thr-based signaling. *Cell Reports, 8*(5), 1583-1594.
