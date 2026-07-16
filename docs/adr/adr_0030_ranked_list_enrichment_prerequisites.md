# ADR-0030: Ranked-List Enrichment Prerequisites

## Status

- **ADR ID:** ADR-0030
- **Title:** Ranked-List Enrichment Prerequisites
- **Status:** Accepted
- **Date:** 2026-06-23
- **Decision Type:** Scientific Architecture and Roadmap

## Context

PhosPy currently supports native enrichment through offline
over-representation analysis (ORA). ORA consumes a selected foreground
identifier set, an explicit background universe, and caller-supplied local
gene-set or PTM-set collections. It tests overlap between selected identifiers
and each set under the explicit background universe.

Phosphoproteomics users may also want a future enrichment mode that uses a full
ranked list from differential or other scoring results. Such a mode can be
useful when there is no defensible hard threshold for selected identifiers, or
when modest but coordinated changes across a pathway or phosphosite set are
scientifically relevant. A ranked-list mode would ask whether members of a
caller-supplied set tend to concentrate toward one end of an ordered statistic
or otherwise show rank-associated structure.

This future mode would be a different statistical contract from ORA. It must
not be described as current support, as online resource integration, or as
GSEA, ssGSEA, PTM-SEA, or package-level parity.

## Problem Statement

A small option such as `method="ranked"` on the existing ORA request would be
unsafe before PhosPy defines the required input shape, identifier semantics,
ranking rules, universe semantics, tie handling, diagnostics, result model, and
documentation obligations.

Ranked-list enrichment also needs ORA diagnostics and set-size filtering to be
stable first. Those existing ORA features establish the identifier, universe,
filtered-set, multiple-testing, and diagnostic patterns that any future ranked
mode should preserve or deliberately extend.

## Decision

PhosPy defers implementation of ranked-list enrichment.

No ranked enrichment workflow, configuration option, result object, set
fetching layer, or numerical kernel should be added until ORA diagnostics and
set-size filters are stable and the prerequisites in this ADR are satisfied.

Future implementation, if added, must remain separate from current offline ORA
semantics. It may share enrichment collection models and multiple-testing
helpers only when the resulting behavior stays explicit, documented, and
method-specific.

This ADR is documentation and design guidance only. It does not implement
ranked enrichment.

## Difference From Current ORA

Current ORA:

- consumes selected identifiers as a foreground set;
- requires an explicit background universe;
- applies the configured ORA universe policy for foreground and set members,
  as specified by
  [ADR-0036](adr_0036_enrichment_universe_and_attrition_policy.md);
- tests each set using overlap counts; and
- reports ORA p-values, adjusted p-values, overlap counts, overlap
  identifiers, enrichment ratio, diagnostics, and provenance.

Future ranked-list enrichment would instead:

- consume an ordered table of eligible identifiers and one ranking statistic;
- treat the ranked table, after explicit filtering, as the tested universe;
- evaluate set members by their positions or statistics in that ranked
  universe;
- report a method-specific rank statistic or score, p-value if the selected
  method provides one, adjusted p-value when p-values are produced, ranking
  diagnostics, set-size diagnostics, and provenance; and
- avoid ORA-specific overlap ratios unless a future method explicitly defines a
  binary subset summary.

The two modes must not be interchangeable in documentation, request
validation, result interpretation, or provenance.

## Required Input Shape

A future ranked-list request must require a tabular input with at least:

- one identifier column;
- one numeric ranking column;
- an explicit `identifier_kind`;
- a homogeneous caller-supplied `GeneSetCollection`, `PtmSetCollection`, or
  `EnrichmentSetCollection`;
- a method-specific ranked enrichment configuration; and
- explicit ranking direction or a method-specific rule that defines whether
  large, small, positive, or negative values represent the top of the list.

Each eligible row must resolve to exactly one tested identifier in the same
namespace as the set collection. Duplicate identifier rows require an explicit
future policy before execution is allowed. Silent de-duplication is not
acceptable.

The ranked input should represent all identifiers eligible for the ranked test,
not only significant or manually selected identifiers. If callers want to test a
pre-filtered ranked universe, that filter must be explicit and recorded in
provenance.

## Ranking Column Requirements

The ranking column must be numeric, finite, and fully populated for all tested
rows after filtering. A future validator must reject missing, non-numeric,
non-finite, or otherwise unordered values unless a method-specific missing-rank
policy is designed, documented, and tested.

The ranking statistic must have documented interpretation. Examples could
include a differential statistic, signed effect size, or other caller-owned
score. PhosPy must not infer ranking direction from column names.

Ranking configuration must record:

- ranking column name;
- sort direction and signed-tail interpretation;
- whether both tails, one tail, or direction-specific tests are run;
- duplicate identifier policy;
- tie policy;
- filtering applied before ranking; and
- the final ranked universe size.

## Gene-Level and Phosphosite-Level Identifier Handling

Gene-level and phosphosite-level ranked enrichment must stay separate.

For gene-level analysis, ranked identifiers and set members must use a
gene-level namespace such as `gene_symbol` or `protein_id`. If the source table
is site-level, future support must require an explicit gene-level aggregation
contract before enrichment. The workflow must not automatically collapse
phosphosites to genes or choose one site per gene without a documented policy.

For phosphosite-level analysis, ranked identifiers and set members must use a
site-level namespace such as `site_key`, `display_id`, or `phosphosite`, with
the existing site-identity policies respected. A PTM set must not be
reinterpreted as a gene set, and a gene set must not be expanded into sites
inside ranked enrichment without a separate explicit mapping contract.

One ranked enrichment run must have one identifier level and one identifier
kind.

## Background and Universe Semantics

For ranked-list enrichment, the tested universe is the eligible ranked
identifier list after validation, optional explicit universe intersection,
duplicate handling, and ranking filters.

A future API may allow an explicit universe argument, but if it does, the
semantics must be different from ORA foreground/background semantics:

- the ranked table supplies the ordered tested identifiers;
- the explicit universe limits which ranked identifiers and set members are
  eligible;
- set sizes are measured after intersection with the final ranked universe;
- identifiers outside the final universe are reported in diagnostics; and
- the set collection alone must never define the universe.

PhosPy must not infer a ranked universe from a reference bundle, gene-set
collection, PTM-set collection, or online database.

## Set-Size Filtering

Future ranked-list enrichment must support set-size filtering before testing.

Set size must be measured after each set is intersected with the final ranked
universe. `min_set_size` and `max_set_size`, if supported, must be applied
before statistic calculation and before multiple-testing correction. Dropped
sets must be recorded with reasons and pre-filter or post-intersection sizes in
diagnostics.

The filtering contract should reuse the current ORA diagnostics pattern where
possible, but it must label ranked-universe sizes rather than ORA
background-overlap sizes.

## Tie Handling

Rank ties are expected in phosphoproteomics outputs, especially when scores are
rounded, censored, or produced by thresholded methods. Future implementation
must define a deterministic tie policy before any ranked test runs.

Acceptable future policies could include average ranks, dense ranks, a
documented secondary sort key, or method-specific grouped-tie handling. Stable
input row order is not enough by itself because it can make results depend on
file ordering.

Diagnostics must report tie counts, largest tie block size, and the selected
tie policy. If the method cannot handle ties safely, validation must reject the
request with a clear error.

## Multiple-Testing Correction

Multiple-testing correction must apply across the sets that are actually tested
after ranked-universe intersection and set-size filtering.

A future ranked method that returns finite p-values should use the same
documented correction vocabulary as current enrichment where suitable:
`benjamini_hochberg`, `bonferroni`, `holm`, `benjamini_yekutieli`, and `none`.
If a future method returns scores without p-values, the result must not invent
adjusted p-values.

The denominator, finite-p-value handling, and excluded-set behavior must be
visible in diagnostics and provenance.

## Result Object Shape

A future ranked enrichment result should not silently reuse ORA-only table
fields whose names imply overlap testing. The result object may share general
result-container conventions with current workflow results, but it must expose
ranked semantics clearly.

The ranked result table should include one row per tested term with fields such
as:

- `term_id` and optional term name/source metadata;
- `identifier_kind` and analysis level;
- ranked method name and ranking column;
- final ranked universe size;
- set size after ranked-universe intersection;
- method-specific enrichment statistic or score;
- p-value when the method produces one;
- adjusted p-value when p-values are adjusted;
- correction method;
- direction or tail tested;
- optional contributing identifiers, if the method defines them; and
- warnings or status flags for edge cases.

The result object should also provide:

- defensive table snapshots;
- typed records if the project keeps that pattern;
- ranking diagnostics;
- set-size filter diagnostics;
- universe diagnostics;
- method metadata;
- set collection summary; and
- provenance sufficient to reproduce ranking, filtering, testing, and
  correction decisions.

## Documentation Requirements

Before ranked-list enrichment can be described as supported, documentation must
be updated in the same change as the executable implementation.

Required documentation includes:

- an API page section that distinguishes ORA from ranked enrichment;
- request, configuration, result, diagnostics, and provenance field
  descriptions;
- examples for gene-level and phosphosite-level usage that do not mix
  identifier levels;
- ranking-column requirements and tie policy;
- background and ranked-universe semantics;
- set-size filtering behavior;
- multiple-testing correction behavior;
- limitations and non-claims;
- scientific coverage status; and
- workflow contract updates that name the owning public entrypoint.

Docs must continue to state that PhosPy does not bundle or fetch GO, KEGG,
Reactome, PTM-SEA, PTMsigDB, Enrichr, gseapy, clusterProfiler, or similar
resources for enrichment unless a future feature explicitly adds local,
license-reviewed resources with provenance. A ranked-list method must not be
documented as GSEA, ssGSEA, PTM-SEA, or package-level parity unless a separate
evidence and scope decision exists.

## Future Ownership

Likely future ownership areas are:

```text
src/phospy/contracts/configs/
src/phospy/science/enrichment/
src/phospy/workflows/enrichment/
src/phospy/provenance/
phospy.contracts.results (future package candidate)
phospy.validation.workflows.enrichment (future package candidate)
tests/
docs/api/enrichment.md
docs/workflow_contracts.md
docs/scientific-coverage.md
```

These paths are future-facing. This ADR does not require new modules today.

Validation modules should own input shape, identifier semantics, ranking
column, duplicate, tie, universe, and set-size rejections. Science-layer code
should own numerical ranked statistics. Workflow orchestration should preserve
the existing validator -> interpreter -> executor style where that pattern
applies. Result and provenance modules should own output contracts and
reproducibility metadata.

## Prerequisites Before Implementation

Implementation should not start until:

- ORA foreground/background diagnostics are stable;
- ORA set-size filtering behavior and diagnostics are stable;
- the future ranked method is selected and scoped by name;
- gene-level and site-level identifier policies are written as separate
  contracts;
- duplicate identifier and tie policies are decided;
- ranked-universe semantics are documented;
- result table and diagnostics fields are specified;
- provenance requirements are specified;
- tests are planned for validation, statistics, diagnostics, result shape,
  provenance, and docs; and
- scientific coverage language is ready to classify the feature without
  overclaiming.

## Consequences

Positive consequences:

- Current ORA behavior remains precise and understandable.
- Future ranked enrichment has a concrete checklist instead of an ambiguous
  method flag.
- Gene-level and site-level enrichment semantics remain separate.
- Documentation can acknowledge user interest without implying support.

Negative consequences:

- Users who need ranked-list pathway or set analysis must use external tools
  for now.
- Future support will require more than a statistic implementation because the
  request, validation, diagnostics, result, provenance, and documentation
  contracts must be added together.

Neutral consequences:

- Existing `EnrichmentWorkflow` remains offline ORA only.
- No online resource fetching or bundled curated enrichment resources are added
  by this decision.

## Non-Goals

This ADR has these non-goals:

- do not implement ranked enrichment;
- do not add public API flags for ranked enrichment;
- do not add GO, KEGG, Reactome, PTM-SEA, PTMsigDB, Enrichr, gseapy,
  clusterProfiler, or other resource fetching;
- do not add curated enrichment resources;
- do not mix gene-level and phosphosite-level semantics;
- do not make GSEA, ssGSEA, PTM-SEA, or package-level parity claims;
- do not change ORA foreground/background behavior;
- do not change current ORA result fields; and
- do not move enrichment ownership into differential, kinase, signalome, importer, or
  visualization workflows.

## Related Records

- [ADR-0005: Result Model Design for PhosPy](adr_0005__result_model_design.md)
- [ADR-0013: Scientific Parity Strategy and Parity-Testing Policy for PhosPy](adr_0013_scientific_parity_strategy_and_parity_testing_policy.md)
- [ADR-0022: Python API as the Supported Scientific Workflow Interface](adr_0022_python_api_as_supported_scientific_workflow_interface.md)
- [ADR-0024: Protein-Scoped Phosphosite Row Identity](adr_0024_protein_scoped_phosphosite_row_identity.md)
- [ADR-0025: Competitive Phosphoproteomics Workflow Coverage Roadmap](adr_0025_competitive_phosphoproteomics_workflow_coverage.md)
- [Workflow Contracts](../workflow_contracts.md)
- [Enrichment Workflow](../api/enrichment.md)
- [Scientific Coverage](../scientific-coverage.md)
