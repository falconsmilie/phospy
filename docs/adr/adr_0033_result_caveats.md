# ADR-0033: Result Caveats and Scientific Warning Ownership

## Status

- **ADR ID:** ADR-0033
- **Title:** Result Caveats and Scientific Warning Ownership
- **Status:** Accepted
- **Date:** 2026-07-06
- **Decision Type:** Scientific Architecture and Result Contract

## Context

PhosPy workflows increasingly expose scientific limitations, assumptions,
policy overrides, and attrition summaries alongside numerical results. These
facts affect interpretation, but they are not all validation failures. Some are
warnings or scope statements that must travel with a successful result.

Logs are not a sufficient contract for these scientific warnings because users
can miss them, saved bundles may not include them, and automated consumers need
stable machine-readable identifiers. Provenance is also not the right place to
hide user-facing caveats: provenance records how a result was produced, while
caveats summarize interpretation risks and limitations that a user should see
without traversing full provenance payloads.

The project therefore needs a clear ownership boundary for where scientific
caveat facts are generated, where user-facing caveats are assembled, and where
they are exposed.

## Decision

Validators and interpreters produce structured facts. Result assemblers expose
user-facing caveats.

Validators, resolved validators, interpreters, and numerical stages may produce
structured diagnostics, eligibility facts, policy outcomes, counters, and
provenance inputs. They should not become the primary owners of final
user-facing caveat prose.

Workflow result assemblers and workflow-owned caveat builders own the final
translation from structured facts into `ResultCaveat` values. Public result
objects expose those caveats as top-level result fields named `caveats`.

Each caveat must have:

- a stable machine-readable `code`;
- a `severity` from the shared result-caveat vocabulary;
- a concise user-facing `message`; and
- compact structured `details` containing only the context needed to interpret
  the caveat.

Severe caveats must not be log-only. If a successful public result has a
scientifically important warning or error-level limitation, that condition must
be available through the result's top-level `caveats` field. When a condition
invalidates execution, validators should still reject the request instead of
returning a caveated result.

## Ownership Boundary

Validation modules own request, dataset, configuration, eligibility, and
scientific invariant checks at their boundary. They raise errors for
non-executable inputs and may emit structured diagnostic facts for executable
but notable conditions.

Interpreters own resolved scientific facts needed for execution, such as
resolved references, effective policy choices, design summaries, eligibility
counters, score-source summaries, and interpreted identifier semantics.

Executors own numerical outputs and method-local diagnostics. They should return
structured metrics or diagnostics when those facts affect interpretation.

Result assemblers own user-facing caveat construction. They consume structured
facts from validators, interpreters, executors, and provenance and attach the
resulting caveats to public result objects.

Public result models own storage and validation of the top-level `caveats`
field. They should validate that caveats use the shared `ResultCaveat` contract,
but should not recalculate scientific caveat facts.

## Caveats and Provenance

Caveats do not replace provenance.

Provenance remains the durable record of how a result was produced: inputs,
configuration, policy choices, resolved resources, hashes, environment details,
stage order, and reproducibility metadata. Caveats are compact interpretation
signals for result consumers.

Caveat `details` may reference provenance concepts, counts, selected policy
names, or caveat-relevant identifiers, but must not duplicate full provenance
payloads. If a user needs full reproducibility context, they should inspect the
result provenance.

## Stable Codes

Caveat codes are part of the machine-readable result contract.

Codes should be stable, descriptive, workflow-scoped where appropriate, and
specific enough for tests and downstream tooling to match. Renaming or removing
a code should be treated as a public result-contract change unless the caveat
was never exposed on a supported public result.

Messages may be clarified over time, but consumers should not need to parse
messages to identify caveat meaning. Machine-readable meaning belongs in
`code`, `severity`, and structured `details`.

## Examples

Differential analysis caveats include direct trusted dataset construction,
declared-scale override use, imputation-aware feature withholding, withheld
features in result tables, and the narrow supported fixed-effect parity
envelope. Differential validators and interpreters own the design, contrast,
scale, imputation, and eligibility facts. The differential caveat builder owns
the public `DifferentialAnalysisResult.caveats` entries.

Kinase workflow caveats include attrition policy warnings, permissive
localisation policy, non-default reference sources, automatic reference
resolution, score-source fallback, and scoring scope that is not
PhosR-equivalent.
Kinase validators and interpreters own reference compatibility, resolved
eligibility, and score-source facts. Kinase result assembly owns the final
top-level caveats on `KinaseWorkflowResult`.

Signalome workflow caveats include inherited upstream kinase attrition,
prediction/reference limitations, permissive localisation policy,
`protein_id` grouping assumptions, the descriptive-only meaning of signalome
network correlations, and dropped fully missing clustering dimensions when
present. Signalome result assembly may reference upstream kinase caveat codes
and structured attrition facts, but must not copy full upstream provenance into
signalome caveat details. Dropped-dimension caveats should use bounded previews
for user-facing messages while preserving full deterministic label lists in
structured provenance when payload size remains reasonable.

Enrichment workflow caveats include offline ORA-only scope, caller-supplied
background-universe assumptions, declared identifier-kind assumptions, and the
absence of rank-based GSEA, ssGSEA, PTM-SEA, ranking, leading-edge, or
enrichment-score semantics. Enrichment interpretation owns the identifier,
background, set-collection, and method facts. Enrichment result assembly owns
the public `EnrichmentWorkflowResult.caveats` entries.

## Non-Goals

Caveats must not become a dumping ground for arbitrary messages, debug logs,
stack traces, documentation excerpts, or broad narrative notes.

Caveats must not duplicate full provenance.

Caveats must not be used to excuse execution on invalid inputs. Validation
failures should remain validation failures.

Caveats must not become a second diagnostics system for every internal counter.
Only facts that materially affect result interpretation, supported-scope
claims, scientific assumptions, or user-visible limitations should be promoted
to caveats.

## Consequences

Positive consequences:

- Scientific warnings remain visible on successful results.
- Downstream tooling can match caveats by stable codes instead of parsing logs
  or prose.
- Validators and interpreters stay focused on structured facts and executable
  eligibility.
- Provenance remains available for reproducibility without becoming the only
  place users can discover caveats.

Negative consequences:

- Workflow assembly code must maintain caveat builders as public result
  contracts evolve.
- Adding a new caveat requires a stable code, tests, and documentation review.

Neutral consequences:

- Existing diagnostics and provenance remain valid and may continue to carry
  detailed machine-readable state.
- Public result objects keep caveats as data fields rather than adding warning
  renderers, logging hooks, or presentation behavior.

## Implementation Guidance

Use the shared `phospy.contracts.result_caveats.ResultCaveat` model for public
workflow result caveats.

Keep workflow-specific caveat builders near workflow assembly code, for
example:

- `phospy.workflows.differential.caveats`
- `phospy.workflows.kinase.caveats`
- `phospy.workflows.signalome.caveats`
- `phospy.workflows.enrichment.caveats`

Result constructors should validate top-level caveats with the common caveat
contract. They should not infer caveats from result tables, provenance, or logs.

`ResultCaveat.details` is JSON-like scientific state. Constructors freeze it
recursively with the shared immutable JSON primitive, reject non-string keys,
duplicate keys, unsupported objects, and non-finite floats, and preserve
`ContractValidationError` at the caveat boundary. Direct detail access exposes
immutable mapping/tuple nodes; `to_payload()` and dataclass deep-copy helpers
return fresh ordinary `dict`/`list` JSON payloads.

Tests for new caveats should assert code, severity, and important structured
details. Tests should not rely on exact message text unless the public wording
itself is the contract under review.

## Related Records

- [ADR-0005: Result Model Design for PhosPy](adr_0005__result_model_design.md)
- [ADR-0007: Validation Domain Architecture for PhosPy](adr_0007_validation_domain_architecture.md)
- [ADR-0013: Scientific Parity Strategy and Parity-Testing Policy for PhosPy](adr_0013_scientific_parity_strategy_and_parity_testing_policy.md)
- [ADR-0022: Python API as the Supported Scientific Workflow Interface](adr_0022_python_api_as_supported_scientific_workflow_interface.md)
- [ADR-0030: Ranked-List Enrichment Prerequisites](adr_0030_ranked_list_enrichment_prerequisites.md)
- [Validation Ownership Map](../validation-ownership.md)
