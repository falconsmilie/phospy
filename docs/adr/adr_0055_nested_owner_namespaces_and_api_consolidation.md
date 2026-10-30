# ADR-0055: Nested Owner Namespaces and API Consolidation

## Status

- **ADR ID:** ADR-0055
- **Title:** Nested Owner Namespaces and API Consolidation
- **Status:** Accepted
- **Date:** 2026-10-30
- **Decision type:** Architecture namespace and public import refinement
- **Amends:** ADR-0053 and ADR-0054
- **Pre-migration revision:** `f180248b2f2c2f0da3323d54549def8128005d72`

The CORE-002 implementation applies the supported-API portion of this
decision: the specialist facade and inventory move beneath `phospy.api`, and
the root conveniences retire. The reusable owner-package moves beneath
`phospy.core` remain sequenced in later tickets.

## Context

ADR-0053 established separate logical owners and ADR-0054 established the
ordinary and specialist aggregate APIs. Their accepted implementation placed
each owner directly beneath `phospy` and placed the specialist facade at
`phospy.advanced`.

The next directory-only refinement groups reusable owner packages beneath an
organisational `phospy.core` namespace and groups both supported API tiers
beneath `phospy.api`. This must not collapse the existing owners into one Core
dependency layer. Architecture checks that infer ownership from the first
component below `phospy` would lose sibling boundaries after that move.

The committed pre-migration candidate was inspected at
`f180248b2f2c2f0da3323d54549def8128005d72`. The authoritative comparison is
`architecture/public-api-baseline.json`, schema version 2, captured with Python
3.11.15. It records three current root
conveniences, 66 ordinary `phospy.api` exports and 126 specialist
`phospy.advanced` exports. Because the starting candidate was committed and
clean, the revision and maintained baseline are sufficient starting-state
evidence; no duplicate uncommitted evidence snapshot is required.

## Decision

### Final physical layout

The following packages remain directly beneath `phospy`:

- `datasets`, `preprocessing`, `differential`, `enrichment`, `kinase` and
  `signalome` as capabilities;
- `io` as the inward adapter; and
- `data` as packaged resources.

The following complete packages move without ownership changes:

| Current package | Final package | Logical owner retained |
| --- | --- | --- |
| `phospy.errors` | `phospy.core.errors` | `errors` |
| `phospy.frames` | `phospy.core.frames` | `frames` |
| `phospy.policies` | `phospy.core.policies` | `policies` |
| `phospy.validation` | `phospy.core.validation` | `validation` |
| `phospy.activities` | `phospy.core.activities` | `activities` |
| `phospy.evidence` | `phospy.core.evidence` | `evidence` |
| `phospy.phosphosites` | `phospy.core.phosphosites` | `phosphosites` |
| `phospy.quantitative` | `phospy.core.quantitative` | `quantitative` |
| `phospy.references` | `phospy.core.references` | `references` |
| `phospy.statistics` | `phospy.core.statistics` | `statistics` |
| `phospy.taxonomy` | `phospy.core.taxonomy` | `taxonomy` |
| `phospy.provenance` | `phospy.core.provenance` | `provenance` |
| `phospy.preparation_records` | `phospy.core.preparation_records` | `preparation_records` |
| `phospy.run_lineage` | `phospy.core.run_lineage` | `run_lineage` |

`phospy.core` is an organisational namespace only. Its future `__init__.py`
is inert: it exports nothing, performs no imports and grants no dependency
permission. When that file is introduced, the owner map will account for that
initializer alone with a narrow `namespace_containers` owner. The owner and
physical-root classifications are added atomically with the initializer; they
are not unused live policy before then. Every child remains assigned to its
existing logical owner.

`phospy.core.run_lineage` retains only the reviewed Run Lineage imports,
including its exact Preprocessing contract dependencies. A physical
`core -> preprocessing -> core` shape is not an owner dependency cycle and is
not used for the owner-DAG gate.

### Supported API layout

The final supported external surfaces are:

| Surface | Policy |
| --- | --- |
| `phospy.api` | Ordinary supported imports; retains the current ordinary symbol set and exact owner-defined identities. |
| `phospy.api.advanced` | Specialist supported imports; replaces `phospy.advanced` and retains its current specialist symbol set and exact identities. |
| `phospy` | Version metadata only; the three convenience exports are retired from the root. |

`AnalysisReadyDatasetBuilder`, `KinaseWorkflow` and `SignalomeWorkflow` remain
available through `phospy.api`. The private inventory moves from
`phospy._api_inventory` to `phospy.api._inventory`. There are no forwarding
packages, `sys.modules` aliases, duplicate definitions or compatibility
registries for the retired locations.

This amends ADR-0054's physical specialist route and root convenience list. It
retains ADR-0054's aggregate-facade model, exact identity rule, ordinary and
specialist symbol classifications, and clean-break compatibility policy.

### Namespace-aware architecture enforcement

Architecture enforcement resolves the shortest namespace that is pure for a
logical owner according to `architecture/module-owners.csv`, rather than
assuming the first component below `phospy` is the package unit. Thus
`phospy.core.frames` and `phospy.core.quantitative` remain separate units at
any depth. The API tiers share the `public_facades` owner, so
`phospy.api.advanced` is the one explicit nested facade boundary needed by
package-level policy.

Physical top-level roots remain useful inventory information and are reported
separately. They do not define owner permissions or owner-level cycles.
Runtime forbidden-import probes derive narrow prefixes from the same owner
map; they may reject `phospy.core.quantitative` without rejecting the entire
`phospy.core` namespace.

All existing substantive gates remain independent and in force:

1. the complete statically observed owner graph and exact owner-edge rules;
2. the policy-filtered owner DAG;
3. the three exact permanent condensation classifications;
4. the complete eager module graph;
5. exact target-path and Signalome-to-Dataset rules;
6. private cross-owner target checks; and
7. public facade identity and removed-route checks.

The owner map and exact path policies change atomically with each later move.
This decision does not temporarily permit both old and new production routes.

## Consequences

- Directory navigation groups reusable packages without weakening their
  independent dependency boundaries.
- Architecture checks can exercise the final depth before production files
  move and will still expose forbidden Core sibling imports.
- Public users ultimately have one ordinary namespace and one nested
  specialist namespace.
- Later migration tickets must update owner rows, exact paths, strict-type
  coverage, facade classifications and generated inventories together.

## Non-goals

This decision does not move production packages in CORE-001, change public
exports before their migration ticket, alter scientific algorithms or payload
schemas, rename workflow `public.py` modules, move `io/bundles/_shared`, or
create a universal Core dependency permission.

## Acceptance criteria

This amendment is implemented when nested synthetic fixtures preserve owner
siblings, forbidden directions, private-target checks and eager-cycle
detection; runtime prefixes remain owner-specific; the current live gates and
inventory still pass; and no production package or public export has moved in
the architecture-preparation ticket.

