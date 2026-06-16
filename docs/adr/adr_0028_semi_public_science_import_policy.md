# ADR-0028: Semi-Public Science Import Policy

## Status

- **ADR ID:** ADR-0028
- **Title:** Semi-Public Science Import Policy
- **Status:** Accepted
- **Date:** 2026-06-15
- **Decision Type:** API Governance

## Context

ADR-0001 defines `phospy` and `phospy.api` as the public API contract. Some
lower-level `phospy.science.*` modules are nevertheless used by extension,
parity, performance, and backend-contract tests. Treating every science module
as either fully public or freely removable is too blunt: several compatibility
routes cannot be removed safely until their status is explicit.

The current ambiguity affects:

- `PreprocessingStageMetadata`
- `phospy.science.signalomes.clustering.protocol`
- `phospy.science.signalomes.clustering.exact_python`
- the rank-fusion helper
  `phospy.science.prediction.scoring.fuse_profile_and_motif_scores_by_rank_weight`

## Decision

PhosPy uses three import-route statuses:

1. **Public** routes are the user-facing API governed by ADR-0001:
   `phospy` and `phospy.api`.
2. **Semi-public** science routes are documented lower-level compatibility
   routes. Their listed non-underscored exports are supported import contracts,
   but they are not promoted into `phospy.api`.
3. **Internal** routes are implementation details. They may be used by PhosPy's
   own tests and modules, but they are not supported external import contracts.

The following routes are semi-public:

| Route | Supported names |
| --- | --- |
| `phospy.science.datasets.preprocessing.stage_registry` | `PreprocessingStageMetadata` and the registry helper functions exported in `__all__`. |
| `phospy.science.signalomes.clustering.protocol` | `ClusterTreeEngine`, `SignalomeClusteringEngine`. |
| `phospy.science.signalomes.clustering.exact_python` | Non-underscored compatibility facade names exported in `__all__`. |
| `phospy.science.prediction.scoring` | `fuse_profile_and_motif_scores_by_rank_weight`. |

Semi-public support is route-specific. It does not make sibling modules,
private implementation modules, or underscored helper functions public.

## Current Classification

`PreprocessingStageMetadata` is semi-public only from
`phospy.science.datasets.preprocessing.stage_registry`. It remains a
compatibility alias for the preprocessing stage contract used by advanced
preprocessing-stage registration. It is intentionally not exported from
`phospy` or `phospy.api`.

`phospy.science.signalomes.clustering.protocol` is semi-public as the clustering
backend protocol route. The supported names are `ClusterTreeEngine` and
`SignalomeClusteringEngine`.

`phospy.science.signalomes.clustering.exact_python` is semi-public as an
exact-Python clustering compatibility facade. Supported external imports are
the non-underscored names listed in its `__all__`.

`fuse_profile_and_motif_scores_by_rank_weight` is semi-public only from
`phospy.science.prediction.scoring`. It is kept for parity checks and advanced
scoring diagnostics. The underscored implementation helper behind it remains
internal.

## Unsupported Routes

Unsupported routes include:

- root or API imports for semi-public science helpers, such as
  `from phospy.api import PreprocessingStageMetadata`
- direct use of underscored helper symbols as public imports
- private implementation modules unless separately documented in the API guide
- wildcard imports that expect underscored names to appear

Existing internal PhosPy tests may still import internal modules when they are
testing implementation details. That does not create an external compatibility
contract.

## Consequences

Positive consequences:

- Import-route support is explicit and testable.
- Follow-up removal or deprecation tickets can distinguish public,
  semi-public, and internal paths.
- `phospy.api` remains the authoritative public API surface.
- Private helpers are not exposed by documenting or testing them as public API.

Negative consequences:

- Some lower-level compatibility routes remain supported until a normal
  deprecation/removal process is followed.
- Contributors must check both public and semi-public route tests before
  moving science modules.

## Validation and Review Criteria

Changes affecting these routes must satisfy:

1. Public import contract tests for `phospy` and `phospy.api`.
2. Semi-public import contract tests for the routes listed in this ADR.
3. Negative tests showing semi-public helpers are not available from
   unsupported public routes.
4. Negative tests showing underscored private helpers are not exported through
   wildcard/public-export routes.
5. Documentation updates when a new semi-public science route is added or
   removed.

## Non-Goals

This ADR does not:

- promote the listed science helpers into `phospy.api`
- make all `phospy.science.*` modules public
- remove any compatibility module
- freeze private helper names
- change workflow behavior or scientific algorithms

## Related Records

- [ADR-0001: Public API Contract for PhosPy](adr_0001_public_api_contract.md)
- [API Guide](../api/guide.md)
