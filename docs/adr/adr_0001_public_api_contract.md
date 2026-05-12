# ADR: Public API Contract for PhosPy

## Document Control

- **ADR ID:** ADR-0001
- **Title:** Public API Contract for PhosPy
- **Status:** Accepted
- **Date:** 2026-05-02
- **Decision Type:** Architecture Decision Record

## Abstract

This ADR defines the supported public API contract and public namespace
ownership for PhosPy. The implementation now includes a broader `phospy.api`
contract and grouped public configuration models. Those decisions are stable
product governance, not incidental documentation.

## Status

Accepted.

This ADR supersedes earlier narrower descriptions of the public contract.

## Context and Problem Statement

The codebase now exposes:

- a governed top-level public workflow namespace (`phospy`)
- a broader stable contract namespace (`phospy.api`)
- grouped configuration objects by user intent rather than one flat config
  surface

Without explicit ADR governance, future changes can accidentally expand or
redefine the public contract through internal modules and compatibility shims.

## Decision Drivers

1. Preserve a clear public surface while allowing internal refactors.
2. Keep workflow-science configuration explicit and reviewable.
3. Prevent hidden API growth through internal module imports.
4. Keep default usage simple while keeping advanced configuration explicit.
5. Keep validation and scientific trade-offs visible in public configuration.

## Decision

### Public Namespace Ownership

1. `phospy` is a first-class public workflow namespace.
2. `phospy.api` is a broader public contract namespace for request/result/config
   models, enums, references, and exceptions.
3. `phospy` exposes core workflow entrypoints:
   `AnalysisReadyDatasetBuilder`, `AnalysisReadyPhosphoDataset`,
   `KinaseWorkflow`, `SignalomeWorkflow`, and `DifferentialAnalysisWorkflow`.
4. `DifferentialAnalysisWorkflow` is intentionally exported from top-level
   `phospy`.
5. Internal modules remain non-public unless explicitly re-exported from an
   approved public namespace.

### Public API Structure

The public contract is workflow-oriented and includes representative stable
types such as:

- Dataset/building: `DatasetBuildRequest`, `AnalysisReadyDatasetBuilder`,
  `AnalysisReadyPhosphoDataset`
- Kinase lane: `KinaseWorkflowRequest`, `KinaseWorkflowResult`,
  `KinaseScoringConfig`, `KinasePredictionConfig`, `KinaseActivityConfig`,
  `KinaseWorkflow`
- Differential lane: `DifferentialAnalysisRequest`,
  `DifferentialAnalysisResult`, `EmpiricalBayesConfig`,
  `MultipleTestingConfig`,
  `DifferentialAnalysisWorkflow`
- Signalome lane: `SignalomeWorkflowRequest`, `SignalomeWorkflowResult`,
  `SignalomeConfig`, `SignalomeScientificConfig`,
  `SignalomeClusteringConfig`, `SignalomeValidationConfig`,
  `SignalomeOutputConfig`, `SignalomePerformanceConfig`, `SignalomeWorkflow`
- Shared contracts: `ReferenceBundle`, `ReferencePreset`, `Organism`,
  public exception families re-exported through `phospy.api`

This list is intentionally representative rather than an exhaustive import dump.

### DifferentialAnalysisWorkflow Compatibility Expectations

For `DifferentialAnalysisWorkflow`, public compatibility includes:

1. top-level import path stability
   (`from phospy import DifferentialAnalysisWorkflow`)
2. stable request/config contract shape (`run(request)`)
3. stable documented result object behavior and fields
4. stable documented error behavior for invalid design/contrast inputs
5. public README/API examples kept in sync with tests

`from phospy.api import DifferentialAnalysisWorkflow` is a supported public
route.

### Grouped Public Configuration Governance

Public config is grouped by user intent, not exposed as one large flat surface.

Required public config groupings include:

- scientific
- validation
- output
- clustering
- performance
- prediction
- preprocessing
- related user-intent groupings where applicable

Public presets are allowed only when they describe concrete behavior.
Presets must not bypass validation.
Presets must not hide scientific or performance trade-offs.
Advanced configuration remains available through explicit config objects.

### Workflow Contract Shape

1. Each public workflow accepts one typed request object and returns one typed
   result object.
2. Public result models remain data contracts, not service facades.
3. Validation remains mandatory at public boundaries.

## Consequences

### Positive Consequences

- Public import behavior is explicit and reviewable.
- Grouped config keeps scientific choices visible.
- Internal module splits can continue without changing public contract.
- Public API growth requires explicit namespace promotion.

### Negative Consequences

- Internal convenience imports are not public by default.
- Contributors must re-export intentionally through `phospy.api` when promoting
  new public contract types.

## Affected Modules

- `src/phospy/__init__.py`
- `src/phospy/api/__init__.py`
- `src/phospy/api/configs/__init__.py`
- `src/phospy/api/configs/common.py`
- `src/phospy/api/configs/dataset.py`
- `src/phospy/api/configs/kinase.py`
- `src/phospy/api/configs/prediction.py`
- `src/phospy/api/configs/signalome.py`

## Scope Boundaries

This ADR defines public API governance and namespace ownership. It does not
define internal package splitting, validation-domain ownership, or test policy
details (covered by ADR-0007, ADR-0010, ADR-0014, ADR-0016, and ADR-0017).

## Validation and Review Criteria

Future changes must satisfy all of the following:

1. `phospy` remains the governed top-level public workflow namespace.
2. `phospy` keeps the governed top-level workflow entrypoints, including
   `DifferentialAnalysisWorkflow`.
3. `phospy.api` remains the stable broader public namespace and is not required
   to mirror every top-level entrypoint.
4. New public presets encode concrete behavior and still pass validation.
5. Internal modules are not treated as public unless re-exported intentionally.
6. Public API changes are reviewed as contract changes, not incidental refactors.

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi,
R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative
phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14),
i349-i356.

YangLab. (n.d.). *PhosR* [Computer software]. GitHub.
https://github.com/PYangLab/PhosR
