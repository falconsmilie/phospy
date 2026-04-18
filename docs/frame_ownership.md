# Frame Ownership Policy

PhosPy uses one package-wide rule for pandas `DataFrame` / `Series` ownership:

1. Public boundary objects own copies of caller-provided frames.
2. Internal stage DTOs may alias already-owned frames.
3. Internal assembly paths may transfer ownership without re-copying.

## Public Boundaries

By default, public boundary models defensively copy incoming pandas objects:

- `AnalysisReadyPhosphoDataset`
- `ReferenceBundle`
- `KinaseScoringResult`
- `KinasePredictionResult`
- `KinaseActivityResult`
- `SignalomeAssignments`
- `SignalomeModules`
- `KinaseNetwork`
- `SignalomeWorkflowResult.expanded_signalome`

This preserves mutation isolation for callers: mutating input frames after
construction must not mutate already-constructed public objects.

## Internal Transfers

Internal workflow/builder assembly paths may call internal `_from_owned(...)`
constructors on boundary models to transfer ownership of frames that are already
owned within PhosPy.

These transfer constructors skip an additional defensive copy and reduce copy
churn in performance-sensitive flow assembly.

## Internal DTO Rule

Internal stage DTOs should not add defensive deep copies by default when they
receive already-owned data from prior stages.

Ownership-transfer safety requirement:

- once a frame is transferred as owned to a downstream stage, upstream code must
  treat it as immutable.
