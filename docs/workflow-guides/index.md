# Workflow Guides

PhosPy has one supported workflow chain:

1. Build `AnalysisReadyPhosphoDataset`
2. Run `KinaseWorkflow`
3. Run `SignalomeWorkflow` (optional, downstream of kinase result)

## Signalome Protein-Identity Contract

The supported signalome lane requires explicit protein identity:

- `kinase_result.dataset.site_metadata.protein_id` must exist and be non-empty
- gene-symbol site-ID prefixes are not a protein-identity fallback
- this is an intentional scientific boundary for protein-aware grouping/module
  assignment, not runtime strictness by accident
- builder input flexibility does not weaken this downstream workflow contract

## Stage-by-Stage Guide

- Beginner walkthrough: [Quickstart: first workflow](../getting-started/quickstart-first-workflow.md)
- API shape and request/result contracts: [API Guide](../api.md)
- Validation rules at each boundary: [Validation Guide](../validation.md)
- CLI workflow execution: [CLI Guide](../cli.md)
- Persisting results to bundles: [Output Bundles](../output_bundles.md)

## Common Failure Points

- Input table shape/required columns:
  [Troubleshooting](../troubleshooting.md#dataset-build-and-input-shape-issues)
- Reference and organism compatibility:
  [Troubleshooting](../troubleshooting.md#reference-and-organism-resolution-issues)
- Signalome prerequisites (`protein_id`):
  [Troubleshooting](../troubleshooting.md#signalome-preconditions-and-runtime-issues)
