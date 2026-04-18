# API Guide

PhosPy currently exposes one supported end-to-end rewrite route:

`DatasetBuildRequest -> AnalysisReadyDatasetBuilder.run(request) -> AnalysisReadyPhosphoDataset -> KinaseWorkflow.run(request)`

## Support Boundary

- `src/phospy/`: supported public rewrite package
- `legacy_archive/phospy_legacy/`: migration reference only; not an installed package target
- Dataset builder: supported
- Kinase workflow: supported
- Signalome workflow: first real vertical slice implemented

## Public Types

Import from top-level `phospy`:

- Dataset and references:
`AnalysisReadyPhosphoDataset`, `Organism`, `ReferencePreset`, `ReferenceBundle`
- Builder:
`DatasetBuildRequest`, `AnalysisReadyDatasetBuilder`
- Workflows and requests:
`KinaseWorkflow`, `KinaseWorkflowRequest`,
`SignalomeWorkflow`, `SignalomeWorkflowRequest`
- Config models:
`KinaseScoringConfig`, `KinasePredictionConfig`, `KinaseActivityConfig`,
`SignalomeConfig`
- Result models:
`KinaseWorkflowResult`, `SignalomeWorkflowResult`,
`KinaseScoringResult`, `KinasePredictionResult`, `KinaseActivityResult`

All public executors use `run(request)`.

## Builder Contract

`DatasetBuildRequest` accepts in-memory pandas `DataFrame` values or file paths.

- Required: `phospho`, `site_metadata`
- Optional: `sample_metadata`, `total`, `organism`, `transformation_state`
- `ReferencePreset.AUTO` requires `dataset.organism` at workflow execution time
- Bundled references are currently available for rat only
- `ReferencePreset.HUMAN` and `ReferencePreset.MOUSE` are part of the public enum,
  but bundled runtime resolution for those presets is intentionally unsupported in
  this release and fails with `UnsupportedOrganismError`
- For non-rat workflows, supply an explicit `ReferenceBundle` instead of relying
  on bundled presets

## Result Contract

`KinaseWorkflowResult` keeps stage outputs nested:

- `result.scoring_result`
- `result.prediction_result`
- `result.activity_result`

No mirrored top-level convenience aliases are part of the rewrite contract.
Reusable save/load output bundles are provided as external services in
`phospy.io` (`save_kinase_workflow_bundle`,
`load_kinase_workflow_bundle`) so result models remain plain containers.

`SignalomeWorkflowResult` currently returns real tables (`module_assignments`,
`signalome_modules`, `kinase_network`) from the first vertical slice, while
`expanded_signalome` remains optional and currently `None`.

See [`output_bundles.md`](output_bundles.md) for manifest format and table
inventory.

## Example

```python
import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    Organism,
    ReferencePreset,
    KinaseWorkflow,
    KinaseWorkflowRequest,
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=pd.DataFrame(
            {"sample_a": [1.0], "sample_b": [1.2]},
            index=["MAPK14;Y182;"],
        ),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=["MAPK14;Y182;"],
        ),
        organism=Organism.RAT,
    )
)

result = KinaseWorkflow().run(
    KinaseWorkflowRequest(dataset=dataset, references=ReferencePreset.AUTO)
)

pred_mat = result.prediction_result.pred_mat
```

## CLI

The rewrite CLI now supports the same narrow lane from files:

1. `dataset-build` (files -> analysis-ready dataset outputs)
2. `kinase` (files -> dataset build + kinase workflow outputs)

Examples:

```bash
phospy dataset-build \
  --phospho ./input/phospho.csv \
  --site-metadata ./input/site_metadata.csv \
  --organism rat \
  --outdir ./out

phospy kinase \
  --phospho ./input/phospho.csv \
  --site-metadata ./input/site_metadata.csv \
  --organism rat \
  --reference auto \
  --outdir ./out
```

Input files currently support `.csv`, `.tsv`/`.txt`, and `.parquet`.
