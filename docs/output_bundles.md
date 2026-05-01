# Output Bundles

PhosPy has two related output paths.

1. The CLI and publisher helpers write simple workflow output directories.
2. Bundle services write reloadable output bundles with a manifest and config
   snapshot.

Use the simple CLI layout when you only need files. Use bundle services when you
want to save and reload a workflow result object.

## CLI and Publisher Output Layout

`publish_dataset(...)` writes:

```text
dataset/
  phospho.*
  site_metadata.*
  sample_metadata.*   # optional
  total.*             # optional
  manifest.json
```

`publish_kinase_workflow(...)` also writes:

```text
kinase/
  scoring/
    profile_scores.*
    motif_scores.*                    # optional
    rank_weighted_fusion_scores.*     # optional
    score_fusion_weights.*            # optional
  prediction/
    pred_mat.*
    substrate_list.*                  # optional
  activity/                           # optional
  references/
    kinase_substrate_map.*
    site_sequences.*
  manifest.json
```

`publish_signalome_workflow(...)` also writes:

```text
signalome/
  module_assignments.*
  signalome_modules.*
  kinase_network_edges.*
  kinase_network_nodes.*                  # optional
  kinase_network_candidate_correlations.* # optional
  expanded_signalome.*                    # optional
  manifest.json
```

Supported table formats are `csv`, `tsv`, and `parquet`.

## Reloadable Bundle Services

Public bundle services live in:

```python
from phospy.io.bundles.kinase import load_kinase_workflow_bundle
from phospy.io.bundles.kinase import save_kinase_workflow_bundle
from phospy.io.bundles.signalome import load_signalome_workflow_bundle
from phospy.io.bundles.signalome import save_signalome_workflow_bundle
```

Bundle writers require a config snapshot:

```python
from pathlib import Path

from phospy.io.bundles.kinase import KinaseWorkflowConfigSnapshot
from phospy.io.bundles.kinase import save_kinase_workflow_bundle

snapshot = KinaseWorkflowConfigSnapshot.from_request(request)
written = save_kinase_workflow_bundle(
    result,
    Path("bundle-out"),
    config_snapshot=snapshot,
    output_format="csv",
)
```

A saved bundle contains:

- data tables
- reference tables
- workflow output tables
- `manifest.json`
- `config/snapshot.json`

Loaders return `LoadedKinaseWorkflowBundle` or `LoadedSignalomeWorkflowBundle`.
Each object contains the reconstructed result, the config snapshot, and the
manifest version.

## Bundle Safety Rules

Bundle-relative manifest paths must stay inside the bundle root. Absolute paths
and paths that escape the root are rejected.

Optional tables remain optional. Missing optional activity, motif, node,
site-membership, or context tables are represented as `None` when loaded.

The manifest records output format, bundle kind, manifest version, table paths,
reference organism, and provenance where available. Dataset metadata in the
manifest includes `processing_state`, so explicit quantitative meaning
(`phosphosite_log_abundance`, `phospho_total_log_ratio`, or mixed state) is
preserved in published outputs.

Bundle loaders reconstruct `dataset.processing_state` and
`dataset.intensity_scale_state` from saved payloads, including mixed
corrected/uncorrected total-protein quantitative meaning and row-level
total-protein correction diagnostics.

## Choosing the Right Output Path

| Need | Use |
| --- | --- |
| Run from terminal and inspect CSV files | CLI |
| Write result tables from Python | `publish_*` helpers |
| Save and reload a result object | bundle services |
| Keep config snapshot beside outputs | bundle services |
