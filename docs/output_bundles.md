# Output Bundles

Workflow bundle persistence is implemented as external services in `phospy.io`,
not methods on result DTOs.

This keeps public result models as nested typed containers and keeps persistence
as an explicit I/O concern.

## Supported Services

```python
from pathlib import Path

from phospy.io import (
    KinaseWorkflowConfigSnapshot,
    SignalomeWorkflowConfigSnapshot,
    load_kinase_workflow_bundle,
    load_signalome_workflow_bundle,
    save_kinase_workflow_bundle,
    save_signalome_workflow_bundle,
)

kinase_snapshot = KinaseWorkflowConfigSnapshot.from_request(kinase_request)
save_kinase_workflow_bundle(
    kinase_result,
    Path("./kinase_bundle"),
    config_snapshot=kinase_snapshot,
)
loaded_kinase = load_kinase_workflow_bundle(Path("./kinase_bundle"))

signalome_snapshot = SignalomeWorkflowConfigSnapshot.from_request(signalome_request)
save_signalome_workflow_bundle(
    signalome_result,
    Path("./signalome_bundle"),
    config_snapshot=signalome_snapshot,
)
loaded_signalome = load_signalome_workflow_bundle(Path("./signalome_bundle"))
```

Loaded bundle objects include:

- reconstructed nested workflow result DTO
- typed config snapshot
- `manifest_version`

## Manifest Contract (v1)

Kinase manifest:

- `bundle_type == "kinase_workflow_result"`
- `manifest_version == 1`
- top-level sections:
  `dataset`, `resolved_references`, `outputs`, `config_snapshot`

Signalome manifest:

- `bundle_type == "signalome_workflow_result"`
- `manifest_version == 1`
- top-level sections:
  `dataset`, `resolved_references`, `upstream_kinase_outputs`,
  `signalome_outputs`, `config_snapshot`

Both manifests store dataset organism and full transformation-state payload.

## Bundle Contents (Default CSV Layout)

Kinase:

```text
manifest.json
config/snapshot.json
dataset/phospho.csv
dataset/site_metadata.csv
dataset/sample_metadata.csv          # optional
dataset/total.csv                    # optional
references/kinase_substrate_map.csv
references/site_sequences.csv
scoring/profile_scores.csv
scoring/motif_scores.csv             # optional
scoring/combined_scores.csv          # optional
scoring/weights.csv                  # optional
prediction/pred_mat.csv
prediction/substrate_list.csv        # optional
activity/weighted_activity.csv       # optional
activity/ksea_scores.csv             # optional
activity/ksea_counts.csv             # optional
activity/target_counts.csv           # optional
activity/target_table.csv            # optional
```

Signalome:

```text
manifest.json
config/snapshot.json
dataset/phospho.csv
dataset/site_metadata.csv
dataset/sample_metadata.csv          # optional
dataset/total.csv                    # optional
references/kinase_substrate_map.csv
references/site_sequences.csv
scoring/profile_scores.csv
scoring/motif_scores.csv             # optional
scoring/combined_scores.csv          # optional
scoring/weights.csv                  # optional
prediction/pred_mat.csv
prediction/substrate_list.csv        # optional
activity/weighted_activity.csv       # optional
activity/ksea_scores.csv             # optional
activity/ksea_counts.csv             # optional
activity/target_counts.csv           # optional
activity/target_table.csv            # optional
signalome/module_assignments.csv
signalome/signalome_modules.csv
signalome/kinase_network_edges.csv
signalome/kinase_network_nodes.csv   # optional
signalome/expanded_signalome.csv     # optional
```

Optional means contract-optional, not always absent.
In the default supported kinase lane, scoring currently populates profile/motif/combined/weight tables.

## Optional Output Semantics

- `activity/*` tables are present only when `kinase_result.activity_result` is present.
- `prediction/substrate_list` is optional.
- `signalome/kinase_network_nodes` is optional.
- `signalome/expanded_signalome` is optional and currently absent (`None`) in the default signalome route.

## Config Snapshot Fields

Kinase config snapshot:

- `scoring_config.min_substrates`
- `prediction_config.top_k`
- `prediction_config.ensemble_size`
- `activity_config` fields when activity is configured

Signalome config snapshot:

- `signalome_config.substrate_support_cutoff`
- `signalome_config.network_correlation_threshold`

Manifest versioning starts at v1 so future format evolution is explicit.
