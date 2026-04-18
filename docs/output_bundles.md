# Output Bundles

`KinaseWorkflowResult` and `SignalomeWorkflowResult` bundle persistence
are implemented as external services in `phospy.io`, not as methods on the
result models.

This follows ADR-005 (result objects stay typed containers) and ADR-015
(fixtures and manifest assets stay in test/docs areas, not runtime reference
resources).

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

`loaded_kinase` returns:

- `result`: reconstructed `KinaseWorkflowResult`
- `config_snapshot`: typed config snapshot for reproducibility
- `manifest_version`

`loaded_signalome` returns:

- `result`: reconstructed `SignalomeWorkflowResult`
- `config_snapshot`: typed config snapshot for reproducibility
- `manifest_version`

## Manifest Contracts (v1)

Kinase:

- File: `manifest.json`
- Required:
`bundle_type == "kinase_workflow_result"`
- Required:
`manifest_version == 1`
- Required sections:
`dataset`, `resolved_references`, `outputs`, `config_snapshot`

Signalome:

- File: `manifest.json`
- Required:
`bundle_type == "signalome_workflow_result"`
- Required:
`manifest_version == 1`
- Required sections:
`dataset`, `resolved_references`, `upstream_kinase_outputs`,
`signalome_outputs`, `config_snapshot`

## Bundle Contents (v1)

Kinase (default `csv` layout):

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
activity/activity_scores.csv         # optional
```

Signalome (default `csv` layout):

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
activity/activity_scores.csv         # optional
signalome/module_assignments.csv
signalome/signalome_modules.csv
signalome/kinase_network_edges.csv
signalome/kinase_network_nodes.csv   # optional
signalome/expanded_signalome.csv     # optional
```

Signalome manifest metadata explicitly captures:

- dataset metadata:
`organism`, full `transformation_state`
- resolved reference metadata:
resolved reference `organism`
- upstream kinase stage output table paths:
scoring, prediction, activity
- signalome output table paths:
module assignments, module matrix, kinase network, expanded signalome
- config snapshot path:
`config/snapshot.json`

The manifest is versioned from day one so future formats can evolve without
guesswork.
