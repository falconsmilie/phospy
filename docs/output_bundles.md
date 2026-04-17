# Output Bundles

`SimpleKinaseWorkflowResult` bundle persistence is implemented as an external
service in `phospy.io`, not as methods on the result models.

This follows ADR-005 (result objects stay typed containers) and ADR-015
(fixtures and manifest assets stay in test/docs areas, not runtime reference
resources).

## Supported Service

```python
from pathlib import Path

from phospy.io import (
    SimpleKinaseWorkflowConfigSnapshot,
    load_simple_kinase_workflow_bundle,
    save_simple_kinase_workflow_bundle,
)

snapshot = SimpleKinaseWorkflowConfigSnapshot.from_request(workflow_request)
save_simple_kinase_workflow_bundle(
    workflow_result,
    Path("./bundle"),
    config_snapshot=snapshot,
)
loaded = load_simple_kinase_workflow_bundle(Path("./bundle"))
```

`loaded` returns:

- `result`: reconstructed `SimpleKinaseWorkflowResult`
- `config_snapshot`: typed config snapshot for reproducibility
- `manifest_version`

## Manifest Contract (v1)

- File: `manifest.json`
- Required:
`bundle_type == "simple_kinase_workflow_result"`
- Required:
`manifest_version == 1`
- Required sections:
`dataset`, `resolved_references`, `outputs`, `config_snapshot`

## Bundle Contents (v1)

Default `csv` layout:

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

Manifest metadata explicitly captures:

- dataset metadata:
`organism`, full `transformation_state`
- resolved reference metadata:
resolved reference `organism`
- stage output table paths:
scoring, prediction, activity
- config snapshot path:
`config/snapshot.json`

The manifest is versioned from day one so future formats can evolve without
guesswork.
