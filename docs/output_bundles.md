# Output Bundles

Workflow bundle persistence is implemented as external services in `phospy.io`,
not methods on result DTOs.

> Audience: users who need reproducible saved outputs and maintainers owning I/O contracts.
> New users can start with [Quickstart](getting-started/quickstart-first-workflow.md)
> and return here when persisting workflow results.

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
  `dataset`, `resolved_references`, `outputs`, `provenance`, `config_snapshot`

Signalome manifest:

- `bundle_type == "signalome_workflow_result"`
- `manifest_version == 1`
- top-level sections:
  `dataset`, `resolved_references`, `upstream_kinase_outputs`,
  `signalome_outputs`, `provenance`, `config_snapshot`
- `signalome_outputs.metadata` includes:
  - `kinase_network_nodes_present`
  - `expanded_signalome_present`
  - `module_selection_diagnostics` payload (strategy, selected count,
    threshold/candidate diagnostics, and degeneracy counters)
  - `score_preconditioning_diagnostics` payload
  - `network_correlation_diagnostics` payload (finite vs undefined candidate
    counts, status-specific counts, and skipped-edge counts)

Both manifests store dataset organism, full `intensity_scale_state` payload,
and full `processing_state` payload.

`provenance` is machine-readable run metadata (`RunProvenance`) and includes:

- environment versions (`phospy`, Python, dependency versions)
- input/output table fingerprints
- preprocessing stage execution trace (shape/hash/drop/imputation metadata)
- reference identity/fingerprints
- workflow parameters and random-state policy

Legacy bundles without top-level `provenance` remain loadable; loaders reconstruct
results with `result.provenance=None` for those manifests.

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
scoring/rank_weighted_fusion_scores.csv          # optional
scoring/score_fusion_weights.csv     # optional
prediction/pred_mat.csv
prediction/substrate_list.csv        # optional
activity/weighted_activity.csv       # optional
activity/thresholded_substrate_mean_activity.csv             # optional
activity/thresholded_substrate_counts.csv             # optional
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
scoring/rank_weighted_fusion_scores.csv          # optional
scoring/score_fusion_weights.csv     # optional
prediction/pred_mat.csv
prediction/substrate_list.csv        # optional
activity/weighted_activity.csv       # optional
activity/thresholded_substrate_mean_activity.csv             # optional
activity/thresholded_substrate_counts.csv             # optional
activity/target_counts.csv           # optional
activity/target_table.csv            # optional
signalome/module_assignments.csv
signalome/signalome_modules.csv
signalome/kinase_network_edges.csv
signalome/kinase_network_nodes.csv   # optional
signalome/kinase_network_candidate_correlations.csv   # optional
signalome/expanded_signalome.csv     # optional
```

When present, `signalome/expanded_signalome.csv` uses the same flattened schema
as `SignalomeWorkflowResult.expanded_signalome` (see `docs/api.md`): focal
kinase, row kind, assignment policy, linked/regulatory module metadata, and
site-level membership rows with stable `site_order`.

Signalome manifests also persist:

- `score_preconditioning_diagnostics` metadata (`policy`, `input_row_count`,
  `dropped_all_missing_row_count`, `retained_row_count`) so dropped
  all-missing downstream-score rows are explicit in published outputs and
  reloads.
- `network_correlation_diagnostics` metadata and optional
  `signalome/kinase_network_candidate_correlations.csv` so undefined
  correlations remain distinguishable from true `0.0` correlations across
  bundle round-trips.

Optional means contract-optional, not always absent.
In the default supported kinase lane, scoring populates `profile_scores` and
`rank_weighted_fusion_scores`; diagnostic `motif_scores` and `score_fusion_weights` are written only when
`scoring_config.include_diagnostic_scoring_tables=True`.
Scoring semantics are upstream-stage stable: they are determined by dataset +
resolved references + scoring config, and are not redefined by prediction mode
or reference input provenance (preset vs equivalent explicit bundle).

## Optional Output Semantics

- `activity/*` tables are present only when `kinase_result.activity_result` is present.
- `prediction/substrate_list` is optional.
- `signalome/kinase_network_nodes` is optional.
- `signalome/kinase_network_candidate_correlations` is optional.
- `signalome/expanded_signalome` is optional by contract for compatibility, but
  is populated in the supported signalome executor lane when the workflow
  completes successfully.
- `scoring/motif_scores` and `scoring/score_fusion_weights` are optional diagnostic tables and
  are absent in the default scoring lane.

## Config Snapshot Coverage and Reload Semantics

Bundle loaders always parse `config/snapshot.json` into typed snapshot DTOs:

- kinase: `LoadedKinaseWorkflowBundle.config_snapshot` (`KinaseWorkflowConfigSnapshot`)
- signalome: `LoadedSignalomeWorkflowBundle.config_snapshot` (`SignalomeWorkflowConfigSnapshot`)

`loaded.result` is reconstructed directly from persisted output tables.
`loaded.config_snapshot` is the scientific configuration record for replay and
interpretation (what policies/thresholds generated those outputs).

Kinase config snapshot (all persisted fields):

- `scoring_config.min_substrates`: substrate-support floor used in scoring.
- `scoring_config.include_diagnostic_scoring_tables`: controls whether
  diagnostic `motif_scores`/`score_fusion_weights` tables are expected to be produced.
- `scoring_config.profile_missing_value_strategy`: profile median behavior
  (`"strict"` vs `"median_skipna"`), which can change downstream scores.
- `prediction_config.top_k`: final rank cutoff in prediction output.
- `prediction_config.deterministic_max_selected_kinases`: deterministic lane
  kinase-selection cap.
- `prediction_config.adaptive_ensemble_runs`: adaptive lane ensemble execution
  count.
- `prediction_config.mode`: prediction lane (`"deterministic_ranking"` or
  `"adaptive_ensemble"`).
- `prediction_config.adaptive_policy`: adaptive sampling policy selection.
- `prediction_config.n_iterations`: adaptive sampling iteration budget.
- `prediction_config.random_state`: optional deterministic seed for adaptive
  sampling.
- `activity_config.enabled`
- `activity_config.threshold`
- `activity_config.min_substrates`
- `activity_config.top_n_substrates`

Signalome config snapshot (all persisted fields):

- `signalome_config.substrate_support_cutoff`: cutoff used in module support.
- `signalome_config.network_correlation_threshold`: correlation threshold used
  for kinase-network edge inclusion.
- `signalome_config.network_policy`: network edge policy
  (`"positive_only"`, `"absolute_threshold"`, `"signed"`).
- `signalome_config.assignment_policy`: module assignment policy
  (`"cutoff_binary"`, `"weighted_top"`).
- `signalome_config.score_preconditioning_policy`: downstream score
  preconditioning policy (`"allow_and_report"` or `"error_on_drop"`).
- `signalome_config.module_count`: explicit module count request when set.
- `signalome_config.module_selection_primary_correlation_threshold`
- `signalome_config.module_selection_fallback_correlation_threshold`
- `signalome_config.module_selection_max_clusters`

Deliberate omissions:

- Signalome bundle snapshots do not include upstream kinase `scoring_config` or
  `prediction_config`; those configs belong to the upstream kinase workflow
  request and are not part of `SignalomeWorkflowRequest`.
- Runtime-derived diagnostics (for example module-selection diagnostics) are
  stored under manifest metadata/output tables, not inside `config_snapshot`.

Legacy reload compatibility (normalization on load):

- Kinase legacy snapshots missing newly added fields are normalized to defaults:
  `include_diagnostic_scoring_tables=True`,
  `profile_missing_value_strategy="strict"`,
  `deterministic_max_selected_kinases=10`,
  `adaptive_ensemble_runs=10`,
  `mode="deterministic_ranking"`, `adaptive_policy="stable"`,
  `n_iterations=5`, and `random_state=None`.
- Legacy `prediction_config.ensemble_size` is accepted and mapped to both
  `deterministic_max_selected_kinases` and `adaptive_ensemble_runs`.
- Signalome legacy `signalome_cutoff` is accepted and mapped to both
  `substrate_support_cutoff` and `network_correlation_threshold`.
  Missing `network_policy`, `assignment_policy`,
  `score_preconditioning_policy`, and module-selection fields are normalized to
  current default values.

Manifest versioning starts at v1 so future format evolution is explicit.

## Internal Ownership (Rewrite Boundary)

Public entrypoints stay at:

- `phospy.io.bundles.kinase`
- `phospy.io.bundles.signalome`

Internal implementation is split by concern ownership under:

- `phospy.io.bundles._kinase`: snapshots, manifest contract, writer orchestration,
  loader orchestration, and result reconstruction
- `phospy.io.bundles._signalome`: same concern split for signalome bundles
- `phospy.io.bundles._shared`: low-level JSON/path/table/coercion helpers used by
  both bundle domains

Compatibility behavior is isolated to signalome compatibility helpers:

- legacy config snapshot cutoff field parsing
- legacy module-assignment candidate/weight normalization

## Where Next

- CLI output paths and flags: [CLI Guide](cli.md)
- Request/result contract details: [API Guide](api.md)
- Validation and boundary guarantees: [Validation Guide](validation.md)
