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

Bundle loading is strict for v1 manifests:

- required manifest markers must be present (`bundle_type`,
  `manifest_version`, `table_format`, required sections)
- `provenance` must be an explicit object (not omitted and not `null`)
- optional table slots must still be declared in manifest table maps as a path
  or explicit `null`
- unsupported or incomplete manifest shapes fail with an explicit
  regenerate-bundle instruction
Published dataset export manifests (`dataset/manifest.json`) also keep this pair
explicit at top level:

- `intensity_scale` (numeric scale, for example `log2`)
- `quantitative_meaning` (scientific meaning, for example
  `phosphosite_log_abundance` or `phospho_total_log_ratio`)

`intensity_scale_state` persists both:

- numeric intensity scale (`linear`/`log2`)
- quantitative meaning (`phosphosite_abundance`,
  `phosphosite_log_abundance`, `phospho_total_log_ratio`, `unknown`)

Example: after subtractive total-protein correction, bundles keep
`intensity_scale_state.label == "log2"` and persist
`intensity_scale_state.quantity == "phospho_total_log_ratio"`.

For `processing_state.total_protein_correction`, persisted fields include:

- `policy`
- `applied`
- `formula`
- `requires_log_scale`
- `input_scale`
- `output_scale`
- `quantitative_meaning`
- `diagnostics` (runtime correction diagnostics, stored in a versioned schema)

`input_scale`/`output_scale` describe numeric scale only.
`quantitative_meaning` describes scientific interpretation.
For example, `log2` phosphosite abundance and `log2` phospho/total ratio are
different quantitative meanings.

`processing_state.total_protein_correction.diagnostics` now persists as a typed
versioned payload:

- `diagnostics_schema_version` (currently `1`)
- typed known fields when present:
  `policy`, `requested_policy`, `resolved_policy`, `formula`,
  `requires_log_scale`, `input_scale`, `output_scale`,
  `quantitative_meaning`, `matched_rows`,
  `identity_mode`, `phosphosite_key`, `total_protein_key`,
  `mapping_phosphosite_key`, `mapping_total_protein_key`,
  `mapping_table_fingerprint`, `duplicate_policy`, `unmatched_policy`,
  `phosphosite_row_count`, `total_protein_row_count`,
  `corrected_row_count`, `uncorrected_row_count`,
  `unused_total_protein_row_count`,
  `total_rows_used_by_multiple_phosphosites`,
  `unmatched_phosphosite_row_ids`, `unused_total_protein_row_ids`,
  `gene_symbol_matching_used`, `gene_symbol_identity_warning`,
  `total_table_hash`, `input_phospho_hash`, `output_phospho_hash`
- unknown top-level diagnostics fields are rejected at load/validation time
- unversioned diagnostics payloads are rejected at load/validation time with a
  schema-version error

Diagnostics values are validated on bundle save/load against the typed schema.
Malformed values (for example wrong scalar types or negative `matched_rows`) and
unknown top-level fields are rejected with explicit bundle field-path errors.

`provenance` is machine-readable run metadata (`RunProvenance`) and includes:

- environment versions (`phospy`, Python, dependency versions)
- input/output table fingerprints
- preprocessing stage execution trace (shape/hash/drop/imputation metadata)
- reference identity/fingerprints
- workflow parameters and random-state policy
  - signalome runs include `workflow_parameters.scale_guard` with
    `site_count`, `cluster_tree_backend`, `candidate_scoring_backend`,
    `candidate_scoring_requested_backend`,
    `max_exact_cluster_tree_sites`, `max_full_correlation_sites`,
    `exact_cluster_tree_built`, `candidate_scoring_mode`,
    `candidate_scoring_evaluated`, `candidate_scoring_skip_reason`,
    `candidate_scoring_sampling`, `candidate_scoring_applies_to`,
    `final_module_assignment_backend`,
    `final_module_assignment_uses_candidate_scoring`, and
    `scale_guard_passed`
  - `candidate_scoring_requested_backend` records the configured backend, while
    `candidate_scoring_mode` records what was actually evaluated
  - when sampled candidate scoring actually runs,
    `workflow_parameters.scale_guard.candidate_scoring_sampling` records
    `sampling_cap`, `sampling_method`, `deterministic_seed_policy`,
    `actual_sampled_pair_count`, and
    `per_cluster_sample_count_summary` (`min`, `max`, `mean`, `total`)
  - `candidate_scoring_applies_to` is
    `candidate_module_count_evaluation_only`; sampled candidate scoring does
    not make cluster-tree construction or final module assignment approximate
  - explicit `module_count` runs record `candidate_scoring_evaluated=false`
    with `candidate_scoring_skip_reason="explicit_module_count"` and
    `candidate_scoring_sampling=null`
  - signalome runs also include
    `workflow_parameters.signalome_score_semantics` with explicit scientific
    meaning fields (downstream score source/meaning, module-selection score
    meaning, candidate-scoring mode/scope, network policy and
    negative-correlation handling, missing/constant profile handling,
    thresholds/limits, clustering backend, and scientific interpretation limits)

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
- `signalome/expanded_signalome` is optional by contract, but
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
- `signalome_config.cluster_tree_backend`: cluster-tree construction backend
  (currently `"exact"`).
- `signalome_config.candidate_scoring_backend`: candidate-scoring backend
  (`"full"` or `"sampled"`).
- `signalome_config.max_exact_cluster_tree_sites`: hard guard limit for exact
  cluster-tree construction.
- `signalome_config.max_full_correlation_sites`: hard guard limit for full
  candidate-correlation scoring.
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

Manifest versioning starts at v1 so future format evolution is explicit.
## Internal Ownership

Public entrypoints stay at:

- `phospy.io.bundles.kinase`
- `phospy.io.bundles.signalome`

Internal implementation is split by concern ownership under:

- `phospy.io.bundles._kinase`: snapshots, manifest contract, writer orchestration,
  loader orchestration, and result reconstruction
- `phospy.io.bundles._signalome`: same concern split for signalome bundles
- `phospy.io.bundles._shared`: low-level JSON/path/table/coercion helpers used by
  both bundle domains

## Where Next

- CLI output paths and flags: [CLI Guide](cli.md)
- Request/result contract details: [API Guide](api.md)
- Validation and boundary guarantees: [Validation Guide](validation.md)
