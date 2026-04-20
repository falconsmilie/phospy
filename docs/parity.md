# Parity to PhosR

PhosPy parity is intentionally narrow and fixture-backed. The rewrite does not
claim full package equivalence with PhosR.

## What Parity Means Here

Parity in this repository is:

- seam-level
- selective
- tied to committed fixtures

Parity here does not mean:

- every PhosR feature is implemented
- every Python path must numerically match PhosR

## Active Parity Coverage

The parity suite currently protects rewrite-era parity families for:

- prediction-science parity on the fragile-support rewrite fixture lane
- kinase workflow parity on the supported L6 rewrite lane
- full promoted L6 downstream prediction/scoring parity against rewrite-owned
  promoted donor references (`profile_scores`, `combined_scores`, `weights`,
  candidate substrates, ranking/top-k summaries)
- adaptive prediction parity from promoted adaptive-sampling fixtures, executed
  in both supported rewrite policy lanes:
  `adaptive_policy="stable"` (default lane) and
  `adaptive_policy="r_parity"`
- adaptive replay-trace parity from promoted replay fixtures, including:
  initial negative-set replay surfaces, per-iteration sample-membership replay
  surfaces, final ensemble probabilities, top-k replay summaries, and
  stable-vs-r_parity comparison metrics under fixed seed
- public end-to-end `predMat` benchmark parity on the rewrite workflow path for
  both supported adaptive policies:
  `adaptive_policy="stable"` and `adaptive_policy="r_parity"`
- public end-to-end `predMat` order-invariance parity on the stable/default lane
  (normalized equality under reference-map order perturbation)
- activity-stage outputs from fixed `predMat` + phospho inputs
- selected signalome regression contracts on the supported L6 lane:
  `module_assignments`, `signalome_modules`, `kinase_network.nodes`,
  `kinase_network.edges`, `expanded_signalome` (selected AKT1 slice)

## Rewrite-owned parity reporting

Parity chatter is emitted by default from the active rewrite suite under
`tests/parity/`. The reporting layer is rewrite-owned (`tests/support/` +
`tests/conftest.py`) and is not routed through `tests_legacy/`.

When parity tests run, terminal output includes grouped scientific summaries for:

- prediction-science parity
- kinase workflow parity
- L6 full scoring/prediction parity
- adaptive prediction parity
- adaptive replay-trace parity
- public end-to-end predMat parity
- public predMat order-invariance parity
- activity-stage parity
- signalome workflow parity

No `PHOSPY_SHOW_*` environment variables are required.

Adaptive policy comparison is part of the active rewrite parity output:

- both supported rewrite policies execute in parity:
  `adaptive_policy="stable"` and `adaptive_policy="r_parity"`
- terminal chatter prints both lanes with clear policy labeling for review:
  `stable (default)` and `r_parity`
- side-by-side comparison metrics are printed in the adaptive parity section
- `svm_mode` remains archival naming and is not a rewrite public API field

Activity parity checks cover:

- `weighted_activity`
- `ksea_scores`
- `ksea_counts`
- `target_counts`
- `target_table`

Activity parity is a hard regression gate in rewrite CI:

- dedicated job: `activity-parity-gate`
- required marker selection: `parity and activity_parity`
- fixture source pinned to `tests/fixtures/rewrite_parity/r_reference_l6/`
  with provenance in
  `tests/fixtures/rewrite_parity/r_reference_l6/PROVENANCE.md`
- active parity assertions compare rewrite runtime outputs to committed
  rewrite-owned fixture expectations; no live `legacy_archive` execution is part
  of this gate

This lane is supported and parity-backed, not provisional.

## Legacy Donor Promotion Inventory

The rewrite suite now carries donor coverage visibility for the strongest legacy
science scenarios. `tests_legacy/` remains archival/provenance, not the primary
regression gate.

| Donor area | Science-gap ticket | Rewrite-owned blockers | Archival donors (historical only) |
| --- | --- | --- | --- |
| profile policy behavior | `SCI-GAP-01` | `tests/unit/test_legacy_donor_science.py::test_profile_policy_donor_locks_strict_median_behavior_and_contract_surface` | `tests_legacy/test_profiles.py::test_build_kinase_substrate_profiles_can_skip_missing_values_when_requested` |
| adaptive sampling / svm_mode | `SCI-GAP-05` | `tests/unit/test_legacy_donor_science.py::test_adaptive_sampling_donor_is_archival_and_svm_mode_is_not_rewrite_contract`; `tests/parity/test_adaptive_prediction_parity.py::test_adaptive_ensemble_outputs_match_promoted_fixture_tolerances` (active rewrite two-policy comparison via `adaptive_policy`) | `tests_legacy/test_prediction.py::test_predict_accepts_explicit_r_parity_mode`; `tests_legacy/test_prediction.py::test_resolve_prediction_sampling_policy_maps_public_modes` |
| signalome clustering/module selection | `SCI-GAP-06` | `tests/unit/test_legacy_donor_science.py::test_signalome_clustering_donor_locks_rewrite_dominant_module_assignment_behavior`; `tests/parity/test_signalome_workflow_parity.py::test_signalome_module_assignments_match_selected_l6_regression_points` | `tests_legacy/test_signalomes.py::test_select_module_count_builds_one_cluster_tree_for_candidate_scoring` |
| weighted-top assignment behavior | `SCI-GAP-08` | `tests/unit/test_legacy_donor_science.py::test_weighted_top_assignment_donor_locks_fractional_metadata_and_non_fractional_module_selection`; `tests/parity/test_signalome_workflow_parity.py::test_signalome_module_assignments_match_selected_l6_regression_points` | `tests_legacy/test_signalomes.py::test_weighted_top_assignment_policy_propagates_fractional_module_shares` |
| network policy variants | `SCI-GAP-09` | `tests/unit/test_legacy_donor_science.py::test_network_policy_variant_donor_locks_signed_edges_and_narrow_config_surface`; `tests/parity/test_signalome_workflow_parity.py::test_signalome_network_edges_match_l6_fixture_pairs_and_sign_counts`; `tests/parity/test_signalome_workflow_parity.py::test_signalome_network_policy_variants_match_fixed_matrix_expectations` | `tests_legacy/test_signalomes.py::test_build_kinase_network_policies_apply_expected_thresholding` |
| expanded signalome outputs | `SCI-GAP-10` | `tests/unit/test_legacy_donor_science.py::test_expanded_signalome_donor_locks_supported_lane_to_materialized_output`; `tests/parity/test_signalome_workflow_parity.py::test_signalome_expanded_slice_matches_l6_selected_akt1_fixture`; `tests/integration/test_signalome_workflow_integration.py::test_signalome_workflow_runs_dataset_to_kinase_to_signalome_path` | `tests_legacy/test_signalomes.py::test_build_expanded_signalomes_uses_neighbor_map_and_preserves_site_order`; `tests_legacy/test_signalomes.py::test_signalome_result_expanded_signalomes_materialize_with_parity` |
| activity parity lock | `SCI-GAP-11` | `tests/unit/test_legacy_donor_science.py::test_activity_parity_lock_donor_uses_rewrite_owned_fixture_path`; `tests/parity/test_activity_stage_parity.py::test_weighted_activity_matches_rewrite_reference_fixture` | `tests_legacy/test_activities.py::test_weighted_activity_matches_reference_on_large_sparse_input` |

Rewrite-side visibility check:

- `tests/unit/test_legacy_donor_inventory.py`

## Fixture Locations

### Rewrite-owned parity inputs and expectations

- `tests/fixtures/rewrite_parity/r_reference_l6/`
- provenance and promotion history:
  `tests/fixtures/rewrite_parity/r_reference_l6/PROVENANCE.md`
- `tests/fixtures/rewrite_parity/r_reference_l6_prediction/`
- provenance and promotion history:
  `tests/fixtures/rewrite_parity/r_reference_l6_prediction/PROVENANCE.md`
- `tests/fixtures/rewrite_parity/fragile_support_reference/`
- provenance and promotion history:
  `tests/fixtures/rewrite_parity/fragile_support_reference/PROVENANCE.md`
- `tests/fixtures/rewrite_parity/adaptive_sampling_edge/`
- provenance and promotion history:
  `tests/fixtures/rewrite_parity/adaptive_sampling_edge/PROVENANCE.md`
- `tests/fixtures/rewrite_parity/adaptive_sampling_replay/`
- provenance and promotion history:
  `tests/fixtures/rewrite_parity/adaptive_sampling_replay/PROVENANCE.md`

These files are the normal source for active parity tests in `tests/parity/`
and helpers in `tests/support/rewrite_fixture_data.py`.

### Rewrite workflow regression expectations

- public predMat benchmark inputs and committed rewrite outputs:
  `tests/fixtures/public_workflow_reference/predmat_input_*.{csv,json}`,
  `tests/fixtures/public_workflow_reference/predmat_rewrite_*.csv`,
  `tests/fixtures/public_workflow_reference/predmat_rewrite_contract.json`
- provenance and promotion history:
  `tests/fixtures/public_workflow_reference/PROVENANCE.md`
- `tests/fixtures/public_workflow_reference/signalome_rewrite_l6_*.csv`
- `tests/fixtures/public_workflow_reference/signalome_rewrite_l6_contract.json`

### Historical reference archive

- `tests_legacy/fixtures/` is retained for provenance and archival material.
- Active rewrite parity tests should not resolve fixtures from this path as their
  normal source.
- `tests_legacy/test_parity-with_metrics.py` is archival/provenance only and is
  not active reporting infrastructure for rewrite parity runs.

## Run the Parity Suite

```bash
pytest tests/parity -m parity -s
```

or:

```bash
make test-parity
```

Useful variants:

```bash
pytest tests/parity -m parity -rs -s
pytest -m parity -vv
pytest tests/parity/test_activity_stage_parity.py -m "parity and activity_parity" -vv
pytest tests/parity/test_signalome_workflow_parity.py -vv
```
