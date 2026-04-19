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

The parity suite currently protects three rewrite-era slices:

- activity-stage outputs from fixed `predMat` + phospho inputs
- selected kinase-scoring/prediction points on the supported L6 lane
- selected signalome regression contracts on the supported L6 lane:
  `module_assignments`, `signalome_modules`, `kinase_network.nodes`,
  `kinase_network.edges`

Activity parity checks cover:

- `weighted_activity`
- `ksea_scores`
- `ksea_counts`
- `target_counts`
- `target_table`

## Legacy Donor Promotion Inventory

The rewrite suite now carries donor coverage visibility for the strongest legacy
science scenarios. `tests_legacy/` remains archival/provenance, not the primary
regression gate.

| Donor area | Science-gap ticket | Rewrite-owned blockers | Archival donors (historical only) |
| --- | --- | --- | --- |
| profile policy behavior | `SCI-GAP-01` | `tests/unit/test_legacy_donor_science.py::test_profile_policy_donor_locks_strict_median_behavior_and_contract_surface` | `tests_legacy/test_profiles.py::test_build_kinase_substrate_profiles_can_skip_missing_values_when_requested` |
| adaptive sampling / svm_mode | `SCI-GAP-05` | `tests/unit/test_legacy_donor_science.py::test_adaptive_sampling_donor_is_archival_and_svm_mode_is_not_rewrite_contract` | `tests_legacy/test_prediction.py::test_predict_accepts_explicit_r_parity_mode` |
| signalome clustering/module selection | `SCI-GAP-06` | `tests/unit/test_legacy_donor_science.py::test_signalome_clustering_donor_locks_rewrite_dominant_module_assignment_behavior`; `tests/parity/test_signalome_workflow_parity.py::test_signalome_module_assignments_match_selected_l6_regression_points` | `tests_legacy/test_signalomes.py::test_select_module_count_builds_one_cluster_tree_for_candidate_scoring` |
| weighted-top assignment behavior | `SCI-GAP-08` | `tests/unit/test_legacy_donor_science.py::test_weighted_top_assignment_donor_locks_fractional_metadata_and_non_fractional_module_selection`; `tests/parity/test_signalome_workflow_parity.py::test_signalome_module_assignments_match_selected_l6_regression_points` | `tests_legacy/test_signalomes.py::test_weighted_top_assignment_policy_propagates_fractional_module_shares` |
| network policy variants | `SCI-GAP-09` | `tests/unit/test_legacy_donor_science.py::test_network_policy_variant_donor_locks_signed_edges_and_narrow_config_surface`; `tests/parity/test_signalome_workflow_parity.py::test_signalome_network_edges_match_l6_fixture_pairs_and_sign_counts` | `tests_legacy/test_signalomes.py::test_build_kinase_network_policies_apply_expected_thresholding` |
| expanded signalome outputs | `SCI-GAP-10` | `tests/unit/test_legacy_donor_science.py::test_expanded_signalome_donor_locks_supported_lane_to_none_output`; `tests/integration/test_signalome_workflow_integration.py::test_signalome_workflow_runs_dataset_to_kinase_to_signalome_path` | `tests_legacy/test_signalomes.py::test_build_expanded_signalomes_uses_neighbor_map_and_preserves_site_order` |
| activity parity lock | `SCI-GAP-11` | `tests/unit/test_legacy_donor_science.py::test_activity_parity_lock_donor_uses_rewrite_owned_fixture_path`; `tests/parity/test_activity_stage_parity.py::test_weighted_activity_matches_legacy_reference_fixture` | `tests_legacy/test_activities.py::test_weighted_activity_matches_reference_on_large_sparse_input` |

Rewrite-side visibility check:

- `tests/unit/test_legacy_donor_inventory.py`

## Fixture Locations

### Rewrite-owned parity inputs and expectations

- `tests/fixtures/rewrite_parity/r_reference_l6/`
- provenance and promotion history:
  `tests/fixtures/rewrite_parity/r_reference_l6/PROVENANCE.md`
- `tests/fixtures/rewrite_parity/fragile_support_reference/`
- provenance and promotion history:
  `tests/fixtures/rewrite_parity/fragile_support_reference/PROVENANCE.md`
- `tests/fixtures/rewrite_parity/adaptive_sampling_edge/`
- provenance and promotion history:
  `tests/fixtures/rewrite_parity/adaptive_sampling_edge/PROVENANCE.md`

These files are the normal source for active parity tests in `tests/parity/`
and helpers in `tests/support/rewrite_fixture_data.py`.

### Rewrite workflow regression expectations

- `tests/fixtures/public_workflow_reference/signalome_rewrite_l6_*.csv`
- `tests/fixtures/public_workflow_reference/signalome_rewrite_l6_contract.json`

### Historical reference archive

- `tests_legacy/fixtures/` is retained for provenance and archival material.
- Active rewrite parity tests should not resolve fixtures from this path as their
  normal source.

## Run the Parity Suite

```bash
pytest -m parity
```

Useful variants:

```bash
pytest -m parity -rs
pytest -m parity -vv
pytest tests/parity/test_signalome_workflow_parity.py -vv
```
