from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.signalomes.clustering import (
    cluster_sites_with_diagnostics,
    fit_cluster_labels,
    select_module_count_with_diagnostics,
)


def test_module_count_selection_clamps_explicit_request_to_site_count() -> None:
    scoring_values = np.asarray(
        [
            [1.0, 2.0, 3.0],
            [1.1, 2.1, 3.1],
            [0.9, 1.9, 2.9],
        ],
        dtype=float,
    )

    diagnostics = select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=5,
    )

    assert diagnostics.strategy == "explicit_module_count"
    assert diagnostics.selected_module_count == 3
    assert diagnostics.requested_module_count == 5
    assert diagnostics.reason == "module_count was provided explicitly by the caller"
    assert diagnostics.candidate_scores == {}


def test_module_count_automatic_selection_surfaces_stable_diagnostics() -> None:
    scoring_values = np.asarray(
        [
            [1.0, 1.1, 1.2],
            [1.1, 1.0, 1.2],
            [0.9, 1.2, 1.0],
            [-1.0, -1.1, -1.2],
            [-1.1, -1.0, -1.2],
            [-0.9, -1.2, -1.0],
        ],
        dtype=float,
    )

    diagnostics = select_module_count_with_diagnostics(scoring_values=scoring_values)

    assert diagnostics.used_automatic_selection
    assert diagnostics.strategy == "correlation_thresholds"
    assert diagnostics.selected_module_count == 2
    assert diagnostics.threshold_used in {0.5, 0.1}
    assert diagnostics.max_clusters_evaluated >= 2
    assert diagnostics.reason
    assert 2 in diagnostics.candidate_scores


def test_cluster_sites_matches_two_pass_partition_for_selected_module_count() -> None:
    scoring_matrix = pd.DataFrame(
        {
            "K1": [1.0, 2.0, 3.0, -1.0, -2.0, -3.0],
            "K2": [2.0, 4.0, 6.0, -2.0, -4.0, -6.0],
            "K3": [3.0, 6.0, 9.0, -3.0, -6.0, -9.0],
        },
        index=[f"P{index};S{index};" for index in range(1, 7)],
    )
    scoring_values = scoring_matrix.to_numpy(dtype=float)
    diagnostics = select_module_count_with_diagnostics(scoring_values=scoring_values)
    module_count = max(
        1, min(diagnostics.selected_module_count, scoring_values.shape[0])
    )
    if module_count == 1:
        baseline_labels = np.ones(scoring_values.shape[0], dtype=int)
    else:
        baseline_labels = fit_cluster_labels(scoring_values, module_count) + 1

    clustered = cluster_sites_with_diagnostics(
        scoring_matrix=scoring_matrix,
        requested_module_count=None,
    )
    observed_labels = clustered.site_clusters.to_numpy(dtype=int, copy=False)
    observed_partition = observed_labels[:, None] == observed_labels[None, :]
    baseline_partition = baseline_labels[:, None] == baseline_labels[None, :]

    assert clustered.module_selection_diagnostics == diagnostics
    assert np.array_equal(observed_partition, baseline_partition)
