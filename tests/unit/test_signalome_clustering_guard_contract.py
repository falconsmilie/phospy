from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import phospy.signalomes.clustering as clustering_module
from phospy.errors.workflows import SignalomeScaleError
from phospy.signalomes.clustering import (
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    cluster_sites,
    cluster_sites_with_diagnostics,
    fit_cluster_labels,
    select_module_count_with_diagnostics,
)
from phospy.signalomes.clustering import exact_python as exact_clustering


def _over_limit_scoring_values() -> np.ndarray:
    site_count = _over_limit_site_count()
    return np.column_stack(
        (
            np.linspace(0.0, 1.0, site_count, dtype=float),
            np.linspace(1.0, 0.0, site_count, dtype=float),
        )
    )


def _over_limit_site_count() -> int:
    return clustering_module.MAX_FULL_CORRELATION_SITE_COUNT + 1


def _over_limit_scoring_matrix() -> pd.DataFrame:
    values = _over_limit_scoring_values()
    return pd.DataFrame(
        values,
        index=[f"P{index};S{index};" for index in range(values.shape[0])],
        columns=["K1", "K2"],
    )


def _small_scoring_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        [[0.95, 0.1], [0.1, 0.95], [0.8, 0.7]],
        index=["P1;S1;", "P2;S2;", "P3;S3;"],
        columns=["K1", "K2"],
        dtype=float,
    )


def _patch_tree_builder(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    tree_calls: list[str] = []

    def _build_tree_should_not_run(scoring_values: object) -> object:
        del scoring_values
        tree_calls.append("called")
        raise AssertionError(
            "_build_cluster_tree should not run when exact guard fails"
        )

    monkeypatch.setattr(
        exact_clustering,
        "_build_cluster_tree",
        _build_tree_should_not_run,
    )
    return tree_calls


def _assert_exact_tree_guard_message(
    message: str,
    *,
    expected_site_count: int,
    expected_limit: int | None = None,
) -> None:
    resolved_limit = (
        clustering_module.MAX_FULL_CORRELATION_SITE_COUNT
        if expected_limit is None
        else int(expected_limit)
    )
    resolved_limit_token = f"{resolved_limit:,}"
    message_lower = message.lower()
    assert "exact cluster-tree construction" in message_lower
    assert f"{expected_site_count:,} sites" in message_lower
    assert f"max_exact_tree_sites={resolved_limit_token}" in message_lower
    assert "tree_engine='exact'" in message_lower


def test_cluster_sites_missing_exact_guard_arg_fails_over_limit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tree_calls = _patch_tree_builder(monkeypatch)
    site_count = _over_limit_site_count()

    with pytest.raises(SignalomeScaleError) as exc_info:
        cluster_sites(
            scoring_matrix=_over_limit_scoring_matrix(),
            requested_module_count=None,
        )

    _assert_exact_tree_guard_message(
        str(exc_info.value),
        expected_site_count=site_count,
    )
    assert tree_calls == []


def test_cluster_sites_with_diagnostics_explicit_none_guard_fails_over_limit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tree_calls = _patch_tree_builder(monkeypatch)
    site_count = _over_limit_site_count()

    with pytest.raises(SignalomeScaleError) as exc_info:
        cluster_sites_with_diagnostics(
            scoring_matrix=_over_limit_scoring_matrix(),
            requested_module_count=2,
            candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
            max_exact_tree_sites=None,
        )

    message = str(exc_info.value).lower()
    _assert_exact_tree_guard_message(
        message,
        expected_site_count=site_count,
    )
    assert "candidate_scoring_policy='sampled'" in message
    assert tree_calls == []


def test_full_candidate_scoring_cannot_bypass_exact_tree_guard(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tree_calls = _patch_tree_builder(monkeypatch)
    site_count = _over_limit_site_count()

    with pytest.raises(SignalomeScaleError) as exc_info:
        cluster_sites_with_diagnostics(
            scoring_matrix=_over_limit_scoring_matrix(),
            requested_module_count=None,
            candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
            max_exact_tree_sites=None,
        )

    message = str(exc_info.value).lower()
    _assert_exact_tree_guard_message(
        message,
        expected_site_count=site_count,
    )
    assert "candidate_scoring_policy='full'" in message
    assert tree_calls == []


def test_full_candidate_scoring_over_full_limit_fails_before_tree_construction(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    site_count = int(_small_scoring_matrix().shape[0])
    tree_calls: list[str] = []

    def _build_tree_should_not_run(scoring_values: object) -> object:
        del scoring_values
        tree_calls.append("called")
        raise AssertionError(
            "_build_cluster_tree should not run when full-correlation guard fails"
        )

    monkeypatch.setattr(
        exact_clustering,
        "_build_cluster_tree",
        _build_tree_should_not_run,
    )

    with pytest.raises(SignalomeScaleError) as exc_info:
        cluster_sites_with_diagnostics(
            scoring_matrix=_small_scoring_matrix(),
            requested_module_count=None,
            candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
            max_exact_tree_sites=10,
            max_full_candidate_scoring_sites=2,
        )

    message = str(exc_info.value).lower()
    assert "full candidate-correlation scoring would evaluate" in message
    assert f"{site_count:,} sites" in message
    assert "max_full_candidate_scoring_sites=2" in message
    assert "exact cluster-tree construction has not been attempted" in message
    assert "use candidate_scoring_policy='sampled'" in message
    assert "candidate module-count evaluation" in message
    assert tree_calls == []


def test_both_limits_exceeded_uses_exact_tree_guard_first_without_building_tree(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tree_calls: list[str] = []

    def _build_tree_should_not_run(scoring_values: object) -> object:
        del scoring_values
        tree_calls.append("called")
        raise AssertionError(
            "_build_cluster_tree should not run when both scale guards fail"
        )

    monkeypatch.setattr(
        exact_clustering,
        "_build_cluster_tree",
        _build_tree_should_not_run,
    )

    scoring_matrix = _small_scoring_matrix()
    site_count = int(scoring_matrix.shape[0])
    with pytest.raises(SignalomeScaleError) as exc_info:
        cluster_sites_with_diagnostics(
            scoring_matrix=scoring_matrix,
            requested_module_count=None,
            candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
            max_exact_tree_sites=2,
            max_full_candidate_scoring_sites=2,
        )

    message = str(exc_info.value).lower()
    _assert_exact_tree_guard_message(
        message,
        expected_site_count=site_count,
        expected_limit=2,
    )
    assert "candidate_scoring_policy='full'" in message
    assert "full candidate-correlation scoring would evaluate" not in message
    assert tree_calls == []


def test_sampled_candidate_scoring_over_full_limit_does_not_use_full_guard() -> None:
    clustered = cluster_sites_with_diagnostics(
        scoring_matrix=_small_scoring_matrix(),
        requested_module_count=None,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        max_exact_tree_sites=10,
        max_full_candidate_scoring_sites=2,
    )

    assert (
        clustered.candidate_scoring_mode == SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
    )


def test_select_module_count_with_diagnostics_missing_guard_arg_fails_over_limit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tree_calls = _patch_tree_builder(monkeypatch)
    site_count = _over_limit_site_count()

    with pytest.raises(SignalomeScaleError) as exc_info:
        select_module_count_with_diagnostics(
            scoring_values=_over_limit_scoring_values()
        )

    _assert_exact_tree_guard_message(
        str(exc_info.value),
        expected_site_count=site_count,
    )
    assert tree_calls == []


def test_fit_cluster_labels_explicit_none_guard_fails_over_limit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    tree_calls = _patch_tree_builder(monkeypatch)
    site_count = _over_limit_site_count()

    with pytest.raises(SignalomeScaleError) as exc_info:
        fit_cluster_labels(
            _over_limit_scoring_values(),
            cluster_count=2,
            max_exact_tree_sites=None,
        )

    _assert_exact_tree_guard_message(
        str(exc_info.value),
        expected_site_count=site_count,
    )
    assert tree_calls == []
