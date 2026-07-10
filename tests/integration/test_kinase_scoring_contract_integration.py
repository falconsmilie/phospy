from __future__ import annotations

import pandas as pd
import pytest

import phospy.workflows.kinase.executor as kinase_executor
from phospy import (
    KinaseWorkflow,
)
from phospy.api import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    ReferenceContextCompatibilityPolicy,
    ReferencePreset,
)
from phospy.api.configs import (
    KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
    KINASE_PREDICTION_MODES,
)
from phospy.science.references.resolution import ReferenceResolver
from tests.support.rewrite_fixture_data import build_rat_l6_dataset

pytestmark = pytest.mark.integration


def _resolved_bundle_for_dataset(dataset):
    return ReferenceResolver().run(
        ReferencePreset.AUTO,
        dataset_organism=dataset.organism,
    )


def _run_workflow(*, dataset, references, mode: str):
    return KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(
                min_substrates=2,
                include_diagnostic_scoring_tables=True,
                reference_context_compatibility_policy=(
                    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                ),
            ),
            prediction_config=KinasePredictionConfig(
                top_k=6,
                deterministic_max_selected_kinases=8,
                adaptive_ensemble_runs=8,
                mode=mode,
                n_iterations=2,
                random_state=19,
            ),
            activity_config=None,
        )
    )


def test_scoring_stage_is_prediction_mode_invariant_for_supported_lane() -> None:
    dataset = build_rat_l6_dataset(n_sites=220)
    references = _resolved_bundle_for_dataset(dataset)
    baseline = _run_workflow(
        dataset=dataset,
        references=references,
        mode=KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
    )

    for mode in sorted(KINASE_PREDICTION_MODES):
        result = _run_workflow(
            dataset=dataset,
            references=references,
            mode=mode,
        )

        pd.testing.assert_frame_equal(
            result.scoring_result.profile_scores,
            baseline.scoring_result.profile_scores,
        )
        assert result.scoring_result.rank_weighted_fusion_scores is not None
        assert baseline.scoring_result.rank_weighted_fusion_scores is not None
        pd.testing.assert_frame_equal(
            result.scoring_result.rank_weighted_fusion_scores,
            baseline.scoring_result.rank_weighted_fusion_scores,
        )
        assert result.scoring_result.motif_scores is not None
        assert baseline.scoring_result.motif_scores is not None
        pd.testing.assert_frame_equal(
            result.scoring_result.motif_scores,
            baseline.scoring_result.motif_scores,
        )
        assert result.scoring_result.score_fusion_weights is not None
        assert baseline.scoring_result.score_fusion_weights is not None
        pd.testing.assert_frame_equal(
            result.scoring_result.score_fusion_weights,
            baseline.scoring_result.score_fusion_weights,
        )


def test_scoring_stage_is_reference_input_form_invariant_for_equivalent_content() -> (
    None
):
    dataset = build_rat_l6_dataset(n_sites=220)
    explicit_bundle = _resolved_bundle_for_dataset(dataset)
    from_preset = _run_workflow(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        mode="deterministic_ranking",
    )
    from_bundle = _run_workflow(
        dataset=dataset,
        references=explicit_bundle,
        mode="deterministic_ranking",
    )

    pd.testing.assert_frame_equal(
        from_preset.scoring_result.profile_scores,
        from_bundle.scoring_result.profile_scores,
        check_dtype=False,
    )
    assert from_preset.scoring_result.rank_weighted_fusion_scores is not None
    assert from_bundle.scoring_result.rank_weighted_fusion_scores is not None
    pd.testing.assert_frame_equal(
        from_preset.scoring_result.rank_weighted_fusion_scores,
        from_bundle.scoring_result.rank_weighted_fusion_scores,
        check_dtype=False,
    )
    assert from_preset.scoring_result.motif_scores is not None
    assert from_bundle.scoring_result.motif_scores is not None
    pd.testing.assert_frame_equal(
        from_preset.scoring_result.motif_scores,
        from_bundle.scoring_result.motif_scores,
        check_dtype=False,
    )
    assert from_preset.scoring_result.score_fusion_weights is not None
    assert from_bundle.scoring_result.score_fusion_weights is not None
    pd.testing.assert_frame_equal(
        from_preset.scoring_result.score_fusion_weights,
        from_bundle.scoring_result.score_fusion_weights,
        check_dtype=False,
    )


def test_supported_lane_is_reference_input_form_invariant_for_equivalent_content(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset = build_rat_l6_dataset(n_sites=220)
    explicit_bundle = _resolved_bundle_for_dataset(dataset)
    captured_sources: list[str] = []
    original_select_downstream = kinase_executor.select_downstream_score_matrix

    def _capture_selected(*, profile_scores, rank_weighted_fusion_scores):
        selected, source = original_select_downstream(
            profile_scores=profile_scores,
            rank_weighted_fusion_scores=rank_weighted_fusion_scores,
        )
        captured_sources.append(source)
        return selected, source

    monkeypatch.setattr(
        kinase_executor,
        "select_downstream_score_matrix",
        _capture_selected,
    )

    from_preset = _run_workflow(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        mode="deterministic_ranking",
    )
    from_bundle = _run_workflow(
        dataset=dataset,
        references=explicit_bundle,
        mode="deterministic_ranking",
    )

    pd.testing.assert_frame_equal(
        from_preset.scoring_result.profile_scores,
        from_bundle.scoring_result.profile_scores,
        check_dtype=False,
    )
    assert from_preset.scoring_result.rank_weighted_fusion_scores is not None
    assert from_bundle.scoring_result.rank_weighted_fusion_scores is not None
    pd.testing.assert_frame_equal(
        from_preset.scoring_result.rank_weighted_fusion_scores,
        from_bundle.scoring_result.rank_weighted_fusion_scores,
        check_dtype=False,
    )
    assert captured_sources == [
        "rank_weighted_fusion_scores",
        "rank_weighted_fusion_scores",
    ]
    pd.testing.assert_frame_equal(
        from_preset.prediction_result.pred_mat,
        from_bundle.prediction_result.pred_mat,
        check_dtype=False,
    )


def test_prediction_stage_consumes_authoritative_rank_weighted_fusion_outputs_across_modes_and_reference_forms(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset = build_rat_l6_dataset(n_sites=220)
    explicit_bundle = _resolved_bundle_for_dataset(dataset)
    captured_sources: list[str] = []
    captured_selected_scores: list[pd.DataFrame] = []
    captured_prediction_inputs: list[pd.DataFrame] = []
    captured_candidate_filters: list[tuple[int, float, int]] = []

    original_select_downstream = kinase_executor.select_downstream_score_matrix
    original_build_candidates = kinase_executor.build_candidate_substrate_list

    def _capture_selected(*, profile_scores, rank_weighted_fusion_scores):
        selected, source = original_select_downstream(
            profile_scores=profile_scores,
            rank_weighted_fusion_scores=rank_weighted_fusion_scores,
        )
        captured_selected_scores.append(selected.copy(deep=True))
        captured_sources.append(source)
        return selected, source

    def _capture_prediction_input(*, scores, top, score_threshold, inclusion):
        captured_prediction_inputs.append(scores.copy(deep=True))
        captured_candidate_filters.append(
            (int(top), float(score_threshold), int(inclusion))
        )
        return original_build_candidates(
            scores=scores,
            top=top,
            score_threshold=score_threshold,
            inclusion=inclusion,
        )

    monkeypatch.setattr(
        kinase_executor,
        "select_downstream_score_matrix",
        _capture_selected,
    )
    monkeypatch.setattr(
        kinase_executor,
        "build_candidate_substrate_list",
        _capture_prediction_input,
    )

    runs = [
        ("deterministic_ranking", ReferencePreset.AUTO),
        ("adaptive_ensemble", ReferencePreset.AUTO),
        ("deterministic_ranking", explicit_bundle),
        ("adaptive_ensemble", explicit_bundle),
    ]
    for mode, references in runs:
        result = _run_workflow(
            dataset=dataset,
            references=references,
            mode=mode,
        )
        assert captured_sources[-1] == "rank_weighted_fusion_scores"
        assert result.scoring_result.rank_weighted_fusion_scores is not None
        pd.testing.assert_frame_equal(
            captured_selected_scores[-1],
            result.scoring_result.rank_weighted_fusion_scores,
            check_dtype=False,
        )
        pd.testing.assert_frame_equal(
            captured_prediction_inputs[-1],
            result.scoring_result.rank_weighted_fusion_scores,
            check_dtype=False,
        )
        assert captured_candidate_filters[-1] == (6, 0.0, 1)

    assert len(captured_sources) == len(runs)
    assert set(captured_sources) == {"rank_weighted_fusion_scores"}
    assert len(captured_candidate_filters) == len(runs)
    assert {
        (score_threshold, inclusion)
        for _top, score_threshold, inclusion in captured_candidate_filters
    } == {(0.0, 1)}


def test_supported_lane_rank_weighted_fusion_scoring_always_enables_profile_only_fallback(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset = build_rat_l6_dataset(n_sites=220)
    references = _resolved_bundle_for_dataset(dataset)
    captured_fallback_flags: list[bool] = []
    original_combine = kinase_executor.fuse_profile_and_motif_scores_by_rank_weight

    def _capture_combine(**kwargs):
        captured_fallback_flags.append(bool(kwargs["allow_profile_only_fallback"]))
        return original_combine(**kwargs)

    monkeypatch.setattr(
        kinase_executor,
        "fuse_profile_and_motif_scores_by_rank_weight",
        _capture_combine,
    )

    _run_workflow(
        dataset=dataset,
        references=references,
        mode="deterministic_ranking",
    )

    assert captured_fallback_flags == [True]
