from __future__ import annotations

import pandas as pd
import pytest

import phospy.workflows.kinase.executor as kinase_executor
from phospy import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    ReferencePreset,
)
from phospy.references.resolution import ReferenceResolver
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
            ),
            prediction_config=KinasePredictionConfig(
                top_k=6,
                ensemble_size=8,
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
    deterministic = _run_workflow(
        dataset=dataset,
        references=references,
        mode="deterministic_ranking",
    )
    adaptive = _run_workflow(
        dataset=dataset,
        references=references,
        mode="adaptive_ensemble",
    )

    pd.testing.assert_frame_equal(
        deterministic.scoring_result.profile_scores,
        adaptive.scoring_result.profile_scores,
        check_dtype=False,
    )
    assert deterministic.scoring_result.combined_scores is not None
    assert adaptive.scoring_result.combined_scores is not None
    pd.testing.assert_frame_equal(
        deterministic.scoring_result.combined_scores,
        adaptive.scoring_result.combined_scores,
        check_dtype=False,
    )
    assert deterministic.scoring_result.motif_scores is not None
    assert adaptive.scoring_result.motif_scores is not None
    pd.testing.assert_frame_equal(
        deterministic.scoring_result.motif_scores,
        adaptive.scoring_result.motif_scores,
        check_dtype=False,
    )
    assert deterministic.scoring_result.weights is not None
    assert adaptive.scoring_result.weights is not None
    pd.testing.assert_frame_equal(
        deterministic.scoring_result.weights,
        adaptive.scoring_result.weights,
        check_dtype=False,
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
    assert from_preset.scoring_result.combined_scores is not None
    assert from_bundle.scoring_result.combined_scores is not None
    pd.testing.assert_frame_equal(
        from_preset.scoring_result.combined_scores,
        from_bundle.scoring_result.combined_scores,
        check_dtype=False,
    )
    assert from_preset.scoring_result.motif_scores is not None
    assert from_bundle.scoring_result.motif_scores is not None
    pd.testing.assert_frame_equal(
        from_preset.scoring_result.motif_scores,
        from_bundle.scoring_result.motif_scores,
        check_dtype=False,
    )
    assert from_preset.scoring_result.weights is not None
    assert from_bundle.scoring_result.weights is not None
    pd.testing.assert_frame_equal(
        from_preset.scoring_result.weights,
        from_bundle.scoring_result.weights,
        check_dtype=False,
    )


def test_prediction_stage_consumes_authoritative_combined_scoring_outputs_across_modes_and_reference_forms(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset = build_rat_l6_dataset(n_sites=220)
    explicit_bundle = _resolved_bundle_for_dataset(dataset)
    captured_sources: list[str] = []
    captured_selected_scores: list[pd.DataFrame] = []
    captured_prediction_inputs: list[pd.DataFrame] = []

    original_select_downstream = kinase_executor.select_downstream_score_matrix
    original_build_candidates = kinase_executor.build_candidate_substrate_list

    def _capture_selected(*, profile_scores, combined_scores):
        selected, source = original_select_downstream(
            profile_scores=profile_scores,
            combined_scores=combined_scores,
        )
        captured_selected_scores.append(selected.copy(deep=True))
        captured_sources.append(source)
        return selected, source

    def _capture_prediction_input(*, scores, top, score_threshold, inclusion):
        captured_prediction_inputs.append(scores.copy(deep=True))
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
        assert captured_sources[-1] == "combined_scores"
        assert result.scoring_result.combined_scores is not None
        pd.testing.assert_frame_equal(
            captured_selected_scores[-1],
            result.scoring_result.combined_scores,
            check_dtype=False,
        )
        pd.testing.assert_frame_equal(
            captured_prediction_inputs[-1],
            result.scoring_result.combined_scores,
            check_dtype=False,
        )

    assert len(captured_sources) == len(runs)
    assert set(captured_sources) == {"combined_scores"}
