from __future__ import annotations

import pandas as pd
import pytest

from phospy.motifs import KinaseMotifScorer
from phospy.validation.errors import InputCompatibilityError, RequestValidationError
from phospy.validation.requests import KinaseWorkflowRequest
from phospy.validation.tables import SiteMatrixSchema
from phospy.workflow import KinaseWorkflow, WorkflowExecutionPlanner


def make_workflow_inputs() -> tuple[
    pd.DataFrame,
    dict[str, list[str]],
    dict[str, str],
    dict[str, list[str]],
]:
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 1.1, 0.9, 1.2, 3.0, 2.9, 3.1, 2.8],
            "sample_2": [2.0, 2.1, 1.9, 2.2, 2.0, 2.1, 1.9, 2.2],
            "sample_3": [3.0, 3.1, 2.9, 3.2, 1.0, 1.1, 0.9, 1.2],
        },
        index=[f"SITE_{i}" for i in range(1, 9)],
    )
    substrate_map = {
        "KINASE_A": ["SITE_1", "SITE_2", "SITE_3", "SITE_4"],
        "KINASE_B": ["SITE_5", "SITE_6", "SITE_7", "SITE_8"],
    }
    site_sequences = {
        "SITE_1": "QQAAAAAYY",
        "SITE_2": "QQAAAAAYY",
        "SITE_3": "QQAAAAAYY",
        "SITE_4": "QQAAAAAYY",
        "SITE_5": "QQTTTTTYY",
        "SITE_6": "QQTTTTTYY",
        "SITE_7": "QQTTTTTYY",
        "SITE_8": "QQTTTTTYY",
    }
    motif_sequences = {
        "KINASE_A": ["QQAAAAAYY", "QQAAAAAYY", "QQAAAAAYY"],
        "KINASE_B": ["QQTTTTTYY", "QQTTTTTYY", "QQTTTTTYY"],
    }
    return phospho_matrix, substrate_map, site_sequences, motif_sequences


def test_kinase_workflow_runs_native_end_to_end_path() -> None:
    phospho_matrix, substrate_map, site_sequences, motif_sequences = (
        make_workflow_inputs()
    )

    workflow = KinaseWorkflow(flank_size=2)

    result = workflow.run(
        phospho_matrix=phospho_matrix,
        substrate_map=substrate_map,
        site_sequences=site_sequences,
        motif_sequences=motif_sequences,
        min_substrates=2,
        min_motif_size=2,
        ensemble_size=3,
        top=4,
        score_threshold=0.75,
        inclusion=3,
        n_iterations=2,
        random_state=17,
    )

    assert list(result.profile_result.profile_matrix.index) == ["KINASE_A", "KINASE_B"]
    assert list(result.motif_result.motif_scores.columns) == ["KINASE_A", "KINASE_B"]
    assert list(result.scoring_result.combined_scores.columns) == [
        "KINASE_A",
        "KINASE_B",
    ]
    assert list(result.prediction_result.pred_matrix.columns) == [
        "KINASE_A",
        "KINASE_B",
    ]
    assert set(result.prediction_result.substrate_list) == {"KINASE_A", "KINASE_B"}
    assert (
        result.prediction_result.pred_matrix.loc[
            ["SITE_1", "SITE_2", "SITE_3", "SITE_4"], "KINASE_A"
        ].mean()
        > result.prediction_result.pred_matrix.loc[
            ["SITE_5", "SITE_6", "SITE_7", "SITE_8"], "KINASE_A"
        ].mean()
    )


def test_kinase_workflow_can_fall_back_to_profile_only_prediction() -> None:
    phospho_matrix, substrate_map, _, _ = make_workflow_inputs()
    workflow = KinaseWorkflow()

    result = workflow.run(
        phospho_matrix=phospho_matrix,
        substrate_map=substrate_map,
        motif_sequences=None,
        allow_profile_only_fallback=True,
        ensemble_size=3,
        top=4,
        score_threshold=0.75,
        inclusion=3,
        n_iterations=2,
        random_state=23,
    )

    assert result.motif_result is None
    assert result.scoring_result.combined_scores is None
    assert list(result.prediction_result.pred_matrix.columns) == [
        "KINASE_A",
        "KINASE_B",
    ]


def test_kinase_workflow_requires_motif_sequences_without_profile_fallback() -> None:
    phospho_matrix, substrate_map, _, _ = make_workflow_inputs()
    workflow = KinaseWorkflow()

    with pytest.raises(RequestValidationError, match="motif_sequences are required"):
        workflow.run(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
        )


def test_kinase_workflow_rejects_empty_substrate_map() -> None:
    phospho_matrix, _, _, _ = make_workflow_inputs()
    workflow = KinaseWorkflow()

    with pytest.raises(RequestValidationError, match="substrate_map must not be empty"):
        workflow.run(
            phospho_matrix=phospho_matrix,
            substrate_map={},
            motif_sequences=None,
            allow_profile_only_fallback=True,
        )


def test_kinase_workflow_rejects_empty_motif_sequences_mapping() -> None:
    phospho_matrix, substrate_map, site_sequences, _ = make_workflow_inputs()
    workflow = KinaseWorkflow()

    with pytest.raises(
        RequestValidationError, match="motif_sequences must not be empty"
    ):
        workflow.run(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            motif_sequences={},
        )


def test_kinase_workflow_requires_site_sequences_when_motifs_are_provided() -> None:
    phospho_matrix, substrate_map, _, motif_sequences = make_workflow_inputs()
    workflow = KinaseWorkflow()

    with pytest.raises(RequestValidationError, match="site_sequences are required"):
        workflow.run(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=None,
            motif_sequences=motif_sequences,
        )


def test_kinase_workflow_accepts_explicit_svm_mode() -> None:
    phospho_matrix, substrate_map, site_sequences, motif_sequences = (
        make_workflow_inputs()
    )

    workflow = KinaseWorkflow(svm_mode="r_parity")

    result = workflow.run(
        phospho_matrix=phospho_matrix,
        substrate_map=substrate_map,
        site_sequences=site_sequences,
        motif_sequences=motif_sequences,
        min_substrates=2,
        min_motif_size=2,
        ensemble_size=2,
        top=4,
        score_threshold=0.75,
        inclusion=3,
        n_iterations=2,
        random_state=17,
        svm_mode="r_parity",
    )

    assert list(result.prediction_result.pred_matrix.columns) == [
        "KINASE_A",
        "KINASE_B",
    ]


def test_kinase_workflow_run_request_validates_boundary_once(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    phospho_matrix, substrate_map, site_sequences, motif_sequences = (
        make_workflow_inputs()
    )
    workflow = KinaseWorkflow(flank_size=2)
    request = KinaseWorkflowRequest.validate_request(
        phospho_matrix=phospho_matrix,
        substrate_map=substrate_map,
        site_sequences=site_sequences,
        motif_sequences=motif_sequences,
        min_substrates=2,
        min_motif_size=2,
        ensemble_size=3,
        top=4,
        score_threshold=0.75,
        inclusion=3,
        n_iterations=2,
        random_state=17,
    )

    matrix_calls: list[str] = []
    motif_calls: list[int] = []
    original_matrix_validate = SiteMatrixSchema.validate
    original_from_sequences = KinaseMotifScorer.from_substrate_sequences

    def counting_matrix_validate(df: pd.DataFrame, *, context: str) -> pd.DataFrame:
        matrix_calls.append(context)
        return original_matrix_validate(df, context=context)

    def counting_from_sequences(
        cls,
        motif_sequences: dict[str, list[str]],
        flank_size: int = 7,
    ) -> KinaseMotifScorer:
        motif_calls.append(flank_size)
        return original_from_sequences(
            motif_sequences=motif_sequences, flank_size=flank_size
        )

    monkeypatch.setattr(
        SiteMatrixSchema,
        "validate",
        staticmethod(counting_matrix_validate),
    )
    monkeypatch.setattr(
        KinaseMotifScorer,
        "from_substrate_sequences",
        classmethod(counting_from_sequences),
    )

    result = workflow.run_request(request)

    assert list(result.prediction_result.pred_matrix.columns) == [
        "KINASE_A",
        "KINASE_B",
    ]
    assert matrix_calls == ["phospho_matrix"]
    assert motif_calls == [2]


def test_kinase_workflow_rejects_inconsistent_motif_widths_at_boundary() -> None:
    phospho_matrix, substrate_map, site_sequences, motif_sequences = (
        make_workflow_inputs()
    )
    workflow = KinaseWorkflow(flank_size=2)
    motif_sequences["KINASE_B"] = ["QTTY", "QTTY", "QTTY"]

    with pytest.raises(
        InputCompatibilityError,
        match="motif_sequences must use the same sequence width across kinases",
    ):
        workflow.run(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            motif_sequences=motif_sequences,
            min_substrates=2,
            min_motif_size=2,
            ensemble_size=3,
            top=4,
            score_threshold=0.75,
            inclusion=3,
            n_iterations=2,
            random_state=17,
        )


def test_workflow_execution_planner_resolves_predictor_mode_from_request() -> None:
    phospho_matrix, substrate_map, site_sequences, motif_sequences = (
        make_workflow_inputs()
    )
    request = KinaseWorkflowRequest.validate_request(
        phospho_matrix=phospho_matrix,
        substrate_map=substrate_map,
        site_sequences=site_sequences,
        motif_sequences=motif_sequences,
        min_substrates=2,
        min_motif_size=2,
        ensemble_size=3,
        top=4,
        score_threshold=0.75,
        inclusion=3,
        n_iterations=2,
        random_state=17,
        svm_mode="r_parity",
    )

    plan = WorkflowExecutionPlanner(
        flank_size=2,
        kernel="rbf",
        default_svm_mode="default",
    ).plan(request)

    assert plan.predictor_svm_mode == "r_parity"
    assert list(plan.validated_inputs.phospho_matrix.index) == list(
        phospho_matrix.index
    )
