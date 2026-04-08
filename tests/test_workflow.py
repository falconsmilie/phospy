from __future__ import annotations

import pandas as pd
import pytest

from phospy.motifs import KinaseMotifScorer
from phospy.prediction import PredMatResult
from phospy.validation.errors import (
    InputCompatibilityError,
    NoCandidateKinasesError,
    RequestValidationError,
)
from phospy.validation.tables import SiteMatrixSchema
from phospy.workflow import KinaseWorkflow, PredMatWorkflow


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


def test_pred_mat_workflow_runs_native_end_to_end_path() -> None:
    phospho_matrix, substrate_map, site_sequences, motif_sequences = (
        make_workflow_inputs()
    )

    workflow = PredMatWorkflow(flank_size=2)

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

    assert list(result.scoring_result.combined_scores.columns) == [
        "KINASE_A",
        "KINASE_B",
    ]
    assert list(result.prediction_result.pred_matrix.columns) == [
        "KINASE_A",
        "KINASE_B",
    ]
    assert result.pred_mat_result is result.prediction_result.pred_mat_result
    assert set(result.prediction_result.substrate_list) == {"KINASE_A", "KINASE_B"}


def test_pred_mat_workflow_raises_domain_error_when_no_candidate_kinases_qualify() -> (
    None
):
    phospho_matrix, substrate_map, site_sequences, motif_sequences = (
        make_workflow_inputs()
    )

    with pytest.raises(
        NoCandidateKinasesError,
        match=(
            r"No candidate kinases qualified for prediction from combined_scores "
            r"using top=2, score_threshold=0\.0, and inclusion=3"
        ),
    ):
        PredMatWorkflow(flank_size=2).run(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            motif_sequences=motif_sequences,
            min_substrates=2,
            min_motif_size=2,
            ensemble_size=2,
            top=2,
            score_threshold=0.0,
            inclusion=3,
            n_iterations=2,
            random_state=17,
        )


def test_pred_mat_workflow_result_has_canonical_pred_mat_result() -> None:
    phospho_matrix, substrate_map, site_sequences, motif_sequences = (
        make_workflow_inputs()
    )

    workflow = PredMatWorkflow(flank_size=2)
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

    assert isinstance(result.pred_mat_result, PredMatResult)
    assert result.pred_mat_result is result.prediction_result.pred_mat_result
    pd.testing.assert_frame_equal(
        result.pred_mat_result.to_frame(copy=False),
        result.prediction_result.pred_mat_result.to_frame(copy=False),
    )
    assert not hasattr(result, "pred_mat")


def test_kinase_workflow_result_tables_are_detached_from_input_matrix() -> None:
    phospho_matrix, substrate_map, site_sequences, motif_sequences = (
        make_workflow_inputs()
    )
    workflow = KinaseWorkflow(flank_size=2)
    original = phospho_matrix.copy(deep=True)

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

    result.profile_result.profile_matrix.loc[
        "KINASE_A", phospho_matrix.columns[0]
    ] = -999.0
    if result.motif_result is not None:
        result.motif_result.motif_scores.loc[
            phospho_matrix.index[0], "KINASE_A"
        ] = -999.0
    result.scoring_result.profile_scores.loc[
        phospho_matrix.index[0], "KINASE_A"
    ] = -999.0
    result.prediction_result.pred_matrix.loc[
        phospho_matrix.index[0], "KINASE_A"
    ] = -999.0

    pd.testing.assert_frame_equal(phospho_matrix, original)


def test_kinase_workflow_limits_motif_and_prediction_outputs_to_sites_with_sequences() -> (
    None
):
    phospho_matrix, substrate_map, site_sequences, motif_sequences = (
        make_workflow_inputs()
    )
    workflow = KinaseWorkflow(flank_size=2)
    partial_site_sequences = {
        site_id: site_sequences[site_id]
        for site_id in ["SITE_1", "SITE_2", "SITE_5", "SITE_6"]
    }

    result = workflow.run(
        phospho_matrix=phospho_matrix,
        substrate_map=substrate_map,
        site_sequences=partial_site_sequences,
        motif_sequences=motif_sequences,
        min_substrates=2,
        min_motif_size=2,
        ensemble_size=2,
        top=2,
        score_threshold=0.0,
        inclusion=1,
        n_iterations=1,
        random_state=17,
    )

    assert list(result.motif_result.motif_scores.index) == [
        "SITE_1",
        "SITE_2",
        "SITE_5",
        "SITE_6",
    ]
    assert list(result.scoring_result.profile_scores.index) == [
        "SITE_1",
        "SITE_2",
        "SITE_5",
        "SITE_6",
    ]
    assert list(result.scoring_result.combined_scores.index) == [
        "SITE_1",
        "SITE_2",
        "SITE_5",
        "SITE_6",
    ]
    assert list(result.prediction_result.pred_matrix.index) == [
        "SITE_1",
        "SITE_2",
        "SITE_5",
        "SITE_6",
    ]


def test_pred_mat_workflow_accepts_partial_site_sequence_coverage() -> None:
    phospho_matrix, substrate_map, site_sequences, motif_sequences = (
        make_workflow_inputs()
    )
    workflow = PredMatWorkflow(flank_size=2)
    partial_site_sequences = {
        site_id: site_sequences[site_id]
        for site_id in ["SITE_1", "SITE_2", "SITE_5", "SITE_6"]
    }

    result = workflow.run(
        phospho_matrix=phospho_matrix,
        substrate_map=substrate_map,
        site_sequences=partial_site_sequences,
        motif_sequences=motif_sequences,
        min_substrates=2,
        min_motif_size=2,
        ensemble_size=2,
        top=2,
        score_threshold=0.0,
        inclusion=1,
        n_iterations=1,
        random_state=17,
    )

    expected_index = ["SITE_1", "SITE_2", "SITE_5", "SITE_6"]
    assert list(result.scoring_result.profile_scores.index) == expected_index
    assert list(result.scoring_result.combined_scores.index) == expected_index
    assert list(result.prediction_result.pred_matrix.index) == expected_index
    assert list(result.pred_mat_result.to_frame(copy=False).index) == expected_index


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


def test_kinase_workflow_run_validated_uses_validated_boundary_request(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    phospho_matrix, substrate_map, site_sequences, motif_sequences = (
        make_workflow_inputs()
    )
    workflow = KinaseWorkflow(flank_size=2)

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

    request = workflow._validate_request(
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
    result = workflow.run_validated(request)

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


def test_validated_workflow_request_resolves_predictor_mode_from_request() -> None:
    phospho_matrix, substrate_map, site_sequences, motif_sequences = (
        make_workflow_inputs()
    )
    request = KinaseWorkflow(flank_size=2)._validate_request(
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

    assert request.predictor_svm_mode == "r_parity"
    assert list(request.phospho_matrix.index) == list(phospho_matrix.index)
