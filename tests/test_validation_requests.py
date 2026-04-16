from __future__ import annotations

import pandas as pd
import pytest

from phospy.datasets import DatasetSchema, PhosphoDataset
from phospy.errors import (
    InputCompatibilityError,
    NoCandidateKinasesError,
    RequestValidationError,
)
from phospy.prediction import PredMatResult
from phospy.prediction.traces import TraceSink
from phospy.preprocessing.matrices import SiteMatrixPolicy
from phospy.validation.requests import (
    CorePipelineRequest,
    KinaseActivityRequest,
    KinaseWorkflowRequest,
    PredictionRequest,
    validate_analysis_request,
    validate_pipeline_construction_request,
    validate_pipeline_runtime_compatibility,
    validate_workflow_request,
)
from phospy.validation.requests.analysis import AnalysisInputs
from phospy.validation.requests.dataset import DatasetInputs
from phospy.validation.requests.pipeline import (
    PipelineInputs,
    build_pipeline_inputs,
)
from phospy.validation.requests.workflow import (
    WorkflowInputs,
    build_workflow_inputs,
)


class _DummyTraceSink(TraceSink):
    def write_rows(self, table_name: str, rows: list[dict[str, object]]) -> None:
        return None

    def read_table(self, table_name: str) -> pd.DataFrame:
        return pd.DataFrame()


def test_core_pipeline_request_requires_existing_paths(tmp_path) -> None:
    total_path = tmp_path / "missing_total.tsv"
    phospho_path = tmp_path / "missing_phospho.tsv"

    with pytest.raises(RequestValidationError, match="Path does not exist"):
        CorePipelineRequest.validate_request(
            total_path=total_path,
            phospho_path=phospho_path,
        )


def test_core_pipeline_request_rejects_duplicate_comparisons(tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    total_path.write_text("genes\tgroup1\nPRKACA\t1\n")
    phospho_path.write_text("uid\tgene_names\n1\tPRKACA\n")

    with pytest.raises(RequestValidationError, match="Duplicate comparison pair"):
        CorePipelineRequest.validate_request(
            total_path=total_path,
            phospho_path=phospho_path,
            comparisons=(("group1", "group4"), ("group1", "group4")),
        )


def test_core_pipeline_request_rejects_unknown_comparison_group(tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    total_path.write_text("genes\tgroup1\nPRKACA\t1\n")
    phospho_path.write_text("uid\tgene_names\n1\tPRKACA\n")

    with pytest.raises(RequestValidationError, match="Unknown comparison group"):
        CorePipelineRequest.validate_request(
            total_path=total_path,
            phospho_path=phospho_path,
            comparisons=(("group1", "group9"),),
        )


def test_core_pipeline_request_rejects_self_comparisons(tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    total_path.write_text("genes\tgroup1\nPRKACA\t1\n")
    phospho_path.write_text("uid\tgene_names\n1\tPRKACA\n")

    with pytest.raises(RequestValidationError, match="Self comparison pair"):
        CorePipelineRequest.validate_request(
            total_path=total_path,
            phospho_path=phospho_path,
            comparisons=(("group1", "group1"),),
        )


def test_core_pipeline_request_rejects_reverse_duplicate_comparisons(tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    total_path.write_text("genes\tgroup1\nPRKACA\t1\n")
    phospho_path.write_text("uid\tgene_names\n1\tPRKACA\n")

    with pytest.raises(
        RequestValidationError,
        match="Duplicate comparison pair regardless of direction",
    ):
        CorePipelineRequest.validate_request(
            total_path=total_path,
            phospho_path=phospho_path,
            comparisons=(("group1", "group4"), ("group4", "group1")),
        )


def test_kinase_workflow_request_rejects_invalid_threshold() -> None:
    phospho_matrix = pd.DataFrame({"sample_1": [1.0]}, index=["SITE_1"])

    with pytest.raises(RequestValidationError, match="score_threshold"):
        KinaseWorkflowRequest.validate_request(
            phospho_matrix=phospho_matrix,
            substrate_map={"KINASE_A": ["SITE_1"]},
            motif_sequences=None,
            allow_profile_only_fallback=True,
            score_threshold=1.5,
        )


@pytest.mark.parametrize(
    ("field_name", "field_value"),
    [("min_substrates", 0), ("top_n_substrates", 0)],
)
def test_kinase_activity_request_rejects_invalid_positive_int_fields(
    field_name: str,
    field_value: int,
) -> None:
    with pytest.raises(RequestValidationError, match=field_name):
        KinaseActivityRequest.validate_request(**{field_name: field_value})


@pytest.mark.parametrize("threshold", [-0.1, 1.1])
def test_kinase_activity_request_rejects_invalid_threshold(
    threshold: float,
) -> None:
    with pytest.raises(RequestValidationError, match="threshold"):
        KinaseActivityRequest.validate_request(threshold=threshold)


def test_kinase_workflow_request_requires_site_sequences_with_motifs() -> None:
    phospho_matrix = pd.DataFrame({"sample_1": [1.0]}, index=["SITE_1"])

    with pytest.raises(RequestValidationError, match="site_sequences are required"):
        KinaseWorkflowRequest.validate_request(
            phospho_matrix=phospho_matrix,
            substrate_map={"KINASE_A": ["SITE_1"]},
            motif_sequences={"KINASE_A": ["QQAAAAAYY"]},
        )


def test_core_pipeline_request_rejects_invalid_max_unmatched_fraction(tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    total_path.write_text("genes\tgroup1\nPRKACA\t1\n")
    phospho_path.write_text("uid\tgene_names\n1\tPRKACA\n")

    with pytest.raises(RequestValidationError, match="max_unmatched_fraction"):
        CorePipelineRequest.validate_request(
            total_path=total_path,
            phospho_path=phospho_path,
            max_unmatched_fraction=1.5,
        )


def test_core_pipeline_request_rejects_directory_paths(tmp_path) -> None:
    total_dir = tmp_path / "total_dir"
    phospho_path = tmp_path / "phospho.tsv"
    total_dir.mkdir()
    phospho_path.write_text(
        "uid	gene_names	gene_p_site	localization_prob	centralized_sequence	p_group1	p_group2	p_group3	p_group4	p_group5	p_group6\n1	PRKACA	PRKACA_S339	0.95	AAAAAA	1	1	1	1	1	1\n"
    )

    with pytest.raises(RequestValidationError, match="Path is not a file"):
        CorePipelineRequest.validate_request(
            total_path=total_dir,
            phospho_path=phospho_path,
        )


def test_kinase_workflow_request_rejects_plain_site_sequence_lists() -> None:
    phospho_matrix = pd.DataFrame({"sample_1": [1.0]}, index=["SITE_1"])

    with pytest.raises(RequestValidationError, match="site_sequences must be provided"):
        KinaseWorkflowRequest.validate_request(
            phospho_matrix=phospho_matrix,
            substrate_map={"KINASE_A": ["SITE_1"]},
            site_sequences=["QQAAAAAYY"],
            motif_sequences={"KINASE_A": ["QQAAAAAYY"]},
        )


def test_kinase_workflow_request_rejects_unordered_substrate_map_iterables() -> None:
    phospho_matrix = pd.DataFrame({"sample_1": [1.0]}, index=["SITE_1"])

    with pytest.raises(
        RequestValidationError,
        match="ordered sequence of values",
    ):
        KinaseWorkflowRequest.validate_request(
            phospho_matrix=phospho_matrix,
            substrate_map={"KINASE_A": {"SITE_1"}},
            motif_sequences=None,
            allow_profile_only_fallback=True,
        )


def test_kinase_workflow_request_rejects_site_sequences_with_null_mapping_value() -> (
    None
):
    phospho_matrix = pd.DataFrame({"sample_1": [1.0]}, index=["SITE_1"])

    with pytest.raises(
        RequestValidationError,
        match="site_sequences values must be non-null sequence strings",
    ):
        KinaseWorkflowRequest.validate_request(
            phospho_matrix=phospho_matrix,
            substrate_map={"KINASE_A": ["SITE_1"]},
            site_sequences={"SITE_1": float("nan")},
            motif_sequences={"KINASE_A": ["QQAAAAAYY"]},
        )


def test_kinase_workflow_request_rejects_site_sequences_with_null_series_value() -> (
    None
):
    phospho_matrix = pd.DataFrame({"sample_1": [1.0]}, index=["SITE_1"])
    site_sequences = pd.Series({"SITE_1": pd.NA}, dtype=object)

    with pytest.raises(
        RequestValidationError,
        match="site_sequences values must be non-null sequence strings",
    ):
        KinaseWorkflowRequest.validate_request(
            phospho_matrix=phospho_matrix,
            substrate_map={"KINASE_A": ["SITE_1"]},
            site_sequences=site_sequences,
            motif_sequences={"KINASE_A": ["QQAAAAAYY"]},
        )


def test_core_pipeline_request_accepts_explicit_sentinel_configuration(
    tmp_path,
) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    total_path.write_text("genes\tgroup1\nPRKACA\t1\n")
    phospho_path.write_text("uid\tgene_names\n1\tPRKACA\n")

    request = CorePipelineRequest.validate_request(
        total_path=total_path,
        phospho_path=phospho_path,
        total_sentinel=99.0,
        phospho_sentinel=88.0,
    )

    assert request.total_sentinel == 99.0
    assert request.phospho_sentinel == 88.0


def test_core_pipeline_request_accepts_site_matrix_policy_mapping(tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    total_path.write_text("genes\tgroup1\nPRKACA\t1\n")
    phospho_path.write_text("uid\tgene_names\n1\tPRKACA\n")

    request = CorePipelineRequest.validate_request(
        total_path=total_path,
        phospho_path=phospho_path,
        site_matrix_policy={
            "duplicate_site_strategy": "first",
            "missing_data_policy": "retain_missing",
        },
    )

    assert request.site_matrix_policy == SiteMatrixPolicy(
        duplicate_site_strategy="first",
        missing_data_policy="retain_missing",
    )


def test_prediction_request_rejects_invalid_trace_level() -> None:
    with pytest.raises(RequestValidationError, match="trace_level"):
        PredictionRequest.validate_request(
            combined_scores=pd.DataFrame({"KINASE_A": [0.9]}, index=["SITE_1"]),
            ensemble_size=2,
            top=2,
            score_threshold=0.8,
            inclusion=1,
            n_iterations=2,
            debug_top_n=1,
            trace_level="broken",
            default_svm_mode="default",
        )


def test_prediction_request_rejects_invalid_trace_sink_format() -> None:
    with pytest.raises(RequestValidationError, match="trace_sink_format"):
        PredictionRequest.validate_request(
            combined_scores=pd.DataFrame({"KINASE_A": [0.9]}, index=["SITE_1"]),
            ensemble_size=2,
            top=2,
            score_threshold=0.8,
            inclusion=1,
            n_iterations=2,
            debug_top_n=1,
            trace_level="full",
            trace_sink_format="json",
            default_svm_mode="default",
        )


def test_prediction_request_uses_boundary_defaults() -> None:
    request = PredictionRequest.validate_request(
        combined_scores=pd.DataFrame({"KINASE_A": [0.9]}, index=["SITE_1"]),
        ensemble_size=2,
        top=2,
        score_threshold=0.8,
        inclusion=1,
        n_iterations=2,
        debug_top_n=1,
        capture_debug_trace=True,
        default_svm_mode="r_parity",
    )

    assert request.svm_mode == "r_parity"
    assert request.trace_level == "summary"


def test_prediction_request_validation_does_not_allocate_trace_sink() -> None:
    request = PredictionRequest.validate_request(
        combined_scores=pd.DataFrame({"KINASE_A": [0.9]}, index=["SITE_1"]),
        ensemble_size=2,
        top=2,
        score_threshold=0.8,
        inclusion=1,
        n_iterations=2,
        debug_top_n=1,
        trace_level="full",
        default_svm_mode="default",
    )

    assert request.trace_level == "full"
    assert request.trace_sink is None


def test_prediction_request_validation_preserves_explicit_trace_sink_spec() -> None:
    trace_sink = _DummyTraceSink()
    request = PredictionRequest.validate_request(
        combined_scores=pd.DataFrame({"KINASE_A": [0.9]}, index=["SITE_1"]),
        ensemble_size=2,
        top=2,
        score_threshold=0.8,
        inclusion=1,
        n_iterations=2,
        debug_top_n=1,
        trace_level="full",
        trace_sink=trace_sink,
        default_svm_mode="default",
    )

    assert isinstance(request.trace_sink, TraceSink)
    assert request.trace_sink is trace_sink


def test_prediction_request_rejects_trace_sink_without_full_trace_level() -> None:
    with pytest.raises(
        RequestValidationError,
        match="trace_sink may only be provided when trace_level='full'",
    ):
        PredictionRequest.validate_request(
            combined_scores=pd.DataFrame({"KINASE_A": [0.9]}, index=["SITE_1"]),
            ensemble_size=2,
            top=2,
            score_threshold=0.8,
            inclusion=1,
            n_iterations=2,
            debug_top_n=1,
            trace_level="summary",
            trace_sink=_DummyTraceSink(),
            default_svm_mode="default",
        )


def test_prediction_request_rejects_non_positive_integer_fields() -> None:
    with pytest.raises(RequestValidationError, match="ensemble_size"):
        PredictionRequest.validate_request(
            combined_scores=pd.DataFrame({"KINASE_A": [0.9]}, index=["SITE_1"]),
            ensemble_size=0,
            top=2,
            score_threshold=0.8,
            inclusion=1,
            n_iterations=2,
            debug_top_n=1,
            default_svm_mode="default",
        )


def test_prediction_request_rejects_out_of_range_score_threshold() -> None:
    with pytest.raises(RequestValidationError, match="score_threshold"):
        PredictionRequest.validate_request(
            combined_scores=pd.DataFrame({"KINASE_A": [0.9]}, index=["SITE_1"]),
            ensemble_size=2,
            top=2,
            score_threshold=1.1,
            inclusion=1,
            n_iterations=2,
            debug_top_n=1,
            default_svm_mode="default",
        )


def test_prediction_request_rejects_empty_combined_scores() -> None:
    with pytest.raises(RequestValidationError, match="combined_scores"):
        PredictionRequest.validate_request(
            combined_scores=pd.DataFrame(columns=["KINASE_A"]),
            ensemble_size=2,
            top=2,
            score_threshold=0.8,
            inclusion=1,
            n_iterations=2,
            debug_top_n=1,
            default_svm_mode="default",
        )


def test_prediction_request_rejects_non_numeric_combined_scores() -> None:
    with pytest.raises(RequestValidationError, match="combined_scores"):
        PredictionRequest.validate_request(
            combined_scores=pd.DataFrame({"KINASE_A": ["bad"]}, index=["SITE_1"]),
            ensemble_size=2,
            top=2,
            score_threshold=0.8,
            inclusion=1,
            n_iterations=2,
            debug_top_n=1,
            default_svm_mode="default",
        )


def test_prediction_request_rejects_duplicate_combined_score_index() -> None:
    with pytest.raises(RequestValidationError, match="duplicate index"):
        PredictionRequest.validate_request(
            combined_scores=pd.DataFrame(
                {"KINASE_A": [0.8, 0.7]},
                index=["SITE_1", "SITE_1"],
            ),
            ensemble_size=2,
            top=2,
            score_threshold=0.8,
            inclusion=1,
            n_iterations=2,
            debug_top_n=1,
            default_svm_mode="default",
        )


def test_core_pipeline_request_validates_comparisons_against_explicit_schema(
    tmp_path,
) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    total_path.write_text("genes\tsample_a\tsample_b\nPRKACA\t1\t1\n")
    phospho_path.write_text(
        "uid\tgene_names\tgene_p_site\tlocalization_prob\tcentralized_sequence\tp_sample_a\tp_sample_b\n"
        "1\tPRKACA\tPRKACA_S339\t0.95\tAAAAAA\t1\t1\n"
    )
    schema = DatasetSchema(
        total_cols=("sample_a", "sample_b"),
        phospho_cols=("p_sample_a", "p_sample_b"),
        corrected_cols=("corrected_a", "corrected_b"),
    )

    request = CorePipelineRequest.validate_request(
        total_path=total_path,
        phospho_path=phospho_path,
        schema=schema,
        comparisons=(("sample_a", "sample_b"),),
    )

    assert request.dataset_schema is schema
    assert request.comparisons == (("sample_a", "sample_b"),)


def test_core_pipeline_request_rejects_unknown_comparison_group_for_schema(
    tmp_path,
) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    total_path.write_text("genes\tsample_a\tsample_b\nPRKACA\t1\t1\n")
    phospho_path.write_text(
        "uid\tgene_names\tgene_p_site\tlocalization_prob\tcentralized_sequence\tp_sample_a\tp_sample_b\n"
        "1\tPRKACA\tPRKACA_S339\t0.95\tAAAAAA\t1\t1\n"
    )
    schema = DatasetSchema(
        total_cols=("sample_a", "sample_b"),
        phospho_cols=("p_sample_a", "p_sample_b"),
        corrected_cols=("corrected_a", "corrected_b"),
    )

    with pytest.raises(RequestValidationError, match="Unknown comparison group"):
        CorePipelineRequest.validate_request(
            total_path=total_path,
            phospho_path=phospho_path,
            schema=schema,
            comparisons=(("group1", "sample_b"),),
        )


def test_core_pipeline_request_does_not_mask_unexpected_comparison_validation_errors(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    total_path.write_text("genes\tgroup1\nPRKACA\t1\n")
    phospho_path.write_text("uid\tgene_names\n1\tPRKACA\n")

    def blow_up(*, schema, comparisons, context: str):
        raise RuntimeError("comparison validator exploded")

    monkeypatch.setattr(
        "phospy.validation.requests.pipeline.validate_dataset_comparisons",
        blow_up,
    )

    with pytest.raises(RuntimeError, match="comparison validator exploded"):
        CorePipelineRequest.validate_request(
            total_path=total_path,
            phospho_path=phospho_path,
            comparisons=(("group1", "group1"),),
        )


def test_validation_modules_expose_use_case_boundaries() -> None:
    assert CorePipelineRequest.__module__ == "phospy.validation.requests.pipeline"
    assert KinaseActivityRequest.__module__ == "phospy.validation.requests.analysis"
    assert KinaseWorkflowRequest.__module__ == "phospy.validation.requests.workflow"
    assert PredictionRequest.__module__ == "phospy.validation.requests.prediction"


def test_pipeline_request_can_be_created_from_public_dataset_boundary() -> None:
    from phospy.datasets import PhosphoDataset

    total_df = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
            "group2": [1.0],
            "group3": [1.0],
            "group4": [1.0],
            "group5": [1.0],
            "group6": [1.0],
        }
    )
    phospho_df = pd.DataFrame(
        {
            "uid": ["u1"],
            "gene_names": ["PRKACA"],
            "gene_p_site": ["PRKACA_S339"],
            "localization_prob": [0.95],
            "centralized_sequence": ["AAAAAA"],
            "p_group1": [1.0],
            "p_group2": [1.0],
            "p_group3": [1.0],
            "p_group4": [1.0],
            "p_group5": [1.0],
            "p_group6": [1.0],
        }
    )

    pipeline_request = validate_pipeline_construction_request(
        dataset=PhosphoDataset(total_df=total_df, phospho_df=phospho_df),
        pred_mat=pd.DataFrame({"PRKACA": [0.9]}, index=["PRKACA;S339;"]),
    )
    assert isinstance(pipeline_request, PipelineInputs)


def test_dataset_validation_internals_are_not_exported_from_validation_package() -> (
    None
):
    import phospy.validation as validation
    import phospy.validation.requests as request_exports

    assert not hasattr(validation, "validate_dataset_request")
    assert not hasattr(validation, "DatasetInputs")
    assert not hasattr(validation, "AnalysisInputs")
    assert not hasattr(validation, "PipelineInputs")
    assert not hasattr(validation, "WorkflowInputs")
    assert not hasattr(request_exports, "DatasetInputs")
    assert not hasattr(request_exports, "AnalysisInputs")
    assert not hasattr(request_exports, "PipelineInputs")
    assert not hasattr(request_exports, "WorkflowInputs")
    assert not hasattr(validation, "validate_workflow_inputs")
    assert not hasattr(request_exports, "PhospyRequestModel")
    assert not hasattr(request_exports, "build_pipeline_inputs")
    assert not hasattr(request_exports, "build_workflow_inputs")
    assert not hasattr(request_exports, "normalize_pred_mat_input")
    assert not hasattr(request_exports, "validate_adapter_value")
    assert not hasattr(request_exports, "validate_existing_file_path")
    assert not hasattr(request_exports, "validate_workflow_inputs")


def test_validation_request_package_exports_match_supported_surface() -> None:
    import phospy.validation.requests as request_exports

    assert set(request_exports.__all__) == {
        "CorePipelineRequest",
        "KinaseActivityRequest",
        "KinaseWorkflowRequest",
        "PredictionRequest",
        "SignalomeRequest",
        "validate_analysis_request",
        "validate_pipeline_construction_request",
        "validate_pipeline_runtime_compatibility",
        "validate_signalome_request",
        "validate_workflow_request",
    }


def test_validated_workflow_and_analysis_requests_can_be_created() -> None:
    workflow_request = validate_workflow_request(
        phospho_matrix=pd.DataFrame(
            {"sample_1": [1.0], "sample_2": [2.0]},
            index=["SITE_1"],
        ),
        substrate_map={"KINASE_A": ["SITE_1"]},
        site_sequences={"SITE_1": "QQAAAAAYY"},
        motif_sequences={"KINASE_A": ["QQAAAAAYY"]},
        allow_profile_only_fallback=False,
        flank_size=2,
        default_svm_mode="default",
    )
    assert isinstance(workflow_request, WorkflowInputs)

    analysis_request = validate_analysis_request(
        pred_mat=pd.DataFrame({"PRKACA": [0.9]}, index=["PRKACA;S339;"]),
        phospho_matrix=pd.DataFrame({"sample_1": [1.0]}, index=["PRKACA;S339;"]),
    )
    assert isinstance(analysis_request, AnalysisInputs)


def test_validate_pipeline_construction_request_rejects_non_dataset_inputs() -> None:
    with pytest.raises(
        RequestValidationError, match="dataset must be a PhosphoDataset"
    ):
        validate_pipeline_construction_request(dataset=object())


def test_validate_pipeline_construction_request_rejects_mixed_preprocessing_config_styles() -> (
    None
):
    from phospy.datasets import PhosphoDataset
    from phospy.preprocessing import CorePreprocessingConfig

    total_df = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
            "group2": [1.0],
            "group3": [1.0],
            "group4": [1.0],
            "group5": [1.0],
            "group6": [1.0],
        }
    )
    phospho_df = pd.DataFrame(
        {
            "uid": ["u1"],
            "gene_names": ["PRKACA"],
            "gene_p_site": ["PRKACA_S339"],
            "localization_prob": [0.95],
            "centralized_sequence": ["AAAAAA"],
            "p_group1": [1.0],
            "p_group2": [1.0],
            "p_group3": [1.0],
            "p_group4": [1.0],
            "p_group5": [1.0],
            "p_group6": [1.0],
        }
    )

    with pytest.raises(
        RequestValidationError,
        match=(
            r"Invalid pipeline construction request: pass either preprocessing_config or scalar "
            r"preprocessing options, not both\."
        ),
    ):
        validate_pipeline_construction_request(
            dataset=PhosphoDataset(total_df=total_df, phospho_df=phospho_df),
            preprocessing_config=CorePreprocessingConfig(),
            min_observed=1,
        )


def test_validate_pipeline_construction_request_rejects_mixed_site_matrix_policy_styles() -> (
    None
):
    from phospy.datasets import PhosphoDataset
    from phospy.preprocessing import CorePreprocessingConfig

    total_df = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
            "group2": [1.0],
            "group3": [1.0],
            "group4": [1.0],
            "group5": [1.0],
            "group6": [1.0],
        }
    )
    phospho_df = pd.DataFrame(
        {
            "uid": ["u1"],
            "gene_names": ["PRKACA"],
            "gene_p_site": ["PRKACA_S339"],
            "localization_prob": [0.95],
            "centralized_sequence": ["AAAAAA"],
            "p_group1": [1.0],
            "p_group2": [1.0],
            "p_group3": [1.0],
            "p_group4": [1.0],
            "p_group5": [1.0],
            "p_group6": [1.0],
        }
    )

    with pytest.raises(
        RequestValidationError,
        match=(
            r"Invalid pipeline construction request: pass either preprocessing_config or scalar "
            r"preprocessing options, not both\."
        ),
    ):
        validate_pipeline_construction_request(
            dataset=PhosphoDataset(total_df=total_df, phospho_df=phospho_df),
            preprocessing_config=CorePreprocessingConfig(),
            site_matrix_policy=SiteMatrixPolicy(missing_data_policy="retain_missing"),
        )


def test_validate_analysis_request_takes_ownership_of_raw_dataframe_inputs() -> None:
    pred_mat = pd.DataFrame({"PRKACA": [0.9]}, index=["PRKACA;S339;"])
    phospho_matrix = pd.DataFrame({"sample_1": [1.0]}, index=["PRKACA;S339;"])

    request = validate_analysis_request(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
    )

    pred_mat.loc["PRKACA;S339;", "PRKACA"] = 0.1
    phospho_matrix.loc["PRKACA;S339;", "sample_1"] = 9.0

    assert request.pred_mat.loc["PRKACA;S339;", "PRKACA"] == 0.9
    assert request.phospho_matrix.loc["PRKACA;S339;", "sample_1"] == 1.0


def test_validate_analysis_request_normalizes_pred_mat_result_at_boundary() -> None:
    pred_mat = pd.DataFrame({"PRKACA": [0.9]}, index=["PRKACA;S339;"])
    phospho_matrix = pd.DataFrame({"sample_1": [1.0]}, index=["PRKACA;S339;"])
    pred_mat_result = PredMatResult(pred_mat)

    request = validate_analysis_request(
        pred_mat=pred_mat_result,
        phospho_matrix=phospho_matrix,
    )

    pred_mat.loc["PRKACA;S339;", "PRKACA"] = 0.1

    assert request.pred_mat.loc["PRKACA;S339;", "PRKACA"] == 0.9


def test_validate_workflow_request_takes_ownership_of_raw_dataframe_inputs() -> None:
    phospho_matrix = pd.DataFrame(
        {"sample_1": [1.0], "sample_2": [2.0]},
        index=["SITE_1"],
    )

    request = validate_workflow_request(
        phospho_matrix=phospho_matrix,
        substrate_map={"KINASE_A": ["SITE_1"]},
        site_sequences={"SITE_1": "QQAAAAAYY"},
        motif_sequences={"KINASE_A": ["QQAAAAAYY"]},
        allow_profile_only_fallback=False,
        flank_size=2,
        default_svm_mode="default",
    )

    phospho_matrix.loc["SITE_1", "sample_1"] = 99.0

    assert request.phospho_matrix.loc["SITE_1", "sample_1"] == 1.0
    assert request.phospho_matrix is not phospho_matrix


def test_build_workflow_inputs_reuses_owned_validated_matrix() -> None:
    raw_request = KinaseWorkflowRequest.validate_request(
        phospho_matrix=pd.DataFrame(
            {"sample_1": [1.0], "sample_2": [2.0]},
            index=["SITE_1"],
        ),
        substrate_map={"KINASE_A": ["SITE_1"]},
        site_sequences={"SITE_1": "QQAAAAAYY"},
        motif_sequences={"KINASE_A": ["QQAAAAAYY"]},
        allow_profile_only_fallback=False,
    )

    request = build_workflow_inputs(
        raw_request,
        flank_size=2,
        default_svm_mode="default",
    )

    pd.testing.assert_frame_equal(request.phospho_matrix, raw_request.phospho_matrix)


def test_kinase_workflow_request_detaches_mapping_backed_sequence_inputs() -> None:
    substrate_map = {"KINASE_A": ["SITE_1"]}
    site_sequences = {"SITE_1": "QQAAAAAYY"}
    motif_sequences = {"KINASE_A": ["QQAAAAAYY"]}

    request = KinaseWorkflowRequest.validate_request(
        phospho_matrix=pd.DataFrame(
            {"sample_1": [1.0], "sample_2": [2.0]},
            index=["SITE_1"],
        ),
        substrate_map=substrate_map,
        site_sequences=site_sequences,
        motif_sequences=motif_sequences,
        allow_profile_only_fallback=False,
    )

    substrate_map["KINASE_A"].append("SITE_2")
    site_sequences["SITE_1"] = "CHANGED"
    motif_sequences["KINASE_A"].append("CHANGED")

    assert request.substrate_map == {"KINASE_A": ("SITE_1",)}
    assert request.site_sequences is not None
    assert request.site_sequences.to_dict() == {"SITE_1": "QQAAAAAYY"}
    assert request.motif_sequences == {"KINASE_A": ("QQAAAAAYY",)}


def test_validate_workflow_request_detaches_mapping_backed_runtime_inputs() -> None:
    substrate_map = {"KINASE_A": ["SITE_1"]}
    site_sequences = {"SITE_1": "QQAAAAAYY"}
    motif_sequences = {"KINASE_A": ["QQAAAAAYY"]}

    request = validate_workflow_request(
        phospho_matrix=pd.DataFrame(
            {"sample_1": [1.0], "sample_2": [2.0]},
            index=["SITE_1"],
        ),
        substrate_map=substrate_map,
        site_sequences=site_sequences,
        motif_sequences=motif_sequences,
        allow_profile_only_fallback=False,
        flank_size=2,
        default_svm_mode="default",
    )

    substrate_map["KINASE_A"].append("SITE_2")
    site_sequences["SITE_1"] = "CHANGED"
    motif_sequences["KINASE_A"].append("CHANGED")

    assert request.substrate_map == {"KINASE_A": ("SITE_1",)}
    assert request.site_sequences is not None
    assert request.site_sequences.to_dict() == {"SITE_1": "QQAAAAAYY"}
    assert request.motif_scorer is not None


def test_validate_workflow_request_detaches_series_backed_site_sequences() -> None:
    site_sequences = pd.Series({"SITE_1": "QQAAAAAYY"}, dtype=object)

    request = validate_workflow_request(
        phospho_matrix=pd.DataFrame(
            {"sample_1": [1.0], "sample_2": [2.0]},
            index=["SITE_1"],
        ),
        substrate_map={"KINASE_A": ["SITE_1"]},
        site_sequences=site_sequences,
        motif_sequences={"KINASE_A": ["QQAAAAAYY"]},
        allow_profile_only_fallback=False,
        flank_size=2,
        default_svm_mode="default",
    )

    site_sequences.loc["SITE_1"] = "CHANGED"

    assert request.site_sequences is not None
    assert request.site_sequences.loc["SITE_1"] == "QQAAAAAYY"


def test_validate_pipeline_construction_request_takes_ownership_of_raw_pred_mat_input() -> (
    None
):
    total_df = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
            "group2": [1.0],
            "group3": [1.0],
            "group4": [1.0],
            "group5": [1.0],
            "group6": [1.0],
        }
    )
    phospho_df = pd.DataFrame(
        {
            "uid": ["u1"],
            "gene_names": ["PRKACA"],
            "gene_p_site": ["PRKACA_S339"],
            "localization_prob": [0.95],
            "centralized_sequence": ["AAAAAA"],
            "p_group1": [1.0],
            "p_group2": [1.0],
            "p_group3": [1.0],
            "p_group4": [1.0],
            "p_group5": [1.0],
            "p_group6": [1.0],
        }
    )
    pred_mat = pd.DataFrame({"PRKACA": [0.9]}, index=["PRKACA;S339;"])

    request = validate_pipeline_construction_request(
        dataset=PhosphoDataset(total_df=total_df, phospho_df=phospho_df),
        pred_mat=pred_mat,
    )

    pred_mat.loc["PRKACA;S339;", "PRKACA"] = 0.1

    assert request.pred_mat is not None
    assert request.pred_mat.loc["PRKACA;S339;", "PRKACA"] == 0.9


def test_validate_pipeline_construction_request_normalizes_pred_mat_result_at_boundary() -> (
    None
):
    total_df = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
            "group2": [1.0],
            "group3": [1.0],
            "group4": [1.0],
            "group5": [1.0],
            "group6": [1.0],
        }
    )
    phospho_df = pd.DataFrame(
        {
            "uid": ["u1"],
            "gene_names": ["PRKACA"],
            "gene_p_site": ["PRKACA_S339"],
            "localization_prob": [0.95],
            "centralized_sequence": ["AAAAAA"],
            "p_group1": [1.0],
            "p_group2": [1.0],
            "p_group3": [1.0],
            "p_group4": [1.0],
            "p_group5": [1.0],
            "p_group6": [1.0],
        }
    )
    pred_mat = pd.DataFrame({"PRKACA": [0.9]}, index=["PRKACA;S339;"])
    pred_mat_result = PredMatResult(pred_mat)

    request = validate_pipeline_construction_request(
        dataset=PhosphoDataset(total_df=total_df, phospho_df=phospho_df),
        pred_mat=pred_mat_result,
    )

    pred_mat.loc["PRKACA;S339;", "PRKACA"] = 0.1

    assert request.pred_mat is not None
    assert request.pred_mat.loc["PRKACA;S339;", "PRKACA"] == 0.9


def test_validate_pipeline_construction_request_builds_explicit_kinase_activity_config() -> (
    None
):
    total_df = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
            "group2": [1.0],
            "group3": [1.0],
            "group4": [1.0],
            "group5": [1.0],
            "group6": [1.0],
        }
    )
    phospho_df = pd.DataFrame(
        {
            "uid": ["u1"],
            "gene_names": ["PRKACA"],
            "gene_p_site": ["PRKACA_S339"],
            "localization_prob": [0.95],
            "centralized_sequence": ["AAAAAA"],
            "p_group1": [1.0],
            "p_group2": [1.0],
            "p_group3": [1.0],
            "p_group4": [1.0],
            "p_group5": [1.0],
            "p_group6": [1.0],
        }
    )

    request = validate_pipeline_construction_request(
        dataset=PhosphoDataset(total_df=total_df, phospho_df=phospho_df),
        pred_mat=pd.DataFrame({"PRKACA": [0.9]}, index=["PRKACA;S339;"]),
        kinase_activity_threshold=0.8,
        kinase_activity_min_substrates=5,
        kinase_activity_top_n_substrates=7,
    )

    assert request.kinase_activity_threshold == 0.8
    assert request.kinase_activity_min_substrates == 5
    assert request.kinase_activity_top_n_substrates == 7


def test_validate_pipeline_construction_request_skips_kinase_activity_config_without_pred_mat() -> (
    None
):
    total_df = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
            "group2": [1.0],
            "group3": [1.0],
            "group4": [1.0],
            "group5": [1.0],
            "group6": [1.0],
        }
    )
    phospho_df = pd.DataFrame(
        {
            "uid": ["u1"],
            "gene_names": ["PRKACA"],
            "gene_p_site": ["PRKACA_S339"],
            "localization_prob": [0.95],
            "centralized_sequence": ["AAAAAA"],
            "p_group1": [1.0],
            "p_group2": [1.0],
            "p_group3": [1.0],
            "p_group4": [1.0],
            "p_group5": [1.0],
            "p_group6": [1.0],
        }
    )

    request = validate_pipeline_construction_request(
        dataset=PhosphoDataset(total_df=total_df, phospho_df=phospho_df),
        pred_mat=None,
        kinase_activity_threshold=0.8,
        kinase_activity_min_substrates=5,
        kinase_activity_top_n_substrates=7,
    )

    assert request.pred_mat is None
    assert request.kinase_activity_threshold is None
    assert request.kinase_activity_min_substrates is None
    assert request.kinase_activity_top_n_substrates is None


def test_build_pipeline_inputs_reuses_owned_dataset_and_pred_mat() -> None:
    total_df = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
            "group2": [1.0],
            "group3": [1.0],
            "group4": [1.0],
            "group5": [1.0],
            "group6": [1.0],
        }
    )
    phospho_df = pd.DataFrame(
        {
            "uid": ["u1"],
            "gene_names": ["PRKACA"],
            "gene_p_site": ["PRKACA_S339"],
            "localization_prob": [0.95],
            "centralized_sequence": ["AAAAAA"],
            "p_group1": [1.0],
            "p_group2": [1.0],
            "p_group3": [1.0],
            "p_group4": [1.0],
            "p_group5": [1.0],
            "p_group6": [1.0],
        }
    )
    dataset = PhosphoDataset(total_df=total_df, phospho_df=phospho_df)
    pred_mat = pd.DataFrame({"PRKACA": [0.9]}, index=["PRKACA;S339;"])

    request = build_pipeline_inputs(
        dataset=dataset,
        pred_mat=pred_mat,
    )

    assert request.dataset is dataset
    assert request.pred_mat is pred_mat


def test_build_pipeline_inputs_applies_site_matrix_policy_override() -> None:
    total_df = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
            "group2": [1.0],
            "group3": [1.0],
            "group4": [1.0],
            "group5": [1.0],
            "group6": [1.0],
        }
    )
    phospho_df = pd.DataFrame(
        {
            "uid": ["u1"],
            "gene_names": ["PRKACA"],
            "gene_p_site": ["PRKACA_S339"],
            "localization_prob": [0.95],
            "centralized_sequence": ["AAAAAA"],
            "p_group1": [1.0],
            "p_group2": [1.0],
            "p_group3": [1.0],
            "p_group4": [1.0],
            "p_group5": [1.0],
            "p_group6": [1.0],
        }
    )
    dataset = PhosphoDataset(total_df=total_df, phospho_df=phospho_df)

    request = build_pipeline_inputs(
        dataset=dataset,
        site_matrix_policy=SiteMatrixPolicy(missing_data_policy="retain_missing"),
    )

    assert (
        request.preprocessing_config.site_matrix_policy.missing_data_policy
        == "retain_missing"
    )


def test_validate_pipeline_runtime_compatibility_builds_analysis_request_after_preprocessing() -> (
    None
):
    total_df = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
            "group2": [1.0],
            "group3": [1.0],
            "group4": [1.0],
            "group5": [1.0],
            "group6": [1.0],
        }
    )
    phospho_df = pd.DataFrame(
        {
            "uid": ["u1"],
            "gene_names": ["PRKACA"],
            "gene_p_site": ["PRKACA_S339"],
            "localization_prob": [0.95],
            "centralized_sequence": ["AAAAAA"],
            "p_group1": [1.0],
            "p_group2": [1.0],
            "p_group3": [1.0],
            "p_group4": [1.0],
            "p_group5": [1.0],
            "p_group6": [1.0],
        }
    )

    request = validate_pipeline_construction_request(
        dataset=PhosphoDataset(total_df=total_df, phospho_df=phospho_df),
        pred_mat=pd.DataFrame({"PRKACA": [0.9]}, index=["PRKACA;S339;"]),
    )

    analysis_request = validate_pipeline_runtime_compatibility(
        request=request,
        site_matrix=pd.DataFrame({"sample_1": [1.0]}, index=["PRKACA;S339;"]),
    )

    assert analysis_request is not None
    assert analysis_request.pred_mat.index.tolist() == ["PRKACA;S339;"]
    assert analysis_request.phospho_matrix.index.tolist() == ["PRKACA;S339;"]


def test_validate_pipeline_runtime_compatibility_reports_post_preprocessing_overlap_failures() -> (
    None
):
    total_df = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
            "group2": [1.0],
            "group3": [1.0],
            "group4": [1.0],
            "group5": [1.0],
            "group6": [1.0],
        }
    )
    phospho_df = pd.DataFrame(
        {
            "uid": ["u1"],
            "gene_names": ["PRKACA"],
            "gene_p_site": ["PRKACA_S339"],
            "localization_prob": [0.95],
            "centralized_sequence": ["AAAAAA"],
            "p_group1": [1.0],
            "p_group2": [1.0],
            "p_group3": [1.0],
            "p_group4": [1.0],
            "p_group5": [1.0],
            "p_group6": [1.0],
        }
    )

    request = validate_pipeline_construction_request(
        dataset=PhosphoDataset(total_df=total_df, phospho_df=phospho_df),
        pred_mat=pd.DataFrame({"PRKACA": [0.9]}, index=["PRKACA;S339;"]),
    )

    with pytest.raises(
        InputCompatibilityError,
        match=(
            r"Pipeline runtime compatibility failed after preprocessing: "
            r"pipeline pred_mat and preprocessed site matrix have no overlapping "
            r"phosphosite IDs"
        ),
    ) as exc_info:
        validate_pipeline_runtime_compatibility(
            request=request,
            site_matrix=pd.DataFrame({"sample_1": [1.0]}, index=["OTHER;S1;"]),
        )
    message = str(exc_info.value)
    assert (
        "shared=0, preprocessed site matrix rows=1, pipeline pred_mat rows=1" in message
    )
    assert "min_overlap=1" in message
    assert "min_fraction=0.500" in message
    assert "Next step: confirm both inputs use the same phosphosite ID space" in message


def test_validate_pipeline_runtime_compatibility_skips_when_pipeline_has_no_pred_mat() -> (
    None
):
    total_df = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
            "group2": [1.0],
            "group3": [1.0],
            "group4": [1.0],
            "group5": [1.0],
            "group6": [1.0],
        }
    )
    phospho_df = pd.DataFrame(
        {
            "uid": ["u1"],
            "gene_names": ["PRKACA"],
            "gene_p_site": ["PRKACA_S339"],
            "localization_prob": [0.95],
            "centralized_sequence": ["AAAAAA"],
            "p_group1": [1.0],
            "p_group2": [1.0],
            "p_group3": [1.0],
            "p_group4": [1.0],
            "p_group5": [1.0],
            "p_group6": [1.0],
        }
    )

    request = validate_pipeline_construction_request(
        dataset=PhosphoDataset(total_df=total_df, phospho_df=phospho_df),
        pred_mat=None,
    )

    assert (
        validate_pipeline_runtime_compatibility(
            request=request,
            site_matrix=pd.DataFrame({"sample_1": [1.0]}, index=["PRKACA;S339;"]),
        )
        is None
    )


def test_prediction_request_takes_ownership_of_raw_combined_scores() -> None:
    combined_scores = pd.DataFrame({"KINASE_A": [0.9]}, index=["SITE_1"])

    request = PredictionRequest.validate_request(
        combined_scores=combined_scores,
        ensemble_size=2,
        top=2,
        score_threshold=0.8,
        inclusion=1,
        n_iterations=2,
        debug_top_n=1,
        default_svm_mode="default",
    )

    combined_scores.loc["SITE_1", "KINASE_A"] = 0.1

    assert request.combined_scores.loc["SITE_1", "KINASE_A"] == 0.9


def test_dataset_inputs_bundles_with_pandas_state_are_not_frozen_dataclasses() -> None:
    assert AnalysisInputs.__dataclass_params__.frozen is False
    assert PipelineInputs.__dataclass_params__.frozen is False
    assert WorkflowInputs.__dataclass_params__.frozen is False
    assert DatasetInputs.__dataclass_params__.frozen is False


def test_validate_analysis_request_rejects_zero_overlap() -> None:
    pred_mat = pd.DataFrame({"PRKACA": [0.9]}, index=["SITE_A"])
    phospho_matrix = pd.DataFrame({"sample_1": [1.0]}, index=["SITE_B"])

    with pytest.raises(
        InputCompatibilityError,
        match="no overlapping phosphosite IDs",
    ) as exc_info:
        validate_analysis_request(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
        )
    message = str(exc_info.value)
    assert "shared=0, phospho_matrix rows=1, pred_mat rows=1" in message
    assert "min_overlap=1" in message
    assert "min_fraction=0.500" in message


def test_validate_analysis_request_reports_empty_pred_mat_with_candidate_context() -> (
    None
):
    pred_mat = pd.DataFrame(index=["SITE_1"])
    phospho_matrix = pd.DataFrame({"sample_1": [1.0]}, index=["SITE_1"])

    with pytest.raises(
        NoCandidateKinasesError,
        match="does not contain any kinase columns",
    ) as exc_info:
        validate_analysis_request(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
        )
    message = str(exc_info.value)
    assert "pred_mat columns=0, phosphosite rows=1" in message
    assert "source=analysis input" in message
    assert "less restrictive top, score_threshold, or inclusion settings" in message


def test_validate_analysis_request_records_partial_overlap_above_default_fraction() -> (
    None
):
    pred_mat = pd.DataFrame(
        {"PRKACA": [0.9, 0.8, 0.7, 0.6, 0.5]},
        index=[f"SITE_{idx}" for idx in range(1, 6)],
    )
    phospho_matrix = pd.DataFrame(
        {"sample_1": [1.0] * 8},
        index=[f"SITE_{idx}" for idx in range(1, 9)],
    )

    inputs = validate_analysis_request(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
    )

    assert inputs.pred_mat.index.tolist() == [f"SITE_{idx}" for idx in range(1, 6)]
    assert inputs.overlap_summary.is_partial is True
    assert inputs.overlap_summary.overlap_count == 5
    assert inputs.overlap_summary.message == (
        "pred_mat and phospho_matrix partially overlap: using 5 shared "
        "phosphosite row(s) (62.5% of phospho_matrix; 100.0% of pred_mat)"
    )


def test_validate_analysis_request_allows_nan_pred_mat_values() -> None:
    pred_mat = pd.DataFrame({"PRKACA": [0.9, float("nan")]}, index=["SITE_1", "SITE_2"])
    phospho_matrix = pd.DataFrame(
        {"sample_1": [1.0, 2.0]},
        index=["SITE_1", "SITE_2"],
    )

    inputs = validate_analysis_request(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
    )

    assert pd.isna(inputs.pred_mat.loc["SITE_2", "PRKACA"])


def test_validate_analysis_request_accepts_full_overlap_without_partial_flag() -> None:
    pred_mat = pd.DataFrame({"PRKACA": [0.9, 0.8]}, index=["SITE_1", "SITE_2"])
    phospho_matrix = pd.DataFrame(
        {"sample_1": [1.0, 2.0]},
        index=["SITE_1", "SITE_2"],
    )

    inputs = validate_analysis_request(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
    )

    assert inputs.overlap_summary.is_partial is False
    assert inputs.phospho_matrix.index.tolist() == ["SITE_1", "SITE_2"]


def test_validate_analysis_request_rejects_insufficient_overlap_fraction() -> None:
    pred_mat = pd.DataFrame({"PRKACA": [0.9]}, index=["SITE_1"])
    phospho_matrix = pd.DataFrame(
        {"sample_1": [1.0] * 20},
        index=[f"SITE_{idx}" for idx in range(1, 21)],
    )

    with pytest.raises(
        InputCompatibilityError,
        match="insufficient overlapping phosphosite IDs",
    ) as exc_info:
        validate_analysis_request(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
        )
    message = str(exc_info.value)
    assert "shared=1, phospho_matrix rows=20, pred_mat rows=1" in message
    assert "phospho_matrix coverage=5.0%" in message
    assert "pred_mat coverage=100.0%" in message
