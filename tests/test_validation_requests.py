from __future__ import annotations

import pandas as pd
import pytest

from phospy import DatasetSchema
from phospy.prediction.traces import TraceSink
from phospy.validation.errors import RequestValidationError
from phospy.validation.requests import (
    CorePipelineRequest,
    KinaseWorkflowRequest,
    PredictionRequest,
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
    total_path.write_text("genes	group1\nPRKACA	1\n")
    phospho_path.write_text("uid	gene_names\n1	PRKACA\n")

    def blow_up(self: DatasetSchema, comparisons, *, context: str = "Dataset schema"):
        raise RuntimeError("comparison validator exploded")

    monkeypatch.setattr(DatasetSchema, "validate_comparisons", blow_up)

    with pytest.raises(RuntimeError, match="comparison validator exploded"):
        CorePipelineRequest.validate_request(
            total_path=total_path,
            phospho_path=phospho_path,
            comparisons=(("group1", "group1"),),
        )
