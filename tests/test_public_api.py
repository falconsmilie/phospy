from __future__ import annotations

import json

import pandas as pd
import pytest

from phospy import (
    DatasetSchema,
    KinaseActivityAnalyzer,
    KinaseWorkflow,
    PhosphoDataset,
    PhosRPipeline,
)
from phospy.core_processing import CorePreprocessingConfig, CoreProcessor
from phospy.dataset_loader import DatasetLoader
from phospy.site_matrix_builder import SiteMatrixBuilder

EXAMPLE_COMPARISONS = [("group1", "group4"), ("group2", "group5"), ("group3", "group6")]


def make_total_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "genes": ["Prkaca", "Prkaca", "Btk", "Lyn"],
            "group1": [1.0, 5.0, 2.0, 3.0],
            "group2": [1.0, 5.0, 2.0, 3.0],
            "group3": [1.0, 5.0, 2.0, 3.0],
            "group4": [1.0, 5.0, 2.0, 3.0],
            "group5": [1.0, 5.0, 2.0, 3.0],
            "group6": [1.0, 5.0, 2.0, 3.0],
        }
    )


def make_phospho_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "uid": ["u1", "u2", "u3", "u4"],
            "gene_names": ["PRKACA", "BTK", "LYN", "PRKACA"],
            "gene_p_site": ["PRKACA_S339", "BTK_Y551", "LYN_Y397", "PRKACA_S339"],
            "localization_prob": [0.95, 0.95, 0.95, 0.95],
            "centralized_sequence": ["AAAAAA", "BBBBBB", "CCCCCC", "DDDDDD"],
            "p_group1": [8.0, 6.0, 7.0, 9.0],
            "p_group2": [7.0, 5.0, 6.0, 8.0],
            "p_group3": [6.0, 4.0, 5.0, 7.0],
            "p_group4": [5.0, 3.0, 4.0, 6.0],
            "p_group5": [4.0, 2.0, 3.0, 5.0],
            "p_group6": [3.0, 1.0, 2.0, 4.0],
        }
    )


def make_pred_mat() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "PRKACA": [0.9, 0.8, 0.7],
            "BTK": [0.2, 0.85, 0.75],
        },
        index=["PRKACA;S339;", "BTK;Y551;", "LYN;Y397;"],
    )


def test_public_root_exports() -> None:
    import phospy

    expected = {
        "CoreOutputs",
        "CoreProcessingResult",
        "DatasetSchema",
        "KinaseActivityAnalyzer",
        "KinaseActivityResult",
        "KinasePredictionResult",
        "KinaseWorkflow",
        "KinaseWorkflowResult",
        "PhosphoDataset",
        "PhosRPipeline",
        "SiteMatrixResult",
    }
    assert set(phospy.__all__) == expected


def test_phospho_dataset_process_core() -> None:
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
        comparisons=EXAMPLE_COMPARISONS,
    )
    result = dataset.process_core()

    assert sorted(result.total_unique["genes"].tolist()) == ["BTK", "LYN", "PRKACA"]
    assert "p_group1_group4" in result.phospho_corrected.columns
    assert "PRKACA;S339;" in result.site_matrix.matrix.index
    assert result.site_matrix.row_drop_stats["retained_rows"] == len(
        result.site_matrix.matrix
    )


def test_phospho_dataset_from_validated_inputs_builds_without_revalidation() -> None:
    loader = DatasetLoader()
    total_df, phospho_df = loader.validate(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
    )

    dataset = PhosphoDataset.from_validated_inputs(
        total_df=total_df,
        phospho_df=phospho_df,
        comparisons=EXAMPLE_COMPARISONS,
    )

    assert dataset.total_df.equals(total_df)
    assert dataset.phospho_df.equals(phospho_df)
    assert dataset.comparisons == tuple(EXAMPLE_COMPARISONS)


def test_kinase_activity_analyzer_analyze() -> None:
    phospho_matrix = pd.DataFrame(
        {
            "phospho_corrected_1": [4.0, 4.0, 4.0],
            "phospho_corrected_2": [5.0, 5.0, 5.0],
        },
        index=["PRKACA;S339;", "BTK;Y551;", "LYN;Y397;"],
    )
    result = KinaseActivityAnalyzer().analyze(
        pred_mat=make_pred_mat(),
        phospho_matrix=phospho_matrix,
        threshold=0.6,
        min_substrates=2,
    )

    assert set(result.weighted_activity.index) == {"PRKACA", "BTK"}


def test_kinase_activity_analyzer_load_and_analyze(tmp_path) -> None:
    pred_mat_path = tmp_path / "predMat.csv"
    make_pred_mat().to_csv(pred_mat_path)
    phospho_matrix = pd.DataFrame(
        {
            "phospho_corrected_1": [4.0, 4.0, 4.0],
            "phospho_corrected_2": [5.0, 5.0, 5.0],
        },
        index=["PRKACA;S339;", "BTK;Y551;", "LYN;Y397;"],
    )

    result = KinaseActivityAnalyzer().load_and_analyze(
        pred_mat_path=pred_mat_path,
        phospho_matrix=phospho_matrix,
        threshold=0.6,
        min_substrates=2,
    )

    assert set(result.weighted_activity.index) == {"PRKACA", "BTK"}


def test_kinase_activity_analyzer_write_outputs(tmp_path) -> None:
    phospho_matrix = pd.DataFrame(
        {
            "phospho_corrected_1": [4.0, 4.0, 4.0],
            "phospho_corrected_2": [5.0, 5.0, 5.0],
        },
        index=["PRKACA;S339;", "BTK;Y551;", "LYN;Y397;"],
    )
    result = KinaseActivityAnalyzer().analyze(
        pred_mat=make_pred_mat(),
        phospho_matrix=phospho_matrix,
        threshold=0.6,
        min_substrates=2,
    )

    outdir = tmp_path / "kinase-output"
    KinaseActivityAnalyzer().write_outputs(result, outdir)

    expected_files = {
        "kinase_activity_matrix.csv",
        "kinase_target_counts.csv",
        "kinase_target_table.csv",
        "ksea_counts.csv",
        "ksea_scores.csv",
    }
    assert expected_files.issubset({path.name for path in outdir.iterdir()})


def test_pipeline_runs_with_class_api(tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    pred_path = tmp_path / "predMat.csv"
    outdir = tmp_path / "out"

    make_total_df().to_csv(total_path, sep="\t", index=False)
    make_phospho_df().to_csv(phospho_path, sep="\t", index=False)
    make_pred_mat().to_csv(pred_path)

    pipeline = PhosRPipeline.from_files(
        total_path=total_path,
        phospho_path=phospho_path,
        pred_mat_path=pred_path,
    )
    outputs = pipeline.run(outdir=outdir)

    assert outputs.kinase_activity is not None


def test_kinase_workflow_runs_with_class_api() -> None:
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 2.0, 10.0, 11.0],
            "sample_2": [1.5, 2.5, 10.5, 11.5],
        },
        index=["SITE_1", "SITE_2", "SITE_3", "SITE_4"],
    )

    workflow = KinaseWorkflow()
    result = workflow.run(
        phospho_matrix=phospho_matrix,
        substrate_map={
            "KINASE_A": ["SITE_1", "SITE_2"],
            "KINASE_B": ["SITE_3", "SITE_4"],
        },
        site_sequences={
            "SITE_1": "QQAAAAAYY",
            "SITE_2": "QQAAAAAYY",
            "SITE_3": "QQTTTTTYY",
            "SITE_4": "QQTTTTTYY",
        },
        motif_sequences={
            "KINASE_A": ["QQAAAAAYY", "QQAAAAAYY"],
            "KINASE_B": ["QQTTTTTYY", "QQTTTTTYY"],
        },
        min_substrates=2,
        min_motif_size=1,
        top=2,
        score_threshold=0.5,
        inclusion=1,
        ensemble_size=2,
        n_iterations=1,
        random_state=7,
    )

    assert result.prediction_result.pred_matrix.shape[1] == 2


def test_pipeline_propagates_max_unmatched_fraction(tmp_path) -> None:
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
            "uid": ["u1", "u2"],
            "gene_names": ["PRKACA", "BTK"],
            "gene_p_site": ["PRKACA_S339", "BTK_Y551"],
            "localization_prob": [0.95, 0.95],
            "centralized_sequence": ["AAAAAA", "BBBBBB"],
            "p_group1": [8.0, 6.0],
            "p_group2": [7.0, 5.0],
            "p_group3": [6.0, 4.0],
            "p_group4": [5.0, 3.0],
            "p_group5": [4.0, 2.0],
            "p_group6": [3.0, 1.0],
        }
    )

    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    total_df.to_csv(total_path, sep="	", index=False)
    phospho_df.to_csv(phospho_path, sep="	", index=False)

    pipeline = PhosRPipeline.from_files(
        total_path=total_path,
        phospho_path=phospho_path,
        max_unmatched_fraction=0.5,
    )
    outputs = pipeline.run()

    assert outputs.core.phospho_corrected.shape[0] == 1
    assert pipeline.preprocessing_config.max_unmatched_fraction == 0.5


def test_phospho_dataset_honours_custom_corrected_columns() -> None:
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
        schema=DatasetSchema(
            corrected_cols=(
                "sample_a",
                "sample_b",
                "sample_c",
                "sample_d",
                "sample_e",
                "sample_f",
            ),
        ),
    )

    total_unique, total_filtered = dataset.prepare_total()
    phospho_filtered = dataset.prepare_phospho()
    corrected = dataset.correct_to_protein(phospho_filtered, total_filtered)

    expected = {"sample_a", "sample_b", "sample_c", "sample_d", "sample_e", "sample_f"}
    assert expected.issubset(corrected.columns)
    assert "phospho_corrected_1" not in corrected.columns


def test_pipeline_propagates_sentinel_configuration(tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"

    total_df = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [99.0],
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
            "p_group1": [99.0],
            "p_group2": [7.0],
            "p_group3": [6.0],
            "p_group4": [5.0],
            "p_group5": [4.0],
            "p_group6": [3.0],
        }
    )
    total_df.to_csv(total_path, sep="	", index=False)
    phospho_df.to_csv(phospho_path, sep="	", index=False)

    pipeline = PhosRPipeline.from_files(
        total_path=total_path,
        phospho_path=phospho_path,
        min_observed=5,
        total_sentinel=99.0,
        phospho_sentinel=99.0,
    )

    total_unique, total_filtered = pipeline.dataset.prepare_total(
        sentinel=pipeline.preprocessing_config.total_sentinel,
        min_observed=pipeline.preprocessing_config.min_observed,
    )
    phospho_filtered = pipeline.dataset.prepare_phospho(
        localization_threshold=pipeline.preprocessing_config.localization_threshold,
        sentinel=pipeline.preprocessing_config.phospho_sentinel,
        min_observed=pipeline.preprocessing_config.min_observed,
    )

    assert pipeline.preprocessing_config.total_sentinel == 99.0
    assert pipeline.preprocessing_config.phospho_sentinel == 99.0
    assert total_unique.shape[0] == 1
    assert pd.isna(total_filtered.iloc[0]["group1"])
    assert pd.isna(phospho_filtered.iloc[0]["p_group1"])


def test_dataset_components_work_together() -> None:
    loader = DatasetLoader()
    total_df, phospho_df = loader.validate(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
    )
    processor = CoreProcessor(
        schema=DatasetSchema(),
        site_matrix_builder=SiteMatrixBuilder(
            value_cols=[
                "phospho_corrected_1",
                "phospho_corrected_2",
                "phospho_corrected_3",
                "phospho_corrected_4",
                "phospho_corrected_5",
                "phospho_corrected_6",
            ]
        ),
    )

    result = processor.process(total_df, phospho_df, config=pipeline_config())

    assert "PRKACA;S339;" in result.site_matrix.matrix.index


def test_dataset_from_files_validates_inputs_once(monkeypatch, tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    make_total_df().to_csv(total_path, sep="	", index=False)
    make_phospho_df().to_csv(phospho_path, sep="	", index=False)

    from phospy.validation.tables import PhosphoInputSchema, TotalInputSchema

    total_calls = 0
    phospho_calls = 0
    original_total_validate = TotalInputSchema.validate
    original_phospho_validate = PhosphoInputSchema.validate

    def counting_total_validate(*args, **kwargs):
        nonlocal total_calls
        total_calls += 1
        return original_total_validate(*args, **kwargs)

    def counting_phospho_validate(*args, **kwargs):
        nonlocal phospho_calls
        phospho_calls += 1
        return original_phospho_validate(*args, **kwargs)

    monkeypatch.setattr(TotalInputSchema, "validate", counting_total_validate)
    monkeypatch.setattr(PhosphoInputSchema, "validate", counting_phospho_validate)

    dataset = PhosphoDataset.from_files(
        total_path=total_path, phospho_path=phospho_path
    )

    assert dataset.total_df.shape[0] == 4
    assert total_calls == 1
    assert phospho_calls == 1


def test_pipeline_from_request_validates_inputs_once(monkeypatch, tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    make_total_df().to_csv(total_path, sep="	", index=False)
    make_phospho_df().to_csv(phospho_path, sep="	", index=False)

    from phospy.validation.requests import CorePipelineRequest
    from phospy.validation.tables import PhosphoInputSchema, TotalInputSchema

    total_calls = 0
    phospho_calls = 0
    original_total_validate = TotalInputSchema.validate
    original_phospho_validate = PhosphoInputSchema.validate

    def counting_total_validate(*args, **kwargs):
        nonlocal total_calls
        total_calls += 1
        return original_total_validate(*args, **kwargs)

    def counting_phospho_validate(*args, **kwargs):
        nonlocal phospho_calls
        phospho_calls += 1
        return original_phospho_validate(*args, **kwargs)

    monkeypatch.setattr(TotalInputSchema, "validate", counting_total_validate)
    monkeypatch.setattr(PhosphoInputSchema, "validate", counting_phospho_validate)

    request = CorePipelineRequest.validate_request(
        total_path=total_path,
        phospho_path=phospho_path,
    )
    pipeline = PhosRPipeline.from_request(request)

    assert pipeline.dataset.total_df.shape[0] == 4
    assert total_calls == 1
    assert phospho_calls == 1


def test_pipeline_writes_run_manifest(tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    pred_path = tmp_path / "predMat.csv"
    outdir = tmp_path / "out"

    make_total_df().to_csv(total_path, sep="	", index=False)
    make_phospho_df().to_csv(phospho_path, sep="	", index=False)
    make_pred_mat().to_csv(pred_path)

    pipeline = PhosRPipeline.from_files(
        total_path=total_path,
        phospho_path=phospho_path,
        pred_mat_path=pred_path,
    )

    pipeline.run(outdir=outdir)

    manifest = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "success"
    assert manifest["has_kinase_activity"] is True
    assert manifest["core_rows"]["site_matrix"] == 3
    assert manifest["preprocessing_config"]["min_observed"] == 4


def test_pipeline_run_does_not_publish_partial_outputs_on_failure(
    monkeypatch, tmp_path
) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    pred_path = tmp_path / "predMat.csv"
    outdir = tmp_path / "out"

    make_total_df().to_csv(total_path, sep="	", index=False)
    make_phospho_df().to_csv(phospho_path, sep="	", index=False)
    make_pred_mat().to_csv(pred_path)

    pipeline = PhosRPipeline.from_files(
        total_path=total_path,
        phospho_path=phospho_path,
        pred_mat_path=pred_path,
    )

    from phospy import pipeline as pipeline_module

    def blow_up(*args, **kwargs) -> None:
        raise RuntimeError("boom")

    monkeypatch.setattr(
        pipeline_module.KinaseActivityWriter,
        "write",
        staticmethod(blow_up),
    )

    with pytest.raises(RuntimeError, match="boom"):
        pipeline.run(outdir=outdir)

    assert not outdir.exists()
    assert not any(path.name.startswith(".out.tmp-") for path in tmp_path.iterdir())


def pipeline_config() -> CorePreprocessingConfig:
    return CorePreprocessingConfig()
