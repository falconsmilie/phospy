from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from phospy import (
    CoreOutputWriter,
    DatasetPreprocessing,
    DatasetSchema,
    DatasetSiteMatrix,
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
        "CoreOutputWriter",
        "CoreOutputs",
        "CoreProcessingResult",
        "DatasetPreprocessing",
        "DatasetSchema",
        "DatasetSiteMatrix",
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


def test_phospho_dataset_preprocessing_run() -> None:
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
        comparisons=EXAMPLE_COMPARISONS,
    )
    result = dataset.preprocessing.run()

    assert sorted(result.total_unique["genes"].tolist()) == ["BTK", "LYN", "PRKACA"]
    assert "p_group1_group4" in result.phospho_corrected.columns
    assert "PRKACA;S339;" in result.site_matrix.matrix.index
    assert result.site_matrix.row_drop_stats["retained_rows"] == len(
        result.site_matrix.matrix
    )


def test_phospho_dataset_preprocessing_facade_is_bound_and_typed() -> None:
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
        comparisons=EXAMPLE_COMPARISONS,
    )

    preprocessing = dataset.preprocessing

    assert isinstance(preprocessing, DatasetPreprocessing)
    assert preprocessing.total_df.equals(dataset.total_df)
    assert preprocessing.phospho_df.equals(dataset.phospho_df)
    assert preprocessing.schema == dataset.schema
    assert preprocessing.comparisons == tuple(EXAMPLE_COMPARISONS)


def test_phospho_dataset_site_matrix_facade_is_bound_and_typed() -> None:
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
    )

    site_matrix = dataset.site_matrix

    assert isinstance(site_matrix, DatasetSiteMatrix)
    assert site_matrix.schema == dataset.schema


def test_phospho_dataset_site_matrix_build_matches_run_output() -> None:
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
        comparisons=EXAMPLE_COMPARISONS,
    )

    core = dataset.preprocessing.run(max_unmatched_fraction=0.1)
    via_service = dataset.site_matrix.build(core.phospho_corrected)

    pd.testing.assert_frame_equal(via_service.phosr_input, core.site_matrix.phosr_input)
    pd.testing.assert_frame_equal(via_service.matrix, core.site_matrix.matrix)
    pd.testing.assert_series_equal(via_service.sequences, core.site_matrix.sequences)
    assert via_service.row_drop_stats == core.site_matrix.row_drop_stats


def test_phospho_dataset_preprocessing_run_matches_core_processor() -> None:
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
        comparisons=EXAMPLE_COMPARISONS,
    )

    via_facade = dataset.preprocessing.run(max_unmatched_fraction=0.1)
    via_processor = CoreProcessor(
        schema=dataset.schema,
        comparisons=dataset.comparisons,
    ).process(
        dataset.total_df,
        dataset.phospho_df,
        config=CorePreprocessingConfig(max_unmatched_fraction=0.1),
    )

    pd.testing.assert_frame_equal(via_facade.total_unique, via_processor.total_unique)
    pd.testing.assert_frame_equal(
        via_facade.total_filtered,
        via_processor.total_filtered,
    )
    pd.testing.assert_frame_equal(
        via_facade.phospho_filtered,
        via_processor.phospho_filtered,
    )
    pd.testing.assert_frame_equal(
        via_facade.phospho_corrected,
        via_processor.phospho_corrected,
    )
    pd.testing.assert_frame_equal(
        via_facade.site_matrix.phosr_input,
        via_processor.site_matrix.phosr_input,
    )
    pd.testing.assert_frame_equal(
        via_facade.site_matrix.matrix,
        via_processor.site_matrix.matrix,
    )
    pd.testing.assert_series_equal(
        via_facade.site_matrix.sequences,
        via_processor.site_matrix.sequences,
    )
    assert (
        via_facade.site_matrix.row_drop_stats
        == via_processor.site_matrix.row_drop_stats
    )


def test_phospho_dataset_does_not_expose_legacy_direct_preprocessing_methods() -> None:
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
    )

    assert not hasattr(dataset, "prepare_total")
    assert not hasattr(dataset, "prepare_phospho")
    assert not hasattr(dataset, "correct_to_protein")
    assert not hasattr(dataset, "add_pairwise_comparisons")
    assert not hasattr(dataset, "process_core")
    assert not hasattr(dataset, "build_site_matrix")
    assert not hasattr(dataset, "write_core_outputs")
    assert not hasattr(dataset.preprocessing, "build_site_matrix")


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


def test_phospho_dataset_defensively_copies_constructor_inputs() -> None:
    total_df = make_total_df()
    phospho_df = make_phospho_df()

    dataset = PhosphoDataset(
        total_df=total_df,
        phospho_df=phospho_df,
        comparisons=EXAMPLE_COMPARISONS,
    )

    total_df.loc[0, "group1"] = 999.0
    phospho_df.loc[0, "p_group1"] = 999.0

    assert dataset.total_df.loc[0, "group1"] != 999.0
    assert dataset.phospho_df.loc[0, "p_group1"] != 999.0


def test_phospho_dataset_accessors_return_deep_copies() -> None:
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
        comparisons=EXAMPLE_COMPARISONS,
    )

    total_view = dataset.total_df
    phospho_view = dataset.phospho_df

    total_view.loc[0, "group1"] = 999.0
    phospho_view.loc[0, "p_group1"] = 999.0

    assert dataset.total_df.loc[0, "group1"] != 999.0
    assert dataset.phospho_df.loc[0, "p_group1"] != 999.0
    assert dataset.inputs.total_df.loc[0, "group1"] != 999.0
    assert dataset.inputs.phospho_df.loc[0, "p_group1"] != 999.0


def test_core_output_writer_writes_csv_outputs(tmp_path) -> None:
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
        comparisons=EXAMPLE_COMPARISONS,
    )
    core = dataset.preprocessing.run(max_unmatched_fraction=0.1)

    outdir = tmp_path / "core-output-csv"
    CoreOutputWriter().write(core, outdir)

    expected_files = {
        "df_phospho_corrected.csv",
        "df_phospho_filtered.csv",
        "df_total_filtered.csv",
        "df_total_unique.csv",
        "mat_phospho_corrected.csv",
        "phosr_input.csv",
        "site_sequences.csv",
    }
    assert expected_files == {path.name for path in outdir.iterdir()}


def test_core_output_writer_writes_tsv_outputs(tmp_path) -> None:
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
        comparisons=EXAMPLE_COMPARISONS,
    )
    core = dataset.preprocessing.run(max_unmatched_fraction=0.1)

    outdir = tmp_path / "core-output-tsv"
    CoreOutputWriter().write(core, outdir, format="tsv")

    expected_files = {
        "df_phospho_corrected.tsv",
        "df_phospho_filtered.tsv",
        "df_total_filtered.tsv",
        "df_total_unique.tsv",
        "mat_phospho_corrected.tsv",
        "phosr_input.tsv",
        "site_sequences.tsv",
    }
    assert expected_files == {path.name for path in outdir.iterdir()}

    total_unique = pd.read_csv(outdir / "df_total_unique.tsv", sep="\t")
    assert total_unique["genes"].tolist() == ["BTK", "LYN", "PRKACA"]


def test_core_output_writer_writes_parquet_outputs(monkeypatch, tmp_path) -> None:
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
        comparisons=EXAMPLE_COMPARISONS,
    )
    core = dataset.preprocessing.run(max_unmatched_fraction=0.1)
    written_paths: list[tuple[Path, bool]] = []

    def fake_to_parquet(
        self: pd.DataFrame,
        path: str | Path,
        *,
        index: bool = True,
        **_: object,
    ) -> None:
        destination = Path(path)
        written_paths.append((destination, index))
        destination.write_text("parquet stub", encoding="utf-8")

    monkeypatch.setattr(pd.DataFrame, "to_parquet", fake_to_parquet)

    outdir = tmp_path / "core-output-parquet"
    CoreOutputWriter().write(core, outdir, format="parquet")

    assert {path.name for path, _ in written_paths} == {
        "df_phospho_corrected.parquet",
        "df_phospho_filtered.parquet",
        "df_total_filtered.parquet",
        "df_total_unique.parquet",
        "mat_phospho_corrected.parquet",
        "phosr_input.parquet",
        "site_sequences.parquet",
    }
    index_flags = {path.name: index for path, index in written_paths}
    assert index_flags["df_total_unique.parquet"] is False
    assert index_flags["mat_phospho_corrected.parquet"] is True
    assert index_flags["site_sequences.parquet"] is True


def test_core_output_writer_rejects_unknown_format(tmp_path) -> None:
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
        comparisons=EXAMPLE_COMPARISONS,
    )
    core = dataset.preprocessing.run(max_unmatched_fraction=0.1)

    with pytest.raises(ValueError, match="Unsupported core output format"):
        CoreOutputWriter().write(core, tmp_path / "core-output-unknown", format="json")


def test_core_output_writer_raises_clear_error_without_parquet_engine(
    monkeypatch,
    tmp_path,
) -> None:
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
        comparisons=EXAMPLE_COMPARISONS,
    )
    core = dataset.preprocessing.run(max_unmatched_fraction=0.1)

    def blow_up(
        self: pd.DataFrame,
        path: str | Path,
        *,
        index: bool = True,
        **_: object,
    ) -> None:
        raise ImportError("missing engine")

    monkeypatch.setattr(pd.DataFrame, "to_parquet", blow_up)

    with pytest.raises(RuntimeError, match="optional pandas parquet engine"):
        CoreOutputWriter().write(
            core, tmp_path / "core-output-parquet", format="parquet"
        )


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

    total_unique, total_filtered = dataset.preprocessing.prepare_total()
    phospho_filtered = dataset.preprocessing.prepare_phospho()
    corrected = dataset.preprocessing.correct_to_protein(
        phospho_filtered,
        total_filtered,
    )

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

    total_unique, total_filtered = pipeline.dataset.preprocessing.prepare_total(
        sentinel=pipeline.preprocessing_config.total_sentinel,
        min_observed=pipeline.preprocessing_config.min_observed,
    )
    phospho_filtered = pipeline.dataset.preprocessing.prepare_phospho(
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
