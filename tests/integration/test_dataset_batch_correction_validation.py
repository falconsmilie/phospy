from __future__ import annotations

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import DifferentialAnalysisWorkflow
from phospy.advanced import (
    DatasetBatchCorrectionConfig,
    DatasetIntensityTransformConfig,
)
from phospy.api import (
    AnalysisReadyDatasetBuilder,
    Contrast,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.errors import PhosPyInputError

pytestmark = pytest.mark.integration


def test_dataset_batch_correction_noop_path_leaves_matrix_unchanged() -> None:
    phospho = _log2_batch_effect_phospho()

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho),
            sample_metadata=_valid_sample_metadata(phospho),
            organism=Organism.RAT,
            input_intensity_scale="log2",
        )
    )

    np.testing.assert_allclose(built.phospho.to_numpy(), phospho.to_numpy())
    assert built.preprocessing_report is not None
    report = built.preprocessing_report.batch_correction
    assert report is not None
    assert report.status == "disabled"
    assert report.confounding_check_status == "not_applicable"
    assert report.batch_levels == ("run_1", "run_2")
    assert report.condition_levels == ("control", "treated")
    assert "batch_correction" not in set(
        built.preprocessing_report.operations.loc[:, "stage"].astype(str).tolist()
    )


def test_dataset_batch_correction_successful_build_applies_corrected_matrix() -> None:
    log2_matrix = _log2_batch_effect_phospho()
    phospho = _linear_from_log2(log2_matrix)

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho),
            sample_metadata=_valid_sample_metadata(phospho),
            organism=Organism.RAT,
            preprocessing_config=_log2_batch_correction_config(),
        )
    )

    expected = _expected_batch_corrected_log2_matrix(
        index=built.phospho.index,
        columns=built.phospho.columns,
    )
    pdt.assert_frame_equal(
        built.phospho,
        expected,
        check_exact=False,
        atol=1e-10,
        rtol=0.0,
    )
    assert built.preprocessing_report is not None
    report = built.preprocessing_report.batch_correction
    assert report is not None
    assert report.status == "applied"
    assert report.confounding_check_status == "passed"
    assert report.batch_levels == ("run_1", "run_2")
    assert report.condition_levels == ("control", "treated")
    assert "batch_correction" in set(
        built.preprocessing_report.operations.loc[:, "stage"].astype(str).tolist()
    )


def test_dataset_batch_correction_confounding_fails_clearly_during_build() -> None:
    phospho = _linear_from_log2(_log2_batch_effect_phospho())

    with pytest.raises(
        PhosPyInputError,
        match="batch and condition are perfectly confounded",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_confounded_sample_metadata(phospho),
                organism=Organism.RAT,
                preprocessing_config=_log2_batch_correction_config(),
            )
        )


def test_dataset_batch_correction_missing_metadata_fails_during_build() -> None:
    phospho = _linear_from_log2(_log2_batch_effect_phospho())

    with pytest.raises(
        PhosPyInputError,
        match="requires sample_metadata input data",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=None,
                organism=Organism.RAT,
                preprocessing_config=_log2_batch_correction_config(),
            )
        )


def test_dataset_batch_correction_report_includes_applied_stage_report() -> None:
    phospho = _linear_from_log2(_log2_batch_effect_phospho())

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho),
            sample_metadata=_valid_sample_metadata(phospho),
            organism=Organism.RAT,
            preprocessing_config=_log2_batch_correction_config(),
        )
    )

    assert built.preprocessing_report is not None
    report = built.preprocessing_report.batch_correction
    assert report is not None
    assert report.to_payload()["status"] == "applied"
    assert report.matrix_shape_before == (2, 4)
    assert report.matrix_shape_after == (2, 4)
    operation = built.preprocessing_report.operations.loc[
        built.preprocessing_report.operations.loc[:, "stage"] == "batch_correction"
    ].iloc[0]
    assert operation["operation"] == "linear_residualize_batch"


def test_dataset_batch_correction_downstream_differential_uses_corrected_matrix() -> (
    None
):
    phospho = _linear_from_log2(_log2_batch_effect_phospho())
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho),
            sample_metadata=_valid_sample_metadata(phospho),
            organism=Organism.RAT,
            preprocessing_config=_log2_batch_correction_config(),
        )
    )

    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=_differential_design(),
            contrasts=(
                Contrast(
                    name="treated_vs_control",
                    numerator_condition="treated",
                    denominator_condition="control",
                ),
            ),
        )
    )

    table = result.table_for("treated_vs_control")
    site_one = dataset.phospho.index[0]
    site_two = dataset.phospho.index[1]
    assert float(table.loc[site_one, "logFC"]) == pytest.approx(4.0)
    assert float(table.loc[site_two, "logFC"]) == pytest.approx(-1.0)


def _log2_batch_effect_phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_1": [9.9, 1.9],
            "sample_2": [14.1, 1.1],
            "sample_3": [15.1, 7.1],
            "sample_4": [18.9, 5.9],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )


def _linear_from_log2(log2_matrix: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        np.power(2.0, log2_matrix.to_numpy(dtype="float64")) - 1.0,
        index=log2_matrix.index.copy(),
        columns=log2_matrix.columns.copy(),
    )


def _expected_batch_corrected_log2_matrix(
    *,
    index: pd.Index,
    columns: pd.Index,
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_1": [9.9, 1.9],
            "sample_2": [14.1, 1.1],
            "sample_3": [10.1, 2.1],
            "sample_4": [13.9, 0.9],
        },
        index=index.copy(),
        columns=columns.copy(),
    )


def _site_metadata(phospho: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "protein_id": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [
                ("A" * 15) + "Y" + ("A" * 15),
                ("A" * 15) + "T" + ("A" * 15),
            ],
            "localisation_confidence": [0.95, 0.92],
        },
        index=phospho.index.copy(),
    )


def _valid_sample_metadata(phospho: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": ["run_1", "run_1", "run_2", "run_2"],
            "condition": ["control", "treated", "control", "treated"],
        },
        index=phospho.columns.copy(),
    )


def _confounded_sample_metadata(phospho: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": ["run_1", "run_1", "run_2", "run_2"],
            "condition": ["control", "control", "treated", "treated"],
        },
        index=phospho.columns.copy(),
    )


def _log2_batch_correction_config() -> DatasetPreprocessingConfig:
    return DatasetPreprocessingConfig(
        intensity_transform=DatasetIntensityTransformConfig(
            policy="log2",
            pseudocount=1.0,
        ),
        batch_correction=DatasetBatchCorrectionConfig(
            method="linear_residualize_batch"
        ),
    )


def _differential_design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="sample_1",
                condition="control",
                biological_replicate_id="control_1",
            ),
            SampleDesignRecord(
                sample_id="sample_2",
                condition="treated",
                biological_replicate_id="treated_1",
            ),
            SampleDesignRecord(
                sample_id="sample_3",
                condition="control",
                biological_replicate_id="control_2",
            ),
            SampleDesignRecord(
                sample_id="sample_4",
                condition="treated",
                biological_replicate_id="treated_2",
            ),
        )
    )
