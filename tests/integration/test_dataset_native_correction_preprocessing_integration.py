from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import DifferentialAnalysisWorkflow
from phospy.api import (
    AnalysisReadyDatasetBuilder,
    Contrast,
    DatasetBatchCorrectionConfig,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.errors import PhosPyInputError
from phospy.science.datasets.builders.interpreter import DatasetBuildRequestInterpreter
from phospy.science.datasets.preprocessing import (
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
    CorrectedPreprocessingOutput,
)

pytestmark = pytest.mark.integration


def test_resolved_native_correction_output_builds_analysis_ready_dataset() -> None:
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho),
            sample_metadata=_sample_metadata(phospho),
            organism=Organism.RAT,
            input_intensity_scale="linear",
            corrected_preprocessing_output=_correction_output(corrected),
        )
    )

    pdt.assert_frame_equal(dataset.phospho, corrected)
    assert dataset.preprocessing_report is not None
    report = dataset.preprocessing_report.batch_correction
    assert report is not None
    assert report.status == "applied"
    assert report.method == "sps_ruv_style"
    assert "batch_correction" in set(
        dataset.preprocessing_report.operations.loc[:, "stage"].astype(str)
    )
    assert dataset.provenance is not None
    assert "batch_correction" in {
        stage.stage for stage in dataset.provenance.preprocessing_stages
    }
    observed_mask = dataset.imputation_observed_mask_dataframe()
    assert observed_mask is not None
    assert observed_mask.to_numpy(dtype=bool).all()


def test_invalid_resolved_correction_output_is_rejected_at_dataset_build() -> None:
    phospho = _phospho()
    misaligned = phospho.copy(deep=True)
    misaligned.index = pd.Index(["other_a", "other_b"], name=phospho.index.name)

    with pytest.raises(
        PhosPyInputError,
        match="corrected_matrix.index must match dataset.phospho.index",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(misaligned),
            )
        )


def test_resolved_correction_output_cannot_be_applied_twice() -> None:
    phospho = _phospho()

    with pytest.raises(
        PhosPyInputError, match="correction must be applied exactly once"
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                preprocessing_config=DatasetPreprocessingConfig(
                    batch_correction=DatasetBatchCorrectionConfig(
                        method="linear_residualize_batch"
                    )
                ),
                corrected_preprocessing_output=_correction_output(phospho + 1.0),
            )
        )


def test_resolved_correction_after_downstream_consumption_is_rejected() -> None:
    phospho = _phospho()

    with pytest.raises(PhosPyInputError, match="already been consumed"):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    phospho + 1.0,
                    consumed_by_downstream=True,
                ),
            )
        )


def test_downstream_workflow_consumes_corrected_analysis_ready_dataset_only() -> None:
    phospho = _phospho()
    corrected = pd.DataFrame(
        {
            "sample_1": [10.0, 2.0],
            "sample_2": [14.0, 1.0],
            "sample_3": [10.0, 2.0],
            "sample_4": [14.0, 1.0],
        },
        index=_resolved_correction_matrix(phospho).index.copy(),
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho),
            sample_metadata=_sample_metadata(phospho),
            organism=Organism.RAT,
            input_intensity_scale="log2",
            corrected_preprocessing_output=_correction_output(corrected),
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
    site_one, site_two = dataset.phospho.index.tolist()
    assert float(table.loc[site_one, "logFC"]) == pytest.approx(4.0)
    assert float(table.loc[site_two, "logFC"]) == pytest.approx(-1.0)


def _correction_output(
    corrected: pd.DataFrame,
    *,
    consumed_by_downstream: bool = False,
) -> CorrectedPreprocessingOutput:
    return CorrectedPreprocessingOutput(
        corrected_matrix=corrected,
        output_observation_mask=pd.DataFrame(
            True,
            index=corrected.index.copy(),
            columns=corrected.columns.copy(),
        ),
        corrected_cell_status=pd.DataFrame(
            "corrected_observed",
            index=corrected.index.copy(),
            columns=corrected.columns.copy(),
        ),
        batch_correction_report=BatchCorrectionReport(
            status="applied",
            policy=BatchCorrectionPolicy(
                method="sps_ruv_style",
                batch_column="batch",
                condition_column="condition",
            ),
            diagnostics=BatchCorrectionDiagnostics(
                number_of_batches=2,
                batch_levels=("run_1", "run_2"),
                condition_levels=("control", "treated"),
                confounding_check_status="passed",
                matrix_shape_before=(2, 4),
                matrix_shape_after=(2, 4),
            ),
        ),
        diagnostics={"executor": {"status": "applied", "method": "sps_ruv_style"}},
        consumed_by_downstream=consumed_by_downstream,
    )


def _resolved_correction_matrix(matrix: pd.DataFrame) -> pd.DataFrame:
    resolved = DatasetBuildRequestInterpreter().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata(_phospho()),
            sample_metadata=_sample_metadata(_phospho()),
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )
    corrected = matrix.copy(deep=True)
    corrected.index = resolved.phospho.index.copy()
    return corrected


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_1": [10.0, 2.0],
            "sample_2": [14.0, 1.0],
            "sample_3": [15.0, 7.0],
            "sample_4": [19.0, 6.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
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


def _sample_metadata(phospho: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": ["run_1", "run_1", "run_2", "run_2"],
            "condition": ["control", "treated", "control", "treated"],
        },
        index=phospho.columns.copy(),
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
