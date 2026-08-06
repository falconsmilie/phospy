from __future__ import annotations

import pandas as pd
import pytest

from phospy.api import (
    AnalysisReadyDatasetBuilder,
    DatasetBatchCorrectionConfig,
    DatasetBuildRequest,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    Organism,
    SpsRuvBatchCorrectionConfig,
)
from phospy.errors import PhosPyInputError
from phospy.science.configs.preprocessing.correction_missingness import (
    CorrectionMissingnessPolicy,
)
from phospy.science.datasets.builders.contracts import InterpretedDatasetBuildRequest
from phospy.science.datasets.builders.executor import DatasetBuildExecutor
from phospy.science.datasets.preprocessing.batch_correction import (
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
)
from phospy.science.datasets.preprocessing.control_sites import (
    ControlSiteSet,
    ControlSiteSourceMetadata,
)
from phospy.science.datasets.preprocessing.correction_output import (
    CorrectedPreprocessingOutput,
)
from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from phospy.science.transformations.models import IntensityScaleKind
from tests.support.dataset_preprocessor_fakes import (
    ConformingDatasetPreprocessorFake,
)


def test_public_builder_rejects_linear_median_centering() -> None:
    phospho = _sps_phospho()

    with pytest.raises(PhosPyInputError) as exc_info:
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_sps_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                preprocessing_config=DatasetPreprocessingConfig(
                    normalisation=DatasetNormalisationConfig(policy="median_center"),
                ),
            )
        )

    _assert_additive_scale_error(
        str(exc_info.value),
        operation="median_center",
    )


def test_public_builder_rejects_linear_fixed_effect_residualisation() -> None:
    phospho = _sps_phospho()

    with pytest.raises(PhosPyInputError) as exc_info:
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_sps_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                preprocessing_config=DatasetPreprocessingConfig(
                    batch_correction=DatasetBatchCorrectionConfig(
                        method="linear_residualize_batch",
                    ),
                ),
            )
        )

    _assert_additive_scale_error(
        str(exc_info.value),
        operation="linear_residualize_batch",
    )


def test_public_builder_rejects_linear_sps_ruv_residualisation() -> None:
    phospho = _sps_phospho()

    with pytest.raises(PhosPyInputError) as exc_info:
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_sps_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                preprocessing_config=DatasetPreprocessingConfig(
                    batch_correction=_sps_ruv_config(),
                ),
            )
        )

    _assert_additive_scale_error(
        str(exc_info.value),
        operation="sps_ruv_style",
    )


def test_public_builder_keeps_valid_log2_median_centering_path() -> None:
    phospho = _sps_phospho()

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_sps_site_metadata(phospho),
            sample_metadata=_sample_metadata(phospho),
            organism=Organism.RAT,
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                normalisation=DatasetNormalisationConfig(policy="median_center"),
            ),
        )
    )

    assert built.processing_state.normalisation.policy == "median_center"
    assert built.intensity_scale_state.kind is IntensityScaleKind.LOG2
    assert built.preprocessing_report is not None
    assert "normalisation" in set(
        built.preprocessing_report.operations.loc[:, "stage"].astype(str)
    )


def test_public_builder_keeps_valid_log2_fixed_effect_residualisation_path() -> None:
    phospho = _sps_phospho()

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_sps_site_metadata(phospho),
            sample_metadata=_sample_metadata(phospho),
            organism=Organism.RAT,
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                batch_correction=DatasetBatchCorrectionConfig(
                    method="linear_residualize_batch",
                ),
            ),
        )
    )

    assert built.intensity_scale_state.kind is IntensityScaleKind.LOG2
    assert built.preprocessing_report is not None
    report = built.preprocessing_report.batch_correction
    assert report is not None
    assert report.status == "applied"
    assert report.method == "linear_residualize_batch"


def test_public_builder_keeps_valid_log2_sps_ruv_residualisation_path() -> None:
    phospho = _sps_phospho()

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_sps_site_metadata(phospho),
            sample_metadata=_sample_metadata(phospho),
            organism=Organism.RAT,
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                batch_correction=_sps_ruv_config(),
            ),
        )
    )

    assert built.intensity_scale_state.kind is IntensityScaleKind.LOG2
    assert built.preprocessing_report is not None
    report = built.preprocessing_report.batch_correction
    assert report is not None
    assert report.status == "applied"
    assert report.method == "sps_ruv_style"


@pytest.mark.parametrize(
    ("operation_key", "operation"),
    (
        ("median_center", "median_center"),
        ("linear_residualize_batch", "linear_residualize_batch"),
        ("sps_ruv_style", "sps_ruv_style"),
    ),
)
def test_rejected_linear_additive_preprocessing_does_not_run_preprocessor(
    operation_key: str,
    operation: str,
) -> None:
    phospho = _sps_phospho()
    plan = PreprocessingPlan.from_config(
        _preprocessing_config_for_operation(operation_key)
    )
    spy = ConformingDatasetPreprocessorFake()

    with pytest.raises(PhosPyInputError, match=f"operation='{operation}'"):
        DatasetBuildExecutor(preprocessor=spy).run(
            InterpretedDatasetBuildRequest(
                phospho=phospho,
                site_metadata=_sps_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                total=None,
                organism=Organism.RAT,
                preprocessing_plan=plan,
                declared_input_intensity_scale_kind=IntensityScaleKind.LINEAR,
                declared_input_intensity_scale_source=(
                    "dataset_build_request.input_intensity_scale"
                ),
            )
        )

    assert spy.preflight_calls == []
    assert spy.run_calls == []


def test_rejected_linear_external_corrected_output_does_not_run_preprocessor() -> None:
    phospho = _sps_phospho()
    spy = ConformingDatasetPreprocessorFake()

    with pytest.raises(PhosPyInputError, match="operation='sps_ruv_style'"):
        DatasetBuildExecutor(preprocessor=spy).run(
            InterpretedDatasetBuildRequest(
                phospho=phospho,
                site_metadata=_sps_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                total=None,
                organism=Organism.RAT,
                preprocessing_plan=PreprocessingPlan.from_config(
                    DatasetPreprocessingConfig()
                ),
                corrected_preprocessing_output=_external_sps_ruv_output(phospho),
                declared_input_intensity_scale_kind=IntensityScaleKind.LINEAR,
                declared_input_intensity_scale_source=(
                    "dataset_build_request.input_intensity_scale"
                ),
            )
        )

    assert spy.preflight_calls == []
    assert spy.run_calls == []


def _assert_additive_scale_error(message: str, *, operation: str) -> None:
    assert f"operation='{operation}'" in message
    assert "current_scale='linear'" in message
    assert "required_scale='log2'" in message
    assert "established log2 phosphosite abundance" in message
    assert "preprocessing_config.intensity_transform.policy='log2'" in message
    assert "input_intensity_scale='log2'" in message
    assert "multiplicative median-scaling" in message


def _sps_phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_1": [10.0, 5.0, 20.0],
            "sample_2": [10.0, 9.0, 20.0],
            "sample_3": [14.0, 8.0, 28.0],
            "sample_4": [14.0, 12.0, 28.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;", "SRC;Y416;"], name="site_id"),
    )


def _sps_site_metadata(phospho: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "SRC"],
            "protein_id": ["MAPK14", "AKT1", "SRC"],
            "site": ["Y182", "T308", "Y416"],
            "site_sequence": [
                ("A" * 15) + "Y" + ("A" * 15),
                ("A" * 15) + "T" + ("A" * 15),
                ("A" * 15) + "Y" + ("A" * 15),
            ],
            "localisation_confidence": [0.95, 0.92, 0.98],
        },
        index=phospho.index.copy(),
    )


def _sample_metadata(phospho: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": ("run_1", "run_1", "run_2", "run_2"),
            "condition": ("control", "treated", "control", "treated"),
            "replicate": ("r1", "r2", "r2", "r1"),
        },
        index=phospho.columns.copy(),
    )


def _sps_ruv_config() -> SpsRuvBatchCorrectionConfig:
    return SpsRuvBatchCorrectionConfig(
        control_site_set=ControlSiteSet.from_site_keys(
            (
                _sps_site_key("MAPK14", "Y", "182"),
                _sps_site_key("SRC", "Y", "416"),
            ),
            source_metadata=ControlSiteSourceMetadata(
                organism="rat",
                identifier_namespace="site_key",
                source_name="manual-curated-controls",
                source_version="manual-v1",
                license="caller local use",
                redistribution="not redistributed",
            ),
        ),
        batch_column="batch",
        condition_columns=("condition",),
        replicate_column="replicate",
        missingness_policy=CorrectionMissingnessPolicy(),
        n_unwanted_factors=1,
        diagnostics_enabled=True,
        provenance_enabled=True,
    )


def _preprocessing_config_for_operation(
    operation_key: str,
) -> DatasetPreprocessingConfig:
    if operation_key == "median_center":
        return DatasetPreprocessingConfig(
            normalisation=DatasetNormalisationConfig(policy="median_center"),
        )
    if operation_key == "linear_residualize_batch":
        return DatasetPreprocessingConfig(
            batch_correction=DatasetBatchCorrectionConfig(
                method="linear_residualize_batch",
            ),
        )
    if operation_key == "sps_ruv_style":
        return DatasetPreprocessingConfig(batch_correction=_sps_ruv_config())
    raise AssertionError(f"unsupported operation key: {operation_key}")


def _external_sps_ruv_output(phospho: pd.DataFrame) -> CorrectedPreprocessingOutput:
    return CorrectedPreprocessingOutput(
        corrected_matrix=phospho + 1.0,
        output_observation_mask=pd.DataFrame(
            True,
            index=phospho.index.copy(),
            columns=phospho.columns.copy(),
        ),
        corrected_cell_status=pd.DataFrame(
            "corrected_observed",
            index=phospho.index.copy(),
            columns=phospho.columns.copy(),
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
                matrix_shape_before=(3, 4),
                matrix_shape_after=(3, 4),
            ),
        ),
        diagnostics={"source": "unit-test"},
    )


def _sps_site_key(protein_identifier: str, residue: str, position: str) -> str:
    return (
        "phospy:v1|organism=rat|protein_namespace=protein_id|"
        f"protein_identifier={protein_identifier}|residue={residue}|"
        f"position={position}"
    )
