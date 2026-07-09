from __future__ import annotations

from collections.abc import Mapping
from dataclasses import replace
from typing import Any, cast

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import DifferentialAnalysisWorkflow
from phospy.api import (
    AnalysisReadyDatasetBuilder,
    Contrast,
    ControlSiteSet,
    CorrectionMissingnessPolicy,
    DatasetBatchCorrectionConfig,
    DatasetBuildRequest,
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DatasetTotalProteinCorrectionConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
    SpsRuvBatchCorrectionConfig,
)
from phospy.contracts.configs.preprocessing import InternalBatchCorrectionStageOrder
from phospy.errors import PhosPyInputError
from phospy.provenance import BatchCorrectionProvenance, fingerprint_matrix
from phospy.science.datasets.builders.interpreter import DatasetBuildRequestInterpreter
from phospy.science.datasets.preprocessing import (
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
    CorrectedPreprocessingOutput,
)
from phospy.science.datasets.preprocessing.control_sites import (
    ControlSiteSourceMetadata,
)
from phospy.validation.datasets.batch_correction import (
    validate_applied_native_sps_ruv_correction_provenance,
)

pytestmark = pytest.mark.integration

_DEFAULT_PROVENANCE: Any = object()
_COMPLETE_EXTERNAL_DEPENDENCY_VERSIONS = {
    "numpy": "test-numpy",
    "pandas": "test-pandas",
    "scipy": "test-scipy",
    "scikit-learn": "test-scikit-learn",
}


def test_resolved_native_correction_output_builds_analysis_ready_dataset() -> None:
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    resolved_input = _resolved_correction_matrix(phospho)

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
    correction_stage = next(
        stage
        for stage in dataset.provenance.preprocessing_stages
        if stage.stage == "batch_correction"
    )
    provenance = correction_stage.batch_correction_provenance
    assert provenance is not None
    assert provenance.requested_method == "sps_ruv_style"
    assert provenance.control_site_source["source_type"] == "caller_supplied"
    assert provenance.selected_site_key_rows
    assert provenance.observation_masks
    assert provenance.input_matrix_fingerprint.name == "batch_correction.native.input"
    assert provenance.input_matrix_fingerprint == fingerprint_matrix(
        resolved_input,
        name="batch_correction.native.input",
    )
    assert provenance.output_matrix_fingerprint is not None
    assert (
        provenance.output_matrix_fingerprint.name == "batch_correction.native.corrected"
    )
    assert provenance.output_matrix_fingerprint == fingerprint_matrix(
        corrected,
        name="batch_correction.native.corrected",
    )
    executor_diagnostics = cast(
        Mapping[str, object], provenance.diagnostics["executor"]
    )
    assert executor_diagnostics["status"] == "applied"
    assert provenance.warnings == ()
    observed_mask = dataset.imputation_observed_mask_dataframe()
    assert observed_mask is not None
    assert observed_mask.to_numpy(dtype=bool).all()
    assert dataset.phospho.notna().to_numpy(dtype=bool).all()


def test_external_corrected_output_allowed_without_downstream_matrix_consuming_stages() -> (
    None
):
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


def test_external_corrected_output_rejected_with_normalization_stage() -> None:
    _assert_external_corrected_output_rejected_with_preprocessing_config(
        DatasetPreprocessingConfig(
            normalisation=DatasetNormalisationConfig(policy="median_center")
        ),
        "normalisation",
    )


def test_external_corrected_output_rejected_with_total_protein_correction_stage() -> (
    None
):
    _assert_external_corrected_output_rejected_with_preprocessing_config(
        DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(policy="log2"),
            total_protein_correction=DatasetTotalProteinCorrectionConfig(
                policy="subtract_log_total"
            ),
        ),
        "total_protein_correction",
        input_intensity_scale="linear",
        total=_total(_phospho()),
    )


def test_external_corrected_output_rejected_with_site_matrix_construction_stage() -> (
    None
):
    _assert_external_corrected_output_rejected_with_preprocessing_config(
        DatasetPreprocessingConfig(
            site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
        ),
        "site_matrix",
    )


def test_external_corrected_output_rejected_with_comparisons_stage() -> None:
    _assert_external_corrected_output_rejected_with_preprocessing_config(
        DatasetPreprocessingConfig(
            comparisons=DatasetComparisonBuildingConfig(
                policy="sample_metadata_pairs",
                sample_group_column="condition",
            )
        ),
        "comparisons",
    )


def test_sps_ruv_corrected_output_without_provenance_is_rejected() -> None:
    phospho = _phospho()

    with pytest.raises(
        PhosPyInputError,
        match="typed BatchCorrectionProvenance",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    _resolved_correction_matrix(phospho + 1.0),
                    provenance=None,
                ),
            )
        )


@pytest.mark.parametrize("status", ("disabled", "rejected"))
def test_external_corrected_output_requires_applied_status(status: str) -> None:
    phospho = _phospho()

    with pytest.raises(
        PhosPyInputError,
        match="applied correction status.*method-specific provenance",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    _resolved_correction_matrix(phospho + 1.0),
                    status=status,
                ),
            )
        )


def test_external_corrected_output_rejects_none_method() -> None:
    phospho = _phospho()

    with pytest.raises(
        PhosPyInputError,
        match="method='none'.*external corrected matrix",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    _resolved_correction_matrix(phospho + 1.0),
                    method="none",
                ),
            )
        )


def test_external_corrected_output_rejects_unsupported_method_label() -> None:
    phospho = _phospho()

    with pytest.raises(
        PhosPyInputError,
        match="supported applied correction method.*method-specific provenance",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    _resolved_correction_matrix(phospho + 1.0),
                    method="linear_residualize_batch",
                ),
            )
        )


def test_sps_ruv_corrected_output_with_empty_dependency_versions_is_rejected() -> None:
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    provenance = replace(
        _complete_sps_ruv_provenance(corrected),
        dependency_versions={},
    )

    with pytest.raises(
        PhosPyInputError,
        match="dependency_versions",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    corrected,
                    provenance=provenance,
                ),
            )
        )


@pytest.mark.parametrize("phospy_version", ("", "unknown"))
def test_sps_ruv_corrected_output_with_missing_phospy_version_is_rejected(
    phospy_version: str,
) -> None:
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    provenance = replace(
        _complete_sps_ruv_provenance(corrected),
        phospy_version=phospy_version,
    )

    with pytest.raises(
        PhosPyInputError,
        match="phospy_version",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    corrected,
                    provenance=provenance,
                ),
            )
        )


def test_sps_ruv_corrected_output_with_missing_python_version_is_rejected() -> None:
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    provenance = replace(
        _complete_sps_ruv_provenance(corrected),
        python_version="",
    )

    with pytest.raises(
        PhosPyInputError,
        match="python_version",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    corrected,
                    provenance=provenance,
                ),
            )
        )


def test_sps_ruv_corrected_output_accepts_complete_external_environment_provenance() -> (
    None
):
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    provenance = _complete_sps_ruv_provenance(corrected)

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho),
            sample_metadata=_sample_metadata(phospho),
            organism=Organism.RAT,
            input_intensity_scale="linear",
            corrected_preprocessing_output=_correction_output(
                corrected,
                provenance=provenance,
            ),
        )
    )

    attached = _attached_batch_correction_provenance(dataset)
    assert attached.phospy_version == "test"
    assert attached.python_version == "3.test"
    assert attached.dependency_versions == _COMPLETE_EXTERNAL_DEPENDENCY_VERSIONS


def test_external_corrected_output_accepts_applied_sps_ruv_style_with_provenance() -> (
    None
):
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho),
            sample_metadata=_sample_metadata(phospho),
            organism=Organism.RAT,
            input_intensity_scale="linear",
            corrected_preprocessing_output=_correction_output(
                corrected,
                method="sps_ruv_style",
                provenance=_complete_sps_ruv_provenance(corrected),
            ),
        )
    )

    pdt.assert_frame_equal(dataset.phospho, corrected)
    assert _attached_batch_correction_provenance(dataset).requested_method == (
        "sps_ruv_style"
    )


def test_sps_ruv_corrected_output_with_not_provided_controls_is_rejected() -> None:
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    provenance = replace(
        _complete_sps_ruv_provenance(corrected),
        selected_site_key_rows=(),
        control_site_source={
            "source_type": "not_provided",
            "reason": "synthesized fallback",
        },
    )

    with pytest.raises(
        PhosPyInputError,
        match="selected controls|control provenance",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    corrected,
                    provenance=provenance,
                ),
            )
        )


def test_sps_ruv_corrected_output_with_one_control_for_one_factor_is_rejected() -> None:
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    provenance = replace(
        _complete_sps_ruv_provenance(corrected),
        selected_site_key_rows=(str(corrected.index[0]),),
    )

    with pytest.raises(
        PhosPyInputError,
        match="selected controls.*unwanted-factor count",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    corrected,
                    provenance=provenance,
                ),
            )
        )


def test_sps_ruv_corrected_output_accepts_two_controls_for_one_factor() -> None:
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    provenance = replace(
        _complete_sps_ruv_provenance(corrected),
        selected_site_key_rows=tuple(str(row) for row in corrected.index[:2]),
    )

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho),
            sample_metadata=_sample_metadata(phospho),
            organism=Organism.RAT,
            input_intensity_scale="linear",
            corrected_preprocessing_output=_correction_output(
                corrected,
                provenance=provenance,
            ),
        )
    )

    assert _attached_batch_correction_provenance(dataset).selected_site_key_rows == (
        str(corrected.index[0]),
        str(corrected.index[1]),
    )


def test_corrected_output_integrator_rejects_selected_controls_absent_from_matrix_index() -> (
    None
):
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    provenance = replace(
        _complete_sps_ruv_provenance(corrected),
        selected_site_key_rows=(str(corrected.index[0]), "absent_site_key_row"),
    )

    with pytest.raises(
        PhosPyInputError,
        match="selected_site_key_rows.*corrected matrix index",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    corrected,
                    provenance=provenance,
                ),
            )
        )


@pytest.mark.parametrize(
    "resolved_parameters",
    (
        {
            "method": "sps_ruv_style",
            "source": "external_corrected_preprocessing_output",
        },
        {"method": "sps_ruv_style", "n_unwanted_factors": None},
        {"method": "sps_ruv_style", "n_unwanted_factors": True},
        {"method": "sps_ruv_style", "n_unwanted_factors": "1"},
        {"method": "sps_ruv_style", "n_unwanted_factors": 1.0},
        {"method": "sps_ruv_style", "n_unwanted_factors": 0},
        {"method": "sps_ruv_style", "n_unwanted_factors": -1},
        {"config": {"method": "sps_ruv_style", "n_unwanted_factors": 0}},
    ),
)
def test_sps_ruv_corrected_output_rejects_missing_or_invalid_unwanted_factors(
    resolved_parameters: Mapping[str, object],
) -> None:
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    provenance = replace(
        _complete_sps_ruv_provenance(corrected),
        resolved_parameters=resolved_parameters,
    )

    with pytest.raises(
        PhosPyInputError,
        match="selected controls.*unwanted-factor count|n_unwanted_factors",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    corrected,
                    provenance=provenance,
                ),
            )
        )


def test_sps_ruv_corrected_output_missing_fingerprints_is_rejected() -> None:
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    provenance = replace(
        _complete_sps_ruv_provenance(corrected),
        output_matrix_fingerprint=None,
    )

    with pytest.raises(
        PhosPyInputError,
        match="input/output matrix fingerprints",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    corrected,
                    provenance=provenance,
                ),
            )
        )


def test_sps_ruv_corrected_output_missing_observation_mask_is_rejected() -> None:
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    provenance = replace(
        _complete_sps_ruv_provenance(corrected),
        observation_masks=(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="observation mask|missingness provenance",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    corrected,
                    provenance=provenance,
                ),
            )
        )


def test_sps_ruv_corrected_output_input_fingerprint_mismatch_is_rejected() -> None:
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    provenance = replace(
        _complete_sps_ruv_provenance(
            corrected,
            input_matrix=_resolved_correction_matrix(phospho),
        ),
        input_matrix_fingerprint=fingerprint_matrix(
            corrected,
            name="batch_correction.native.input",
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="input_matrix_fingerprint.*pre-correction dataset.phospho",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    corrected,
                    input_matrix=_resolved_correction_matrix(phospho),
                    provenance=provenance,
                ),
            )
        )


def test_sps_ruv_corrected_output_output_fingerprint_mismatch_is_rejected() -> None:
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    provenance = replace(
        _complete_sps_ruv_provenance(
            corrected,
            input_matrix=_resolved_correction_matrix(phospho),
        ),
        output_matrix_fingerprint=fingerprint_matrix(
            corrected + 10.0,
            name="batch_correction.native.corrected",
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="output_matrix_fingerprint.*corrected_matrix",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    corrected,
                    input_matrix=_resolved_correction_matrix(phospho),
                    provenance=provenance,
                ),
            )
        )


def test_sps_ruv_corrected_output_observation_mask_fingerprint_mismatch_is_rejected() -> (
    None
):
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    mismatched_mask = pd.DataFrame(
        False,
        index=corrected.index.copy(),
        columns=corrected.columns.copy(),
    )
    provenance = replace(
        _complete_sps_ruv_provenance(
            corrected,
            input_matrix=_resolved_correction_matrix(phospho),
        ),
        observation_masks=(
            fingerprint_matrix(
                mismatched_mask.astype("int8"),
                name="batch_correction.native.observation_mask",
            ),
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="observation_masks.*output_observation_mask",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    corrected,
                    input_matrix=_resolved_correction_matrix(phospho),
                    provenance=provenance,
                ),
            )
        )


def test_sps_ruv_corrected_output_missing_control_metadata_requires_rationale() -> None:
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    provenance = replace(
        _complete_sps_ruv_provenance(corrected),
        control_site_source={
            "source_type": "caller_supplied",
            "source_version_unavailable_reason": "local controls",
        },
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "control source.*organism|control source.*identifier_namespace|"
            "control source.*license|control source.*redistribution"
        ),
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    corrected,
                    provenance=provenance,
                ),
            )
        )


@pytest.mark.parametrize("field_name", ("license", "redistribution"))
def test_sps_ruv_corrected_output_missing_control_usage_metadata_is_rejected(
    field_name: str,
) -> None:
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    control_site_source = dict(
        _complete_sps_ruv_provenance(corrected).control_site_source
    )
    control_site_source.pop(field_name)
    provenance = replace(
        _complete_sps_ruv_provenance(corrected),
        control_site_source=control_site_source,
    )

    with pytest.raises(
        PhosPyInputError,
        match=f"control source.*{field_name}.*explicit rationale",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                corrected_preprocessing_output=_correction_output(
                    corrected,
                    provenance=provenance,
                ),
            )
        )


def test_sps_ruv_corrected_output_accepts_missing_control_usage_metadata_rationale() -> (
    None
):
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    control_site_source = dict(
        _complete_sps_ruv_provenance(corrected).control_site_source
    )
    control_site_source.pop("license")
    control_site_source.pop("redistribution")
    control_site_source["license_missing_reason"] = (
        "caller-local controls are not licensed data"
    )
    control_site_source["redistribution_missing_reason"] = (
        "caller-local controls are not redistributed"
    )
    provenance = replace(
        _complete_sps_ruv_provenance(corrected),
        control_site_source=control_site_source,
    )

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho),
            sample_metadata=_sample_metadata(phospho),
            organism=Organism.RAT,
            input_intensity_scale="linear",
            corrected_preprocessing_output=_correction_output(
                corrected,
                provenance=provenance,
            ),
        )
    )

    attached = _attached_batch_correction_provenance(dataset)
    assert attached.control_site_source["license_missing_reason"]
    assert attached.control_site_source["redistribution_missing_reason"]


def test_sps_ruv_corrected_output_accepts_control_metadata_missing_rationale() -> None:
    phospho = _phospho()
    corrected = _resolved_correction_matrix(phospho + 1.0)
    provenance = replace(
        _complete_sps_ruv_provenance(corrected),
        control_site_source={
            "source_type": "caller_supplied",
            "metadata_missing_reason": {
                "organism": "caller-local controls did not declare organism",
                "identifier_namespace": (
                    "caller-local controls did not declare namespace"
                ),
                "source_version": "caller-local controls have no source version",
                "license": "caller-local controls are not licensed data",
                "redistribution": "caller-local controls are not redistributed",
            },
            "source_version_unavailable_reason": (
                "caller-local controls have no source version"
            ),
        },
    )

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho),
            sample_metadata=_sample_metadata(phospho),
            organism=Organism.RAT,
            input_intensity_scale="linear",
            corrected_preprocessing_output=_correction_output(
                corrected,
                provenance=provenance,
            ),
        )
    )

    assert _attached_batch_correction_provenance(dataset).control_site_source[
        "metadata_missing_reason"
    ]


def test_public_sps_ruv_preprocessing_config_builds_corrected_dataset_with_provenance() -> (
    None
):
    phospho = _sps_phospho()

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_sps_site_metadata(phospho),
            sample_metadata=_sps_sample_metadata(phospho),
            organism=Organism.RAT,
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                batch_correction=SpsRuvBatchCorrectionConfig(
                    control_site_set=ControlSiteSet.from_site_keys(
                        (
                            _sps_site_key("MAPK14", "Y", "182"),
                            _sps_site_key("SRC", "Y", "416"),
                        ),
                        source_metadata=_control_source_metadata(),
                    ),
                    batch_column="batch",
                    condition_columns=("condition",),
                    replicate_column="replicate",
                    missingness_policy=CorrectionMissingnessPolicy(),
                    n_unwanted_factors=1,
                    diagnostics_enabled=True,
                    provenance_enabled=True,
                )
            ),
        )
    )

    assert dataset.preprocessing_report is not None
    report = dataset.preprocessing_report.batch_correction
    assert report is not None
    assert report.status == "applied"
    assert report.method == "sps_ruv_style"
    assert report.batch_column == "batch"
    assert report.condition_column == "condition"
    assert dataset.provenance is not None
    correction_stage = next(
        stage
        for stage in dataset.provenance.preprocessing_stages
        if stage.stage == "batch_correction"
    )
    provenance = correction_stage.batch_correction_provenance
    assert provenance is not None
    assert provenance.requested_method == "sps_ruv_style"
    assert provenance.phospy_version
    assert provenance.phospy_version != "unknown"
    assert provenance.python_version
    assert provenance.python_version != "unknown"
    assert {"numpy", "pandas", "scipy", "scikit-learn"}.issubset(
        set(provenance.dependency_versions)
    )
    assert provenance.selected_site_key_rows
    stage_order = tuple(
        stage.stage for stage in dataset.provenance.preprocessing_stages
    )
    assert stage_order.index("missing_data") < stage_order.index("batch_correction")
    assert correction_stage.parameters["executed_stage_order"] == list(
        provenance.preprocessing_stage_order
    )
    assert provenance.preprocessing_stage_order == (
        "missing_data",
        "batch_correction",
        "downstream_workflows",
    )
    interpreter_plan = cast(
        Mapping[str, object],
        provenance.resolved_parameters["interpreter_plan"],
    )
    assert interpreter_plan["executed_stage_order"] == list(
        provenance.preprocessing_stage_order
    )
    executor_diagnostics = provenance.diagnostics["executor"]
    assert isinstance(executor_diagnostics, Mapping)
    assert executor_diagnostics["status"] == "applied"
    assert provenance.output_matrix_fingerprint is not None
    assert not dataset.phospho.equals(phospho)


def test_workflow_generated_selected_control_provenance_still_passes_validation() -> (
    None
):
    phospho = _sps_phospho()

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_sps_site_metadata(phospho),
            sample_metadata=_sps_sample_metadata(phospho),
            organism=Organism.RAT,
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                batch_correction=SpsRuvBatchCorrectionConfig(
                    control_site_set=ControlSiteSet.from_site_keys(
                        (
                            _sps_site_key("MAPK14", "Y", "182"),
                            _sps_site_key("SRC", "Y", "416"),
                        ),
                        source_metadata=_control_source_metadata(),
                    ),
                    batch_column="batch",
                    condition_columns=("condition",),
                    replicate_column="replicate",
                    missingness_policy=CorrectionMissingnessPolicy(),
                    n_unwanted_factors=1,
                    diagnostics_enabled=True,
                    provenance_enabled=True,
                )
            ),
        )
    )

    assert dataset.preprocessing_report is not None
    report = dataset.preprocessing_report.batch_correction
    assert report is not None
    provenance = _attached_batch_correction_provenance(dataset)

    validate_applied_native_sps_ruv_correction_provenance(
        method=report.method,
        status=report.status,
        provenance=provenance,
    )
    assert len(set(provenance.selected_site_key_rows)) == len(
        provenance.selected_site_key_rows
    )
    assert set(provenance.selected_site_key_rows).issubset(
        set(dataset.phospho.index.astype(str))
    )


def test_internal_sps_ruv_config_still_runs_at_batch_correction_stage() -> None:
    phospho = _sps_phospho()

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_sps_site_metadata(phospho),
            sample_metadata=_sps_sample_metadata(phospho),
            organism=Organism.RAT,
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                batch_correction=SpsRuvBatchCorrectionConfig(
                    control_site_set=ControlSiteSet.from_site_keys(
                        (
                            _sps_site_key("MAPK14", "Y", "182"),
                            _sps_site_key("SRC", "Y", "416"),
                        ),
                        source_metadata=_control_source_metadata(),
                    ),
                    batch_column="batch",
                    condition_columns=("condition",),
                    replicate_column="replicate",
                    missingness_policy=CorrectionMissingnessPolicy(),
                    n_unwanted_factors=1,
                    diagnostics_enabled=True,
                    provenance_enabled=True,
                ),
                normalisation=DatasetNormalisationConfig(policy="median_center"),
            ),
        )
    )

    assert dataset.provenance is not None
    stage_order = tuple(
        stage.stage for stage in dataset.provenance.preprocessing_stages
    )
    assert stage_order.index("batch_correction") < stage_order.index("normalisation")
    correction_stage = next(
        stage
        for stage in dataset.provenance.preprocessing_stages
        if stage.stage == "batch_correction"
    )
    assert correction_stage.batch_correction_provenance is not None
    assert correction_stage.batch_correction_provenance.preprocessing_stage_order == (
        "missing_data",
        "batch_correction",
        "downstream_workflows",
    )


def test_public_sps_ruv_preprocessing_preserves_multiple_condition_columns_in_reports() -> (
    None
):
    phospho = _sps_multi_condition_phospho()
    condition_columns = ("condition", "timepoint")
    expected_joint_levels = [
        "condition=control|timepoint=early",
        "condition=treated|timepoint=early",
        "condition=control|timepoint=late",
        "condition=treated|timepoint=late",
    ]
    expected_condition_terms = [
        "intercept",
        "condition[condition=treated|timepoint=early]",
        "condition[condition=control|timepoint=late]",
        "condition[condition=treated|timepoint=late]",
    ]

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_sps_site_metadata(phospho),
            sample_metadata=_sps_multi_condition_sample_metadata(phospho),
            organism=Organism.RAT,
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                batch_correction=SpsRuvBatchCorrectionConfig(
                    control_site_set=ControlSiteSet.from_site_keys(
                        (
                            _sps_site_key("MAPK14", "Y", "182"),
                            _sps_site_key("SRC", "Y", "416"),
                        ),
                        source_metadata=_control_source_metadata(),
                    ),
                    batch_column="batch",
                    condition_columns=condition_columns,
                    replicate_column="replicate",
                    missingness_policy=CorrectionMissingnessPolicy(),
                    n_unwanted_factors=1,
                    diagnostics_enabled=True,
                    provenance_enabled=True,
                )
            ),
        )
    )

    assert dataset.preprocessing_report is not None
    report = dataset.preprocessing_report.batch_correction
    assert report is not None
    assert report.condition_column == "condition"
    assert report.condition_columns == condition_columns
    report_payload = report.to_payload()
    assert report_payload["condition_columns"] == list(condition_columns)
    policy_payload = report_payload["policy"]
    assert isinstance(policy_payload, Mapping)
    assert policy_payload["condition_columns"] == list(condition_columns)

    assert dataset.provenance is not None
    correction_stage = next(
        stage
        for stage in dataset.provenance.preprocessing_stages
        if stage.stage == "batch_correction"
    )
    assert correction_stage.parameters["condition_columns"] == list(condition_columns)
    stage_diagnostics = cast(Mapping[str, object], correction_stage.diagnostics)
    stage_executor_diagnostics = cast(
        Mapping[str, object], stage_diagnostics["executor"]
    )
    stage_design_summary = cast(
        Mapping[str, object], stage_executor_diagnostics["design_summary"]
    )
    assert stage_design_summary["condition_columns"] == list(condition_columns)
    assert stage_design_summary["condition_levels"] == expected_joint_levels
    assert (
        stage_design_summary["condition_terms_to_preserve"] == expected_condition_terms
    )
    assert stage_design_summary["samples_per_condition"] == {
        level: 2 for level in expected_joint_levels
    }

    provenance = correction_stage.batch_correction_provenance
    assert provenance is not None
    assert provenance.design_metadata["condition_columns"] == list(condition_columns)
    provenance_condition_by_sample = cast(
        Mapping[str, object],
        provenance.batch_metadata["condition_by_sample"],
    )
    assert (
        provenance_condition_by_sample["sample_2"]
        == "condition=treated|timepoint=early"
    )
    assert (
        provenance.design_metadata["condition_terms_to_preserve"]
        == expected_condition_terms
    )
    provenance_config = cast(
        Mapping[str, object], provenance.resolved_parameters["config"]
    )
    assert provenance_config["condition_columns"] == list(condition_columns)
    interpreter_plan = cast(
        Mapping[str, object],
        provenance.resolved_parameters["interpreter_plan"],
    )
    assert interpreter_plan["condition_terms_to_preserve"] == expected_condition_terms
    design_payload = cast(
        Mapping[str, object], interpreter_plan["resolved_design_matrix"]
    )
    assert design_payload["columns"] == [
        *expected_condition_terms,
        "batch[run_2]",
    ]
    provenance_executor_diagnostics = cast(
        Mapping[str, object], provenance.diagnostics["executor"]
    )
    provenance_design_summary = cast(
        Mapping[str, object], provenance_executor_diagnostics["design_summary"]
    )
    assert provenance_design_summary["condition_columns"] == list(condition_columns)
    assert provenance_design_summary["condition_levels"] == expected_joint_levels


def test_public_sps_ruv_preprocessing_invalid_controls_fail_before_execution() -> None:
    phospho = _sps_phospho()

    with pytest.raises(
        PhosPyInputError,
        match="batch correction validation failed before correction execution",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_sps_site_metadata(phospho),
                sample_metadata=_sps_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="log2",
                preprocessing_config=DatasetPreprocessingConfig(
                    batch_correction=SpsRuvBatchCorrectionConfig(
                        control_site_set=ControlSiteSet.from_site_keys(
                            (_sps_site_key("MAPK14", "Y", "182"),),
                            source_metadata=_control_source_metadata(),
                        ),
                        batch_column="batch",
                        condition_columns=("condition",),
                        replicate_column="replicate",
                        missingness_policy=CorrectionMissingnessPolicy(),
                        n_unwanted_factors=1,
                    )
                ),
            )
        )


def test_public_sps_ruv_preprocessing_rejects_unsupported_stage_order() -> None:
    phospho = _sps_phospho()

    with pytest.raises(
        PhosPyInputError,
        match=(
            "after_intensity_transform_before_missing_data.*unsupported.*"
            "supported stage order is missing_data -> batch_correction -> "
            "downstream_workflows.*provenance must match"
        ),
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_sps_site_metadata(phospho),
                sample_metadata=_sps_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="log2",
                preprocessing_config=DatasetPreprocessingConfig(
                    batch_correction=SpsRuvBatchCorrectionConfig(
                        control_site_set=ControlSiteSet.from_site_keys(
                            (
                                _sps_site_key("MAPK14", "Y", "182"),
                                _sps_site_key("SRC", "Y", "416"),
                            ),
                            source_metadata=_control_source_metadata(),
                        ),
                        batch_column="batch",
                        condition_columns=("condition",),
                        replicate_column="replicate",
                        missingness_policy=CorrectionMissingnessPolicy(),
                        n_unwanted_factors=1,
                        stage_order=(
                            InternalBatchCorrectionStageOrder.AFTER_INTENSITY_TRANSFORM_BEFORE_MISSING_DATA
                        ),
                    )
                ),
            )
        )


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
            "sample_1": [9.9, 1.9],
            "sample_2": [14.1, 1.1],
            "sample_3": [10.1, 2.1],
            "sample_4": [13.9, 0.9],
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
    assert float(cast(float, table.loc[site_one, "logFC"])) == pytest.approx(4.0)
    assert float(cast(float, table.loc[site_two, "logFC"])) == pytest.approx(-1.0)


def _assert_external_corrected_output_rejected_with_preprocessing_config(
    preprocessing_config: DatasetPreprocessingConfig,
    expected_stage: str,
    *,
    input_intensity_scale: str = "linear",
    total: pd.DataFrame | None = None,
) -> None:
    phospho = _phospho()

    with pytest.raises(
        PhosPyInputError,
        match=(
            "external corrected output cannot be integrated after downstream "
            f"preprocessing stages.*{expected_stage}.*"
            "only matrix-changing preprocessing input.*"
            "SpsRuvBatchCorrectionConfig"
        ),
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_sample_metadata(phospho),
                total=total,
                organism=Organism.RAT,
                input_intensity_scale=input_intensity_scale,
                preprocessing_config=preprocessing_config,
                corrected_preprocessing_output=_correction_output(phospho + 1.0),
            )
        )


def _correction_output(
    corrected: pd.DataFrame,
    *,
    consumed_by_downstream: bool = False,
    input_matrix: pd.DataFrame | None = None,
    method: str = "sps_ruv_style",
    provenance: BatchCorrectionProvenance | None | Any = _DEFAULT_PROVENANCE,
    status: str = "applied",
) -> CorrectedPreprocessingOutput:
    resolved_provenance = (
        _complete_sps_ruv_provenance(corrected, input_matrix=input_matrix)
        if provenance is _DEFAULT_PROVENANCE
        else provenance
    )
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
            status=status,
            policy=BatchCorrectionPolicy(
                method=method,
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
        diagnostics={"executor": {"status": status, "method": method}},
        provenance=cast(BatchCorrectionProvenance | None, resolved_provenance),
        consumed_by_downstream=consumed_by_downstream,
    )


def _complete_sps_ruv_provenance(
    corrected: pd.DataFrame,
    *,
    input_matrix: pd.DataFrame | None = None,
) -> BatchCorrectionProvenance:
    resolved_input = (
        _resolved_correction_matrix(_phospho())
        if input_matrix is None
        else input_matrix.copy(deep=True)
    )
    mask = pd.DataFrame(
        True,
        index=corrected.index.copy(),
        columns=corrected.columns.copy(),
    )
    return BatchCorrectionProvenance(
        requested_method="sps_ruv_style",
        resolved_parameters={
            "method": "sps_ruv_style",
            "n_unwanted_factors": 1,
            "source": "external_corrected_preprocessing_output",
        },
        preprocessing_stage_order=(
            "missing_data",
            "batch_correction",
            "downstream_workflows",
        ),
        control_site_source={
            "source_type": "caller_supplied",
            "organism": "rat",
            "identifier_namespace": "site_key",
            "source_name": "manual-curated-controls",
            "source_version": "manual-v1",
            "license": "caller local use",
            "redistribution": "not redistributed",
        },
        selected_site_key_rows=tuple(str(row) for row in corrected.index[:2]),
        batch_metadata={
            "column": "batch",
            "levels": ["run_1", "run_2"],
            "sample_order": list(corrected.columns.astype(str)),
        },
        replicate_metadata=None,
        design_metadata={
            "condition_columns": ["condition"],
            "preserve_condition_effects": True,
        },
        missing_value_policy={
            "policy": "reject_missing",
            "imputation_policy": "none",
        },
        observation_masks=(
            fingerprint_matrix(
                mask.astype("int8"),
                name="batch_correction.native.observation_mask",
            ),
        ),
        input_matrix_fingerprint=fingerprint_matrix(
            resolved_input,
            name="batch_correction.native.input",
        ),
        output_matrix_fingerprint=fingerprint_matrix(
            corrected,
            name="batch_correction.native.corrected",
        ),
        diagnostics={"executor": {"status": "applied", "method": "sps_ruv_style"}},
        warnings=(),
        phospy_version="test",
        python_version="3.test",
        dependency_versions=_COMPLETE_EXTERNAL_DEPENDENCY_VERSIONS,
    )


def _attached_batch_correction_provenance(dataset: object) -> BatchCorrectionProvenance:
    provenance = getattr(dataset, "provenance", None)
    assert provenance is not None
    for stage in provenance.preprocessing_stages:
        if stage.stage != "batch_correction":
            continue
        batch_provenance = stage.batch_correction_provenance
        assert batch_provenance is not None
        return batch_provenance
    raise AssertionError("batch_correction provenance was not attached")


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


def _total(phospho: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_1": [3.0, 2.0],
            "sample_2": [3.1, 2.1],
            "sample_3": [3.2, 2.2],
            "sample_4": [3.3, 2.3],
        },
        index=pd.Index(["MAPK14", "AKT1"], name=phospho.index.name),
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


def _sps_multi_condition_phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_1": [10.0, 5.0, 20.0],
            "sample_2": [12.0, 8.0, 21.0],
            "sample_3": [11.0, 6.0, 24.0],
            "sample_4": [13.0, 9.0, 25.0],
            "sample_5": [14.0, 9.0, 28.0],
            "sample_6": [16.0, 12.0, 29.0],
            "sample_7": [15.0, 10.0, 32.0],
            "sample_8": [17.0, 13.0, 33.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;", "SRC;Y416;"], name="site_id"),
    )


def _sps_site_key(protein_identifier: str, residue: str, position: str) -> str:
    return (
        "phospy:v1|organism=rat|protein_namespace=protein_id|"
        f"protein_identifier={protein_identifier}|residue={residue}|"
        f"position={position}"
    )


def _control_source_metadata() -> ControlSiteSourceMetadata:
    return ControlSiteSourceMetadata(
        organism="rat",
        identifier_namespace="site_key",
        source_name="manual-curated-controls",
        source_version="manual-v1",
        license="caller local use",
        redistribution="not redistributed",
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


def _sps_sample_metadata(phospho: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": ("run_1", "run_1", "run_2", "run_2"),
            "condition": ("control", "treated", "control", "treated"),
            "replicate": ("r1", "r2", "r2", "r1"),
        },
        index=phospho.columns.copy(),
    )


def _sps_multi_condition_sample_metadata(phospho: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": (
                "run_1",
                "run_1",
                "run_1",
                "run_1",
                "run_2",
                "run_2",
                "run_2",
                "run_2",
            ),
            "condition": (
                "control",
                "treated",
                "control",
                "treated",
                "control",
                "treated",
                "control",
                "treated",
            ),
            "timepoint": (
                "early",
                "early",
                "late",
                "late",
                "early",
                "early",
                "late",
                "late",
            ),
            "replicate": ("r1", "r2", "r3", "r4", "r2", "r1", "r4", "r3"),
        },
        index=phospho.columns.copy(),
    )
