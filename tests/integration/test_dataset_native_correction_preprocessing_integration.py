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
    DatasetPreprocessingConfig,
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

pytestmark = pytest.mark.integration

_DEFAULT_PROVENANCE: Any = object()


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
    assert provenance.output_matrix_fingerprint is not None
    assert (
        provenance.output_matrix_fingerprint.name == "batch_correction.native.corrected"
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
        match="control source.*organism|control source.*identifier_namespace",
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
                        )
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
                            (_sps_site_key("MAPK14", "Y", "182"),)
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
                            )
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
    assert float(cast(float, table.loc[site_one, "logFC"])) == pytest.approx(4.0)
    assert float(cast(float, table.loc[site_two, "logFC"])) == pytest.approx(-1.0)


def _correction_output(
    corrected: pd.DataFrame,
    *,
    consumed_by_downstream: bool = False,
    provenance: BatchCorrectionProvenance | None | Any = _DEFAULT_PROVENANCE,
) -> CorrectedPreprocessingOutput:
    resolved_provenance = (
        _complete_sps_ruv_provenance(corrected)
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
        provenance=cast(BatchCorrectionProvenance | None, resolved_provenance),
        consumed_by_downstream=consumed_by_downstream,
    )


def _complete_sps_ruv_provenance(corrected: pd.DataFrame) -> BatchCorrectionProvenance:
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
            "source_version_unavailable_reason": "caller-local controls",
        },
        selected_site_key_rows=(str(corrected.index[0]),),
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
            corrected,
            name="batch_correction.native.input",
        ),
        output_matrix_fingerprint=fingerprint_matrix(
            corrected,
            name="batch_correction.native.corrected",
        ),
        diagnostics={"executor": {"status": "applied", "method": "sps_ruv_style"}},
        warnings=(),
        phospy_version="test",
        dependency_versions={},
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


def _sps_site_key(protein_identifier: str, residue: str, position: str) -> str:
    return (
        "phospy:v1|organism=rat|protein_namespace=protein_id|"
        f"protein_identifier={protein_identifier}|residue={residue}|"
        f"position={position}"
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
            "replicate": ("r1", "r1", "r2", "r2"),
        },
        index=phospho.columns.copy(),
    )
