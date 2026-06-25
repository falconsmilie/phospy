from __future__ import annotations

from collections.abc import Mapping
from typing import cast

import pandas as pd

from phospy.api import (
    AnalysisReadyDatasetBuilder,
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
    DatasetBatchCorrectionConfig,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    Organism,
)
from phospy.provenance.serialization import from_payload, to_payload
from phospy.science.datasets.models import DatasetPreprocessingReport


def _batch_correction_policy(method: str = "none") -> BatchCorrectionPolicy:
    return BatchCorrectionPolicy(
        method=method,
        batch_column="batch",
        condition_column="condition",
        design_preservation_policy="preserve_condition_effects",
        preserve_condition_effects=True,
    )


def test_disabled_batch_correction_report_can_represent_noop_status() -> None:
    report = BatchCorrectionReport(
        status="disabled",
        policy=_batch_correction_policy(),
        diagnostics=BatchCorrectionDiagnostics(
            number_of_batches=2,
            batch_levels=("run_1", "run_2"),
            condition_levels=("control", "treated"),
            confounding_check_status="not_applicable",
            matrix_shape_before=(2, 4),
            matrix_shape_after=(2, 4),
            limitations=("batch correction disabled by preprocessing configuration",),
        ),
    )

    assert report.status == "disabled"
    assert report.method == "none"
    assert report.batch_column == "batch"
    assert report.condition_column == "condition"
    assert report.number_of_batches == 2
    assert report.matrix_shape_before == (2, 4)
    assert report.matrix_shape_after == (2, 4)
    assert report.warnings == ()


def test_applied_batch_correction_report_construction_is_typed() -> None:
    report = BatchCorrectionReport(
        status="applied",
        policy=_batch_correction_policy("linear_residualize_batch"),
        diagnostics=BatchCorrectionDiagnostics(
            number_of_batches=2,
            batch_levels=("run_1", "run_2"),
            condition_levels=("control", "treated"),
            confounding_check_status="passed",
            matrix_shape_before=(3, 4),
            matrix_shape_after=(3, 4),
            limitations=("linear residualisation preserves matrix shape",),
        ),
    )

    assert report.status == "applied"
    assert report.method == "linear_residualize_batch"
    assert report.design_preservation_policy == "preserve_condition_effects"
    assert report.confounding_check_status == "passed"


def test_rejected_batch_correction_report_construction_is_explicit() -> None:
    report = BatchCorrectionReport(
        status="rejected",
        policy=_batch_correction_policy("linear_residualize_batch"),
        diagnostics=BatchCorrectionDiagnostics(
            number_of_batches=1,
            batch_levels=("run_1",),
            condition_levels=("control", "treated"),
            confounding_check_status="confounded",
            matrix_shape_before=(3, 4),
            matrix_shape_after=(3, 4),
            warnings=("condition is confounded with batch",),
            limitations=("batch correction was not applied",),
        ),
    )

    assert report.status == "rejected"
    assert report.confounding_check_status == "confounded"
    assert report.warnings == ("condition is confounded with batch",)
    assert "BatchCorrectionReport" in repr(report)


def test_batch_correction_report_payload_is_structured() -> None:
    report = BatchCorrectionReport(
        status="disabled",
        policy=_batch_correction_policy(),
        diagnostics=BatchCorrectionDiagnostics(
            number_of_batches=2,
            batch_levels=("run_1", "run_2"),
            condition_levels=("control", "treated"),
            confounding_check_status="not_applicable",
            matrix_shape_before=(2, 4),
            matrix_shape_after=(2, 4),
        ),
    )

    payload = report.to_payload()

    assert payload["status"] == "disabled"
    assert payload["method"] == "none"
    assert payload["batch_levels"] == ["run_1", "run_2"]
    assert payload["matrix_shape_before"] == [2, 4]
    assert isinstance(payload["policy"], dict)
    assert isinstance(payload["diagnostics"], dict)


def test_batch_correction_report_preserves_plural_condition_columns() -> None:
    report = BatchCorrectionReport(
        status="applied",
        policy=BatchCorrectionPolicy(
            method="sps_ruv_style",
            batch_column="batch",
            condition_column="condition",
            condition_columns=("condition", "timepoint"),
        ),
    )

    payload = report.to_payload()

    assert report.condition_column == "condition"
    assert report.condition_columns == ("condition", "timepoint")
    assert payload["condition_column"] == "condition"
    assert payload["condition_columns"] == ["condition", "timepoint"]
    policy_payload = payload["policy"]
    assert isinstance(policy_payload, dict)
    assert policy_payload["condition_columns"] == ["condition", "timepoint"]


def test_batch_correction_report_integrates_with_preprocessing_report() -> None:
    batch_report = BatchCorrectionReport(
        status="disabled",
        policy=_batch_correction_policy(),
        diagnostics=BatchCorrectionDiagnostics(matrix_shape_before=(1, 2)),
    )

    preprocessing_report = DatasetPreprocessingReport.from_rows(
        batch_correction=batch_report
    )

    assert preprocessing_report.batch_correction is batch_report
    assert preprocessing_report.batch_correction_summary() is batch_report


def test_dataset_builder_records_disabled_batch_correction_report_without_execution() -> (
    None
):
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            sample_metadata=_sample_metadata(),
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    assert built.preprocessing_report is not None
    report = built.preprocessing_report.batch_correction
    assert report is not None
    assert report.status == "disabled"
    assert report.method == "none"
    assert report.batch_levels == ("run_1", "run_2")
    assert report.condition_levels == ("control", "treated")
    assert report.matrix_shape_before == (2, 4)
    assert report.matrix_shape_after == (2, 4)
    assert "batch_correction" not in set(
        built.preprocessing_report.operations.loc[:, "stage"].astype(str).tolist()
    )


def test_dataset_builder_records_applied_declared_batch_correction_execution() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            sample_metadata=_sample_metadata(),
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                batch_correction=DatasetBatchCorrectionConfig(
                    method="linear_residualize_batch"
                )
            ),
        )
    )

    assert built.preprocessing_report is not None
    report = built.preprocessing_report.batch_correction
    assert report is not None
    assert report.status == "applied"
    assert report.method == "linear_residualize_batch"
    assert report.confounding_check_status == "passed"
    assert report.matrix_shape_before == report.matrix_shape_after == (2, 4)
    assert report.warnings == ()
    assert "batch_correction" in set(
        built.preprocessing_report.operations.loc[:, "stage"].astype(str).tolist()
    )


def test_dataset_builder_attaches_full_native_batch_correction_provenance() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            sample_metadata=_sample_metadata(),
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                batch_correction=DatasetBatchCorrectionConfig(
                    method="linear_residualize_batch"
                )
            ),
        )
    )

    assert built.provenance is not None
    stage = next(
        item
        for item in built.provenance.preprocessing_stages
        if item.stage == "batch_correction"
    )
    provenance = stage.batch_correction_provenance
    assert provenance is not None
    assert provenance.requested_method == "linear_residualize_batch"
    assert provenance.resolved_parameters["source"] == "native_preprocessing_stage"
    assert "batch_correction" in provenance.preprocessing_stage_order
    assert provenance.control_site_source["source_type"] == "not_applicable"
    assert provenance.selected_site_key_rows == ()
    assert provenance.batch_metadata["column"] == "batch"
    assert provenance.batch_metadata["sample_order"] == [
        "sample_1",
        "sample_2",
        "sample_3",
        "sample_4",
    ]
    assert provenance.design_metadata["condition_columns"] == ["condition"]
    assert provenance.replicate_metadata is None
    assert provenance.missing_value_policy["policy"] == (
        "reject_missing_at_batch_correction"
    )
    assert provenance.imputation_policy["policy"] == "forbid"
    assert provenance.observation_masks
    assert provenance.input_matrix_fingerprint.name == "batch_correction.native.input"
    assert provenance.output_matrix_fingerprint is not None
    assert (
        provenance.output_matrix_fingerprint.name == "batch_correction.native.corrected"
    )
    stage_diagnostics = cast(
        Mapping[str, object], provenance.diagnostics["stage_diagnostics"]
    )
    assert stage_diagnostics["status"] == "applied"
    assert provenance.warnings == ()
    assert provenance.rejected_entities == ()
    assert provenance.phospy_version
    assert provenance.python_version
    assert "numpy" in provenance.dependency_versions

    payload = to_payload(built.provenance)
    preprocessing_stages = cast(list[object], payload["preprocessing_stages"])
    stage_payload = next(
        item
        for item in preprocessing_stages
        if isinstance(item, dict) and item["stage"] == "batch_correction"
    )
    assert "batch_correction_provenance" in stage_payload
    restored = from_payload(payload)
    restored_stage = next(
        item
        for item in restored.preprocessing_stages
        if item.stage == "batch_correction"
    )
    assert restored_stage.batch_correction_provenance == provenance


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_1": [1.0, 2.0],
            "sample_2": [2.0, 3.0],
            "sample_3": [4.0, 5.0],
            "sample_4": [5.0, 6.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )


def _site_metadata() -> pd.DataFrame:
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
        index=_phospho().index.copy(),
    )


def _sample_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": ["run_1", "run_2", "run_1", "run_2"],
            "condition": ["control", "control", "treated", "treated"],
        },
        index=_phospho().columns.copy(),
    )
