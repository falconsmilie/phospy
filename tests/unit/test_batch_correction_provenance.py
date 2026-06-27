from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import PhosPyInputError
from phospy.provenance import (
    BatchCorrectionProvenance,
    BatchCorrectionRejectedEntity,
    batch_correction_provenance_from_payload,
    batch_correction_provenance_to_payload,
    fingerprint_matrix,
)
from phospy.validation.datasets.batch_correction import (
    validate_applied_native_sps_ruv_correction_provenance,
)

_COMPLETE_EXTERNAL_DEPENDENCY_VERSIONS = {
    "numpy": "test-numpy",
    "pandas": "test-pandas",
    "scipy": "test-scipy",
    "scikit-learn": "test-scikit-learn",
}


def _matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [4.0, 5.0, 6.0],
        },
        index=pd.Index(["AKT1_S473", "MAPK1_T202", "GSK3B_S9"], name="site_key"),
    )


def _mask() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [True, True, False],
            "sample_b": [True, False, True],
        },
        index=_matrix().index.copy(),
    )


def test_batch_correction_provenance_serialization_records_required_fields() -> None:
    input_fingerprint = fingerprint_matrix(_matrix(), name="batch_correction.input")
    output_fingerprint = fingerprint_matrix(_matrix(), name="batch_correction.output")
    mask_fingerprint = fingerprint_matrix(_mask(), name="batch_correction.observed")
    provenance = BatchCorrectionProvenance(
        requested_method="sps_ruv_style",
        resolved_parameters={
            "method": "sps_ruv_style",
            "n_unwanted_factors": 2,
            "temporary_imputation": False,
        },
        preprocessing_stage_order=(
            "intensity_transform",
            "batch_correction",
            "missing_data",
        ),
        control_site_source={
            "source_type": "caller_supplied",
            "identifier_namespace": "site_key",
            "source_version": "manual-v1",
        },
        selected_site_key_rows=("AKT1_S473", "GSK3B_S9"),
        batch_metadata={
            "column": "batch",
            "levels": ["run_1", "run_2"],
            "sample_to_batch": {"sample_a": "run_1", "sample_b": "run_2"},
        },
        replicate_metadata={
            "column": "bio_rep",
            "sample_to_replicate": {"sample_a": "rep_1", "sample_b": "rep_1"},
        },
        design_metadata={
            "condition_columns": ["condition"],
            "rank": 2,
            "preserve_condition_effects": True,
        },
        missing_value_policy={
            "policy": "reject_missing",
            "imputation_policy": "none",
        },
        observation_masks=(mask_fingerprint,),
        input_matrix_fingerprint=input_fingerprint,
        output_matrix_fingerprint=output_fingerprint,
        diagnostics={"condition_number": 3.5, "status": "recorded_only"},
        warnings=("future correction not executed",),
        rejected_entities=(
            BatchCorrectionRejectedEntity(
                entity_type="site",
                identifier="MAPK1_T202",
                reason="not_selected_control",
                details={"row_position": 1},
            ),
            BatchCorrectionRejectedEntity(
                entity_type="sample",
                identifier="sample_b",
                reason="missing_replicate_metadata",
            ),
        ),
        phospy_version="1.2.3",
        python_version="3.12.0",
        dependency_versions={"numpy": "2.0.0", "limma": None},
    )

    payload = batch_correction_provenance_to_payload(provenance)
    restored = batch_correction_provenance_from_payload(payload)

    expected_keys = {
        "requested_method",
        "resolved_parameters",
        "preprocessing_stage_order",
        "control_site_source",
        "selected_site_key_rows",
        "batch_metadata",
        "replicate_metadata",
        "design_metadata",
        "missing_value_policy",
        "observation_masks",
        "input_matrix_fingerprint",
        "output_matrix_fingerprint",
        "diagnostics",
        "warnings",
        "rejected_entities",
        "phospy_version",
        "python_version",
        "dependency_versions",
    }
    assert expected_keys.issubset(payload.keys())
    assert restored == provenance
    assert payload["selected_site_key_rows"] == ["AKT1_S473", "GSK3B_S9"]
    assert payload["warnings"] == ["future correction not executed"]
    assert restored.python_version == "3.12.0"


def test_batch_correction_provenance_can_record_rejected_without_output_matrix() -> (
    None
):
    provenance = BatchCorrectionProvenance(
        requested_method="sps_ruv_style",
        resolved_parameters={"method": "sps_ruv_style"},
        preprocessing_stage_order=("batch_correction",),
        control_site_source={"source_type": "dataset_metadata"},
        selected_site_key_rows=(),
        batch_metadata={"column": "batch"},
        replicate_metadata=None,
        design_metadata={"condition_columns": ["condition"]},
        missing_value_policy={"policy": "reject_missing"},
        observation_masks=(),
        input_matrix_fingerprint=fingerprint_matrix(
            _matrix(),
            name="batch_correction.input",
        ),
        output_matrix_fingerprint=None,
        diagnostics={"status": "rejected"},
        warnings=("correction rejected before execution",),
        rejected_entities=(
            BatchCorrectionRejectedEntity(
                entity_type="row",
                identifier="MAPK1_T202",
                reason="missing_value",
            ),
        ),
        phospy_version="unknown",
        dependency_versions={},
    )

    payload = batch_correction_provenance_to_payload(provenance)
    restored = batch_correction_provenance_from_payload(payload)

    assert payload["output_matrix_fingerprint"] is None
    assert restored.output_matrix_fingerprint is None
    assert restored.rejected_entities[0].entity_type == "row"


def test_matrix_fingerprint_is_deterministic_for_identical_matrices() -> None:
    first = fingerprint_matrix(_matrix(), name="batch_correction.input")
    second = fingerprint_matrix(
        _matrix().copy(deep=True), name="batch_correction.input"
    )

    assert first.exact_hash_value == second.exact_hash_value
    assert first.tolerance_hash_value == second.tolerance_hash_value


def test_matrix_fingerprint_changes_when_values_change() -> None:
    changed = _matrix()
    changed.loc["AKT1_S473", "sample_a"] = 9.0

    assert (
        fingerprint_matrix(
            _matrix(),
            name="batch_correction.input",
        ).tolerance_hash_value
        != fingerprint_matrix(
            changed, name="batch_correction.input"
        ).tolerance_hash_value
    )


def test_matrix_fingerprint_changes_when_row_order_changes() -> None:
    reordered = _matrix().iloc[[2, 1, 0], :]

    assert (
        fingerprint_matrix(
            _matrix(),
            name="batch_correction.input",
        ).tolerance_hash_value
        != fingerprint_matrix(
            reordered,
            name="batch_correction.input",
        ).tolerance_hash_value
    )


def test_matrix_fingerprint_changes_when_column_order_changes() -> None:
    reordered = _matrix().loc[:, ["sample_b", "sample_a"]]

    assert (
        fingerprint_matrix(
            _matrix(),
            name="batch_correction.input",
        ).tolerance_hash_value
        != fingerprint_matrix(
            reordered,
            name="batch_correction.input",
        ).tolerance_hash_value
    )


def test_applied_provenance_rejects_duplicate_selected_site_key_rows() -> None:
    provenance = _complete_applied_sps_ruv_provenance(
        selected_site_key_rows=("AKT1_S473", "GSK3B_S9", "AKT1_S473"),
    )

    with pytest.raises(
        PhosPyInputError,
        match="selected_site_key_rows.*duplicate",
    ):
        validate_applied_native_sps_ruv_correction_provenance(
            method="sps_ruv_style",
            status="applied",
            provenance=provenance,
        )


def test_applied_provenance_rejects_blank_selected_site_key_rows() -> None:
    provenance = _complete_applied_sps_ruv_provenance(
        selected_site_key_rows=("AKT1_S473", " "),
    )

    with pytest.raises(
        PhosPyInputError,
        match="selected_site_key_rows.*blank",
    ):
        validate_applied_native_sps_ruv_correction_provenance(
            method="sps_ruv_style",
            status="applied",
            provenance=provenance,
        )


def test_applied_provenance_rejects_sentinel_selected_site_key_rows() -> None:
    provenance = _complete_applied_sps_ruv_provenance(
        selected_site_key_rows=("AKT1_S473", "NoT_PrOvIdEd"),
    )

    with pytest.raises(
        PhosPyInputError,
        match="selected_site_key_rows.*sentinel",
    ):
        validate_applied_native_sps_ruv_correction_provenance(
            method="sps_ruv_style",
            status="applied",
            provenance=provenance,
        )


def test_applied_provenance_counts_unique_selected_controls_for_factor_count() -> None:
    provenance = _complete_applied_sps_ruv_provenance(
        selected_site_key_rows=("site_a", "site_a"),
        n_unwanted_factors=1,
    )

    with pytest.raises(
        PhosPyInputError,
        match="unique_selected_controls=1.*required_selected_controls=2.*duplicate",
    ):
        validate_applied_native_sps_ruv_correction_provenance(
            method="sps_ruv_style",
            status="applied",
            provenance=provenance,
        )


def _complete_applied_sps_ruv_provenance(
    *,
    selected_site_key_rows: tuple[str, ...] = ("AKT1_S473", "GSK3B_S9"),
    n_unwanted_factors: int = 1,
) -> BatchCorrectionProvenance:
    matrix = _matrix()
    return BatchCorrectionProvenance(
        requested_method="sps_ruv_style",
        resolved_parameters={
            "method": "sps_ruv_style",
            "n_unwanted_factors": n_unwanted_factors,
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
        selected_site_key_rows=selected_site_key_rows,
        batch_metadata={
            "column": "batch",
            "levels": ["run_1", "run_2"],
            "sample_order": list(matrix.columns.astype(str)),
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
                _mask().astype("int8"),
                name="batch_correction.native.observation_mask",
            ),
        ),
        input_matrix_fingerprint=fingerprint_matrix(
            matrix,
            name="batch_correction.native.input",
        ),
        output_matrix_fingerprint=fingerprint_matrix(
            matrix,
            name="batch_correction.native.corrected",
        ),
        diagnostics={"executor": {"status": "applied", "method": "sps_ruv_style"}},
        warnings=(),
        phospy_version="test",
        python_version="3.test",
        dependency_versions=_COMPLETE_EXTERNAL_DEPENDENCY_VERSIONS,
    )
