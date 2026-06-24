from __future__ import annotations

import json

import pytest

from phospy.contracts.configs.preprocessing import (
    CorrectedMissingCellAction,
    CorrectionMaskPolicy,
    CorrectionMissingnessPolicy,
    ObservationMask,
    OriginallyMissingCellTracking,
    RowSampleEligibilityImpact,
    TemporaryImputationMethod,
    TemporaryImputationPolicy,
)
from phospy.errors import PhosPyInputError


def test_observation_mask_distinguishes_original_missing_from_observed_cells() -> None:
    mask = ObservationMask(
        feature_ids=("site_a", "site_b"),
        sample_ids=("sample_1", "sample_2"),
        originally_missing_cells=(("site_a", "sample_2"),),
    )

    assert mask.is_originally_missing("site_a", "sample_2") is True
    assert mask.is_originally_observed("site_a", "sample_1") is True
    assert mask.is_originally_observed("site_b", "sample_2") is True
    assert mask.to_payload()["originally_missing_cells"] == [
        {"feature_id": "site_a", "sample_id": "sample_2"}
    ]


def test_observation_mask_rejects_unknown_or_duplicate_missing_coordinates() -> None:
    mask = ObservationMask(
        feature_ids=("site_a",),
        sample_ids=("sample_1",),
    )

    with pytest.raises(
        PhosPyInputError,
        match="feature_id must be present in feature_ids",
    ):
        mask.is_originally_missing("site_b", "sample_1")

    with pytest.raises(
        PhosPyInputError,
        match="feature_id must be present in feature_ids",
    ):
        ObservationMask(
            feature_ids=("site_a",),
            sample_ids=("sample_1",),
            originally_missing_cells=(("site_b", "sample_1"),),
        )

    with pytest.raises(
        PhosPyInputError,
        match="must not contain duplicate coordinates",
    ):
        ObservationMask(
            feature_ids=("site_a",),
            sample_ids=("sample_1",),
            originally_missing_cells=(
                ("site_a", "sample_1"),
                ("site_a", "sample_1"),
            ),
        )


def test_supported_missingness_policy_serializes_provenance_ready_payload() -> None:
    mask = ObservationMask(
        feature_ids=("site_a", "site_b"),
        sample_ids=("sample_1", "sample_2"),
        originally_missing_cells=(("site_b", "sample_1"),),
    )
    policy = CorrectionMissingnessPolicy(
        temporary_imputation=TemporaryImputationPolicy(
            allowed=True,
            method="row_median_temporary",  # type: ignore[arg-type]
            method_parameters={"min_observed_values": 2},  # type: ignore[arg-type]
            random_seed=123,
        ),
        originally_missing_cells_tracked_by="observation_mask",  # type: ignore[arg-type]
        correction_mask_policy=CorrectionMaskPolicy(
            corrected_missing_cell_action="restore_missing",  # type: ignore[arg-type]
        ),
        row_sample_eligibility_impact=(  # type: ignore[arg-type]
            "exclude_rows_with_insufficient_observed_values"
        ),
        observation_mask=mask,
    )

    payload = policy.to_payload()

    assert payload["temporary_imputation"] == {
        "allowed": True,
        "method": "row_median_temporary",
        "method_parameters": {"min_observed_values": 2},
        "random_seed": 123,
        "supported": True,
        "unsupported_reason": None,
        "imputed_values_are_observed_evidence": False,
    }
    assert payload["originally_missing_cells_tracked_by"] == "observation_mask"
    assert payload["correction_mask_policy"] == {
        "corrected_missing_cell_action": "restore_missing",
        "preserve_observation_mask": True,
        "supported": True,
        "unsupported_reason": None,
    }
    assert (
        payload["row_sample_eligibility_impact"]
        == "exclude_rows_with_insufficient_observed_values"
    )
    json.dumps(payload)


def test_missingness_policy_can_represent_unsupported_states() -> None:
    policy = CorrectionMissingnessPolicy(
        temporary_imputation=TemporaryImputationPolicy(
            allowed=True,
            method=TemporaryImputationMethod.UNSUPPORTED,
            supported=False,
            unsupported_reason="future method is not implemented",
        ),
        originally_missing_cells_tracked_by=OriginallyMissingCellTracking.UNSUPPORTED,
        correction_mask_policy=CorrectionMaskPolicy(
            corrected_missing_cell_action=CorrectedMissingCellAction.UNSUPPORTED,
            supported=False,
            unsupported_reason="withheld-cell export is not implemented",
        ),
        row_sample_eligibility_impact=RowSampleEligibilityImpact.UNSUPPORTED,
        supported=False,
        unsupported_reason="native correction is not implemented",
    )

    payload = policy.to_payload()

    assert payload["supported"] is False
    assert payload["unsupported_reason"] == "native correction is not implemented"
    assert payload["temporary_imputation"]["supported"] is False  # type: ignore[index]
    assert (
        payload["correction_mask_policy"]["corrected_missing_cell_action"]  # type: ignore[index]
        == "unsupported"
    )


def test_temporary_imputation_policy_rejects_observed_evidence_claims() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="must not treat temporary imputed values as observed",
    ):
        TemporaryImputationPolicy(
            allowed=True,
            method=TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY,
            imputed_values_are_observed_evidence=True,
        )


def test_temporary_imputation_requires_original_missing_cell_tracking() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="must track originally missing cells",
    ):
        CorrectionMissingnessPolicy(
            temporary_imputation=TemporaryImputationPolicy(
                allowed=True,
                method=TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY,
            ),
        )


def test_correction_mask_policy_must_preserve_observation_mask() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="must preserve the observation mask",
    ):
        CorrectionMaskPolicy(preserve_observation_mask=False)
