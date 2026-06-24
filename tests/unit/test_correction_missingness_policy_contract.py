from __future__ import annotations

import json

import pandas as pd
import pytest

from phospy.contracts.configs.preprocessing import (
    CorrectedMissingCellAction,
    CorrectionMaskPolicy,
    CorrectionMissingnessCompatibilityValidator,
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
                method_parameters={"min_observed_values": 2},  # type: ignore[arg-type]
            ),
        )


def test_correction_mask_policy_must_preserve_observation_mask() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="must preserve the observation mask",
    ):
        CorrectionMaskPolicy(preserve_observation_mask=False)


def test_complete_matrix_passes_without_missingness_policy_when_allowed() -> None:
    CorrectionMissingnessCompatibilityValidator().run(
        phospho=_complete_matrix(),
        allow_complete_observed_data=True,
    )


def test_missing_matrix_rejects_without_explicit_supported_policy() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="missing values but no explicit supported correction missingness policy",
    ):
        CorrectionMissingnessCompatibilityValidator().run(phospho=_missing_matrix())


def test_missing_matrix_passes_with_observation_mask_preserving_policy() -> None:
    CorrectionMissingnessCompatibilityValidator().run(
        phospho=_missing_matrix(),
        policy=_mask_preserving_policy(),
    )


def test_missing_matrix_rejects_temporary_imputation_without_observation_mask() -> None:
    policy = CorrectionMissingnessPolicy(
        temporary_imputation=TemporaryImputationPolicy(
            allowed=True,
            method=TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY,
            method_parameters={"min_observed_values": 2},  # type: ignore[arg-type]
        ),
        originally_missing_cells_tracked_by=(
            OriginallyMissingCellTracking.EXISTING_IMPUTATION_PROVENANCE
        ),
    )

    with pytest.raises(PhosPyInputError, match="observation mask"):
        CorrectionMissingnessCompatibilityValidator().run(
            phospho=_missing_matrix(),
            policy=policy,
        )


def test_complete_matrix_rejects_temporary_imputation_policy_without_mask() -> None:
    policy = CorrectionMissingnessPolicy(
        temporary_imputation=TemporaryImputationPolicy(
            allowed=True,
            method=TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY,
            method_parameters={"min_observed_values": 2},  # type: ignore[arg-type]
        ),
        originally_missing_cells_tracked_by=(
            OriginallyMissingCellTracking.EXISTING_IMPUTATION_PROVENANCE
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="temporary imputation requires an observation mask",
    ):
        CorrectionMissingnessCompatibilityValidator().run(
            phospho=_complete_matrix(),
            policy=policy,
        )


def test_temporary_imputation_requires_policy_parameters() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="missing required parameter\\(s\\).*min_observed_values",
    ):
        TemporaryImputationPolicy(
            allowed=True,
            method=TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY,
        )


def test_random_temporary_imputation_requires_seed() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="random_seed is required.*method='minprob_temporary'",
    ):
        TemporaryImputationPolicy(
            allowed=True,
            method=TemporaryImputationMethod.MINPROB_TEMPORARY,
            method_parameters={"q": 0.01, "width": 0.3},  # type: ignore[arg-type]
        )


def test_unsupported_missingness_policy_fails_before_execution() -> None:
    unsupported = CorrectionMissingnessPolicy(
        supported=False,
        unsupported_reason="future correction backend is not implemented",
    )

    with pytest.raises(
        PhosPyInputError,
        match="policy is unsupported: future correction backend is not implemented",
    ):
        CorrectionMissingnessCompatibilityValidator().run(
            phospho=_complete_matrix(),
            policy=unsupported,
        )


def test_missing_matrix_rejects_unsupported_mask_action_before_execution() -> None:
    policy = CorrectionMissingnessPolicy(
        originally_missing_cells_tracked_by=OriginallyMissingCellTracking.OBSERVATION_MASK,
        correction_mask_policy=CorrectionMaskPolicy(
            corrected_missing_cell_action=CorrectedMissingCellAction.UNSUPPORTED,
            supported=False,
            unsupported_reason="exporting missing flags is not implemented",
        ),
        observation_mask=_observation_mask(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="correction mask policy is unsupported",
    ):
        CorrectionMissingnessCompatibilityValidator().run(
            phospho=_missing_matrix(),
            policy=policy,
        )


def test_missing_matrix_rejects_unsafe_row_sample_eligibility_state() -> None:
    policy = _mask_preserving_policy(
        row_sample_eligibility_impact=(
            RowSampleEligibilityImpact.REQUIRE_COMPLETE_CASES
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="row/sample eligibility impact 'require_complete_cases' makes correction unsafe",
    ):
        CorrectionMissingnessCompatibilityValidator().run(
            phospho=_missing_matrix(),
            policy=policy,
        )


def test_observation_mask_must_match_matrix_missing_coordinates() -> None:
    policy = CorrectionMissingnessPolicy(
        originally_missing_cells_tracked_by=OriginallyMissingCellTracking.OBSERVATION_MASK,
        observation_mask=ObservationMask(
            feature_ids=("site_a", "site_b"),
            sample_ids=("sample_1", "sample_2"),
            originally_missing_cells=(("site_a", "sample_2"),),
        ),
    )

    with pytest.raises(
        PhosPyInputError,
        match="originally_missing_cells must match",
    ):
        CorrectionMissingnessCompatibilityValidator().run(
            phospho=_missing_matrix(),
            policy=policy,
        )


def _complete_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_1": [1.0, 2.0],
            "sample_2": [3.0, 4.0],
        },
        index=pd.Index(["site_a", "site_b"], name="site_key"),
    )


def _missing_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_1": [1.0, float("nan")],
            "sample_2": [3.0, 4.0],
        },
        index=pd.Index(["site_a", "site_b"], name="site_key"),
    )


def _observation_mask() -> ObservationMask:
    return ObservationMask(
        feature_ids=("site_a", "site_b"),
        sample_ids=("sample_1", "sample_2"),
        originally_missing_cells=(("site_b", "sample_1"),),
    )


def _mask_preserving_policy(
    *,
    row_sample_eligibility_impact: RowSampleEligibilityImpact = (
        RowSampleEligibilityImpact.NO_CHANGE
    ),
) -> CorrectionMissingnessPolicy:
    return CorrectionMissingnessPolicy(
        temporary_imputation=TemporaryImputationPolicy(
            allowed=True,
            method=TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY,
            method_parameters={"min_observed_values": 2},  # type: ignore[arg-type]
        ),
        originally_missing_cells_tracked_by=OriginallyMissingCellTracking.OBSERVATION_MASK,
        correction_mask_policy=CorrectionMaskPolicy(
            corrected_missing_cell_action=CorrectedMissingCellAction.RESTORE_MISSING,
        ),
        row_sample_eligibility_impact=row_sample_eligibility_impact,
        observation_mask=_observation_mask(),
    )
