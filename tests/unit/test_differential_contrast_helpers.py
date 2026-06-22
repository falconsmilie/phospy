from __future__ import annotations

import pytest

from phospy.api import (
    Contrast,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    SampleDesignRecord,
    all_pairwise_contrasts,
    contrasts_vs_control,
)
from phospy.errors import WorkflowValidationError


def _design(condition_order: tuple[str, ...]) -> ExperimentalDesign:
    records: list[SampleDesignRecord] = []
    for condition_index, condition in enumerate(condition_order):
        records.extend(
            (
                SampleDesignRecord(
                    sample_id=f"sample_{condition_index}_1",
                    condition=condition,
                    biological_replicate_id=f"replicate_{condition_index}_1",
                ),
                SampleDesignRecord(
                    sample_id=f"sample_{condition_index}_2",
                    condition=condition,
                    biological_replicate_id=f"replicate_{condition_index}_2",
                ),
            )
        )
    return ExperimentalDesign(samples=tuple(records))


def test_all_pairwise_contrasts_use_design_condition_order() -> None:
    design = _design(("control", "low_dose", "high_dose"))

    contrasts = all_pairwise_contrasts(design)

    assert contrasts == (
        Contrast(
            name="low_dose_vs_control",
            numerator_condition="low_dose",
            denominator_condition="control",
        ),
        Contrast(
            name="high_dose_vs_control",
            numerator_condition="high_dose",
            denominator_condition="control",
        ),
        Contrast(
            name="high_dose_vs_low_dose",
            numerator_condition="high_dose",
            denominator_condition="low_dose",
        ),
    )


def test_control_contrasts_use_treatment_vs_control_direction() -> None:
    design = _design(("control", "low_dose", "high_dose"))

    contrasts = contrasts_vs_control(design, "control")

    assert contrasts == (
        Contrast(
            name="low_dose_vs_control",
            numerator_condition="low_dose",
            denominator_condition="control",
        ),
        Contrast(
            name="high_dose_vs_control",
            numerator_condition="high_dose",
            denominator_condition="control",
        ),
    )


def test_pairwise_contrast_ordering_is_deterministic_for_design_order() -> None:
    design = _design(("treatment_b", "control", "treatment_a"))

    contrasts = all_pairwise_contrasts(design)

    assert tuple(contrast.name for contrast in contrasts) == (
        "control_vs_treatment_b",
        "treatment_a_vs_treatment_b",
        "treatment_a_vs_control",
    )
    assert contrasts == all_pairwise_contrasts(design)


def test_unknown_control_condition_fails_clearly() -> None:
    design = _design(("control", "treatment"))

    with pytest.raises(
        WorkflowValidationError,
        match="control condition 'vehicle' is not present",
    ):
        contrasts_vs_control(design, "vehicle")


def test_duplicate_generated_pairwise_names_fail_clearly() -> None:
    design = _design(("C", "A_vs_B", "B_vs_C", "A"))

    with pytest.raises(
        WorkflowValidationError,
        match="duplicate contrast names: A_vs_B_vs_C",
    ):
        all_pairwise_contrasts(design)


def test_helper_contrasts_are_request_compatible() -> None:
    design = _design(("control", "treatment"))
    contrasts = contrasts_vs_control(design, "control")

    request = DifferentialAnalysisRequest(
        dataset=object(),  # type: ignore[arg-type]
        design=design,
        contrasts=contrasts,
    )

    assert request.contrasts is contrasts
    assert request.contrasts[0].name == "treatment_vs_control"
    assert request.config.minimum_condition_replicates == 2


def test_manual_contrasts_still_work_unchanged() -> None:
    design = _design(("control", "treatment"))
    manual_contrasts = (
        Contrast(
            name="manual_treatment_minus_control",
            numerator_condition="treatment",
            denominator_condition="control",
        ),
    )

    request = DifferentialAnalysisRequest(
        dataset=object(),  # type: ignore[arg-type]
        design=design,
        contrasts=manual_contrasts,
    )

    assert request.contrasts is manual_contrasts
