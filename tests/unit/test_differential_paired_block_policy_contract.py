from __future__ import annotations

from dataclasses import fields
from typing import get_type_hints

import pytest

from phospy.api import (
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    PAIRED_DESIGN_POLICY_REJECT,
    SUPPORTED_PAIRED_DESIGN_POLICIES,
    DifferentialAnalysisConfig,
    ExperimentalDesign,
    SampleDesignRecord,
)
from phospy.errors import WorkflowValidationError


def test_paired_design_policy_defaults_to_reject() -> None:
    config = DifferentialAnalysisConfig()

    assert config.paired_design_policy == PAIRED_DESIGN_POLICY_REJECT


def test_fixed_block_policy_can_be_declared() -> None:
    config = DifferentialAnalysisConfig(
        paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK
    )

    assert config.paired_design_policy == "fixed_block"


def test_block_identifier_can_be_supplied_explicitly() -> None:
    record = SampleDesignRecord(
        sample_id="donor_1_treated",
        condition="treated",
        block_id=" donor_1 ",
    )

    assert record.block_id == "donor_1"


def test_sample_names_do_not_infer_block_identifier() -> None:
    record = SampleDesignRecord(
        sample_id="donor_1_treated",
        condition="treated",
    )

    assert record.block_id is None


def test_condition_only_design_still_uses_narrow_sample_contract() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="control_1", condition="control"),
            SampleDesignRecord(sample_id="control_2", condition="control"),
            SampleDesignRecord(sample_id="treated_1", condition="treated"),
            SampleDesignRecord(sample_id="treated_2", condition="treated"),
        )
    )

    assert design.sample_ids() == (
        "control_1",
        "control_2",
        "treated_1",
        "treated_2",
    )
    assert design.condition_labels() == ("control", "treated")
    assert all(record.block_id is None for record in design.samples)


def test_unsupported_paired_design_policy_value_is_rejected() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="differential.paired_design_policy",
    ):
        DifferentialAnalysisConfig(
            paired_design_policy="mixed_effects"  # type: ignore[arg-type]
        )


def test_public_paired_block_contract_remains_narrow() -> None:
    sample_fields = {field.name for field in fields(SampleDesignRecord)}
    config_hints = get_type_hints(DifferentialAnalysisConfig)

    assert "block_id" in sample_fields
    assert {"block", "pair_id", "subject_id"}.isdisjoint(sample_fields)
    assert "paired_design_policy" in config_hints
    assert set(SUPPORTED_PAIRED_DESIGN_POLICIES) == {"reject", "fixed_block"}
