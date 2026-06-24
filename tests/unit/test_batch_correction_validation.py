from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import PhosPyInputError
from phospy.science.datasets.builders.preprocessing import DatasetPreprocessor
from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from phospy.validation.datasets.batch_correction import (
    BatchCorrectionAdequacyValidator,
    BatchDesignMetadataValidator,
    DesignRankValidator,
    ReplicateStructureValidator,
    SampleMetadataAlignmentValidator,
)


def test_batch_validation_accepts_valid_removable_batch_effect_design() -> None:
    BatchCorrectionAdequacyValidator().run(
        batch_by_sample={
            "sample_1": "run_1",
            "sample_2": "run_1",
            "sample_3": "run_2",
            "sample_4": "run_2",
        },
        condition_by_sample={
            "sample_1": "control",
            "sample_2": "treated",
            "sample_3": "control",
            "sample_4": "treated",
        },
        sample_order=("sample_1", "sample_2", "sample_3", "sample_4"),
        preserve_condition_effects=True,
    )


def test_design_metadata_validator_accepts_valid_batch_condition_replicates() -> None:
    resolved = BatchDesignMetadataValidator().run(
        phospho=_phospho(),
        sample_metadata=_batch_design_sample_metadata(),
        batch_column="batch",
        condition_columns=("condition",),
        replicate_column="replicate",
        require_replicate_column=True,
    )

    assert resolved.sample_order == ("sample_1", "sample_2", "sample_3", "sample_4")
    assert resolved.batch_labels == ("run_1", "run_1", "run_2", "run_2")
    assert resolved.condition_labels == (
        "control",
        "treated",
        "control",
        "treated",
    )
    assert resolved.replicate_labels == ("r1", "r1", "r2", "r2")


def test_sample_metadata_alignment_validator_rejects_misaligned_metadata() -> None:
    sample_metadata = _batch_design_sample_metadata()
    sample_metadata.index = pd.Index(
        ["sample_1", "sample_2", "sample_3", "sample_extra"]
    )

    with pytest.raises(
        PhosPyInputError,
        match="misaligned with matrix columns.*sample_4.*sample_extra",
    ):
        SampleMetadataAlignmentValidator().run(
            phospho=_phospho(),
            sample_metadata=sample_metadata,
            required_columns=("batch", "condition"),
        )


def test_sample_metadata_alignment_validator_rejects_duplicated_metadata() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="sample_metadata.index contains duplicate sample labels.*sample_2",
    ):
        SampleMetadataAlignmentValidator().run(
            phospho=_phospho(),
            sample_metadata=_batch_design_sample_metadata(
                index=("sample_1", "sample_2", "sample_2", "sample_4")
            ),
            required_columns=("batch", "condition"),
        )


def test_design_metadata_validator_rejects_missing_condition_columns() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="missing required column 'condition'",
    ):
        BatchDesignMetadataValidator().run(
            phospho=_phospho(),
            sample_metadata=_batch_design_sample_metadata().drop(columns=["condition"]),
            batch_column="batch",
            condition_columns=("condition",),
        )


def test_replicate_structure_validator_rejects_missing_required_column() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="requires replicate_column metadata",
    ):
        ReplicateStructureValidator().run(
            sample_metadata=_batch_design_sample_metadata().drop(columns=["replicate"]),
            sample_order=("sample_1", "sample_2", "sample_3", "sample_4"),
            replicate_column=None,
            required=True,
        )


def test_design_rank_validator_rejects_rank_deficient_matrix() -> None:
    design = pd.DataFrame(
        {
            "condition_A": [1.0, 1.0, 0.0, 0.0],
            "duplicated_condition_A": [1.0, 1.0, 0.0, 0.0],
        },
        index=pd.Index(["sample_1", "sample_2", "sample_3", "sample_4"]),
    )

    with pytest.raises(
        PhosPyInputError,
        match="rank-deficient.*condition_A.*duplicated_condition_A",
    ):
        DesignRankValidator().run(design, context="test design")


def test_batch_validation_rejects_perfectly_confounding_batch_and_condition() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="perfectly confounded.*condition signal",
    ):
        BatchCorrectionAdequacyValidator().run(
            batch_by_sample={
                "sample_1": "run_1",
                "sample_2": "run_1",
                "sample_3": "run_2",
                "sample_4": "run_2",
            },
            condition_by_sample={
                "sample_1": "control",
                "sample_2": "control",
                "sample_3": "treated",
                "sample_4": "treated",
            },
            sample_order=("sample_1", "sample_2", "sample_3", "sample_4"),
            preserve_condition_effects=True,
        )


def test_batch_validation_rejects_missing_batch_labels() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="requires batch labels for every sample.*'sample_2'",
    ):
        BatchCorrectionAdequacyValidator().run(
            batch_by_sample={
                "sample_1": "run_1",
                "sample_3": "run_2",
                "sample_4": "run_2",
            },
            condition_by_sample={
                "sample_1": "control",
                "sample_2": "treated",
                "sample_3": "control",
                "sample_4": "treated",
            },
            sample_order=("sample_1", "sample_2", "sample_3", "sample_4"),
            preserve_condition_effects=True,
        )


def test_batch_validation_rejects_missing_condition_labels() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="requires condition labels for every sample.*'sample_4'",
    ):
        BatchCorrectionAdequacyValidator().run(
            batch_by_sample={
                "sample_1": "run_1",
                "sample_2": "run_1",
                "sample_3": "run_2",
                "sample_4": "run_2",
            },
            condition_by_sample={
                "sample_1": "control",
                "sample_2": "treated",
                "sample_3": "control",
                "sample_4": pd.NA,
            },
            sample_order=("sample_1", "sample_2", "sample_3", "sample_4"),
            preserve_condition_effects=True,
        )


def test_batch_validation_rejects_one_batch_level() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="at least two batch levels",
    ):
        BatchCorrectionAdequacyValidator().run(
            batch_by_sample={
                "sample_1": "run_1",
                "sample_2": "run_1",
            },
            condition_by_sample={
                "sample_1": "control",
                "sample_2": "treated",
            },
            sample_order=("sample_1", "sample_2"),
            preserve_condition_effects=True,
        )


def test_batch_validation_rejects_single_sample_design() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="at least two samples",
    ):
        BatchCorrectionAdequacyValidator().run(
            batch_by_sample={"sample_1": "run_1"},
            condition_by_sample={"sample_1": "control"},
            sample_order=("sample_1",),
            preserve_condition_effects=True,
        )


def test_batch_validation_rejects_non_estimable_design_without_residual_degrees() -> (
    None
):
    with pytest.raises(
        PhosPyInputError,
        match="requires more samples than estimable condition-plus-batch design",
    ):
        BatchCorrectionAdequacyValidator().run(
            batch_by_sample={
                "sample_1": "run_1",
                "sample_2": "run_1",
                "sample_3": "run_2",
                "sample_4": "run_2",
            },
            condition_by_sample={
                "sample_1": "condition_1",
                "sample_2": "condition_1",
                "sample_3": "condition_2",
                "sample_4": "condition_3",
            },
            sample_order=("sample_1", "sample_2", "sample_3", "sample_4"),
            preserve_condition_effects=True,
        )


def test_batch_validation_rejects_rank_deficient_preservation_design() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="condition preservation design is saturated",
    ):
        BatchCorrectionAdequacyValidator().run(
            batch_by_sample={
                "sample_1": "run_1",
                "sample_2": "run_1",
                "sample_3": "run_2",
                "sample_4": "run_2",
            },
            condition_by_sample={
                "sample_1": "condition_1",
                "sample_2": "condition_2",
                "sample_3": "condition_3",
                "sample_4": "condition_4",
            },
            sample_order=("sample_1", "sample_2", "sample_3", "sample_4"),
            preserve_condition_effects=True,
        )


def test_batch_validation_rejects_when_condition_effects_not_preserved() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="preserve_condition_effects=True",
    ):
        BatchCorrectionAdequacyValidator().run(
            batch_by_sample={
                "sample_1": "run_1",
                "sample_2": "run_1",
                "sample_3": "run_2",
                "sample_4": "run_2",
            },
            condition_by_sample={
                "sample_1": "control",
                "sample_2": "treated",
                "sample_3": "control",
                "sample_4": "treated",
            },
            sample_order=("sample_1", "sample_2", "sample_3", "sample_4"),
            preserve_condition_effects=False,
        )


def test_batch_validation_fails_during_preprocessing_stage_execution() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="perfectly confounded",
    ):
        DatasetPreprocessor().run(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            sample_metadata=_confounded_sample_metadata(),
            total=None,
            plan=PreprocessingPlan(
                localisation_mode="ignore",
                batch_correction_method="linear_residualize_batch",
                stage_order=("batch_correction",),
            ),
        )


def test_batch_correction_plan_rejects_requested_method_without_stage() -> None:
    with pytest.raises(PhosPyInputError, match="stage_order.*batch_correction"):
        PreprocessingPlan(batch_correction_method="linear_residualize_batch")


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
            "site_sequence": ["SEQ_A", "SEQ_R"],
            "localisation_confidence": [0.95, 0.92],
        },
        index=_phospho().index.copy(),
    )


def _confounded_sample_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": ["run_1", "run_1", "run_2", "run_2"],
            "condition": ["control", "control", "treated", "treated"],
        },
        index=_phospho().columns.copy(),
    )


def _batch_design_sample_metadata(
    *,
    index: tuple[str, ...] = ("sample_1", "sample_2", "sample_3", "sample_4"),
) -> pd.DataFrame:
    source = {
        "sample_1": ("run_1", "control", "r1"),
        "sample_2": ("run_1", "treated", "r1"),
        "sample_3": ("run_2", "control", "r2"),
        "sample_4": ("run_2", "treated", "r2"),
    }
    return pd.DataFrame(
        {
            "batch": [source[sample][0] for sample in index],
            "condition": [source[sample][1] for sample in index],
            "replicate": [source[sample][2] for sample in index],
        },
        index=pd.Index(index, name="sample"),
    )
