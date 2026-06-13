from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import PhosPyInputError
from phospy.science.datasets.builders.preprocessing import DatasetPreprocessor
from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from phospy.validation.datasets.batch_correction import (
    BatchCorrectionAdequacyValidator,
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


def test_batch_validation_fails_before_preprocessing_pipeline_execution() -> None:
    preprocessor = DatasetPreprocessor(pipeline=_PipelineShouldNotRun())

    with pytest.raises(
        PhosPyInputError,
        match="perfectly confounded",
    ):
        preprocessor.run(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            sample_metadata=_confounded_sample_metadata(),
            plan=PreprocessingPlan(batch_correction_method="linear_residualize_batch"),
        )


class _PipelineShouldNotRun:
    def run_with_trace(self, state: object) -> object:
        raise AssertionError("preprocessing pipeline ran before batch validation")


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
