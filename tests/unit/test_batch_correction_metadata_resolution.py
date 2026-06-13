from __future__ import annotations

import pandas as pd
import pytest

from phospy.api import (
    AnalysisReadyDatasetBuilder,
    DatasetBatchCorrectionConfig,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    Organism,
)
from phospy.errors import PhosPyInputError
from phospy.science.datasets.preprocessing.batch_correction_metadata import (
    BatchCorrectionMetadataResolver,
)


def test_resolves_batch_and_condition_metadata_in_matrix_sample_order() -> None:
    resolved = BatchCorrectionMetadataResolver().run(
        phospho=_phospho(),
        sample_metadata=_sample_metadata(index=["sample_c", "sample_a", "sample_b"]),
        batch_column="batch",
        condition_column="condition",
    )

    assert resolved.sample_order == ("sample_b", "sample_a", "sample_c")
    assert resolved.batch_by_sample == {
        "sample_b": "run_2",
        "sample_a": "run_1",
        "sample_c": "run_1",
    }
    assert resolved.condition_by_sample == {
        "sample_b": "treated",
        "sample_a": "control",
        "sample_c": "treated",
    }
    assert resolved.batch_labels == ("run_2", "run_1", "run_1")
    assert resolved.condition_labels == ("treated", "control", "treated")


def test_rejects_missing_batch_column() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="batch_column.*missing sample_metadata column 'batch'",
    ):
        BatchCorrectionMetadataResolver().run(
            phospho=_phospho(),
            sample_metadata=_sample_metadata().drop(columns=["batch"]),
            batch_column="batch",
            condition_column="condition",
        )


def test_rejects_missing_condition_column() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="condition_column.*missing sample_metadata column 'condition'",
    ):
        BatchCorrectionMetadataResolver().run(
            phospho=_phospho(),
            sample_metadata=_sample_metadata().drop(columns=["condition"]),
            batch_column="batch",
            condition_column="condition",
        )


def test_rejects_missing_sample_metadata_row() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="missing rows.*'sample_c'",
    ):
        BatchCorrectionMetadataResolver().run(
            phospho=_phospho(),
            sample_metadata=_sample_metadata(index=["sample_a", "sample_b"]),
            batch_column="batch",
            condition_column="condition",
        )


def test_rejects_extra_sample_metadata_row() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="rows not present in phospho columns.*'sample_extra'",
    ):
        BatchCorrectionMetadataResolver().run(
            phospho=_phospho(),
            sample_metadata=_sample_metadata(
                index=["sample_a", "sample_b", "sample_c", "sample_extra"]
            ),
            batch_column="batch",
            condition_column="condition",
        )


def test_rejects_duplicate_sample_metadata_row() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="sample_metadata.index contains duplicate sample labels.*'sample_b'",
    ):
        BatchCorrectionMetadataResolver().run(
            phospho=_phospho(),
            sample_metadata=_sample_metadata(
                index=["sample_a", "sample_b", "sample_b"]
            ),
            batch_column="batch",
            condition_column="condition",
        )


def test_rejects_missing_batch_label_value() -> None:
    sample_metadata = _sample_metadata()
    sample_metadata.loc["sample_b", "batch"] = pd.NA

    with pytest.raises(
        PhosPyInputError,
        match="contains missing batch labels.*'sample_b'",
    ):
        BatchCorrectionMetadataResolver().run(
            phospho=_phospho(),
            sample_metadata=sample_metadata,
            batch_column="batch",
            condition_column="condition",
        )


def test_dataset_builder_resolves_declared_batch_metadata_during_execution() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_adequate_phospho(),
            site_metadata=_adequate_site_metadata(),
            sample_metadata=_sample_metadata(
                index=["sample_d", "sample_b", "sample_a", "sample_c"]
            ),
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
    assert report.confounding_check_status == "passed"
    assert report.batch_levels == ("run_1", "run_2")
    assert report.condition_levels == ("control", "treated")
    assert "batch_correction" in set(
        built.preprocessing_report.operations.loc[:, "stage"].astype(str).tolist()
    )


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_b": [1.0, 2.0],
            "sample_a": [2.0, 3.0],
            "sample_c": [4.0, 5.0],
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


def _adequate_phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [2.0, 3.0],
            "sample_c": [4.0, 5.0],
            "sample_d": [5.0, 6.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )


def _adequate_site_metadata() -> pd.DataFrame:
    site_metadata = _site_metadata()
    site_metadata.index = _adequate_phospho().index.copy()
    return site_metadata


def _sample_metadata(index: list[str] | None = None) -> pd.DataFrame:
    source = {
        "sample_a": ("run_1", "control"),
        "sample_b": ("run_2", "treated"),
        "sample_c": ("run_1", "treated"),
        "sample_d": ("run_2", "control"),
        "sample_extra": ("run_9", "control"),
    }
    resolved_index = ["sample_a", "sample_b", "sample_c"] if index is None else index
    return pd.DataFrame(
        {
            "batch": [source[sample][0] for sample in resolved_index],
            "condition": [source[sample][1] for sample in resolved_index],
        },
        index=pd.Index(resolved_index, name="sample_id"),
    )
