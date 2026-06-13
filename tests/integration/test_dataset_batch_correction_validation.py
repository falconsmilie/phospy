from __future__ import annotations

import numpy as np
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

pytestmark = pytest.mark.integration


def test_dataset_batch_correction_valid_design_is_checked_before_noop_report() -> None:
    phospho = _phospho()

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho),
            sample_metadata=_valid_sample_metadata(phospho),
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                batch_correction=DatasetBatchCorrectionConfig(
                    method="linear_residualize_batch"
                )
            ),
        )
    )

    np.testing.assert_allclose(built.phospho.to_numpy(), phospho.to_numpy())
    assert built.preprocessing_report is not None
    report = built.preprocessing_report.batch_correction
    assert report is not None
    assert report.status == "rejected"
    assert report.confounding_check_status == "passed"
    assert report.batch_levels == ("run_1", "run_2")
    assert report.condition_levels == ("control", "treated")
    assert "batch_correction" not in set(
        built.preprocessing_report.operations.loc[:, "stage"].astype(str).tolist()
    )


def test_dataset_batch_correction_confounding_fails_clearly_before_execution() -> None:
    phospho = _phospho()

    with pytest.raises(
        PhosPyInputError,
        match="batch and condition are perfectly confounded",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho),
                sample_metadata=_confounded_sample_metadata(phospho),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                preprocessing_config=DatasetPreprocessingConfig(
                    batch_correction=DatasetBatchCorrectionConfig(
                        method="linear_residualize_batch"
                    )
                ),
            )
        )


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


def _valid_sample_metadata(phospho: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": ["run_1", "run_1", "run_2", "run_2"],
            "condition": ["control", "treated", "control", "treated"],
        },
        index=phospho.columns.copy(),
    )


def _confounded_sample_metadata(phospho: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": ["run_1", "run_1", "run_2", "run_2"],
            "condition": ["control", "control", "treated", "treated"],
        },
        index=phospho.columns.copy(),
    )
