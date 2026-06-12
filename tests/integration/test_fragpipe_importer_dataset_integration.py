from __future__ import annotations

from pathlib import Path

import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import Organism
from phospy.api.configs import (
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
)
from phospy.api.requests import (
    DATASET_MULTI_SITE_POLICY_SPLIT,
    DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
)
from phospy.io.readers import (
    FragPipePTMProphetImporter,
    FragPipePTMProphetImportRequest,
)

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests" / "fixtures" / "fragpipe"

pytestmark = pytest.mark.integration


def test_fragpipe_ptmprophet_importer_output_feeds_dataset_builder_split_policy() -> (
    None
):
    imported = FragPipePTMProphetImporter().run(
        FragPipePTMProphetImportRequest(
            source=FIXTURES / "ptmprophet_sites.tsv",
        )
    )

    dataset = AnalysisReadyDatasetBuilder().run(
        imported.to_dataset_build_request(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
            organism=Organism.HUMAN,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                localisation=DatasetLocalisationConfig(min_confidence=0.5)
            ),
        )
    )

    display_ids = set(dataset.site_metadata.loc[:, "display_id"].tolist())
    assert display_ids == {
        "MAPK1;S10;",
        "AKT1;S473;",
        "AKT1;T475;",
        "GSK3B;S10;",
        "GSK3B;T11;",
    }
    assert dataset.phospho.shape == (5, 2)
    assert dataset.site_metadata.loc[
        :,
        "localisation_confidence",
    ].min() == pytest.approx(0.5)
    assert dataset.preprocessing_report is not None
    resolution_rows = dataset.preprocessing_report.operations.loc[
        dataset.preprocessing_report.operations.loc[:, "stage"]
        == "peptide_evidence_resolution"
    ]
    assert int(resolution_rows.shape[0]) == 1
    parameters = resolution_rows.iloc[0]["parameters"]
    assert parameters["input_mode"] == DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE
    assert parameters["multi_site_policy"] == DATASET_MULTI_SITE_POLICY_SPLIT
    assert int(parameters["split_observations"]) == 2
