from __future__ import annotations

from pathlib import Path

import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import Organism
from phospy.api.requests import (
    DATASET_MULTI_SITE_POLICY_SPLIT,
    DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
)
from phospy.io.readers import (
    MaxQuantPhosphositeImporter,
    MaxQuantPhosphositeImportRequest,
)

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests" / "fixtures" / "maxquant"

pytestmark = pytest.mark.integration


def test_maxquant_importer_output_feeds_site_level_dataset_builder() -> None:
    imported = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=FIXTURES / "phospho_sty_sites_standard.txt"
        )
    )

    dataset = AnalysisReadyDatasetBuilder().run(
        imported.to_dataset_build_request(
            organism=Organism.HUMAN,
            input_intensity_scale="linear",
        )
    )

    assert dataset.phospho.shape == (2, 2)
    assert dataset.site_metadata.loc[:, "display_id"].tolist() == [
        "MAPK1;S10;",
        "AKT1;S473;",
    ]
    assert dataset.site_metadata.loc[:, "protein_id"].tolist() == ["P28482", "P31749"]
    assert dataset.site_metadata.index.name == "site_key"
    assert dataset.site_metadata.index.astype(str).str.startswith("phospy:v1").all()


def test_maxquant_multisite_peptide_evidence_feeds_dataset_builder_split_policy() -> (
    None
):
    imported = MaxQuantPhosphositeImporter().run(
        MaxQuantPhosphositeImportRequest(
            source=FIXTURES / "phospho_sty_sites_multisite.txt"
        )
    )

    dataset = AnalysisReadyDatasetBuilder().run(
        imported.to_dataset_build_request(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
            organism=Organism.HUMAN,
            input_intensity_scale="linear",
        )
    )

    display_ids = set(dataset.site_metadata.loc[:, "display_id"].tolist())
    assert display_ids == {"MAPK1;S10;", "MAPK1;T12;", "AKT1;S473;"}
    assert dataset.phospho.shape == (3, 2)
    assert dataset.preprocessing_report is not None
    resolution_rows = dataset.preprocessing_report.operations.loc[
        dataset.preprocessing_report.operations.loc[:, "stage"]
        == "peptide_evidence_resolution"
    ]
    assert int(resolution_rows.shape[0]) == 1
    parameters = resolution_rows.iloc[0]["parameters"]
    assert parameters["input_mode"] == DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE
    assert parameters["multi_site_policy"] == DATASET_MULTI_SITE_POLICY_SPLIT
    assert int(parameters["split_observations"]) == 1
