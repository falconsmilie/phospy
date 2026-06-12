from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import Organism, PhosphositeImportRequest
from phospy.io.readers import MappedPhosphositeTableImporter

pytestmark = pytest.mark.integration


def test_phosphosite_importer_output_feeds_dataset_builder() -> None:
    source = pd.DataFrame(
        {
            "gene": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "protein": ["MAPK14", "GSK3B"],
            "sequence_window": [
                ("A" * 15) + "Y" + ("A" * 15),
                ("A" * 15) + "S" + ("A" * 15),
            ],
            "localisation": [0.95, 0.92],
            "intensity_a": [10.0, 20.0],
            "intensity_b": [11.0, 21.0],
        }
    )

    imported = MappedPhosphositeTableImporter().run(
        PhosphositeImportRequest(
            source=source,
            sample_intensity_columns={
                "intensity_a": "sample_a",
                "intensity_b": "sample_b",
            },
            gene_symbol_column="gene",
            site_column="site",
            protein_id_column="protein",
            site_sequence_column="sequence_window",
            localisation_confidence_column="localisation",
            localisation_confidence_scale="probability",
        )
    )

    dataset = AnalysisReadyDatasetBuilder().run(
        imported.to_dataset_build_request(
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    assert dataset.phospho.shape == (2, 2)
    assert dataset.site_metadata.loc[:, "display_id"].tolist() == [
        "MAPK14;Y182;",
        "GSK3B;S9;",
    ]
    assert dataset.site_metadata.index.name == "site_key"
    assert dataset.site_metadata.loc[:, "site_key"].tolist() == (
        dataset.phospho.index.tolist()
    )
