from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import DatasetBuildRequest, Organism
from phospy.errors.input import PhosPyInputError
from phospy.science.sites.site_keys import decode_site_key


def _legacy_display_indexed_frames(
    *,
    protein_ids: tuple[str, ...] = ("P28482",),
    invalid_site_key_values: tuple[str, ...] | None = ("MAPK14;Y182;",),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    display_ids = ["MAPK14;Y182;"] * len(protein_ids)
    phospho = pd.DataFrame(
        {
            "sample_a": [float(position + 1) for position in range(len(protein_ids))],
            "sample_b": [float(position + 2) for position in range(len(protein_ids))],
        },
        index=pd.Index(display_ids, name="site_id"),
    )
    site_metadata_data: dict[str, object] = {
        "gene_symbol": ["MAPK14"] * len(protein_ids),
        "site": ["Y182"] * len(protein_ids),
        "protein_id": list(protein_ids),
        "site_sequence": ["AAAAAYAAAAA"] * len(protein_ids),
        "localisation_confidence": [0.95] * len(protein_ids),
    }
    if invalid_site_key_values is not None:
        site_metadata_data["site_key"] = list(invalid_site_key_values)
    site_metadata = pd.DataFrame(
        site_metadata_data,
        index=phospho.index.copy(),
    )
    return phospho, site_metadata


def test_builder_reindexes_legacy_display_indexed_input_to_encoded_site_key() -> None:
    phospho, site_metadata = _legacy_display_indexed_frames()
    original_phospho = phospho.copy(deep=True)
    original_site_metadata = site_metadata.copy(deep=True)

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    built_phospho = dataset.phospho
    built_site_metadata = dataset.site_metadata
    assert built_phospho.index.name == "site_key"
    assert built_site_metadata.index.name == "site_key"
    assert built_phospho.index.equals(built_site_metadata.index)
    assert built_site_metadata.loc[:, "site_key"].tolist() == (
        built_site_metadata.index.tolist()
    )
    assert "MAPK14;Y182;" not in built_phospho.index.tolist()
    assert built_site_metadata.loc[:, "display_id"].tolist() == ["MAPK14;Y182;"]

    decoded = decode_site_key(
        built_phospho.index[0],
        field_name="test.dataset.phospho.index[0]",
        error_type=ValueError,
    )
    assert decoded.organism == "rat"
    assert decoded.protein_namespace == "protein_id"
    assert decoded.protein_identifier == "P28482"
    assert decoded.residue == "Y"
    assert decoded.position == 182

    pdt.assert_frame_equal(phospho, original_phospho)
    pdt.assert_frame_equal(site_metadata, original_site_metadata)


def test_builder_keeps_duplicate_display_ids_separate_by_encoded_site_key() -> None:
    phospho, site_metadata = _legacy_display_indexed_frames(
        protein_ids=("P28482", "Q5S007"),
        invalid_site_key_values=("MAPK14;Y182;", "MAPK14;Y182;"),
    )

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    built_phospho = dataset.phospho
    built_site_metadata = dataset.site_metadata
    assert built_phospho.shape[0] == 2
    assert built_phospho.index.name == "site_key"
    assert built_phospho.index.is_unique
    assert built_phospho.index.equals(built_site_metadata.index)
    assert built_site_metadata.loc[:, "site_key"].tolist() == (
        built_site_metadata.index.tolist()
    )
    assert built_site_metadata.loc[:, "display_id"].tolist() == [
        "MAPK14;Y182;",
        "MAPK14;Y182;",
    ]
    assert built_site_metadata.loc[:, "display_id"].nunique() == 1
    assert built_site_metadata.loc[:, "site_key"].nunique() == 2


def test_builder_rejects_duplicate_derived_site_key_values_for_as_input() -> None:
    phospho, site_metadata = _legacy_display_indexed_frames(
        protein_ids=("P28482", "P28482"),
        invalid_site_key_values=("MAPK14;Y182;", "MAPK14;Y182;"),
    )

    with pytest.raises(PhosPyInputError, match="site_key must be unique"):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                input_intensity_scale="linear",
            )
        )
