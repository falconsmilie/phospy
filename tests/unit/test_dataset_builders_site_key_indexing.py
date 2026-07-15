from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import DatasetBuildRequest, Organism
from phospy.errors import DatasetValidationError
from phospy.errors.input import PhosPyInputError
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    decode_site_key,
    encode_site_key,
)


def _site_key(
    *,
    protein_identifier: str = "P28482",
    organism: str = "rat",
    protein_namespace: str = "protein_accession",
    residue: str = "Y",
    position: int = 182,
) -> str:
    return encode_site_key(
        build_protein_scoped_site_key(
            organism=organism,
            protein_namespace=protein_namespace,
            protein_identifier=protein_identifier,
            residue=residue,
            position=position,
            field_name="test.site_key",
            error_type=ValueError,
        )
    )


def _explicit_identity_frames(
    *,
    site_keys: tuple[str, ...] | None = None,
    protein_accessions: tuple[str, ...] = ("P28482",),
    display_ids: tuple[str, ...] | None = None,
    include_protein_context: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    resolved_display_ids = display_ids or ("MAPK14;Y182;",) * len(protein_accessions)
    resolved_site_keys = site_keys or tuple(
        _site_key(protein_identifier=protein_accession)
        for protein_accession in protein_accessions
    )
    row_count = len(resolved_display_ids)
    phospho = pd.DataFrame(
        {
            "sample_a": [float(position + 1) for position in range(row_count)],
            "sample_b": [float(position + 2) for position in range(row_count)],
        },
        index=pd.Index(list(resolved_display_ids), name="site_id"),
    )
    site_metadata_data: dict[str, object] = {
        "gene_symbol": ["MAPK14"] * len(resolved_display_ids),
        "site": ["Y182"] * len(resolved_display_ids),
        "organism": ["rat"] * len(resolved_display_ids),
        "display_id": list(resolved_display_ids),
        "site_key": list(resolved_site_keys),
        "site_sequence": ["AAAAAYAAAAA"] * len(resolved_display_ids),
        "localisation_confidence": [0.95] * len(resolved_display_ids),
    }
    if include_protein_context:
        site_metadata_data["protein_accession"] = list(protein_accessions)
    site_metadata = pd.DataFrame(
        site_metadata_data,
        index=phospho.index.copy(),
    )
    return phospho, site_metadata


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


def test_builder_rejects_legacy_display_indexed_input_without_protein_context() -> None:
    phospho, site_metadata = _legacy_display_indexed_frames(
        invalid_site_key_values=("MAPK14;Y182;",),
    )
    site_metadata = site_metadata.drop(columns=["protein_id"])

    with pytest.raises(
        PhosPyInputError,
        match="protein context is required to derive site_key",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                input_intensity_scale="linear",
            )
        )


def test_builder_accepts_explicit_site_key_when_it_matches_metadata() -> None:
    phospho, site_metadata = _explicit_identity_frames()

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    built_site_metadata = dataset.site_metadata
    expected_site_key = _site_key(protein_identifier="P28482")
    assert dataset.phospho.index.tolist() == [expected_site_key]
    assert built_site_metadata.index.tolist() == [expected_site_key]
    assert built_site_metadata.loc[:, "site_key"].tolist() == [expected_site_key]
    assert built_site_metadata.loc[:, "display_id"].tolist() == ["MAPK14;Y182;"]


def test_builder_rejects_human_rows_with_rat_dataset_request() -> None:
    phospho, site_metadata = _explicit_identity_frames()
    human_site_key = _site_key(organism="human", protein_identifier="P28482")
    site_metadata.loc[:, "organism"] = ["human"]
    site_metadata.loc[:, "site_key"] = [human_site_key]

    with pytest.raises(
        DatasetValidationError,
        match="dataset\\.organism must match every .* row_examples",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                input_intensity_scale="linear",
            )
        )


@pytest.mark.parametrize(
    "wrong_site_key",
    (
        pytest.param(_site_key(organism="human"), id="organism"),
        pytest.param(_site_key(protein_namespace="protein_id"), id="protein_namespace"),
        pytest.param(_site_key(protein_identifier="Q5S007"), id="protein_identifier"),
        pytest.param(_site_key(residue="T"), id="residue"),
        pytest.param(_site_key(position=183), id="position"),
    ),
)
def test_builder_rejects_explicit_site_key_that_disagrees_with_metadata(
    wrong_site_key: str,
) -> None:
    phospho, site_metadata = _explicit_identity_frames(site_keys=(wrong_site_key,))

    with pytest.raises(
        PhosPyInputError,
        match="explicit identity fields must match metadata-derived",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                input_intensity_scale="linear",
            )
        )


def test_builder_rejects_explicit_t309_site_key_when_metadata_says_t308() -> None:
    stale_site_key = _site_key(
        organism="human",
        protein_namespace="uniprot",
        protein_identifier="P31749",
        residue="T",
        position=309,
    )
    expected_site_key = _site_key(
        organism="human",
        protein_namespace="uniprot",
        protein_identifier="P31749",
        residue="T",
        position=308,
    )
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["AKT1"],
            "site": ["T308"],
            "organism": ["human"],
            "protein_namespace": ["uniprot"],
            "protein_identifier": ["P31749"],
            "display_id": ["AKT1;T308;"],
            "site_key": [stale_site_key],
            "site_sequence": ["AAAAATAAAAA"],
            "localisation_confidence": [0.95],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(PhosPyInputError) as exc_info:
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.HUMAN,
                input_intensity_scale="linear",
            )
        )

    message = str(exc_info.value)
    assert "explicit identity fields must match metadata-derived" in message
    assert "site_key" in message
    assert stale_site_key in message
    assert expected_site_key in message


def test_builder_rejects_explicit_site_key_without_protein_context() -> None:
    phospho, site_metadata = _explicit_identity_frames(include_protein_context=False)

    with pytest.raises(
        PhosPyInputError,
        match="protein context is required to derive site_key",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                input_intensity_scale="linear",
            )
        )


def test_builder_accepts_duplicate_explicit_display_id_with_distinct_site_key() -> None:
    phospho, site_metadata = _explicit_identity_frames(
        protein_accessions=("P28482", "Q5S007"),
    )
    expected_site_keys = site_metadata.loc[:, "site_key"].astype(str).tolist()

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    built_site_metadata = dataset.site_metadata
    assert dataset.phospho.shape[0] == 2
    assert dataset.phospho.index.astype(str).tolist() == expected_site_keys
    assert built_site_metadata.loc[:, "display_id"].tolist() == [
        "MAPK14;Y182;",
        "MAPK14;Y182;",
    ]
    assert built_site_metadata.loc[:, "site_key"].astype(str).tolist() == (
        expected_site_keys
    )
    assert built_site_metadata.loc[:, "site_key"].nunique() == 2
    assert int(built_site_metadata.loc[:, "display_id"].nunique()) == 1
    assert dataset.phospho.loc[expected_site_keys, "sample_a"].tolist() == [1.0, 2.0]


def test_builder_derives_distinct_site_keys_for_duplicate_display_ids() -> None:
    duplicate_display_id = "MAPK14;Y182;"
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.1, 2.1],
        },
        index=pd.Index([duplicate_display_id, duplicate_display_id], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "organism": ["rat", "rat"],
            "protein_namespace": ["protein_id", "protein_id"],
            "protein_identifier": ["MAPK14_A", "MAPK14_B"],
            "display_id": [duplicate_display_id, duplicate_display_id],
            "site_sequence": [
                "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
            ],
            "protein_id": ["MAPK14_A", "MAPK14_B"],
            "localisation_confidence": [0.95, 0.95],
        },
        index=phospho.index.copy(),
    )

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    built_site_metadata = dataset.site_metadata
    assert dataset.phospho.shape[0] == 2
    assert dataset.phospho.index.name == "site_key"
    assert dataset.phospho.index.is_unique
    assert built_site_metadata.index.equals(dataset.phospho.index)
    assert built_site_metadata.loc[:, "site_key"].astype(str).tolist() == (
        dataset.phospho.index.astype(str).tolist()
    )
    assert built_site_metadata.loc[:, "display_id"].tolist() == [
        duplicate_display_id,
        duplicate_display_id,
    ]
    assert int(built_site_metadata.loc[:, "display_id"].nunique()) == 1
    assert built_site_metadata.loc[:, "protein_identifier"].tolist() == [
        "MAPK14_A",
        "MAPK14_B",
    ]
    assert dataset.phospho.loc[:, "sample_a"].tolist() == [1.0, 2.0]


def test_builder_rejects_duplicate_explicit_site_key_values() -> None:
    duplicated_site_key = _site_key(protein_identifier="P28482")
    phospho, site_metadata = _explicit_identity_frames(
        protein_accessions=("P28482", "P28482"),
        site_keys=(duplicated_site_key, duplicated_site_key),
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


def test_builder_does_not_mutate_original_phospho_dataframe() -> None:
    phospho, site_metadata = _explicit_identity_frames()
    original_phospho = phospho.copy(deep=True)

    AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    pdt.assert_frame_equal(phospho, original_phospho)


def test_builder_does_not_mutate_original_site_metadata_dataframe() -> None:
    phospho, site_metadata = _explicit_identity_frames()
    original_site_metadata = site_metadata.copy(deep=True)

    AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

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
