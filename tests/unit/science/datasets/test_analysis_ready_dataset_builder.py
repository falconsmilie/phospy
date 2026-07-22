from __future__ import annotations

from collections.abc import Mapping

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import DatasetBuildRequest, Organism
from phospy.errors import DatasetValidationError, PhosPyInputError

_CENTRED_Y_SEQUENCE = "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA"
_CENTRED_T_SEQUENCE = "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA"


def _valid_builder_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    site_ids = pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id")
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=site_ids.copy(),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "protein_id": ["MAPK14", "AKT1"],
            "site_sequence": [_CENTRED_Y_SEQUENCE, _CENTRED_T_SEQUENCE],
            "localisation_confidence": [0.95, 0.96],
        },
        index=site_ids.copy(),
    )
    sample_metadata = pd.DataFrame(
        {"condition": ["control", "treated"]},
        index=pd.Index(["sample_a", "sample_b"], name="sample_id"),
    )
    return phospho, site_metadata, sample_metadata


def test_builder_accepts_valid_phospho_site_and_sample_metadata() -> None:
    phospho, site_metadata, sample_metadata = _valid_builder_inputs()

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    assert dataset.phospho.index.name == "site_key"
    assert dataset.site_metadata.index.equals(dataset.phospho.index)
    assert dataset.site_metadata.loc[:, "site_key"].tolist() == (
        dataset.phospho.index.tolist()
    )
    assert dataset.sample_metadata is not None
    assert dataset.sample_metadata.index.tolist() == ["sample_a", "sample_b"]


def test_builder_rejects_missing_site_sequence_when_not_established() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [1.5]},
        index=pd.Index(["FAKE1;S1;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["FAKE1"],
            "site": ["S1"],
            "protein_id": ["FAKE1"],
            "localisation_confidence": [0.95],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(
        PhosPyInputError,
        match="cannot construct AnalysisReadyPhosphoDataset.*site_sequence",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                input_intensity_scale="linear",
            )
        )


def test_builder_rejects_non_site_key_indexed_matrix_that_cannot_align_to_metadata() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=pd.Index(["row_a", "row_b"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "protein_id": ["MAPK14"],
            "site_sequence": [_CENTRED_Y_SEQUENCE],
            "localisation_confidence": [0.95],
        },
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )

    with pytest.raises(
        PhosPyInputError,
        match="phospho and site_metadata row counts must match",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                input_intensity_scale="linear",
            )
        )


def test_builder_rejects_mismatched_sample_labels() -> None:
    phospho, site_metadata, _sample_metadata = _valid_builder_inputs()
    sample_metadata = pd.DataFrame(
        {"condition": ["control", "treated"]},
        index=pd.Index(["sample_a", "sample_c"], name="sample_id"),
    )

    with pytest.raises(
        DatasetValidationError,
        match="dataset.sample_metadata.index must exactly match",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                organism=Organism.RAT,
                input_intensity_scale="linear",
            )
        )


def test_builder_records_construction_provenance() -> None:
    phospho, site_metadata, sample_metadata = _valid_builder_inputs()

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    provenance = dataset.provenance
    assert provenance is not None
    construction = provenance.workflow_parameters.get("construction")
    assert isinstance(construction, Mapping)
    assert construction["method"] == "AnalysisReadyDatasetBuilder.run"
    assert construction["dataset_type"] == "AnalysisReadyPhosphoDataset"
    assert (
        construction["model_constructor"]
        == "AnalysisReadyPhosphoDataset._from_builder_output"
    )

    input_identities = construction["input_table_identities"]
    assert isinstance(input_identities, Mapping)
    for table_name in (
        "dataset.phospho",
        "dataset.site_metadata",
        "dataset.sample_metadata",
    ):
        identity = input_identities[table_name]
        assert isinstance(identity, Mapping)
        assert identity["rows"] > 0
        assert identity["columns"] > 0
        assert isinstance(identity["exact_hash_value"], str)
        assert isinstance(identity["tolerance_hash_value"], str)

    establishment = construction["processing_state_establishment"]
    assert isinstance(establishment, Mapping)
    assert str(establishment["source"]).endswith(
        "DatasetTransformationStateResolver.run"
    )
    assert establishment["analysis_ready_boundary"] == (
        "AnalysisReadyPhosphoDataset._from_builder_output"
    )
    assert isinstance(establishment["intensity_scale_establishment"], Mapping)

    input_names = {fingerprint.name for fingerprint in provenance.input_tables}
    assert {
        "dataset.phospho",
        "dataset.site_metadata",
        "dataset.sample_metadata",
    }.issubset(input_names)
