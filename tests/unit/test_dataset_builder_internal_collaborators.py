from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import DatasetBuildRequest
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.builders.contracts import (
    InterpretedDatasetBuildRequest,
    PreprocessedDatasetBuildTables,
)
from phospy.science.datasets.builders.interpreter_collaborators import (
    DatasetBuildSourceResolver,
)
from phospy.science.datasets.builders.normalizer import NormalizedDatasetInputs
from phospy.science.datasets.builders.provenance_assembler import (
    DatasetRunProvenanceAssembler,
)
from phospy.science.datasets.builders.site_identity_derivation import (
    DatasetSiteIdentityDeriver,
)
from phospy.science.datasets.builders.site_sequence_boundary import (
    AnalysisReadySiteSequenceValidator,
)
from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from phospy.science.references.models import Organism
from phospy.science.sites.site_keys import decode_site_key


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )


def _site_metadata(include_sequence: bool = True) -> pd.DataFrame:
    data: dict[str, object] = {
        "gene_symbol": ["MAPK14"],
        "site": ["Y182"],
        "protein_accession": ["P28482"],
    }
    if include_sequence:
        data["site_sequence"] = ["AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA"]
    return pd.DataFrame(data, index=_phospho().index.copy())


def test_source_resolver_reads_and_normalizes_site_level_inputs() -> None:
    reader_calls: list[str] = []
    normalizer_calls: list[tuple[pd.DataFrame, pd.DataFrame]] = []
    phospho = _phospho()
    site_metadata = _site_metadata()
    sample_metadata = pd.DataFrame(
        {"group": ["A"]},
        index=pd.Index(["sample_a"], name="sample"),
    )

    class ReaderSpy:
        def run(self, value: object, *, field_name: str) -> pd.DataFrame:
            reader_calls.append(field_name)
            assert isinstance(value, pd.DataFrame)
            return value

    class NormalizerSpy:
        def run(
            self,
            *,
            phospho: pd.DataFrame,
            site_metadata: pd.DataFrame,
            sample_metadata: pd.DataFrame | None,
            total: pd.DataFrame | None,
        ) -> NormalizedDatasetInputs:
            assert sample_metadata is not None
            assert total is None
            normalizer_calls.append((phospho, site_metadata))
            return NormalizedDatasetInputs(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                total=None,
                site_identifier_normalisation=None,
            )

    resolved = DatasetBuildSourceResolver(
        reader=ReaderSpy(),
        normalizer=NormalizerSpy(),
    ).run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            organism=Organism.HUMAN,
            input_intensity_scale="linear",
        )
    )

    assert reader_calls == ["sample_metadata", "phospho", "site_metadata"]
    assert len(normalizer_calls) == 1
    assert resolved.site_resolution_mode == "site_level_resolved"
    assert resolved.multi_site_policy is None
    assert resolved.peptide_evidence_resolution is None
    assert resolved.site_metadata.loc[:, "display_id"].tolist() == ["MAPK14;Y182;"]
    assert resolved.site_metadata.loc[:, "site_key"].shape[0] == 1


def test_source_resolver_derives_site_identity_for_display_indexed_rows_with_context() -> (
    None
):
    phospho = _phospho()
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "protein_id": ["P28482"],
            "site_sequence": ["AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA"],
        },
        index=phospho.index.copy(),
    )

    resolved = DatasetBuildSourceResolver().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.HUMAN,
            input_intensity_scale="linear",
        )
    )

    assert resolved.site_metadata.loc[:, "display_id"].tolist() == ["MAPK14;Y182;"]
    key = decode_site_key(
        resolved.site_metadata.loc[:, "site_key"].iloc[0],
        field_name="test.source_resolver.site_key",
        error_type=ValueError,
    )
    assert key.organism == "human"
    assert key.protein_namespace == "protein_id"
    assert key.protein_identifier == "P28482"
    assert key.residue == "Y"
    assert key.position == 182


def test_site_sequence_boundary_validator_rejects_missing_sequence_column() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="site_sequence is missing, blank, or invalid after builder enrichment",
    ):
        AnalysisReadySiteSequenceValidator().run(
            site_metadata=_site_metadata(include_sequence=False),
            preprocessing_trace=None,
        )


def test_run_provenance_assembler_records_opaque_token_mode() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()
    plan = PreprocessingPlan.default()
    request = InterpretedDatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        organism=None,
        preprocessing_plan=plan,
        allow_opaque_site_values=True,
    )
    preprocessed = PreprocessedDatasetBuildTables(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        preprocessing_trace=None,
    )

    provenance = DatasetRunProvenanceAssembler().run(
        request=request,
        preprocessed=preprocessed,
        validated_site_metadata=site_metadata,
        resolved_phospho=phospho,
        resolved_total=None,
        preprocessing_trace=None,
        intensity_scale_label="linear",
        intensity_scale_establishment={"source": "test"},
        quantitative_meaning="phosphosite_abundance",
        allow_opaque_site_values=True,
    )

    assert provenance.workflow_parameters["site_token_validation"] == {
        "mode": "opaque_opt_in"
    }


def _duplicate_identity_frames(
    *,
    site_rows: tuple[str, str] = ("Y182", "Y182"),
    protein_id: tuple[object, object] | None = None,
    protein_accession: tuple[object, object] | None = None,
    isoform_id: tuple[object, object] | None = None,
    include_context_columns: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=pd.Index(["row_a", "row_b"], name="source_row"),
    )
    site_metadata_data: dict[str, object] = {
        "gene_symbol": ["MAPK14", "MAPK14"],
        "site": list(site_rows),
        "site_sequence": ["AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA"] * 2,
    }
    if include_context_columns and protein_id is not None:
        site_metadata_data["protein_id"] = list(protein_id)
    if include_context_columns and protein_accession is not None:
        site_metadata_data["protein_accession"] = list(protein_accession)
    if include_context_columns and isoform_id is not None:
        site_metadata_data["isoform_id"] = list(isoform_id)
    site_metadata = pd.DataFrame(site_metadata_data, index=phospho.index.copy())
    return phospho, site_metadata


def test_source_resolver_rejects_duplicate_site_key_from_same_protein_context() -> None:
    phospho, site_metadata = _duplicate_identity_frames(
        protein_id=("P28482", "P28482"),
        protein_accession=("P28482-1", "P28482-1"),
    )

    with pytest.raises(
        PhosPyInputError,
        match="site_key must be unique",
    ) as exc_info:
        DatasetBuildSourceResolver().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.HUMAN,
                input_intensity_scale="linear",
            )
        )
    message = str(exc_info.value)
    assert "site_matrix.policy='as_input'" in message


def test_source_resolver_allows_duplicate_display_ids_with_distinct_protein_id() -> (
    None
):
    phospho, site_metadata = _duplicate_identity_frames(
        protein_id=("P28482", "Q5S007"),
    )

    resolved = DatasetBuildSourceResolver().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.HUMAN,
            input_intensity_scale="linear",
        )
    )

    assert resolved.site_metadata.loc[:, "display_id"].tolist() == [
        "MAPK14;Y182;",
        "MAPK14;Y182;",
    ]
    assert resolved.site_metadata.loc[:, "site_key"].nunique() == 2


def test_source_resolver_allows_duplicate_display_ids_with_distinct_accession() -> None:
    phospho, site_metadata = _duplicate_identity_frames(
        protein_id=("P28482", "P28482"),
        protein_accession=("P28482-1", "P28482-2"),
    )

    resolved = DatasetBuildSourceResolver().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.HUMAN,
            input_intensity_scale="linear",
        )
    )

    assert resolved.site_metadata.loc[:, "display_id"].nunique() == 1
    assert resolved.site_metadata.loc[:, "site_key"].nunique() == 2


def test_source_resolver_allows_duplicate_display_ids_with_distinct_isoform_id() -> (
    None
):
    phospho, site_metadata = _duplicate_identity_frames(
        protein_id=("P28482", "P28482"),
        protein_accession=("P28482-1", "P28482-1"),
        isoform_id=("isoform-1", "isoform-2"),
    )

    resolved = DatasetBuildSourceResolver().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.HUMAN,
            input_intensity_scale="linear",
        )
    )

    assert resolved.site_metadata.loc[:, "display_id"].nunique() == 1
    assert resolved.site_metadata.loc[:, "site_key"].nunique() == 2


def test_source_resolver_rejects_display_indexed_rows_when_protein_context_is_missing() -> (
    None
):
    phospho, site_metadata = _duplicate_identity_frames(include_context_columns=False)

    with pytest.raises(
        PhosPyInputError,
        match="protein context is required to derive site_key",
    ) as exc_info:
        DatasetBuildSourceResolver().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.HUMAN,
                input_intensity_scale="linear",
            )
        )
    assert "gene_symbol is display metadata and is not protein context" in str(
        exc_info.value
    )


def test_source_resolver_allows_distinct_display_ids_with_different_context() -> None:
    phospho, site_metadata = _duplicate_identity_frames(
        site_rows=("Y182", "T308"),
        protein_id=("P28482", "Q5S007"),
        protein_accession=("P28482-1", "Q5S007-1"),
    )

    resolved = DatasetBuildSourceResolver().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.HUMAN,
            input_intensity_scale="linear",
        )
    )
    assert resolved.site_metadata.loc[:, "site"].tolist() == ["Y182", "T308"]


def test_site_identity_deriver_derives_display_id_and_site_key_without_mutation() -> (
    None
):
    input_site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["mapk14"],
            "site": ["y182"],
            "protein_accession": ["P28482"],
            "isoform_id": [None],
        },
        index=pd.Index(["row_a"], name="source_row"),
    )
    original_columns = input_site_metadata.columns.tolist()
    original = input_site_metadata.copy(deep=True)

    derived = DatasetSiteIdentityDeriver().run(
        site_metadata=input_site_metadata,
        organism=Organism.HUMAN,
    )

    assert input_site_metadata.equals(original)
    assert input_site_metadata.columns.tolist() == original_columns
    assert derived.loc["row_a", "display_id"] == "MAPK14;Y182;"
    key = decode_site_key(
        derived.loc["row_a", "site_key"],
        field_name="test.site_key",
        error_type=ValueError,
    )
    assert key.organism == "human"
    assert key.protein_namespace == "protein_accession"
    assert key.protein_identifier == "P28482"
    assert key.residue == "Y"
    assert key.position == 182
    assert key.isoform_id is None


def test_site_identity_deriver_uses_accession_over_protein_id_and_falls_back_namespace() -> (
    None
):
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["AKT1", "AKT1"],
            "site": ["S473", "S473"],
            "organism": ["human", "human"],
            "protein_accession": ["P31749", ""],
            "protein_id": ["AKT1_HUMAN", "AKT1_HUMAN"],
            "protein_namespace": ["", ""],
            "isoform_id": [None, "2"],
        },
        index=pd.Index(["row_a", "row_b"], name="source_row"),
    )

    derived = DatasetSiteIdentityDeriver().run(
        site_metadata=site_metadata, organism=None
    )

    key_a = decode_site_key(
        derived.loc["row_a", "site_key"],
        field_name="test.site_key.row_a",
        error_type=ValueError,
    )
    key_b = decode_site_key(
        derived.loc["row_b", "site_key"],
        field_name="test.site_key.row_b",
        error_type=ValueError,
    )
    assert key_a.protein_namespace == "protein_accession"
    assert key_a.protein_identifier == "P31749"
    assert key_b.protein_namespace == "protein_id"
    assert key_b.protein_identifier == "AKT1_HUMAN"
    assert key_a.isoform_id is None
    assert key_b.isoform_id == "2"


def test_site_identity_deriver_same_display_id_different_protein_context_changes_site_key() -> (
    None
):
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "organism": ["human", "human"],
            "protein_accession": ["P28482", "Q9Y2J8"],
        },
        index=pd.Index(["row_a", "row_b"], name="source_row"),
    )

    derived = DatasetSiteIdentityDeriver().run(
        site_metadata=site_metadata, organism=None
    )

    assert derived.loc["row_a", "display_id"] == "MAPK14;Y182;"
    assert derived.loc["row_b", "display_id"] == "MAPK14;Y182;"
    assert derived.loc["row_a", "site_key"] != derived.loc["row_b", "site_key"]


def test_site_identity_deriver_same_protein_context_produces_same_site_key() -> None:
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "organism": ["human", "human"],
            "protein_accession": ["P28482", "P28482"],
            "isoform_id": [None, None],
        },
        index=pd.Index(["row_a", "row_b"], name="source_row"),
    )

    derived = DatasetSiteIdentityDeriver().run(
        site_metadata=site_metadata, organism=None
    )
    assert derived.loc["row_a", "site_key"] == derived.loc["row_b", "site_key"]


def test_site_identity_deriver_rejects_missing_organism() -> None:
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "protein_accession": ["P28482"],
        },
        index=pd.Index(["row_a"], name="source_row"),
    )

    with pytest.raises(PhosPyInputError, match="organism is required"):
        DatasetSiteIdentityDeriver().run(site_metadata=site_metadata, organism=None)


def test_site_identity_deriver_rejects_missing_protein_context() -> None:
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "organism": ["human"],
            "protein_accession": [""],
            "protein_id": [None],
        },
        index=pd.Index(["row_a"], name="source_row"),
    )

    with pytest.raises(
        PhosPyInputError, match="protein context is required to derive site_key"
    ):
        DatasetSiteIdentityDeriver().run(site_metadata=site_metadata, organism=None)


def test_site_identity_deriver_rejects_gene_symbol_only_protein_context() -> None:
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "organism": ["human"],
        },
        index=pd.Index(["row_a"], name="source_row"),
    )

    with pytest.raises(PhosPyInputError, match="gene_symbol is display metadata"):
        DatasetSiteIdentityDeriver().run(site_metadata=site_metadata, organism=None)


@pytest.mark.parametrize(
    "site_token",
    [
        pytest.param("A123", id="invalid_residue"),
        pytest.param("S0", id="invalid_position"),
        pytest.param("Y182-T180", id="non_strict_token"),
    ],
)
def test_site_identity_deriver_rejects_invalid_site_token(site_token: str) -> None:
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": [site_token],
            "organism": ["human"],
            "protein_accession": ["P28482"],
        },
        index=pd.Index(["row_a"], name="source_row"),
    )

    with pytest.raises(PhosPyInputError, match="strict residue/position"):
        DatasetSiteIdentityDeriver().run(site_metadata=site_metadata, organism=None)


def test_dataset_builder_does_not_mutate_caller_owned_input_frames() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "protein_id": ["P28482"],
            "site_sequence": ["AAAAAYAAAAA"],
            "localisation_confidence": [0.95],
        },
        index=phospho.index.copy(),
    )
    original_phospho = phospho.copy(deep=True)
    original_site_metadata = site_metadata.copy(deep=True)

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.HUMAN,
            input_intensity_scale="linear",
        )
    )

    pd.testing.assert_frame_equal(phospho, original_phospho)
    pd.testing.assert_frame_equal(site_metadata, original_site_metadata)
    assert built.site_metadata.loc[:, "display_id"].tolist() == ["MAPK14;Y182;"]
    assert built.site_metadata.loc[:, "site_key"].shape[0] == 1
