from __future__ import annotations

import pandas as pd
import pytest

import phospy.io as phospy_io
from phospy.api.requests import DatasetBuildRequest
from phospy.datasets.builders.public import AnalysisReadyDatasetBuilder
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors import (
    DatasetValidationError,
    ReferenceCompatibilityError,
    ReferenceResolutionError,
    ReferenceValidationError,
    TransformationValidationError,
    UnsupportedInputFormatError,
    UnsupportedOrganismError,
)
from phospy.references.models import Organism, ReferenceBundle, ReferencePreset
from phospy.references.resolution import ReferenceResolver
from phospy.transformations.models import (
    MatrixTransformationState,
    TransformationState,
)


def _phospho() -> pd.DataFrame:
    return pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"])


def _site_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
        },
        index=["MAPK14;Y182;"],
    )


def _valid_bundle(organism: Organism = Organism.RAT) -> ReferenceBundle:
    return ReferenceBundle(
        organism=organism,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["MAP2K6"], "substrate_site": ["MAPK14;Y182;"]}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
    )


def test_dataset_requires_site_sequence_column() -> None:
    bad_site_metadata = _site_metadata().drop(columns=["site_sequence"])
    with pytest.raises(DatasetValidationError):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=bad_site_metadata,
            organism=Organism.RAT,
            transformation_state=TransformationState.raw(has_total_matrix=False),
        )


def test_dataset_rejects_blank_gene_symbol_values() -> None:
    bad_site_metadata = _site_metadata().copy(deep=True)
    bad_site_metadata.loc[:, "gene_symbol"] = ["   "]
    with pytest.raises(
        DatasetValidationError,
        match="dataset.site_metadata.gene_symbol must contain non-empty string values",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=bad_site_metadata,
            organism=Organism.RAT,
            transformation_state=TransformationState.raw(has_total_matrix=False),
        )


def test_dataset_rejects_blank_site_values() -> None:
    bad_site_metadata = _site_metadata().copy(deep=True)
    bad_site_metadata.loc[:, "site"] = [""]
    with pytest.raises(
        DatasetValidationError,
        match="dataset.site_metadata.site must contain non-empty string values",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=bad_site_metadata,
            organism=Organism.RAT,
            transformation_state=TransformationState.raw(has_total_matrix=False),
        )


def test_dataset_rejects_blank_site_sequence_values() -> None:
    bad_site_metadata = _site_metadata().copy(deep=True)
    bad_site_metadata.loc[:, "site_sequence"] = [" \t "]
    with pytest.raises(
        DatasetValidationError,
        match=(
            "dataset.site_metadata.site_sequence must contain non-empty string values"
        ),
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=bad_site_metadata,
            organism=Organism.RAT,
            transformation_state=TransformationState.raw(has_total_matrix=False),
        )


def test_dataset_rejects_empty_phospho_matrix() -> None:
    with pytest.raises(DatasetValidationError):
        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame(),
            site_metadata=pd.DataFrame(),
            organism=Organism.RAT,
            transformation_state=TransformationState.raw(has_total_matrix=False),
        )


def test_dataset_rejects_nan_in_phospho_matrix() -> None:
    phospho = _phospho()
    phospho.loc["MAPK14;Y182;", "sample_a"] = float("nan")
    with pytest.raises(DatasetValidationError, match="must not contain missing values"):
        AnalysisReadyPhosphoDataset(
            phospho=phospho,
            site_metadata=_site_metadata(),
            organism=Organism.RAT,
            transformation_state=TransformationState.raw(has_total_matrix=False),
        )


def test_dataset_rejects_inf_in_phospho_matrix() -> None:
    phospho = _phospho()
    phospho.loc["MAPK14;Y182;", "sample_a"] = float("inf")
    with pytest.raises(
        DatasetValidationError, match="must contain finite numeric values"
    ):
        AnalysisReadyPhosphoDataset(
            phospho=phospho,
            site_metadata=_site_metadata(),
            organism=Organism.RAT,
            transformation_state=TransformationState.raw(has_total_matrix=False),
        )


def test_builder_rejects_sparse_missingness_in_phospho_matrix() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan")],
            "sample_b": [2.0, 3.0],
        },
        index=["MAPK14;Y182;", "AKT1;T308;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "RPHFPQFSYSASGTA",
            ],
        },
        index=phospho.index,
    )
    with pytest.raises(DatasetValidationError, match="must not contain missing values"):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
            )
        )


def test_dataset_validates_transformation_state_for_total_matrix() -> None:
    with pytest.raises(TransformationValidationError):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            total=pd.DataFrame({"sample_a": [2.0]}, index=["GENEA"]),
            organism=Organism.RAT,
            transformation_state=TransformationState.raw(has_total_matrix=False),
        )


def test_dataset_rejects_mixed_transformation_kind_between_phospho_and_total() -> None:
    with pytest.raises(TransformationValidationError):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            total=pd.DataFrame({"sample_a": [2.0]}, index=["GENEA"]),
            organism=Organism.RAT,
            transformation_state=TransformationState(
                phospho=MatrixTransformationState.linear(),
                total=MatrixTransformationState.log2(),
            ),
        )


def test_dataset_rejects_inf_in_total_matrix() -> None:
    with pytest.raises(
        DatasetValidationError, match="dataset.total must contain finite numeric values"
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            total=pd.DataFrame({"sample_a": [float("-inf")]}, index=["MAPK14"]),
            organism=Organism.RAT,
            transformation_state=TransformationState.raw(has_total_matrix=True),
        )


def test_dataset_rejects_total_matrix_sample_mismatch() -> None:
    with pytest.raises(
        DatasetValidationError,
        match="dataset.total.columns must exactly match dataset.phospho.columns",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            total=pd.DataFrame({"sample_b": [2.0]}, index=["MAPK14"]),
            organism=Organism.RAT,
            transformation_state=TransformationState.raw(has_total_matrix=True),
        )


def test_reference_bundle_requires_non_empty_resources() -> None:
    with pytest.raises(ReferenceValidationError):
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["AAAAAAA"]},
                index=pd.Index(["GENEA;S1;"]),
            ),
        )


def test_reference_resolver_auto_uses_dataset_organism() -> None:
    resolved = ReferenceResolver().run(
        ReferencePreset.AUTO,
        dataset_organism=Organism.RAT,
    )
    assert resolved.organism is Organism.RAT
    assert not resolved.kinase_substrate_map.empty
    assert not resolved.site_sequences.empty


def test_reference_resolver_rejects_preset_dataset_mismatch() -> None:
    with pytest.raises(ReferenceCompatibilityError):
        ReferenceResolver().run(
            ReferencePreset.HUMAN,
            dataset_organism=Organism.RAT,
        )


def test_reference_resolver_rejects_unsupported_human_bundled_preset() -> None:
    with pytest.raises(UnsupportedOrganismError):
        ReferenceResolver().run(
            ReferencePreset.HUMAN,
            dataset_organism=Organism.HUMAN,
        )


def test_reference_resolver_auto_requires_dataset_organism() -> None:
    with pytest.raises(ReferenceResolutionError):
        ReferenceResolver().run(
            ReferencePreset.AUTO,
            dataset_organism=None,
        )


def test_builder_rejects_non_dataframe_phospho_input() -> None:
    with pytest.raises(UnsupportedInputFormatError):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho="not-a-dataframe",
                site_metadata=_site_metadata(),
                organism=Organism.RAT,
            )
        )


def test_builder_supports_file_path_inputs(tmp_path) -> None:
    phospho_path = tmp_path / "phospho.csv"
    site_metadata_path = tmp_path / "site_metadata.csv"
    _phospho().to_csv(phospho_path)
    _site_metadata().to_csv(site_metadata_path)

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho_path,
            site_metadata=site_metadata_path,
            organism=Organism.RAT,
        )
    )
    assert list(built.phospho.index) == ["MAPK14;Y182;"]
    assert built.site_metadata.loc["MAPK14;Y182;", "site_sequence"]


def test_builder_derives_site_sequence_for_supported_organism() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata().drop(columns=["site_sequence"]),
            organism=Organism.RAT,
        )
    )
    assert built.site_metadata.loc["MAPK14;Y182;", "site_sequence"]


def test_builder_derives_gene_symbol_and_site_from_supported_index_convention() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=pd.DataFrame(
                {"site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"]},
                index=["MAPK14;Y182;"],
            ),
            organism=Organism.RAT,
        )
    )
    assert built.site_metadata.loc["MAPK14;Y182;", "gene_symbol"] == "MAPK14"
    assert built.site_metadata.loc["MAPK14;Y182;", "site"] == "Y182"


def test_builder_rejects_ambiguous_sequence_column_without_explicit_convention() -> (
    None
):
    with pytest.raises(
        UnsupportedInputFormatError,
        match="column 'sequence' is ambiguous and unsupported",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=_phospho(),
                site_metadata=pd.DataFrame(
                    {
                        "gene_symbol": ["MAPK14"],
                        "site": ["Y182"],
                        "sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                    },
                    index=["MAPK14;Y182;"],
                ),
                organism=Organism.RAT,
            )
        )


def test_builder_rejects_ambiguous_multi_alias_gene_symbol_mapping() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match="has ambiguous columns for 'gene_symbol'",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=_phospho(),
                site_metadata=pd.DataFrame(
                    {
                        "gene_symbol": ["MAPK14"],
                        "gene": ["MAPK14"],
                        "site": ["Y182"],
                        "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                    },
                    index=["MAPK14;Y182;"],
                ),
                organism=Organism.RAT,
            )
        )


def test_builder_fails_when_missing_gene_or_site_cannot_be_derived_from_index() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match="site_metadata is missing required metadata columns",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=_phospho(),
                site_metadata=pd.DataFrame(
                    {"site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"]},
                    index=["MAPK14"],
                ),
                organism=Organism.RAT,
            )
        )


def test_builder_fails_fast_when_site_sequence_derivation_is_unsupported() -> None:
    with pytest.raises(UnsupportedInputFormatError):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=_phospho(),
                site_metadata=_site_metadata().drop(columns=["site_sequence"]),
                organism=None,
            )
        )


def test_io_namespace_does_not_expose_parallel_dataset_file_builder() -> None:
    assert not hasattr(phospy_io, "build_dataset_from_files")


def test_reference_resolver_accepts_explicit_bundle() -> None:
    bundle = _valid_bundle(organism=Organism.MOUSE)
    resolved = ReferenceResolver().run(
        bundle,
        dataset_organism=Organism.MOUSE,
    )
    assert resolved is bundle


def test_builder_establishes_transformation_state_with_supported_path() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            organism=Organism.RAT,
        )
    )
    assert built.transformation_state.label == "linear"
    assert (
        built.transformation_state.phospho.established_by
        == "phospy.transformations.transformers.identity"
    )
