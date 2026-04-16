from __future__ import annotations

import pandas as pd
import pytest

from phospy.datasets.models import (
    AnalysisReadyPhosphoDataset,
    AnalysisReadyPreprocessingProvenance,
    AnalysisReadyRowCounts,
    AnalysisReadySiteMatrixStats,
    DatasetSchema,
    _build_site_to_protein_mapping_series,
    _evaluate_candidate_metadata_column_completeness,
    _find_ambiguous_fallback_identifiers,
    _parse_canonical_entities_for_site_index,
    _resolve_site_to_protein_policy,
)


def _make_analysis_ready_dataset(
    *,
    site_ids: list[str],
    site_metadata: pd.DataFrame,
) -> AnalysisReadyPhosphoDataset:
    site_index = pd.Index(site_ids, name="site_id")
    matrix = pd.DataFrame(
        {"sample_1": [float(index) for index in range(len(site_ids))]},
        index=site_index,
    )
    sequences = pd.Series(
        ["SEQUENCE"] * len(site_index),
        index=site_index,
        dtype="string",
        name="centralized_sequence",
    )
    provenance = AnalysisReadyPreprocessingProvenance(
        source="mapping helper tests",
        schema=DatasetSchema(),
        comparisons=None,
        row_counts=AnalysisReadyRowCounts(
            total_unique=0,
            total_filtered=0,
            phospho_filtered=0,
            phospho_corrected=0,
            phospho_matrix_sites=len(site_index),
        ),
        site_matrix_stats=AnalysisReadySiteMatrixStats(
            input_rows=len(site_index),
            dropped_missing_sequence=0,
            dropped_incomplete_values=0,
            deduplicated_site_rows=0,
            retained_rows=len(site_index),
        ),
    )
    return AnalysisReadyPhosphoDataset.from_external(
        phospho_matrix=matrix,
        site_metadata=site_metadata,
        site_sequences=sequences,
        phospho_corrected=pd.DataFrame(index=site_index),
        provenance=provenance,
    )


def test_analysis_ready_dataset_ownership_accessors_are_explicit() -> None:
    site_ids = ["PRKACA;S339;", "BTK;Y551;"]
    site_index = pd.Index(site_ids, name="site_id")
    dataset = AnalysisReadyPhosphoDataset.from_external(
        phospho_matrix=pd.DataFrame(
            {"sample_1": [1.0, 2.0], "sample_2": [3.0, 4.0]},
            index=site_index,
        ),
        site_metadata=pd.DataFrame(
            {"protein_id": ["P17612", "Q06187"]},
            index=site_index,
        ),
        site_sequences=pd.Series(
            ["AAAAAAA", "BBBBBBB"],
            index=site_index,
            dtype="string",
            name="centralized_sequence",
        ),
        phospho_corrected=pd.DataFrame(
            {"p_group1_group4": [0.1, -0.2]},
            index=site_index,
        ),
        provenance=AnalysisReadyPreprocessingProvenance(
            source="ownership test",
            schema=DatasetSchema(),
            comparisons=None,
            row_counts=AnalysisReadyRowCounts(
                total_unique=2,
                total_filtered=2,
                phospho_filtered=2,
                phospho_corrected=2,
                phospho_matrix_sites=2,
            ),
            site_matrix_stats=AnalysisReadySiteMatrixStats(
                input_rows=2,
                dropped_missing_sequence=0,
                dropped_incomplete_values=0,
                deduplicated_site_rows=0,
                retained_rows=2,
            ),
        ),
    )

    detached_matrix = dataset.to_phospho_matrix()
    detached_metadata = dataset.to_site_metadata()
    detached_sequences = dataset.to_site_sequences()
    detached_corrected = dataset.to_phospho_corrected()
    owned_matrix = dataset.to_owned_phospho_matrix()
    owned_metadata = dataset.to_owned_site_metadata()
    owned_sequences = dataset.to_owned_site_sequences()
    owned_corrected = dataset.to_owned_phospho_corrected()

    assert detached_matrix is not owned_matrix
    assert detached_metadata is not owned_metadata
    assert detached_sequences is not owned_sequences
    assert detached_corrected is not owned_corrected
    assert dataset.to_mutable_phospho_matrix_unsafe() is owned_matrix
    assert dataset.to_mutable_site_metadata_unsafe() is owned_metadata
    assert dataset.to_mutable_site_sequences_unsafe() is owned_sequences
    assert dataset.to_mutable_phospho_corrected_unsafe() is owned_corrected

    matrix_original = float(owned_matrix.iloc[0, 0])
    metadata_original = str(owned_metadata.iloc[0, 0])
    sequence_original = str(owned_sequences.iloc[0])
    corrected_original = float(owned_corrected.iloc[0, 0])

    detached_matrix.iloc[0, 0] = matrix_original + 10.0
    detached_metadata.iloc[0, 0] = "DETACHED_ONLY"
    detached_sequences.iloc[0] = "DETACHEDSEQ"
    detached_corrected.iloc[0, 0] = corrected_original + 10.0

    assert float(owned_matrix.iloc[0, 0]) == matrix_original
    assert str(owned_metadata.iloc[0, 0]) == metadata_original
    assert str(owned_sequences.iloc[0]) == sequence_original
    assert float(owned_corrected.iloc[0, 0]) == corrected_original

    mutable_matrix = dataset.to_mutable_phospho_matrix_unsafe()
    mutable_metadata = dataset.to_mutable_site_metadata_unsafe()
    mutable_sequences = dataset.to_mutable_site_sequences_unsafe()
    mutable_corrected = dataset.to_mutable_phospho_corrected_unsafe()

    mutable_matrix.iloc[0, 0] = matrix_original + 20.0
    mutable_metadata.iloc[0, 0] = "MUTATED_SHARED"
    mutable_sequences.iloc[0] = "MUTABLESEQ"
    mutable_corrected.iloc[0, 0] = corrected_original + 20.0

    assert float(dataset.to_owned_phospho_matrix().iloc[0, 0]) == matrix_original + 20.0
    assert str(dataset.to_owned_site_metadata().iloc[0, 0]) == "MUTATED_SHARED"
    assert str(dataset.to_owned_site_sequences().iloc[0]) == "MUTABLESEQ"
    assert (
        float(dataset.to_owned_phospho_corrected().iloc[0, 0])
        == corrected_original + 20.0
    )


def test_resolve_site_to_protein_policy_normalizes_candidate_columns() -> None:
    policy = _resolve_site_to_protein_policy(
        metadata_columns=(" protein ", "  ", "gene"),
        fallback_policy="metadata",
    )

    assert policy.fallback_policy == "metadata"
    assert policy.candidate_columns == ("protein", "gene")


def test_evaluate_candidate_metadata_column_completeness_returns_diagnostic() -> None:
    metadata_index = pd.Index(
        ["SITE_1", "SITE_2", "SITE_3", "SITE_4"],
        name="site_id",
    )
    values = pd.Series(["", " ", pd.NA, ""], index=metadata_index, dtype="string")

    resolved, diagnostic = _evaluate_candidate_metadata_column_completeness(
        values=values,
        metadata_index=metadata_index,
        candidate_column="protein_id",
    )

    assert resolved is None
    assert diagnostic is not None
    assert (
        "column 'protein_id' has missing/empty values for: SITE_1, SITE_2, SITE_3"
        in (diagnostic)
    )
    assert diagnostic.endswith(", ...")


def test_parse_canonical_entities_for_site_index_extracts_only_parseable_ids() -> None:
    site_index = pd.Index(["PROTEIN_1;S1;", "SITE_2", "protein_3;t3;"], name="site_id")

    canonical_entities = _parse_canonical_entities_for_site_index(site_index)

    assert canonical_entities.iloc[0] == "PROTEIN_1"
    assert pd.isna(canonical_entities.iloc[1])
    assert canonical_entities.iloc[2] == "PROTEIN_3"


def test_find_ambiguous_fallback_identifiers_detects_shared_fallback_values() -> None:
    site_index = pd.Index(["SITE_1", "SITE_2", "SITE_3"], name="site_id")
    fallback_values = pd.Series(
        ["SHARED", "SHARED", "UNIQUE"],
        index=site_index,
        dtype="string",
    )
    canonical_entities = pd.Series(
        ["PROTEIN_A", "PROTEIN_B", pd.NA],
        index=site_index,
        dtype="string",
    )

    ambiguous_identifiers = _find_ambiguous_fallback_identifiers(
        fallback_values=fallback_values,
        canonical_entities=canonical_entities,
    )

    assert ambiguous_identifiers == ["SHARED"]


def test_build_site_to_protein_mapping_series_shapes_final_output() -> None:
    site_index = pd.Index(["SITE_1", "SITE_2"], name="site_id")
    resolved_values = pd.Series(["P1", "P2"], index=site_index, dtype="string")

    mapping = _build_site_to_protein_mapping_series(
        site_index=site_index,
        resolved_values=resolved_values,
    )

    assert mapping.name == "protein_id"
    assert mapping.dtype == object
    assert mapping.index.name == "site_id"
    assert mapping.to_dict() == {"SITE_1": "P1", "SITE_2": "P2"}


def test_resolve_site_to_protein_mapping_falls_back_to_next_complete_column() -> None:
    site_ids = ["PROTEIN_1;S1;", "PROTEIN_2;S2;"]
    site_index = pd.Index(site_ids, name="site_id")
    site_metadata = pd.DataFrame(
        {
            "protein_id": ["PROTEIN_1", ""],
            "protein": ["PROTEIN_1", "PROTEIN_2"],
        },
        index=site_index,
    )
    dataset = _make_analysis_ready_dataset(
        site_ids=site_ids,
        site_metadata=site_metadata,
    )

    mapping = dataset.resolve_site_to_protein_mapping(
        fallback_policy="metadata",
        metadata_columns=("protein_id", "protein"),
    )

    assert mapping.to_dict() == {
        "PROTEIN_1;S1;": "PROTEIN_1",
        "PROTEIN_2;S2;": "PROTEIN_2",
    }


def test_resolve_site_to_protein_mapping_reports_strict_diagnostics() -> None:
    site_ids = ["PROTEIN_1;S1;", "PROTEIN_2;S2;"]
    site_index = pd.Index(site_ids, name="site_id")
    dataset = _make_analysis_ready_dataset(
        site_ids=site_ids,
        site_metadata=pd.DataFrame(
            {"protein_id": ["PROTEIN_1", "PROTEIN_2"]},
            index=site_index,
        ),
    )

    result = dataset.resolve_site_to_protein_mapping_with_diagnostics()

    assert result.diagnostics.fallback_policy == "strict"
    assert result.diagnostics.chosen_identifier_column == "protein_id"
    assert result.diagnostics.fallback_mode == "strict_protein_id"
    assert result.diagnostics.ambiguous_identifier_count == 0
    assert result.diagnostics.ambiguous_identifiers == ()
    assert result.diagnostics.checked_columns == ("protein_id",)
    assert result.mapping.attrs["site_to_protein_resolution_diagnostics"] == (
        result.diagnostics
    )
    assert dataset.last_site_to_protein_resolution_diagnostics == result.diagnostics


def test_resolve_site_to_protein_mapping_reports_metadata_fallback_diagnostics() -> (
    None
):
    site_ids = ["PROTEIN_1;S1;", "PROTEIN_2;S2;"]
    site_index = pd.Index(site_ids, name="site_id")
    dataset = _make_analysis_ready_dataset(
        site_ids=site_ids,
        site_metadata=pd.DataFrame(
            {
                "protein_id": ["PROTEIN_1", ""],
                "protein": ["PROTEIN_1", "PROTEIN_2"],
            },
            index=site_index,
        ),
    )

    result = dataset.resolve_site_to_protein_mapping_with_diagnostics(
        fallback_policy="metadata",
        metadata_columns=("protein_id", "protein"),
    )

    assert result.diagnostics.fallback_policy == "metadata"
    assert result.diagnostics.chosen_identifier_column == "protein"
    assert result.diagnostics.fallback_mode == "metadata_identifier_column"
    assert result.diagnostics.checked_columns == ("protein_id", "protein")
    assert result.diagnostics.incomplete_candidate_diagnostics and (
        "column 'protein_id' has missing/empty values"
        in result.diagnostics.incomplete_candidate_diagnostics[0]
    )
    assert result.diagnostics.ambiguous_identifier_count == 0
    assert dataset.last_site_to_protein_resolution_diagnostics == result.diagnostics


def test_resolve_site_to_protein_mapping_reports_ambiguous_gene_diagnostics() -> None:
    site_ids = ["PROTEIN_A;S1;", "PROTEIN_B;S2;"]
    site_index = pd.Index(site_ids, name="site_id")
    dataset = _make_analysis_ready_dataset(
        site_ids=site_ids,
        site_metadata=pd.DataFrame(
            {"gene": ["SHARED_GENE", "SHARED_GENE"]},
            index=site_index,
        ),
    )

    with pytest.warns(UserWarning) as captured_warnings:
        result = dataset.resolve_site_to_protein_mapping_with_diagnostics(
            fallback_policy="metadata",
            metadata_columns=("gene",),
            allow_gene_symbol_fallback=True,
            allow_ambiguous_fallback=True,
        )

    warning_messages = [str(record.message) for record in captured_warnings]
    assert any(
        "Proceeding because allow_ambiguous_fallback=True" in message
        for message in warning_messages
    )

    assert result.diagnostics.chosen_identifier_column == "gene"
    assert result.diagnostics.fallback_mode == "metadata_gene_symbol"
    assert result.diagnostics.gene_symbol_fallback_used
    assert result.diagnostics.ambiguous_fallback_allowed
    assert result.diagnostics.ambiguous_identifier_count == 1
    assert result.diagnostics.ambiguous_identifiers == ("SHARED_GENE",)
