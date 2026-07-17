from __future__ import annotations

import inspect
from dataclasses import dataclass

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.errors.validation import DatasetValidationError
from phospy.provenance.derived_quantitative import (
    DerivedQuantitativeDataProvenance,
    DerivedSampleMapping,
    build_derived_quantitative_run_provenance,
)
from phospy.provenance.hashing import fingerprint_optional_table, fingerprint_table
from phospy.provenance.models import RunProvenance, TableFingerprint
from phospy.science.datasets.derived_quantitative import (
    DerivedAnalysisReadyPhosphoDataset,
)
from phospy.science.references.models import Organism
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.processing_state import imputed_processing_state
from tests.support.site_keys import protein_site_key_index, site_key_context_columns


@dataclass(frozen=True, slots=True)
class _DerivedTables:
    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame
    total: pd.DataFrame
    comparisons: pd.DataFrame
    mask: pd.DataFrame


def _mutate_site_metadata_score(tables: _DerivedTables) -> None:
    tables.site_metadata.loc[:, "metadata_score"] = [9.0, 0.2]


def _mutate_observation_mask(tables: _DerivedTables) -> None:
    tables.mask.loc[:, "bio_a"] = [False, True]


def _mutate_total_matrix(tables: _DerivedTables) -> None:
    tables.total.loc[:, "bio_a"] = [10.0, 2.0]


def _mutate_comparisons(tables: _DerivedTables) -> None:
    tables.comparisons.loc[:, "bio_b_vs_bio_a"] = [40.0, 4.0]


def test_public_derived_constructor_rejects_stale_matrix_hash() -> None:
    tables = _valid_tables()
    lineage = _lineage_for_tables(tables)
    provenance = _provenance_for(lineage)
    stale_actual = _copy_tables(tables)
    stale_actual.phospho.iloc[0, 0] = 99.0

    with pytest.raises(
        DatasetValidationError,
        match="dataset\\.phospho.*exact_hash_value",
    ):
        _build_dataset(
            stale_actual,
            lineage=lineage,
            provenance=provenance,
        )


@pytest.mark.parametrize(
    ("mutator", "expected_table"),
    [
        pytest.param(
            _mutate_site_metadata_score,
            "dataset.site_metadata",
            id="metadata",
        ),
        pytest.param(
            _mutate_observation_mask,
            "dataset.imputation_observation_mask",
            id="mask",
        ),
        pytest.param(
            _mutate_total_matrix,
            "dataset.total",
            id="total",
        ),
        pytest.param(
            _mutate_comparisons,
            "dataset.comparisons",
            id="comparisons",
        ),
    ],
)
def test_public_derived_constructor_rejects_optional_table_mismatches(
    mutator,
    expected_table: str,
) -> None:
    tables = _valid_tables()
    lineage = _lineage_for_tables(tables)
    provenance = _provenance_for(lineage)
    stale_actual = _copy_tables(tables)
    mutator(stale_actual)

    with pytest.raises(DatasetValidationError, match=expected_table):
        _build_dataset(
            stale_actual,
            lineage=lineage,
            provenance=provenance,
        )


def test_public_derived_constructor_rejects_shape_mismatch() -> None:
    tables = _valid_tables()
    lineage = _lineage_for_tables(tables)
    provenance = _provenance_for(lineage)
    one_row = _copy_tables(tables)
    retained_index = one_row.phospho.index[:1]
    one_row.phospho.drop(index=one_row.phospho.index[1:], inplace=True)
    one_row.site_metadata.drop(index=one_row.site_metadata.index[1:], inplace=True)
    one_row.mask.drop(index=one_row.mask.index[1:], inplace=True)
    one_row.comparisons.drop(index=one_row.comparisons.index[1:], inplace=True)
    assert one_row.phospho.index.equals(retained_index)

    with pytest.raises(DatasetValidationError, match="dataset\\.phospho.*rows"):
        _build_dataset(
            one_row,
            lineage=lineage,
            provenance=provenance,
        )


def test_public_derived_constructor_rejects_row_order_mismatch() -> None:
    tables = _valid_tables()
    lineage = _lineage_for_tables(tables)
    provenance = _provenance_for(lineage)
    row_order = list(reversed(tables.phospho.index.tolist()))
    reversed_rows = _DerivedTables(
        phospho=tables.phospho.loc[row_order, :].copy(deep=True),
        site_metadata=tables.site_metadata.loc[row_order, :].copy(deep=True),
        sample_metadata=tables.sample_metadata.copy(deep=True),
        total=tables.total.copy(deep=True),
        comparisons=tables.comparisons.loc[row_order, :].copy(deep=True),
        mask=tables.mask.loc[row_order, :].copy(deep=True),
    )

    with pytest.raises(
        DatasetValidationError,
        match="dataset\\.phospho.*index_structure",
    ):
        _build_dataset(
            reversed_rows,
            lineage=lineage,
            provenance=provenance,
        )


def test_public_derived_constructor_rejects_sample_mapping_mismatch() -> None:
    tables = _valid_tables()
    renamed = _copy_tables(tables)
    rename_map = {"bio_a": "other_a", "bio_b": "other_b"}
    for matrix in (renamed.phospho, renamed.total, renamed.mask):
        matrix.rename(columns=rename_map, inplace=True)
    renamed.sample_metadata.rename(index=rename_map, inplace=True)
    lineage = _lineage_for_tables(renamed)
    provenance = _provenance_for(lineage)

    with pytest.raises(DatasetValidationError, match="sample_mapping"):
        _build_dataset(
            tables,
            lineage=lineage,
            provenance=provenance,
        )


def test_public_derived_constructor_rejects_run_provenance_output_mismatch() -> None:
    tables = _valid_tables()
    lineage = _lineage_for_tables(tables)
    provenance_lineage = _lineage_for_tables(_mutated_phospho_tables(tables))
    provenance = _provenance_for(provenance_lineage)

    with pytest.raises(DatasetValidationError, match="run_provenance.output_tables"):
        _build_dataset(
            tables,
            lineage=lineage,
            provenance=provenance,
        )


def test_public_derived_constructor_rejects_free_text_lineage_substitution() -> None:
    tables = _valid_tables()
    lineage = _lineage_for_tables(tables)
    valid_provenance = _provenance_for(lineage)
    provenance = RunProvenance(
        environment=valid_provenance.environment,
        input_tables=lineage.parent_dataset_fingerprints,
        preprocessing_stages=(),
        reference=None,
        workflow_name=lineage.derivation_type,
        workflow_parameters={"derived_quantitative_data": "lineage is fine"},
        random_state=None,
        random_seed_policy=None,
        output_tables=lineage.derived_dataset_fingerprints,
        reference_context=None,
    )

    with pytest.raises(DatasetValidationError, match="typed derived_quantitative_data"):
        _build_dataset(
            tables,
            lineage=lineage,
            provenance=provenance,
        )


def test_public_derived_constructor_copies_caller_tables() -> None:
    tables = _valid_tables()
    lineage = _lineage_for_tables(tables)
    provenance = _provenance_for(lineage)

    dataset = _build_dataset(tables, lineage=lineage, provenance=provenance)
    expected_phospho = dataset.phospho
    expected_metadata = dataset.site_metadata
    expected_sample_metadata = dataset.sample_metadata
    expected_total = dataset.total
    expected_comparisons = dataset.comparisons
    expected_mask = dataset.imputation_observed_mask_dataframe()

    tables.phospho.iloc[0, 0] = 99.0
    tables.site_metadata.loc[:, "metadata_score"] = [9.0, 9.0]
    tables.sample_metadata.loc[:, "condition"] = ["changed", "changed"]
    tables.total.iloc[0, 0] = 99.0
    tables.comparisons.iloc[0, 0] = 99.0
    tables.mask.iloc[0, 0] = False

    pdt.assert_frame_equal(dataset.phospho, expected_phospho)
    pdt.assert_frame_equal(dataset.site_metadata, expected_metadata)
    assert dataset.sample_metadata is not None
    assert expected_sample_metadata is not None
    pdt.assert_frame_equal(dataset.sample_metadata, expected_sample_metadata)
    assert dataset.total is not None
    assert expected_total is not None
    pdt.assert_frame_equal(dataset.total, expected_total)
    assert dataset.comparisons is not None
    assert expected_comparisons is not None
    pdt.assert_frame_equal(dataset.comparisons, expected_comparisons)
    assert dataset.imputation_observed_mask_dataframe() is not None
    assert expected_mask is not None
    pdt.assert_frame_equal(dataset.imputation_observed_mask_dataframe(), expected_mask)


def test_owned_transfer_path_is_package_private() -> None:
    assert (
        "_assume_owned"
        not in inspect.signature(DerivedAnalysisReadyPhosphoDataset).parameters
    )
    assert not hasattr(
        DerivedAnalysisReadyPhosphoDataset,
        "from_owned_derived_tables",
    )
    assert hasattr(
        DerivedAnalysisReadyPhosphoDataset,
        "_from_owned_derived_tables",
    )


def _valid_tables() -> _DerivedTables:
    display_ids = ["MAPK14;Y182;", "GSK3B;S9;"]
    sites = ["Y182", "S9"]
    site_index = protein_site_key_index(
        protein_identifiers=["MAPK14", "GSK3B"],
        sites=sites,
    )
    phospho = pd.DataFrame(
        {
            "bio_a": [2.0, 4.0],
            "bio_b": [6.0, 8.0],
        },
        index=site_index.copy(),
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": display_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": sites,
            "site_sequence": [("A" * 15) + site[0] + ("A" * 15) for site in sites],
            "protein_id": ["MAPK14", "GSK3B"],
            "metadata_score": [0.1, 0.2],
        },
        index=site_index.copy(),
    )
    sample_metadata = pd.DataFrame(
        {
            "condition": ["A", "B"],
            "biological_replicate_id": ["bio_a", "bio_b"],
        },
        index=pd.Index(["bio_a", "bio_b"], name="sample_id"),
    )
    total = pd.DataFrame(
        {
            "bio_a": [1.0, 2.0],
            "bio_b": [3.0, 4.0],
        },
        index=pd.Index(["MAPK14", "GSK3B"], name="protein_id"),
    )
    comparisons = pd.DataFrame(
        {"bio_b_vs_bio_a": [4.0, 4.0]},
        index=site_index.copy(),
    )
    mask = pd.DataFrame(
        {
            "bio_a": [True, True],
            "bio_b": [False, True],
        },
        index=site_index.copy(),
    )
    return _DerivedTables(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=total,
        comparisons=comparisons,
        mask=mask,
    )


def _copy_tables(tables: _DerivedTables) -> _DerivedTables:
    return _DerivedTables(
        phospho=tables.phospho.copy(deep=True),
        site_metadata=tables.site_metadata.copy(deep=True),
        sample_metadata=tables.sample_metadata.copy(deep=True),
        total=tables.total.copy(deep=True),
        comparisons=tables.comparisons.copy(deep=True),
        mask=tables.mask.copy(deep=True),
    )


def _mutated_phospho_tables(tables: _DerivedTables) -> _DerivedTables:
    mutated = _copy_tables(tables)
    mutated.phospho.iloc[0, 0] = 99.0
    return mutated


def _lineage_for_tables(tables: _DerivedTables) -> DerivedQuantitativeDataProvenance:
    parent_phospho = pd.DataFrame(
        {
            "bio_a_t1": [1.0, 3.0],
            "bio_a_t2": [3.0, 5.0],
            "bio_b_t1": [5.0, 7.0],
            "bio_b_t2": [7.0, 9.0],
        },
        index=tables.phospho.index.copy(),
    )
    return DerivedQuantitativeDataProvenance(
        derivation_type="unit_derived_quantitative_state",
        parent_dataset_type="AnalysisReadyPhosphoDataset",
        derived_dataset_type="DerivedAnalysisReadyPhosphoDataset",
        parent_dataset_fingerprints=(
            fingerprint_table(parent_phospho, name="dataset.phospho"),
        ),
        derived_dataset_fingerprints=_fingerprints_for_tables(tables),
        sample_mapping=tuple(
            DerivedSampleMapping(
                output_sample_id=str(sample_id),
                input_sample_ids=(f"{sample_id}_t1", f"{sample_id}_t2"),
                condition=str(tables.sample_metadata.loc[sample_id, "condition"]),
                biological_replicate_id=str(sample_id),
                technical_replicate_ids=("t1", "t2"),
            )
            for sample_id in tables.phospho.columns.astype(str).tolist()
        ),
        aggregation_method="mean",
        input_intensity_scale="linear",
        output_intensity_scale="linear",
        quantitative_meaning="phosphosite_abundance",
        missingness_policy={
            "policy": "impute_row_median",
            "complete_matrix": True,
            "imputed": True,
        },
        matrices_transformed={
            "phospho": True,
            "sample_metadata": True,
            "total_protein": True,
            "imputation_observation_mask": True,
            "comparisons": False,
        },
        implementation="tests.unit.derived_quantitative_state",
        implementation_version="1",
        parameters={"aggregation_axis": "samples"},
    )


def _fingerprints_for_tables(
    tables: _DerivedTables,
) -> tuple[TableFingerprint, ...]:
    entries = (
        ("dataset.phospho", tables.phospho),
        ("dataset.site_metadata", tables.site_metadata),
        ("dataset.sample_metadata", tables.sample_metadata),
        ("dataset.total", tables.total),
        ("dataset.comparisons", tables.comparisons),
        ("dataset.imputation_observation_mask", tables.mask),
    )
    fingerprints: list[TableFingerprint] = []
    for name, table in entries:
        fingerprint = fingerprint_optional_table(table, name=name)
        if fingerprint is not None:
            fingerprints.append(fingerprint)
    return tuple(fingerprints)


def _provenance_for(lineage: DerivedQuantitativeDataProvenance) -> RunProvenance:
    return build_derived_quantitative_run_provenance(lineage=lineage)


def _build_dataset(
    tables: _DerivedTables,
    *,
    lineage: DerivedQuantitativeDataProvenance,
    provenance: RunProvenance,
) -> DerivedAnalysisReadyPhosphoDataset:
    base_processing_state = supported_linear_processing_state(has_total_matrix=True)
    return DerivedAnalysisReadyPhosphoDataset(
        phospho=tables.phospho,
        site_metadata=tables.site_metadata,
        sample_metadata=tables.sample_metadata,
        total=tables.total,
        comparisons=tables.comparisons,
        imputation_observation_mask=tables.mask,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=True
        ),
        processing_state=imputed_processing_state(base_processing_state),
        derived_lineage=lineage,
        provenance=provenance,
    )
