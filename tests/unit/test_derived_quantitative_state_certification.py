from __future__ import annotations

import inspect
from dataclasses import dataclass, replace

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.errors import PhosPyInputError
from phospy.errors.validation import DatasetValidationError
from phospy.provenance.derived_quantitative import (
    DerivedQuantitativeDataProvenance,
    DerivedSampleMapping,
    build_derived_quantitative_run_provenance,
)
from phospy.provenance.hashing import fingerprint_optional_table_strict
from phospy.provenance.models import RunProvenance, TableFingerprint
from phospy.science.datasets.derived_quantitative import (
    CertifiedDerivedQuantitativeParentState,
    DerivedAnalysisReadyPhosphoDataset,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.references.models import Organism
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
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


@dataclass(frozen=True, slots=True)
class _ParentTables:
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
    parent_state = _parent_state_for_tables(tables)
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
            parent_state=parent_state,
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
    parent_state = _parent_state_for_tables(tables)
    lineage = _lineage_for_tables(tables)
    provenance = _provenance_for(lineage)
    stale_actual = _copy_tables(tables)
    mutator(stale_actual)

    with pytest.raises(DatasetValidationError, match=expected_table):
        _build_dataset(
            stale_actual,
            lineage=lineage,
            provenance=provenance,
            parent_state=parent_state,
        )


def test_public_derived_constructor_rejects_shape_mismatch() -> None:
    tables = _valid_tables()
    parent_state = _parent_state_for_tables(tables)
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
            parent_state=parent_state,
        )


def test_public_derived_constructor_rejects_row_order_mismatch() -> None:
    tables = _valid_tables()
    parent_state = _parent_state_for_tables(tables)
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
            parent_state=parent_state,
        )


def test_derived_lineage_rejects_sample_mapping_mismatch_before_constructor() -> None:
    tables = _valid_tables()
    valid_lineage = _lineage_for_tables(tables)

    with pytest.raises(PhosPyInputError, match="sample_mapping output_sample_id"):
        replace(
            valid_lineage,
            sample_mapping=tuple(
                replace(
                    mapping,
                    output_sample_id=f"other_{mapping.output_sample_id}",
                )
                for mapping in valid_lineage.sample_mapping
            ),
        )


def test_public_derived_constructor_rejects_fabricated_parent_input_samples() -> None:
    tables = _valid_tables()
    lineage = _lineage_for_tables(
        tables,
        input_sample_ids_by_output={
            "bio_a": ("ghost_a_t1", "ghost_a_t2"),
            "bio_b": ("ghost_b_t1", "ghost_b_t2"),
        },
    )
    provenance = _provenance_for(lineage)

    with pytest.raises(DatasetValidationError, match="input_sample_ids"):
        _build_dataset(
            tables,
            lineage=lineage,
            provenance=provenance,
        )


def test_public_derived_constructor_rejects_false_input_intensity_scale() -> None:
    tables = _valid_tables()
    lineage = _lineage_for_tables(tables, input_intensity_scale="log2")
    provenance = _provenance_for(lineage)

    with pytest.raises(DatasetValidationError, match="input_intensity_scale"):
        _build_dataset(
            tables,
            lineage=lineage,
            provenance=provenance,
        )


def test_public_derived_constructor_rejects_false_output_intensity_scale() -> None:
    tables = _valid_tables()
    lineage = _lineage_for_tables(tables, output_intensity_scale="log2")
    provenance = _provenance_for(lineage)

    with pytest.raises(DatasetValidationError, match="output_intensity_scale"):
        _build_dataset(
            tables,
            lineage=lineage,
            provenance=provenance,
        )


def test_public_derived_constructor_rejects_false_quantitative_meaning() -> None:
    tables = _valid_tables()
    lineage = _lineage_for_tables(tables, quantitative_meaning="unknown")
    provenance = _provenance_for(lineage)

    with pytest.raises(DatasetValidationError, match="quantitative_meaning"):
        _build_dataset(
            tables,
            lineage=lineage,
            provenance=provenance,
        )


def test_public_derived_constructor_rejects_false_matrix_transformation_flags() -> None:
    tables = _valid_tables()
    lineage = _lineage_for_tables(
        tables,
        matrices_transformed={
            "phospho": False,
            "sample_metadata": True,
            "total_protein": True,
            "imputation_observation_mask": True,
            "comparisons": False,
        },
    )
    provenance = _provenance_for(lineage)

    with pytest.raises(DatasetValidationError, match="matrices_transformed"):
        _build_dataset(
            tables,
            lineage=lineage,
            provenance=provenance,
        )


@pytest.mark.parametrize(
    ("input_sample_ids_by_output", "expected_message"),
    [
        pytest.param(
            {
                "bio_a": ("bio_a_t1", "bio_a_t1"),
                "bio_b": ("bio_b_t1", "bio_b_t2"),
            },
            "duplicate",
            id="duplicate",
        ),
        pytest.param(
            {
                "bio_a": ("bio_a_t1",),
                "bio_b": ("bio_b_t1", "bio_b_t2"),
            },
            "omitted",
            id="omitted",
        ),
    ],
)
def test_public_derived_constructor_rejects_duplicate_or_omitted_parent_samples(
    input_sample_ids_by_output: dict[str, tuple[str, ...]],
    expected_message: str,
) -> None:
    tables = _valid_tables()
    lineage = _lineage_for_tables(
        tables,
        derivation_type="technical_replicate_aggregation",
        input_sample_ids_by_output=input_sample_ids_by_output,
    )
    provenance = _provenance_for(lineage)

    with pytest.raises(DatasetValidationError, match=expected_message):
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
    actual_mask = dataset.imputation_observed_mask_dataframe()
    assert actual_mask is not None
    assert expected_mask is not None
    pdt.assert_frame_equal(actual_mask, expected_mask)


def test_public_derived_constructor_isolated_from_caller_parent_and_output_frames() -> (
    None
):
    tables = _valid_tables()
    parent_tables = _parent_tables_for(tables)
    parent_dataset = _parent_dataset_from_tables(parent_tables)
    parent_state = CertifiedDerivedQuantitativeParentState.from_dataset(parent_dataset)
    lineage = _lineage_for_tables(tables)
    provenance = _provenance_for(lineage)

    dataset = _build_dataset(
        tables,
        lineage=lineage,
        provenance=provenance,
        parent_state=parent_state,
    )
    expected_phospho = dataset.phospho
    expected_parent_samples = parent_state.phospho_sample_ids
    expected_parent_fingerprints = parent_state.parent_dataset_fingerprints

    parent_tables.phospho.iloc[0, 0] = 99.0
    parent_tables.sample_metadata.loc[:, "condition"] = ["changed"] * 4
    parent_export = parent_dataset.phospho
    parent_export.iloc[0, 0] = 123.0
    tables.phospho.iloc[0, 0] = 456.0

    pdt.assert_frame_equal(dataset.phospho, expected_phospho)
    assert parent_state.phospho_sample_ids == expected_parent_samples
    assert parent_state.parent_dataset_fingerprints == expected_parent_fingerprints


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


def _lineage_for_tables(
    tables: _DerivedTables,
    *,
    derivation_type: str = "unit_derived_quantitative_state",
    input_sample_ids_by_output: dict[str, tuple[str, ...]] | None = None,
    input_intensity_scale: str = "linear",
    output_intensity_scale: str = "linear",
    quantitative_meaning: str = "phosphosite_abundance",
    matrices_transformed: dict[str, bool] | None = None,
) -> DerivedQuantitativeDataProvenance:
    parent_state = _parent_state_for_tables(tables)
    return DerivedQuantitativeDataProvenance(
        derivation_type=derivation_type,
        parent_dataset_type="AnalysisReadyPhosphoDataset",
        derived_dataset_type="DerivedAnalysisReadyPhosphoDataset",
        parent_dataset_fingerprints=parent_state.parent_dataset_fingerprints,
        derived_dataset_fingerprints=_fingerprints_for_tables(tables),
        sample_mapping=tuple(
            DerivedSampleMapping(
                output_sample_id=str(sample_id),
                input_sample_ids=(
                    input_sample_ids_by_output[str(sample_id)]
                    if input_sample_ids_by_output is not None
                    else (f"{sample_id}_t1", f"{sample_id}_t2")
                ),
                condition=str(tables.sample_metadata.loc[sample_id, "condition"]),
                biological_replicate_id=str(sample_id),
                technical_replicate_ids=("t1", "t2"),
            )
            for sample_id in tables.phospho.columns.astype(str).tolist()
        ),
        aggregation_method="mean",
        input_intensity_scale=input_intensity_scale,
        output_intensity_scale=output_intensity_scale,
        quantitative_meaning=quantitative_meaning,
        missingness_policy={
            "policy": "impute_row_median",
            "complete_matrix": True,
            "imputed": True,
        },
        matrices_transformed=(
            {
                "phospho": True,
                "sample_metadata": True,
                "total_protein": True,
                "imputation_observation_mask": True,
                "comparisons": False,
            }
            if matrices_transformed is None
            else matrices_transformed
        ),
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
        fingerprint = fingerprint_optional_table_strict(table, name=name)
        if fingerprint is not None:
            fingerprints.append(fingerprint)
    return tuple(fingerprints)


def _provenance_for(lineage: DerivedQuantitativeDataProvenance) -> RunProvenance:
    return build_derived_quantitative_run_provenance(lineage=lineage)


def _parent_tables_for(tables: _DerivedTables) -> _ParentTables:
    parent_phospho = pd.DataFrame(
        {
            "bio_a_t1": [1.0, 3.0],
            "bio_a_t2": [3.0, 5.0],
            "bio_b_t1": [5.0, 7.0],
            "bio_b_t2": [7.0, 9.0],
        },
        index=tables.phospho.index.copy(),
    )
    parent_sample_metadata = pd.DataFrame(
        {
            "condition": ["A", "A", "B", "B"],
            "biological_replicate_id": ["bio_a", "bio_a", "bio_b", "bio_b"],
            "technical_replicate_id": ["t1", "t2", "t1", "t2"],
        },
        index=pd.Index(parent_phospho.columns.tolist(), name="sample_id"),
    )
    parent_total = pd.DataFrame(
        {
            "bio_a_t1": [0.5, 1.5],
            "bio_a_t2": [1.5, 2.5],
            "bio_b_t1": [2.5, 3.5],
            "bio_b_t2": [3.5, 4.5],
        },
        index=tables.total.index.copy(),
    )
    parent_mask = pd.DataFrame(
        {
            "bio_a_t1": [True, True],
            "bio_a_t2": [True, True],
            "bio_b_t1": [False, True],
            "bio_b_t2": [False, True],
        },
        index=tables.mask.index.copy(),
    )
    return _ParentTables(
        phospho=parent_phospho,
        site_metadata=tables.site_metadata.copy(deep=True),
        sample_metadata=parent_sample_metadata,
        total=parent_total,
        comparisons=tables.comparisons.copy(deep=True),
        mask=parent_mask,
    )


def _parent_dataset_for(tables: _DerivedTables) -> AnalysisReadyPhosphoDataset:
    return _parent_dataset_from_tables(_parent_tables_for(tables))


def _parent_dataset_from_tables(
    parent_tables: _ParentTables,
) -> AnalysisReadyPhosphoDataset:
    base_processing_state = supported_linear_processing_state(has_total_matrix=True)
    return trusted_analysis_ready_dataset_from_tables(
        phospho=parent_tables.phospho,
        site_metadata=parent_tables.site_metadata,
        sample_metadata=parent_tables.sample_metadata,
        total=parent_tables.total,
        comparisons=parent_tables.comparisons,
        imputation_observation_mask=parent_tables.mask,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=True
        ),
        processing_state=imputed_processing_state(base_processing_state),
    )


def _parent_state_for_tables(
    tables: _DerivedTables,
) -> CertifiedDerivedQuantitativeParentState:
    return CertifiedDerivedQuantitativeParentState.from_dataset(
        _parent_dataset_for(tables)
    )


def _build_dataset(
    tables: _DerivedTables,
    *,
    lineage: DerivedQuantitativeDataProvenance,
    provenance: RunProvenance,
    parent_state: CertifiedDerivedQuantitativeParentState | None = None,
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
        parent_state=(
            _parent_state_for_tables(tables) if parent_state is None else parent_state
        ),
        derived_lineage=lineage,
        provenance=provenance,
    )
