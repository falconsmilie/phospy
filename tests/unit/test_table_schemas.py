from __future__ import annotations

import pandas as pd
import pytest
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from phospy.errors.validation import (
    DatasetValidationError,
    PhosPyValidationError,
    ReferenceValidationError,
    WorkflowValidationError,
)
from phospy.science.references.models import Organism, ReferenceBundle
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    decode_site_key,
    encode_site_key,
)
from phospy.science.tables.activity import ActivityMatrix, ActivityTargetTable
from phospy.science.tables.datasets import (
    PhosphoIntensityMatrix,
    SampleMetadataTable,
    SiteMetadataTable,
    TotalProteinMatrix,
)
from phospy.science.tables.kinase import KinasePredictionMatrix, KinaseScoreMatrix
from phospy.science.tables.references import (
    KinaseSubstrateReference,
    SiteSequenceReference,
)
from phospy.science.tables.signalome import (
    SignalomeProteinSiteContext,
    SignalomeSiteContext,
)

_PROPERTY_SETTINGS = settings(
    max_examples=30,
    deadline=None,
    derandomize=True,
)


def _canonical_site_ids(
    *,
    min_size: int = 1,
    max_size: int = 6,
) -> st.SearchStrategy[list[str]]:
    return st.lists(
        st.integers(min_value=1, max_value=9999),
        min_size=min_size,
        max_size=max_size,
        unique=True,
    ).map(lambda ids: [f"G{site_id};S{site_id};" for site_id in ids])


def _parse_site_id(site_id: str) -> tuple[str, str]:
    gene_symbol, site, _ = site_id.split(";")
    return gene_symbol, site


def _site_key_for_display_id(
    display_id: str,
    *,
    protein_identifier: str | None = None,
) -> str:
    gene_symbol, site = _parse_site_id(display_id)
    key = build_protein_scoped_site_key(
        organism="rat",
        protein_namespace="protein_id",
        protein_identifier=protein_identifier or gene_symbol,
        residue=site[0],
        position=int(site[1:]),
        field_name="tests.unit.test_table_schemas.site_key",
        error_type=ValueError,
    )
    return encode_site_key(key)


def _site_keys_for_display_ids(display_ids: list[str]) -> list[str]:
    return [_site_key_for_display_id(display_id) for display_id in display_ids]


_MAPK14_Y182 = _site_key_for_display_id("MAPK14;Y182;", protein_identifier="P28482")
_MAPK14_T185 = _site_key_for_display_id("MAPK14;T185;", protein_identifier="P28482")
_AKT1_T308 = _site_key_for_display_id("AKT1;T308;", protein_identifier="P31749")
_GSK3B_S9 = _site_key_for_display_id("GSK3B;S9;", protein_identifier="GSK3B")
_DISPLAY_ID_BY_SITE_KEY = {
    _MAPK14_Y182: "MAPK14;Y182;",
    _MAPK14_T185: "MAPK14;T185;",
    _AKT1_T308: "AKT1;T308;",
    _GSK3B_S9: "GSK3B;S9;",
}


def _phospho_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=pd.Index([_MAPK14_Y182, _AKT1_T308], name="site_key"),
    )


def _site_metadata_frame(index: pd.Index) -> pd.DataFrame:
    site_keys = index.astype(str).tolist()
    decoded_keys = [
        decode_site_key(
            site_key,
            field_name="tests.unit.test_table_schemas.site_key",
            error_type=ValueError,
        )
        for site_key in site_keys
    ]
    display_ids = [
        _DISPLAY_ID_BY_SITE_KEY.get(
            site_key,
            f"{key.protein_identifier};{key.residue}{key.position};",
        )
        for site_key, key in zip(site_keys, decoded_keys, strict=True)
    ]
    parsed_display = [_parse_site_id(display_id) for display_id in display_ids]
    return pd.DataFrame(
        {
            "site_key": site_keys,
            "display_id": display_ids,
            "organism": [key.organism.value for key in decoded_keys],
            "protein_namespace": [key.protein_namespace for key in decoded_keys],
            "protein_identifier": [key.protein_identifier for key in decoded_keys],
            "gene_symbol": [gene for gene, _ in parsed_display],
            "site": [site for _, site in parsed_display],
            "site_sequence": [
                ("A" * 15) + site[0] + ("A" * 15) for _, site in parsed_display
            ],
            "protein_group_id": [key.protein_identifier for key in decoded_keys],
        },
        index=index.copy(),
    )


def _site_key_context_columns(site_keys: list[str] | pd.Index) -> dict[str, list[str]]:
    decoded_keys = [
        decode_site_key(
            site_key,
            field_name="tests.unit.test_table_schemas.site_key",
            error_type=ValueError,
        )
        for site_key in pd.Index(site_keys).astype(str).tolist()
    ]
    return {
        "organism": [key.organism.value for key in decoded_keys],
        "protein_namespace": [key.protein_namespace for key in decoded_keys],
        "protein_identifier": [key.protein_identifier for key in decoded_keys],
    }


def _site_membership_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_key": [_MAPK14_Y182, _AKT1_T308],
            "display_id": ["MAPK14;Y182;", "AKT1;T308;"],
            "site_id": ["MAPK14;Y182;", "AKT1;T308;"],
            "site": ["Y182", "T308"],
            "protein_group_id": ["P28482", "P31749"],
            "protein_accession": ["P28482", "P31749"],
            "isoform_id": ["", ""],
            "site_cluster": [1, pd.NA],
            "protein_module_id": [1, 0],
            "included_in_module_table": [True, False],
            "excluded_reason": ["", "dropped_all_missing_downstream_scores"],
            "gene_symbol": ["MAPK14", "AKT1"],
            "top_kinase": ["MAP2K6", ""],
            "top_kinase_score": [0.91, float("nan")],
            "top_kinase_weight": [0.83, float("nan")],
            "n_supported_kinases": [2, 0],
        }
    )


def _protein_site_context_frame() -> pd.DataFrame:
    site_keys_json = f'["{_MAPK14_Y182}","{_MAPK14_T185}"]'
    display_ids_json = '["MAPK14;Y182;","MAPK14;T185;"]'
    site_key_to_display_id_json = (
        f'{{"{_MAPK14_Y182}":"MAPK14;Y182;","{_MAPK14_T185}":"MAPK14;T185;"}}'
    )
    return pd.DataFrame(
        {
            "protein_group_id": ["P28482"],
            "n_sites": [2],
            "site_ids": ['["MAPK14;Y182;","MAPK14;T185;"]'],
            "site_keys": [site_keys_json],
            "display_ids": [display_ids_json],
            "site_clusters": ["[1,2]"],
            "n_distinct_site_clusters": [2],
            "protein_module_id": [1],
            "multi_site_protein": [True],
            "ambiguous_module_context": [True],
            "gene_symbol": ["MAPK14"],
            "site": ["Y182,T185"],
            "protein_accession": ["P28482"],
            "isoform_id": [""],
            "top_kinases_by_site": ['{"MAPK14;Y182;":"MAP2K6","MAPK14;T185;":"MAPK1"}'],
            "module_ids_by_site": ['{"MAPK14;Y182;":1,"MAPK14;T185;":2}'],
            "site_key_to_display_id": [site_key_to_display_id_json],
        }
    )


def test_dataset_schema_valid_phospho_matrix_passes() -> None:
    wrapper = PhosphoIntensityMatrix(frame=_phospho_frame())
    assert wrapper.frame.shape == (2, 2)


def test_dataset_schema_non_numeric_phospho_column_fails() -> None:
    bad = _phospho_frame().astype(object)
    bad.loc[:, "sample_a"] = ["x", "y"]
    with pytest.raises(DatasetValidationError, match="numeric columns"):
        PhosphoIntensityMatrix(frame=bad)


def test_dataset_schema_all_boolean_phospho_matrix_fails() -> None:
    bad = pd.DataFrame(
        {
            "sample_a": [True, False],
            "sample_b": [False, True],
        },
        index=_phospho_frame().index.copy(),
    )
    with pytest.raises(DatasetValidationError) as exc_info:
        PhosphoIntensityMatrix(frame=bad)
    message = str(exc_info.value)
    assert "boolean columns are invalid" in message
    assert "sample_a" in message
    assert "sample_b" in message


def test_dataset_schema_mixed_float_and_boolean_columns_fail() -> None:
    index = _phospho_frame().index.copy()
    bad = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": pd.Series([True, False], index=index, dtype="bool"),
        },
        index=index,
    )
    with pytest.raises(DatasetValidationError) as exc_info:
        PhosphoIntensityMatrix(frame=bad)
    message = str(exc_info.value)
    assert "boolean columns are invalid" in message
    assert "sample_b" in message


def test_dataset_schema_nullable_boolean_phospho_column_fails() -> None:
    index = _phospho_frame().index.copy()
    bad = pd.DataFrame(
        {
            "sample_a": pd.Series([True, pd.NA], index=index, dtype="boolean"),
            "sample_b": [1.5, 2.5],
        },
        index=index,
    )
    with pytest.raises(DatasetValidationError) as exc_info:
        PhosphoIntensityMatrix(frame=bad)
    message = str(exc_info.value)
    assert "boolean columns are invalid" in message
    assert "sample_a" in message


def test_dataset_schema_valid_integer_phospho_matrix_passes() -> None:
    frame = pd.DataFrame(
        {
            "sample_a": [1, 2],
            "sample_b": [3, 4],
        },
        index=_phospho_frame().index.copy(),
    )
    wrapper = PhosphoIntensityMatrix(frame=frame)
    assert wrapper.frame.shape == (2, 2)


def test_dataset_schema_numeric_looking_string_column_fails() -> None:
    bad = _phospho_frame().astype(object)
    bad.loc[:, "sample_a"] = ["1.0", "2.0"]
    with pytest.raises(DatasetValidationError, match="non-numeric columns: sample_a"):
        PhosphoIntensityMatrix(frame=bad)


@given(site_ids=_canonical_site_ids(min_size=2))
@_PROPERTY_SETTINGS
def test_dataset_schema_property_non_numeric_string_columns_are_rejected(
    site_ids: list[str],
) -> None:
    site_keys = _site_keys_for_display_ids(site_ids)
    frame = pd.DataFrame(
        {
            "sample_a": [1.0 for _ in site_ids],
            "sample_b": [2.0 for _ in site_ids],
        },
        index=pd.Index(site_keys, name="site_key"),
    ).astype(object)
    frame.loc[:, "sample_a"] = ["x" for _ in site_ids]
    with pytest.raises(DatasetValidationError, match="non-numeric columns: sample_a"):
        PhosphoIntensityMatrix(frame=frame)


@given(site_ids=_canonical_site_ids(min_size=2))
@_PROPERTY_SETTINGS
def test_dataset_schema_property_boolean_columns_are_rejected(
    site_ids: list[str],
) -> None:
    index = pd.Index(_site_keys_for_display_ids(site_ids), name="site_key")
    frame = pd.DataFrame(
        {
            "sample_a": [1.0 for _ in site_ids],
            "sample_b": pd.Series(
                [(idx % 2) == 0 for idx in range(len(site_ids))],
                index=index,
                dtype="bool",
            ),
        },
        index=index,
    )
    with pytest.raises(DatasetValidationError, match="boolean columns are invalid"):
        PhosphoIntensityMatrix(frame=frame)


def test_dataset_schema_missing_phospho_value_fails() -> None:
    bad = _phospho_frame().copy(deep=True)
    bad.loc[_MAPK14_Y182, "sample_a"] = float("nan")
    with pytest.raises(DatasetValidationError, match="must not contain missing values"):
        PhosphoIntensityMatrix(frame=bad)


def test_dataset_schema_duplicate_phospho_index_fails() -> None:
    bad = _phospho_frame().copy(deep=True)
    bad.index = pd.Index([_MAPK14_Y182, _MAPK14_Y182], name="site_key")
    with pytest.raises(DatasetValidationError, match="duplicate_site_key"):
        PhosphoIntensityMatrix(frame=bad)


def test_dataset_schema_duplicate_phospho_columns_fail() -> None:
    bad = _phospho_frame().copy(deep=True)
    bad.columns = pd.Index(["sample_a", "sample_a"])
    with pytest.raises(DatasetValidationError, match="columns must be unique"):
        PhosphoIntensityMatrix(frame=bad)


@given(data=st.data())
@_PROPERTY_SETTINGS
def test_dataset_schema_property_duplicate_index_rejected_for_unique_site_ids(
    data: st.DataObject,
) -> None:
    site_ids = data.draw(_canonical_site_ids(min_size=2), label="site_ids")
    site_keys = _site_keys_for_display_ids(site_ids)
    duplicate_source = data.draw(st.sampled_from(site_keys), label="duplicate_source")
    duplicate_target_idx = data.draw(
        st.integers(min_value=0, max_value=len(site_ids) - 1),
        label="duplicate_target_idx",
    )
    duplicated_index = list(site_keys)
    duplicated_index[duplicate_target_idx] = duplicate_source
    assume(len(set(duplicated_index)) < len(duplicated_index))

    frame = pd.DataFrame(
        {
            "sample_a": [float(i + 1) for i in range(len(duplicated_index))],
            "sample_b": [float(i + 10) for i in range(len(duplicated_index))],
        },
        index=pd.Index(duplicated_index, name="site_key"),
    )
    with pytest.raises(
        DatasetValidationError,
        match="duplicate_site_key",
    ):
        PhosphoIntensityMatrix(frame=frame)


@given(data=st.data())
@_PROPERTY_SETTINGS
def test_prediction_schema_property_duplicate_columns_rejected_for_unique_kinase_labels(
    data: st.DataObject,
) -> None:
    kinase_labels = data.draw(
        st.lists(
            st.integers(min_value=1, max_value=1000).map(lambda i: f"K{i}"),
            min_size=2,
            max_size=6,
            unique=True,
        ),
        label="kinase_labels",
    )
    duplicate_source = data.draw(
        st.sampled_from(kinase_labels),
        label="duplicate_source",
    )
    duplicate_target_idx = data.draw(
        st.integers(min_value=0, max_value=len(kinase_labels) - 1),
        label="duplicate_target_idx",
    )
    duplicated_columns = list(kinase_labels)
    duplicated_columns[duplicate_target_idx] = duplicate_source
    assume(len(set(duplicated_columns)) < len(duplicated_columns))

    frame = pd.DataFrame(
        [list(range(1, len(duplicated_columns) + 1))],
        index=pd.Index([_MAPK14_Y182], name="site_key"),
        columns=pd.Index(duplicated_columns, name="kinase"),
        dtype=float,
    )
    with pytest.raises(
        PhosPyValidationError,
        match="prediction_result.pred_mat.columns must be unique",
    ):
        KinasePredictionMatrix(frame=frame)


def test_dataset_schema_non_canonical_site_index_fails() -> None:
    bad = _phospho_frame().copy(deep=True)
    bad.index = pd.Index([f"{_MAPK14_Y182} ", _AKT1_T308], name="site_key")
    with pytest.raises(DatasetValidationError, match="valid PhosPy site_key"):
        PhosphoIntensityMatrix(frame=bad)


def test_dataset_schema_site_metadata_missing_required_columns_fails() -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    bad = _site_metadata_frame(phospho.frame.index).drop(columns=["gene_symbol"])
    with pytest.raises(DatasetValidationError, match="missing required columns"):
        SiteMetadataTable(frame=bad, expected_index=phospho.frame.index)


def test_site_metadata_table_requires_site_sequence_column() -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    frame = _site_metadata_frame(phospho.frame.index).drop(columns=["site_sequence"])
    with pytest.raises(
        DatasetValidationError,
        match="dataset.site_metadata is missing required columns: site_sequence",
    ):
        SiteMetadataTable(frame=frame, expected_index=phospho.frame.index)


@given(data=st.data())
@_PROPERTY_SETTINGS
def test_site_metadata_property_required_columns_allow_extra_columns(
    data: st.DataObject,
) -> None:
    site_ids = data.draw(_canonical_site_ids(min_size=1), label="site_ids")
    extra_column_names = data.draw(
        st.lists(
            st.integers(min_value=1, max_value=50).map(lambda i: f"extra_{i}"),
            min_size=1,
            max_size=3,
            unique=True,
        ),
        label="extra_column_names",
    )
    parsed = [_parse_site_id(site_id) for site_id in site_ids]
    site_keys = _site_keys_for_display_ids(site_ids)
    frame = pd.DataFrame(
        {
            "site_key": site_keys,
            "display_id": site_ids,
            **_site_key_context_columns(site_keys),
            "gene_symbol": [gene for gene, _ in parsed],
            "site": [site for _, site in parsed],
            "site_sequence": [("A" * 15) + site[0] + ("A" * 15) for _, site in parsed],
        },
        index=pd.Index(site_keys, name="site_key"),
    )
    for idx, column_name in enumerate(extra_column_names):
        frame.loc[:, column_name] = idx

    wrapper = SiteMetadataTable(frame=frame, expected_index=frame.index)
    assert set(extra_column_names).issubset(set(wrapper.frame.columns))


@given(
    site_ids=_canonical_site_ids(min_size=1),
    missing_column=st.sampled_from(
        (
            "site_key",
            "display_id",
            "organism",
            "protein_namespace",
            "protein_identifier",
            "gene_symbol",
            "site",
            "site_sequence",
        )
    ),
)
@_PROPERTY_SETTINGS
def test_site_metadata_property_removing_required_column_fails(
    site_ids: list[str],
    missing_column: str,
) -> None:
    parsed = [_parse_site_id(site_id) for site_id in site_ids]
    site_keys = _site_keys_for_display_ids(site_ids)
    frame = pd.DataFrame(
        {
            "site_key": site_keys,
            "display_id": site_ids,
            **_site_key_context_columns(site_keys),
            "gene_symbol": [gene for gene, _ in parsed],
            "site": [site for _, site in parsed],
            "site_sequence": [("A" * 15) + site[0] + ("A" * 15) for _, site in parsed],
            "extra_col": ["x" for _ in site_ids],
        },
        index=pd.Index(site_keys, name="site_key"),
    ).drop(columns=[missing_column])
    with pytest.raises(DatasetValidationError, match="missing required columns"):
        SiteMetadataTable(frame=frame, expected_index=frame.index)


def test_dataset_schema_site_metadata_index_mismatch_fails() -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    bad = _site_metadata_frame(pd.Index([_MAPK14_Y182, _GSK3B_S9], name="site_key"))
    with pytest.raises(DatasetValidationError) as exc_info:
        SiteMetadataTable(frame=bad, expected_index=phospho.frame.index)
    message = str(exc_info.value)
    assert (
        "dataset.site_metadata.index must exactly match dataset.phospho.index"
        in message
    )
    assert "Only in dataset.site_metadata.index" in message
    assert "Only in dataset.phospho.index" in message
    assert "First positional mismatch: position 1" in message


def test_dataset_schema_site_metadata_identity_mismatch_fails() -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    bad = _site_metadata_frame(phospho.frame.index).copy(deep=True)
    bad.loc[_MAPK14_Y182, "gene_symbol"] = "MAPK1"
    with pytest.raises(DatasetValidationError, match="site-identity coherence failed"):
        SiteMetadataTable(frame=bad, expected_index=phospho.frame.index)


@pytest.mark.parametrize("value", [0.0, 0.75, 1.0, "0.5", pd.NA, None, "  "])
def test_site_metadata_localisation_probability_valid_values_pass(
    value: object,
) -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    frame = _site_metadata_frame(phospho.frame.index).copy(deep=True)
    frame.loc[:, "localisation_probability"] = [value, 0.8]
    wrapped = SiteMetadataTable(frame=frame, expected_index=phospho.frame.index)
    assert "localisation_probability" in wrapped.frame.columns


@pytest.mark.parametrize("value", [-0.1, 1.2, "high", "unknown", True, float("inf")])
def test_site_metadata_localisation_probability_invalid_values_fail(
    value: object,
) -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    frame = _site_metadata_frame(phospho.frame.index).copy(deep=True)
    frame.loc[:, "localisation_probability"] = [value, 0.8]
    with pytest.raises(
        DatasetValidationError,
        match="localisation_probability must contain values in \\[0.0, 1.0\\] or missing",
    ):
        SiteMetadataTable(frame=frame, expected_index=phospho.frame.index)


def test_site_metadata_rejects_malformed_site_tokens_by_default() -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    frame = _site_metadata_frame(phospho.frame.index).copy(deep=True)
    frame.loc[_MAPK14_Y182, "site"] = "Y18X"
    with pytest.raises(
        DatasetValidationError,
        match="site values must use strict 'S/T/Y<position>' tokens",
    ):
        SiteMetadataTable(frame=frame, expected_index=phospho.frame.index)


@pytest.mark.parametrize(
    "invalid_site_value",
    ["FOO", "A123", "S0", "S", "123", "", None],
)
def test_site_metadata_rejects_invalid_site_tokens_by_default(
    invalid_site_value: object,
) -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    frame = _site_metadata_frame(phospho.frame.index).copy(deep=True)
    frame.loc[_MAPK14_Y182, "site"] = invalid_site_value
    with pytest.raises(DatasetValidationError):
        SiteMetadataTable(frame=frame, expected_index=phospho.frame.index)


def test_site_metadata_rejects_opaque_site_tokens_at_site_key_boundary() -> None:
    site_key = _site_key_for_display_id("MAPK14;S1;", protein_identifier="P28482")
    phospho = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index([site_key], name="site_key"),
    )
    frame = pd.DataFrame(
        {
            "site_key": [site_key],
            "display_id": ["MAPK14;FOO;"],
            **_site_key_context_columns([site_key]),
            "gene_symbol": ["MAPK14"],
            "site": ["FOO"],
            "site_sequence": ["AAAAAYAAAAA"],
        },
        index=phospho.index.copy(),
    )
    with pytest.raises(DatasetValidationError, match="strict 'S/T/Y<position>'"):
        SiteMetadataTable(
            frame=frame,
            expected_index=phospho.index.copy(),
            allow_opaque_site_values=True,
        )


@pytest.mark.parametrize("site_token", ["S1", "T45", "Y999"])
def test_site_metadata_accepts_valid_sty_site_tokens(site_token: str) -> None:
    display_id = f"MAPK14;{site_token};"
    site_key = _site_key_for_display_id(display_id, protein_identifier="P28482")
    phospho = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index([site_key], name="site_key"),
    )
    frame = pd.DataFrame(
        {
            "site_key": [site_key],
            "display_id": [display_id],
            **_site_key_context_columns([site_key]),
            "gene_symbol": ["MAPK14"],
            "site": [site_token],
            "site_sequence": [("A" * 5) + site_token[0] + ("A" * 5)],
        },
        index=phospho.index.copy(),
    )
    wrapped = SiteMetadataTable(frame=frame, expected_index=phospho.index.copy())
    assert wrapped.frame.loc[site_key, "site"] == site_token


def test_site_metadata_rejects_residue_site_inconsistency() -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    frame = _site_metadata_frame(phospho.frame.index).copy(deep=True)
    frame.loc[:, "residue"] = ["S", "T"]
    with pytest.raises(
        DatasetValidationError,
        match="residue column must match parsed site residue",
    ):
        SiteMetadataTable(frame=frame, expected_index=phospho.frame.index)


def test_site_metadata_rejects_site_position_site_inconsistency() -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    frame = _site_metadata_frame(phospho.frame.index).copy(deep=True)
    frame.loc[:, "site_position"] = [123, 308]
    with pytest.raises(
        DatasetValidationError,
        match="site position column must match parsed site position",
    ):
        SiteMetadataTable(frame=frame, expected_index=phospho.frame.index)


def test_site_metadata_rejects_site_sequence_centre_residue_mismatch() -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    frame = _site_metadata_frame(phospho.frame.index).copy(deep=True)
    frame.loc[_MAPK14_Y182, "site_sequence"] = "AAAAASAAAAA"
    with pytest.raises(
        DatasetValidationError,
        match="site_sequence central residue must agree with site/residue metadata",
    ):
        SiteMetadataTable(frame=frame, expected_index=phospho.frame.index)


@pytest.mark.parametrize("centre_residue", ["S", "T", "Y"])
def test_site_metadata_accepts_phosphorylatable_site_sequence_centre(
    centre_residue: str,
) -> None:
    display_id = f"MAPK14;{centre_residue}182;"
    site_key = _site_key_for_display_id(display_id, protein_identifier="P28482")
    phospho = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index([site_key], name="site_key"),
    )
    frame = pd.DataFrame(
        {
            "site_key": [site_key],
            "display_id": [display_id],
            **_site_key_context_columns([site_key]),
            "gene_symbol": ["MAPK14"],
            "site": [f"{centre_residue}182"],
            "site_sequence": [f"AAAAA{centre_residue}AAAAA"],
        },
        index=phospho.index.copy(),
    )

    wrapped = SiteMetadataTable(frame=frame, expected_index=phospho.index.copy())
    assert (
        wrapped.frame.loc[phospho.index[0], "site_sequence"]
        == f"AAAAA{centre_residue}AAAAA"
    )


def test_site_metadata_rejects_non_phosphorylatable_site_sequence_centre() -> None:
    site_key = _site_key_for_display_id("GSK3B;S9;", protein_identifier="GSK3B")
    phospho = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index([site_key], name="site_key"),
    )
    frame = pd.DataFrame(
        {
            "site_key": [site_key],
            "display_id": ["GSK3B;S9;"],
            **_site_key_context_columns([site_key]),
            "gene_symbol": ["GSK3B"],
            "site": ["S9"],
            "site_sequence": ["AAAAAKAAAAA"],
        },
        index=phospho.index.copy(),
    )
    with pytest.raises(
        DatasetValidationError,
        match="must contain a centred phosphorylatable residue \\(S/T/Y\\)",
    ):
        SiteMetadataTable(frame=frame, expected_index=phospho.index.copy())


@pytest.mark.parametrize(
    ("site", "site_sequence"),
    [
        ("S10", "AAAAATAAAAA"),
        ("T10", "AAAAASAAAAA"),
        ("Y10", "AAAAASAAAAA"),
        ("Y10", "AAAAATAAAAA"),
    ],
)
def test_site_metadata_rejects_site_sequence_centre_residue_site_mismatch(
    site: str,
    site_sequence: str,
) -> None:
    display_id = f"GENE;{site};"
    site_key = _site_key_for_display_id(display_id, protein_identifier="GENE")
    phospho = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index([site_key], name="site_key"),
    )
    frame = pd.DataFrame(
        {
            "site_key": [site_key],
            "display_id": [display_id],
            **_site_key_context_columns([site_key]),
            "gene_symbol": ["GENE"],
            "site": [site],
            "site_sequence": [site_sequence],
        },
        index=phospho.index.copy(),
    )
    with pytest.raises(
        DatasetValidationError,
        match="site_sequence central residue must agree with site/residue metadata",
    ):
        SiteMetadataTable(frame=frame, expected_index=phospho.index.copy())


def test_site_metadata_accepts_lowercase_site_sequence_when_centre_is_valid() -> None:
    site_key = _site_key_for_display_id("MAPK14;Y182;", protein_identifier="P28482")
    phospho = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index([site_key], name="site_key"),
    )
    frame = pd.DataFrame(
        {
            "site_key": [site_key],
            "display_id": ["MAPK14;Y182;"],
            **_site_key_context_columns([site_key]),
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["aaaaayaaaaa"],
        },
        index=phospho.index.copy(),
    )

    wrapped = SiteMetadataTable(frame=frame, expected_index=phospho.index.copy())
    assert wrapped.frame.loc[site_key, "site_sequence"] == "aaaaayaaaaa"


@pytest.mark.parametrize("invalid_value", ["A1A", "AA AA", "**", "S"])
def test_site_metadata_rejects_implausible_site_sequence_values(
    invalid_value: str,
) -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    frame = _site_metadata_frame(phospho.frame.index).copy(deep=True)
    frame.loc[_MAPK14_Y182, "site_sequence"] = invalid_value
    with pytest.raises(
        DatasetValidationError,
        match="dataset.site_metadata.site_sequence must be plausible amino-acid context strings",
    ):
        SiteMetadataTable(frame=frame, expected_index=phospho.frame.index)


@pytest.mark.parametrize(
    "invalid_letter_sequence",
    ["AAAAABAAAAA", "AAAAAJAAAAA", "AAAAAOAAAAA", "AAAAAUAAAAA", "AAAAAZAAAAA"],
)
def test_site_metadata_rejects_non_policy_amino_acid_letters(
    invalid_letter_sequence: str,
) -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    frame = _site_metadata_frame(phospho.frame.index).copy(deep=True)
    frame.loc[_MAPK14_Y182, "site_sequence"] = invalid_letter_sequence
    with pytest.raises(
        DatasetValidationError,
        match="allowed residues: ACDEFGHIKLMNPQRSTVWY",
    ):
        SiteMetadataTable(frame=frame, expected_index=phospho.frame.index)


@pytest.mark.parametrize("blank_value", ["", "   "])
def test_site_metadata_rejects_blank_site_sequence_values(blank_value: str) -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    frame = _site_metadata_frame(phospho.frame.index).copy(deep=True)
    frame.loc[_MAPK14_Y182, "site_sequence"] = blank_value
    with pytest.raises(
        DatasetValidationError,
        match="site_sequence must contain non-empty string values",
    ):
        SiteMetadataTable(frame=frame, expected_index=phospho.frame.index)


@pytest.mark.parametrize(
    "sequence_value", ["AAAAXYAAAAA", "AAAA-YAAAAA", "AAAA_YAAAAA"]
)
def test_site_metadata_accepts_explicit_unknown_and_gap_sequence_policy(
    sequence_value: str,
) -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    frame = _site_metadata_frame(phospho.frame.index).copy(deep=True)
    frame.loc[_MAPK14_Y182, "site_sequence"] = sequence_value
    wrapped = SiteMetadataTable(frame=frame, expected_index=phospho.frame.index)
    assert wrapped.frame.loc[_MAPK14_Y182, "site_sequence"] == sequence_value


def test_dataset_schema_sample_metadata_index_mismatch_fails() -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    sample_metadata = pd.DataFrame(
        {"group": ["g1", "g2"]},
        index=pd.Index(["sample_a", "sample_x"], name="sample_id"),
    )
    with pytest.raises(DatasetValidationError, match="must exactly match"):
        SampleMetadataTable(frame=sample_metadata, expected_index=phospho.frame.columns)


def test_dataset_schema_total_matrix_column_mismatch_fails() -> None:
    phospho = PhosphoIntensityMatrix(frame=_phospho_frame())
    total = pd.DataFrame(
        {"sample_a": [1.0], "sample_x": [2.0]},
        index=pd.Index(["P28482"], name="protein_id"),
    )
    with pytest.raises(DatasetValidationError) as exc_info:
        TotalProteinMatrix(frame=total, expected_sample_index=phospho.frame.columns)
    message = str(exc_info.value)
    assert "dataset.total.columns must exactly match dataset.phospho.columns" in message
    assert "Only in dataset.total.columns: 'sample_x'" in message
    assert "Only in dataset.phospho.columns: 'sample_b'" in message
    assert "First positional mismatch: position 1" in message


def test_reference_schema_valid_kinase_substrate_reference_passes() -> None:
    frame = pd.DataFrame(
        {
            "kinase": ["MAP2K6", "AKT1"],
            "substrate_site": ["MAPK14;Y182;", "AKT1;T308;"],
        }
    )
    wrapper = KinaseSubstrateReference(frame=frame)
    assert wrapper.frame.shape == (2, 2)


def test_reference_schema_missing_required_column_fails() -> None:
    with pytest.raises(ReferenceValidationError, match="missing required columns"):
        KinaseSubstrateReference(
            frame=pd.DataFrame({"kinase": ["MAP2K6"]}),
        )


def test_reference_schema_duplicate_pairs_fail() -> None:
    frame = pd.DataFrame(
        {
            "kinase": ["MAP2K6", "MAP2K6"],
            "substrate_site": ["MAPK14;Y182;", "MAPK14;Y182;"],
        }
    )
    with pytest.raises(ReferenceValidationError, match="duplicate"):
        KinaseSubstrateReference(frame=frame)


def test_reference_schema_non_canonical_kinase_normalizes() -> None:
    frame = pd.DataFrame(
        {
            "kinase": [" MAP2K6 "],
            "substrate_site": ["MAPK14;Y182;"],
        }
    )
    wrapper = KinaseSubstrateReference(frame=frame)
    assert wrapper.frame.loc[0, "kinase"] == "MAP2K6"


def test_reference_schema_non_canonical_substrate_site_normalizes() -> None:
    frame = pd.DataFrame(
        {
            "kinase": ["MAP2K6"],
            "substrate_site": [" mapk14 ; y182 "],
        }
    )
    wrapper = KinaseSubstrateReference(frame=frame)
    assert wrapper.frame.loc[0, "substrate_site"] == "MAPK14;Y182;"


def test_reference_schema_missing_site_sequence_for_substrate_fails() -> None:
    with pytest.raises(ReferenceValidationError, match="missing sequence entries"):
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {
                    "kinase": ["MAP2K6"],
                    "substrate_site": ["MAPK14;Y182;"],
                }
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["A" * 31]},
                index=pd.Index(["AKT1;T308;"], name="site_id"),
            ),
        )


def test_prediction_schema_valid_prediction_matrix_passes() -> None:
    pred_mat = pd.DataFrame(
        {"MAP2K6": [0.9, 0.8]},
        index=pd.Index([_MAPK14_Y182, _AKT1_T308], name="site_key"),
    )
    wrapper = KinasePredictionMatrix(frame=pred_mat)
    assert wrapper.frame.shape == (2, 1)


def test_prediction_schema_allows_nan_when_missing_policy_allows() -> None:
    pred_mat = pd.DataFrame(
        {"MAP2K6": [0.9, float("nan")]},
        index=pd.Index([_MAPK14_Y182, _AKT1_T308], name="site_key"),
    )
    wrapper = KinasePredictionMatrix(frame=pred_mat)
    assert pd.isna(wrapper.frame.loc[_AKT1_T308, "MAP2K6"])


@pytest.mark.parametrize("invalid_value", [float("inf"), float("-inf")])
def test_prediction_schema_rejects_infinite_values_when_missing_allowed(
    invalid_value: float,
) -> None:
    pred_mat = pd.DataFrame(
        {"MAP2K6": [0.9, invalid_value]},
        index=pd.Index([_MAPK14_Y182, _AKT1_T308], name="site_key"),
    )
    with pytest.raises(PhosPyValidationError, match="finite numeric values"):
        KinasePredictionMatrix(frame=pred_mat)


@given(
    site_ids=_canonical_site_ids(min_size=2),
    invalid_value=st.sampled_from((float("inf"), float("-inf"))),
)
@_PROPERTY_SETTINGS
def test_prediction_schema_property_rejects_non_finite_numeric_entries(
    site_ids: list[str],
    invalid_value: float,
) -> None:
    site_keys = _site_keys_for_display_ids(site_ids)
    frame = pd.DataFrame(
        {
            "K1": [0.9 for _ in site_ids],
            "K2": [0.1 for _ in site_ids],
        },
        index=pd.Index(site_keys, name="site_key"),
        dtype=float,
    )
    frame.iloc[0, 1] = invalid_value
    with pytest.raises(PhosPyValidationError, match="finite numeric values"):
        KinasePredictionMatrix(frame=frame)


def test_prediction_schema_duplicate_kinase_columns_fail() -> None:
    pred_mat = pd.DataFrame(
        [[0.9, 0.8]],
        index=pd.Index([_MAPK14_Y182], name="site_key"),
        columns=pd.Index(["MAP2K6", "MAP2K6"], name="kinase"),
    )
    with pytest.raises(PhosPyValidationError, match="columns must be unique"):
        KinasePredictionMatrix(frame=pred_mat)


def test_prediction_schema_display_site_index_fails() -> None:
    pred_mat = pd.DataFrame(
        {"MAP2K6": [0.9]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    with pytest.raises(PhosPyValidationError, match="must start with 'phospy:v1'"):
        KinasePredictionMatrix(frame=pred_mat)


def test_kinase_score_schema_display_site_index_fails() -> None:
    score_matrix = pd.DataFrame(
        {"MAP2K6": [0.9]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    with pytest.raises(PhosPyValidationError, match="must start with 'phospy:v1'"):
        KinaseScoreMatrix(frame=score_matrix)


def test_prediction_schema_non_canonical_site_key_index_fails() -> None:
    pred_mat = pd.DataFrame(
        {"MAP2K6": [0.9]},
        index=pd.Index([f"{_MAPK14_Y182} "], name="site_key"),
    )
    with pytest.raises(PhosPyValidationError, match="analysis-ready encoded site_key"):
        KinasePredictionMatrix(frame=pred_mat)


def test_prediction_schema_out_of_range_score_fails() -> None:
    pred_mat = pd.DataFrame(
        {"MAP2K6": [1.2]},
        index=pd.Index([_MAPK14_Y182], name="site_key"),
    )
    with pytest.raises(PhosPyValidationError, match="between 0.0 and 1.0"):
        KinasePredictionMatrix(frame=pred_mat)


def test_prediction_schema_boolean_kinase_score_column_fails() -> None:
    pred_mat = pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.8],
            "AKT1": [True, False],
        },
        index=pd.Index([_MAPK14_Y182, _AKT1_T308], name="site_key"),
    )
    with pytest.raises(PhosPyValidationError) as exc_info:
        KinasePredictionMatrix(frame=pred_mat)
    message = str(exc_info.value)
    assert "boolean columns are invalid" in message
    assert "AKT1" in message


def test_activity_schema_non_numeric_matrix_fails() -> None:
    bad = pd.DataFrame(
        {"sample_a": ["x"]},
        index=pd.Index(["MAP2K6"], name="kinase"),
    )
    with pytest.raises(PhosPyValidationError, match="numeric columns"):
        ActivityMatrix(frame=bad)


def test_activity_schema_target_table_missing_required_columns_fails() -> None:
    with pytest.raises(PhosPyValidationError, match="missing required columns"):
        ActivityTargetTable(
            frame=pd.DataFrame({"site_id": ["MAPK14;Y182;"], "score": [0.8]}),
        )


def test_signalome_schema_valid_site_membership_table_passes() -> None:
    wrapper = SignalomeSiteContext(frame=_site_membership_frame())
    assert wrapper.frame.shape[0] == 2


def test_signalome_schema_missing_required_site_membership_column_fails() -> None:
    bad = _site_membership_frame().drop(columns=["top_kinase_weight"])
    with pytest.raises(WorkflowValidationError, match="missing required columns"):
        SignalomeSiteContext(frame=bad)


def test_signalome_schema_invalid_site_id_fails() -> None:
    bad = _site_membership_frame().copy(deep=True)
    bad.loc[0, "display_id"] = "MAPK14;Y182; "
    with pytest.raises(
        WorkflowValidationError,
        match="recommended site identifier format",
    ):
        SignalomeSiteContext(frame=bad)


def test_signalome_schema_invalid_boolean_integer_numeric_columns_fail() -> None:
    bad = _site_membership_frame().copy(deep=True)
    bad = bad.astype({"included_in_module_table": object})
    bad.loc[0, "included_in_module_table"] = "yes"
    with pytest.raises(WorkflowValidationError, match="boolean"):
        SignalomeSiteContext(frame=bad)


def test_signalome_schema_valid_protein_site_context_passes() -> None:
    wrapper = SignalomeProteinSiteContext(frame=_protein_site_context_frame())
    assert wrapper.frame.shape[0] == 1


def test_signalome_schema_malformed_json_columns_fail() -> None:
    bad = _protein_site_context_frame().copy(deep=True)
    bad.loc[0, "site_ids"] = "[invalid"
    with pytest.raises(WorkflowValidationError, match="parseable JSON"):
        SignalomeProteinSiteContext(frame=bad)


def test_reference_site_sequence_wrapper_missing_required_column_fails() -> None:
    with pytest.raises(ReferenceValidationError, match="missing required columns"):
        SiteSequenceReference(
            frame=pd.DataFrame(
                {"sequence": ["A" * 31]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
        )
