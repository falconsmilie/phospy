from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.errors.input import PhosPyInputError
from phospy.provenance.hashing import (
    fingerprint_table,
    hash_table_exact,
    hash_table_tolerance,
)
from phospy.provenance.models import (
    ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2,
    PREPROCESSING_STAGE_DETERMINISM_PURE,
    PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3,
    EnvironmentProvenance,
    PreprocessingStageProvenance,
    RunProvenance,
)
from phospy.provenance.serialization import from_payload, to_payload


def _base_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, np.nan, 3.0],
            "sample_b": [4.0, 5.0, 6.0],
        },
        index=pd.Index(["A;S1;", "B;S2;", "C;S3;"], name="site_id"),
    )


def _current_run_provenance() -> RunProvenance:
    fingerprint = fingerprint_table(_base_table(), name="dataset.phospho")
    stage = PreprocessingStageProvenance(
        stage="missing_data",
        operation="drop_rows",
        parameters={"policy": "drop_any_missing"},
        input_shape=(3, 2),
        output_shape=(2, 2),
        input_hash=fingerprint.tolerance_hash_value,
        output_hash=fingerprint.tolerance_hash_value,
        phospho_input_hash=fingerprint.tolerance_hash_value,
        phospho_output_hash=fingerprint.tolerance_hash_value,
        dropped_row_ids=("B;S2;",),
        dropped_row_count=1,
        schema_version=PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3,
        consumed_input_tables=(fingerprint,),
        produced_output_tables=(fingerprint,),
        determinism=PREPROCESSING_STAGE_DETERMINISM_PURE,
    )
    return RunProvenance(
        environment=EnvironmentProvenance(
            package_name="phospy",
            package_version="0",
            python_version="3.13",
            dependency_versions={"pandas": "2"},
            schema_version=ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2,
        ),
        input_tables=(fingerprint,),
        preprocessing_stages=(stage,),
        reference=None,
        workflow_name="test_workflow",
        workflow_parameters={"mode": "test"},
        random_state=None,
        random_seed_policy=None,
        output_tables=(fingerprint,),
        scientific_policies=(),
    )


def _collect_keys(value: object) -> set[str]:
    if isinstance(value, dict):
        keys = {str(key) for key in value}
        for item in value.values():
            keys.update(_collect_keys(item))
        return keys
    if isinstance(value, list):
        keys: set[str] = set()
        for item in value:
            keys.update(_collect_keys(item))
        return keys
    return set()


def _current_payload() -> dict[str, object]:
    return to_payload(_current_run_provenance())


def test_hash_is_deterministic_for_identical_table() -> None:
    table = _base_table()
    assert hash_table_tolerance(table, name="dataset.phospho") == hash_table_tolerance(
        table.copy(deep=True),
        name="dataset.phospho",
    )


def test_hash_changes_when_row_order_changes() -> None:
    table = _base_table()
    reordered = table.iloc[[2, 1, 0], :]
    assert hash_table_tolerance(table, name="dataset.phospho") != hash_table_tolerance(
        reordered,
        name="dataset.phospho",
    )


def test_hash_changes_when_column_order_changes() -> None:
    table = _base_table()
    reordered = table.loc[:, ["sample_b", "sample_a"]]
    assert hash_table_tolerance(table, name="dataset.phospho") != hash_table_tolerance(
        reordered,
        name="dataset.phospho",
    )


def test_hash_changes_when_value_changes() -> None:
    table = _base_table()
    changed = table.copy(deep=True)
    changed.loc["A;S1;", "sample_a"] = 9.0
    assert hash_table_tolerance(table, name="dataset.phospho") != hash_table_tolerance(
        changed,
        name="dataset.phospho",
    )


def test_hash_changes_when_dtype_changes() -> None:
    table = _base_table()
    changed = table.astype({"sample_b": "int64"})
    assert hash_table_tolerance(table, name="dataset.phospho") != hash_table_tolerance(
        changed,
        name="dataset.phospho",
    )


def test_hash_distinguishes_numeric_and_string_column_labels() -> None:
    numeric_columns = pd.DataFrame(
        {1: [1.0, 2.0]},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    string_columns = pd.DataFrame(
        {"1": [1.0, 2.0]},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    assert hash_table_tolerance(numeric_columns, name="table") != hash_table_tolerance(
        string_columns, name="table"
    )


def test_hash_distinguishes_numeric_and_string_index_labels() -> None:
    numeric_index = pd.DataFrame(
        {"x": [1.0, 2.0]},
        index=pd.Index([1, 2], name="row_id"),
    )
    string_index = pd.DataFrame(
        {"x": [1.0, 2.0]},
        index=pd.Index(["1", "2"], name="row_id"),
    )
    assert hash_table_tolerance(numeric_index, name="table") != hash_table_tolerance(
        string_index, name="table"
    )


def test_hash_changes_when_row_index_name_changes() -> None:
    first = pd.DataFrame(
        {"x": [1.0, 2.0]},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    second = first.copy(deep=True)
    second.index = second.index.rename("site_id")
    assert hash_table_tolerance(first, name="table") != hash_table_tolerance(
        second, name="table"
    )


def test_hash_changes_when_column_index_name_changes() -> None:
    first = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["A"], name="row_id"),
    )
    first.columns = first.columns.rename("sample_id")
    second = first.copy(deep=True)
    second.columns = second.columns.rename("run_id")
    assert hash_table_tolerance(first, name="table") != hash_table_tolerance(
        second, name="table"
    )


def test_hash_distinguishes_range_index_from_equivalent_integer_index() -> None:
    range_index_table = pd.DataFrame({"x": [1.0, 2.0, 3.0]})
    integer_index_table = range_index_table.copy(deep=True)
    integer_index_table.index = pd.Index([0, 1, 2], dtype="int64")
    assert hash_table_tolerance(
        range_index_table, name="table"
    ) != hash_table_tolerance(integer_index_table, name="table")


def test_hash_supports_multiindex_rows() -> None:
    first = pd.DataFrame(
        {"x": [1.0, 2.0]},
        index=pd.MultiIndex.from_tuples(
            [(1, "A"), (2, "B")], names=["batch", "replicate"]
        ),
    )
    second = first.copy(deep=True)
    second.index = pd.MultiIndex.from_tuples(
        [("1", "A"), (2, "B")], names=["batch", "replicate"]
    )
    assert hash_table_tolerance(first, name="table") != hash_table_tolerance(
        second, name="table"
    )


def test_hash_supports_multiindex_columns() -> None:
    columns_first = pd.MultiIndex.from_tuples(
        [(1, "treated"), ("1", "control")], names=["sample_id", "group"]
    )
    columns_second = pd.MultiIndex.from_tuples(
        [("1", "treated"), ("1", "control")], names=["sample_id", "group"]
    )
    first = pd.DataFrame(
        [[1.0, 2.0], [3.0, 4.0]],
        index=pd.Index(["A", "B"], name="row_id"),
        columns=columns_first,
    )
    second = pd.DataFrame(
        [[1.0, 2.0], [3.0, 4.0]],
        index=pd.Index(["A", "B"], name="row_id"),
        columns=columns_second,
    )
    assert hash_table_tolerance(first, name="table") != hash_table_tolerance(
        second, name="table"
    )


def test_hash_changes_when_display_is_identical_but_dtype_differs() -> None:
    int64_table = pd.DataFrame(
        {"x": pd.Series([1, 2], dtype="int64")},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    nullable_int_table = pd.DataFrame(
        {"x": pd.Series([1, 2], dtype="Int64")},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    assert hash_table_tolerance(int64_table, name="table") != hash_table_tolerance(
        nullable_int_table, name="table"
    )


def test_missing_value_representation_is_stable() -> None:
    first = pd.DataFrame(
        {"x": [1.0, np.nan]},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    second = pd.DataFrame(
        {"x": [1.0, None]},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    assert hash_table_tolerance(first, name="table") == hash_table_tolerance(
        second, name="table"
    )


def test_fingerprint_captures_structural_metadata() -> None:
    table = _base_table()
    fingerprint = fingerprint_table(table, name="dataset.phospho")
    assert fingerprint.name == "dataset.phospho"
    assert fingerprint.rows == 3
    assert fingerprint.columns == 2
    assert fingerprint.index_name == "site_id"
    assert fingerprint.column_names == ("sample_a", "sample_b")
    assert fingerprint.index_structure is not None
    assert fingerprint.index_structure["type"] == "index"
    assert fingerprint.column_index_structure is not None
    assert fingerprint.column_index_structure["type"] == "index"
    assert fingerprint.exact_hash_algorithm == "sha256-stable-json-v1"
    assert isinstance(fingerprint.exact_hash_value, str)
    assert len(fingerprint.exact_hash_value) == 64
    assert fingerprint.tolerance_hash_algorithm == "sha256-float-round-8dp-v1"
    assert isinstance(fingerprint.tolerance_hash_value, str)
    assert len(fingerprint.tolerance_hash_value) == 64


def test_current_provenance_payload_omits_legacy_hash_and_determinism_aliases() -> None:
    payload = _current_payload()

    keys = _collect_keys(payload)

    assert "hash_algorithm" not in keys
    assert "hash_value" not in keys
    assert "is_deterministic" not in keys


@pytest.mark.parametrize(
    "mutation",
    [
        pytest.param("legacy_table_hash_fields", id="legacy-table-hash-fields"),
        pytest.param("legacy_stage_determinism", id="legacy-stage-determinism"),
        pytest.param("legacy_environment_schema", id="legacy-environment-schema"),
        pytest.param("legacy_stage_schema", id="legacy-stage-schema"),
    ],
)
def test_legacy_provenance_payloads_are_rejected_clearly(mutation: str) -> None:
    payload = _current_payload()
    if mutation == "legacy_table_hash_fields":
        table_payload = payload["input_tables"][0]
        assert isinstance(table_payload, dict)
        table_payload["hash_algorithm"] = "sha256"
        table_payload["hash_value"] = table_payload["tolerance_hash_value"]
    elif mutation == "legacy_stage_determinism":
        stage_payload = payload["preprocessing_stages"][0]
        assert isinstance(stage_payload, dict)
        stage_payload["is_deterministic"] = True
    elif mutation == "legacy_environment_schema":
        environment_payload = payload["environment"]
        assert isinstance(environment_payload, dict)
        environment_payload["schema_version"] = 1
    elif mutation == "legacy_stage_schema":
        stage_payload = payload["preprocessing_stages"][0]
        assert isinstance(stage_payload, dict)
        stage_payload["schema_version"] = 2
    else:
        raise AssertionError(f"Unknown mutation: {mutation}")

    with pytest.raises(
        PhosPyInputError,
        match=(
            "Legacy provenance schemas are no longer supported. "
            "Regenerate the result with the current PhosPy version."
        ),
    ):
        from_payload(payload)


def test_exact_hash_changes_for_sub_8dp_float_differences() -> None:
    first = pd.DataFrame(
        {"sample_a": [1.123456781]},
        index=pd.Index(["A;S1;"], name="site_id"),
    )
    second = pd.DataFrame(
        {"sample_a": [1.123456784]},
        index=pd.Index(["A;S1;"], name="site_id"),
    )
    assert hash_table_exact(first, name="dataset.phospho") != hash_table_exact(
        second,
        name="dataset.phospho",
    )


def test_tolerance_hash_can_stay_stable_for_sub_8dp_float_differences() -> None:
    first = pd.DataFrame(
        {"sample_a": [1.123456781]},
        index=pd.Index(["A;S1;"], name="site_id"),
    )
    second = pd.DataFrame(
        {"sample_a": [1.123456784]},
        index=pd.Index(["A;S1;"], name="site_id"),
    )
    assert hash_table_tolerance(first, name="dataset.phospho") == hash_table_tolerance(
        second,
        name="dataset.phospho",
    )


def test_exact_hash_changes_when_row_order_changes() -> None:
    table = _base_table()
    reordered = table.iloc[[2, 1, 0], :]
    assert hash_table_exact(table, name="dataset.phospho") != hash_table_exact(
        reordered,
        name="dataset.phospho",
    )


def test_exact_hash_changes_when_column_order_changes() -> None:
    table = _base_table()
    reordered = table.loc[:, ["sample_b", "sample_a"]]
    assert hash_table_exact(table, name="dataset.phospho") != hash_table_exact(
        reordered,
        name="dataset.phospho",
    )


def test_exact_hash_changes_when_index_labels_change() -> None:
    first = pd.DataFrame(
        {"x": [1.0, 2.0]},
        index=pd.Index(["A", "B"], name="row_id"),
    )
    second = first.copy(deep=True)
    second.index = pd.Index(["A", "B_CHANGED"], name="row_id")
    assert hash_table_exact(first, name="table") != hash_table_exact(
        second,
        name="table",
    )
