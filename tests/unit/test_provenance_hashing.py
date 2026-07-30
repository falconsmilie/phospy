from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from datetime import date, datetime, time
from decimal import Decimal

import numpy as np
import pandas as pd
import pytest

from phospy.errors.input import PhosPyInputError
from phospy.errors.provenance import ProvenanceFingerprintError
from phospy.provenance.hashing import (
    _index_structure,
    _normalize_value,
    _update,
    fingerprint_optional_table_normalized_axes,
    fingerprint_table,
    hash_table_exact,
    hash_table_tolerance,
)
from phospy.provenance.models import (
    ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2,
    PREPROCESSING_EXTERNAL_NONDETERMINISM_CAVEAT_CODE,
    PREPROCESSING_STAGE_DETERMINISM_PURE,
    PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3,
    DeterminismKind,
    EnvironmentProvenance,
    PreprocessingStageProvenance,
    ReproducibilityCaveat,
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


def test_normalized_axis_fingerprint_is_stable_for_row_and_column_order() -> None:
    table = pd.DataFrame(
        {
            "sample_b": [4.0, 3.0],
            "sample_a": [2.0, 1.0],
        },
        index=pd.Index(["B;S2;", "A;S1;"], name="site_id"),
    )
    reordered = table.loc[["A;S1;", "B;S2;"], ["sample_a", "sample_b"]]

    first = fingerprint_optional_table_normalized_axes(
        table,
        name="dataset.phospho",
    )
    second = fingerprint_optional_table_normalized_axes(
        reordered,
        name="dataset.phospho",
    )

    assert first is not None
    assert second is not None
    assert first.exact_hash_value == second.exact_hash_value
    assert first.tolerance_hash_value == second.tolerance_hash_value
    assert first.column_names == ("sample_a", "sample_b")
    assert second.column_names == ("sample_a", "sample_b")


def test_normalized_axis_fingerprint_sorts_integer_labels_numerically() -> None:
    table = pd.DataFrame(
        [[40.0, 30.0], [20.0, 10.0]],
        index=pd.Index([10, 2], name="row_id"),
        columns=pd.Index([10, 2], name="sample_id"),
    )
    reordered = table.iloc[[1, 0], [1, 0]]

    first = fingerprint_optional_table_normalized_axes(table, name="integer_labels")
    second = fingerprint_optional_table_normalized_axes(
        reordered, name="integer_labels"
    )

    assert first is not None
    assert second is not None
    assert first.exact_hash_value == second.exact_hash_value
    assert first.tolerance_hash_value == second.tolerance_hash_value
    assert first.column_names == ("2", "10")


def test_normalized_axis_fingerprint_distinguishes_mixed_int_string_labels() -> None:
    table = pd.DataFrame(
        [[4.0, 3.0], [2.0, 1.0]],
        index=pd.Index(["1", 1], dtype=object, name="row_id"),
        columns=pd.Index(["2", 2], dtype=object, name="sample_id"),
    )
    reordered = table.iloc[[1, 0], [1, 0]]

    first = fingerprint_optional_table_normalized_axes(table, name="mixed_labels")
    second = fingerprint_optional_table_normalized_axes(reordered, name="mixed_labels")

    assert first is not None
    assert second is not None
    assert first.exact_hash_value == second.exact_hash_value
    assert first.tolerance_hash_value == second.tolerance_hash_value
    assert first.index_structure is not None
    assert first.index_structure["values"] == (
        {"kind": "int", "value": 1},
        {"kind": "str", "value": "1"},
    )
    assert first.column_index_structure is not None
    assert first.column_index_structure["values"] == (
        {"kind": "int", "value": 2},
        {"kind": "str", "value": "2"},
    )


def test_normalized_axis_fingerprint_supports_tuple_labels_and_multiindex() -> None:
    tuple_index = pd.MultiIndex.from_tuples(
        [("B", 2), ("A", 1)],
        names=["group", "replicate"],
    )
    tuple_columns = pd.MultiIndex.from_tuples(
        [("sample", 2), ("sample", 1)],
        names=["kind", "replicate"],
    )
    table = pd.DataFrame(
        [[4.0, 3.0], [2.0, 1.0]],
        index=tuple_index,
        columns=tuple_columns,
    )
    reordered = table.iloc[[1, 0], [1, 0]]

    first = fingerprint_optional_table_normalized_axes(table, name="tuple_labels")
    second = fingerprint_optional_table_normalized_axes(reordered, name="tuple_labels")

    assert first is not None
    assert second is not None
    assert first.exact_hash_value == second.exact_hash_value
    assert first.tolerance_hash_value == second.tolerance_hash_value
    assert first.index_structure is not None
    assert first.index_structure["type"] == "multi_index"
    assert first.column_index_structure is not None
    assert first.column_index_structure["type"] == "multi_index"


def test_normalized_axis_fingerprint_uses_same_order_for_exact_and_tolerance() -> None:
    table = pd.DataFrame(
        {"b": [1.123456781, 2.0], "a": [3.0, 4.0]},
        index=pd.Index(["row_b", "row_a"], name="row_id"),
    )
    expected_order = table.loc[["row_a", "row_b"], ["a", "b"]]

    fingerprint = fingerprint_optional_table_normalized_axes(
        table,
        name="same_axis_order",
    )

    assert fingerprint is not None
    assert fingerprint.exact_hash_value == hash_table_exact(
        expected_order,
        name="same_axis_order",
    )
    assert fingerprint.tolerance_hash_value == hash_table_tolerance(
        expected_order,
        name="same_axis_order",
    )


def test_normalized_axis_fingerprint_rejects_unsupported_label_type() -> None:
    table = pd.DataFrame(
        {"x": [1.0]},
        index=pd.Index([object()], dtype=object, name="row_id"),
    )

    with pytest.raises(ProvenanceFingerprintError) as exc_info:
        fingerprint_optional_table_normalized_axes(table, name="unsupported")

    assert str(exc_info.value) == (
        "normalized provenance fingerprint for table 'unsupported' cannot sort "
        "row axis label at position 0: unsupported axis label type object. "
        "Supported labels under policy typed-axis-label-sort-v1 are non-missing "
        "strings, integers, and tuple/MultiIndex labels composed only of supported "
        "labels. Convert labels to a supported, collision-safe representation "
        "before fingerprinting."
    )


def test_normalized_axis_fingerprint_rejects_duplicate_canonical_labels() -> None:
    table = pd.DataFrame(
        {"x": [1.0, 2.0]},
        index=pd.Index(["row", "row"], name="row_id"),
    )

    with pytest.raises(
        ProvenanceFingerprintError,
        match="positions 0 and 1 share the same canonical typed key",
    ):
        fingerprint_optional_table_normalized_axes(table, name="duplicate_labels")


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


def test_external_nondeterminism_caveat_round_trips_in_stage_provenance() -> None:
    provenance = _current_run_provenance()
    caveat = ReproducibilityCaveat(
        code=PREPROCESSING_EXTERNAL_NONDETERMINISM_CAVEAT_CODE,
        severity="warning",
        message="External preprocessing execution requires external state.",
        details={
            "stage": "missing_data",
            "determinism_kind": "externally_nondeterministic",
        },
    )
    stage = replace(
        provenance.preprocessing_stages[0],
        determinism=DeterminismKind.EXTERNALLY_NONDETERMINISTIC,
        reproducibility_caveats=(caveat,),
    )
    payload = to_payload(replace(provenance, preprocessing_stages=(stage,)))
    stage_payload = payload["preprocessing_stages"][0]
    assert isinstance(stage_payload, dict)
    assert stage_payload["determinism"] == "externally_nondeterministic"
    caveat_payloads = stage_payload["reproducibility_caveats"]
    assert isinstance(caveat_payloads, list)
    assert (
        caveat_payloads[0]["code"] == PREPROCESSING_EXTERNAL_NONDETERMINISM_CAVEAT_CODE
    )

    restored = from_payload(payload)

    restored_stage = restored.preprocessing_stages[0]
    assert restored_stage.determinism is DeterminismKind.EXTERNALLY_NONDETERMINISTIC
    assert restored_stage.reproducibility_caveats == (caveat,)


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


def test_table_hash_fast_scalar_path_matches_stable_json_reference() -> None:
    table = pd.DataFrame(
        {
            "float": [1.123456789, np.nan, np.inf, -np.inf],
            "integer": np.asarray([1, 2, 3, 4], dtype=np.int64),
            "text": ["plain", 'needs"escaping', "unicode-μ", "slash\\value"],
            "object": [
                Decimal("1.2300"),
                datetime(2026, 7, 27, 12, 30, 15),
                date(2026, 7, 27),
                (time(12, 30), {"nested": [1.123456789, "value"]}),
            ],
        },
        index=pd.Index(["A;S1;", "B;T2;", "C;Y3;", "D;S4;"], name="site_id"),
    )

    assert hash_table_exact(table, name="mixed") == _reference_table_hash(
        table,
        name="mixed",
        round_floats=False,
    )
    assert hash_table_tolerance(table, name="mixed") == _reference_table_hash(
        table,
        name="mixed",
        round_floats=True,
    )


def _reference_table_hash(
    table: pd.DataFrame,
    *,
    name: str,
    round_floats: bool,
) -> str:
    hasher = hashlib.sha256()
    _update(hasher, name)
    _update(hasher, [int(table.shape[0]), int(table.shape[1])])
    _update(hasher, _index_structure(table.index, round_floats=round_floats))
    _update(hasher, _index_structure(table.columns, round_floats=round_floats))
    _update(hasher, [str(dtype) for dtype in table.dtypes.tolist()])
    values = table.to_numpy(dtype=object, copy=False)
    for row in values:
        for value in row:
            encoded = json.dumps(
                _normalize_value(value, round_floats=round_floats),
                sort_keys=True,
                separators=(",", ":"),
                ensure_ascii=True,
                allow_nan=False,
            ).encode("utf-8")
            hasher.update(encoded)
            hasher.update(b"\n")
    return hasher.hexdigest()


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
