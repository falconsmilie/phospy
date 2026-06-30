"""Shared private table parsing helpers for phosphosite readers."""

from __future__ import annotations

import math
import re
from collections.abc import Callable, Mapping, Sequence

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.validation.datasets.importers import normalise_sample_column_mapping

_MULTI_VALUE_SPLIT_PATTERN = re.compile(r"\s*[,;]\s*")


def require_non_empty_unique_columns(
    source: pd.DataFrame,
    *,
    importer_label: str,
) -> None:
    if not isinstance(source, pd.DataFrame):
        raise PhosPyInputError(
            f"{importer_label} import source must be a pandas DataFrame"
        )
    if source.empty:
        raise PhosPyInputError(f"{importer_label} import source must not be empty")
    column_labels = [str(column) for column in source.columns.tolist()]
    if len(set(column_labels)) != len(column_labels):
        raise PhosPyInputError(f"{importer_label} import source columns must be unique")


def resolve_required_column(
    columns: pd.Index,
    *,
    explicit: str | None,
    candidates: tuple[str, ...],
    field_name: str,
    importer_label: str,
    validate_column_name: Callable[..., str | None],
) -> str:
    resolved = resolve_column(
        columns,
        explicit=explicit,
        candidates=candidates,
        field_name=field_name,
        importer_label=importer_label,
        required=True,
        validate_column_name=validate_column_name,
    )
    if resolved is None:  # pragma: no cover - resolve_column raises first.
        raise PhosPyInputError(
            f"{importer_label} importer could not infer {field_name}; configure an "
            "explicit column mapping."
        )
    return resolved


def resolve_column(
    columns: pd.Index,
    *,
    explicit: str | None,
    candidates: tuple[str, ...],
    field_name: str,
    importer_label: str,
    required: bool,
    validate_column_name: Callable[..., str | None],
) -> str | None:
    override = validate_column_name(explicit, field_name=field_name)
    if override is not None:
        if override in columns:
            return override
        raise PhosPyInputError(f"{field_name}={override!r} is not present in source")
    for candidate in candidates:
        match = find_column(columns, candidate, importer_label=importer_label)
        if match is not None:
            return match
    if not required:
        return None
    accepted = ", ".join(repr(candidate) for candidate in candidates)
    raise PhosPyInputError(
        f"{importer_label} importer could not infer {field_name}; configure an "
        f"explicit column mapping. tried: {accepted}"
    )


def find_column(
    columns: pd.Index,
    wanted: str,
    *,
    importer_label: str,
) -> str | None:
    if wanted in columns:
        return str(wanted)
    normalised_wanted = normalise_column_label(wanted)
    matches = [
        str(column)
        for column in columns.tolist()
        if normalise_column_label(str(column)) == normalised_wanted
    ]
    if len(matches) > 1:
        joined = ", ".join(repr(item) for item in matches)
        raise PhosPyInputError(
            f"{importer_label} importer found ambiguous source columns for "
            f"{wanted!r}: {joined}"
        )
    if matches:
        return matches[0]
    return None


def normalise_column_label(value: str) -> str:
    return " ".join(str(value).strip().lower().split())


def resolve_intensity_columns(
    source: pd.DataFrame,
    value: Mapping[str, str] | Sequence[str] | None,
    *,
    intensity_column_prefixes: Sequence[str],
    importer_label: str,
    request_label: str,
    mapping_class_name: str,
    reject_duplicate_inferred_sample_ids: bool,
) -> dict[str, str]:
    if value is not None:
        mapping = normalise_sample_column_mapping(value)
        missing = [column for column in mapping if column not in source.columns]
        if missing:
            joined = ", ".join(missing)
            raise PhosPyInputError(
                f"{importer_label} intensity column mapping includes missing columns: "
                f"{joined}"
            )
        return mapping
    if isinstance(intensity_column_prefixes, str) or not isinstance(
        intensity_column_prefixes,
        Sequence,
    ):
        raise PhosPyInputError(
            f"{request_label} import request intensity_column_prefixes must be a "
            "sequence of strings"
        )
    prefixes = tuple(str(prefix) for prefix in intensity_column_prefixes)
    if not prefixes or any(prefix.strip() == "" for prefix in prefixes):
        raise PhosPyInputError(
            f"{request_label} import request intensity_column_prefixes must contain "
            "non-empty strings"
        )
    mapping: dict[str, str] = {}
    inferred_sample_columns: dict[str, list[str]] = {}
    for column in source.columns.astype(str).tolist():
        for prefix in prefixes:
            if not column.startswith(prefix):
                continue
            sample_id = column[len(prefix) :].strip() or column
            mapping[column] = sample_id
            inferred_sample_columns.setdefault(sample_id, []).append(column)
            break
    if reject_duplicate_inferred_sample_ids:
        duplicate_sample_ids = {
            sample_id: columns
            for sample_id, columns in inferred_sample_columns.items()
            if len(columns) > 1
        }
        if duplicate_sample_ids:
            details = "; ".join(
                f"{sample_id!r}: {', '.join(columns)}"
                for sample_id, columns in duplicate_sample_ids.items()
            )
            raise PhosPyInputError(
                f"{importer_label} importer inferred multiple intensity columns for "
                f"the same inferred sample IDs. Configure "
                f"{mapping_class_name}.intensity_columns to select one quantitative "
                f"family or provide unique sample IDs. ambiguous_sample_ids={details}"
            )
    if not mapping:
        accepted = ", ".join(repr(prefix) for prefix in prefixes)
        raise PhosPyInputError(
            f"{importer_label} importer could not infer intensity columns. Configure "
            f"{mapping_class_name}.intensity_columns or adjust "
            f"intensity_column_prefixes. tried prefixes: {accepted}"
        )
    return normalise_sample_column_mapping(mapping)


def resolve_flag_series(
    source: pd.DataFrame,
    *,
    column: str | None,
    field_name: str,
) -> pd.Series | None:
    if column is None:
        return None
    return pd.Series(
        [
            parse_flag(
                value,
                field_name=field_name,
                row_position=position,
            )
            for position, value in enumerate(source.loc[:, column].tolist())
        ],
        index=source.index.copy(),
        dtype=bool,
    )


def raise_for_forbidden_flags(
    values: pd.Series,
    *,
    policy: str,
    error_policy: str,
    importer_label: str,
    label: str,
) -> None:
    if policy != error_policy:
        return
    count = int(values.astype(bool).sum())
    if not count:
        return
    raise PhosPyInputError(
        f"{importer_label} importer encountered {count} {label} row(s) under "
        "policy='error'"
    )


def parse_flag(
    value: object,
    *,
    field_name: str,
    row_position: int,
) -> bool:
    if is_missing(value):
        return False
    if isinstance(value, bool):
        return bool(value)
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        numeric = float(value)
        if numeric == 0.0:
            return False
        if numeric == 1.0:
            return True
    token = str(value).strip().lower()
    if token in {"", "-", "false", "f", "no", "n", "0"}:
        return False
    if token in {"+", "true", "t", "yes", "y", "1"}:
        return True
    raise PhosPyInputError(
        f"{field_name} row_position={row_position} contains unsupported flag "
        f"value {value!r}; expected '+', blank, boolean, or 0/1 style values"
    )


def build_row_ids(
    *,
    source: pd.DataFrame,
    explicit_column: str | None,
    protein_values: list[str],
    site_values: list[str],
    source_row_numbers: list[int],
    importer_label: str,
    generated_prefix: str,
) -> list[str]:
    if explicit_column is not None:
        return [
            required_text(
                value,
                field_name=f"{importer_label} {explicit_column}",
                row_position=position,
            )
            for position, value in enumerate(source.loc[:, explicit_column])
        ]
    return [
        f"{generated_prefix}:{protein}:{site}:row{row_number}"
        for protein, site, row_number in zip(
            protein_values,
            site_values,
            source_row_numbers,
            strict=True,
        )
    ]


def build_unique_feature_ids(
    *,
    source: pd.DataFrame,
    explicit_column: str | None,
    source_row_numbers: list[int],
    importer_label: str,
    generated_prefix: str,
) -> list[str]:
    if explicit_column is not None:
        return [
            required_text(
                value,
                field_name=f"{importer_label} {explicit_column}",
                row_position=position,
            )
            for position, value in enumerate(source.loc[:, explicit_column])
        ]
    return [
        f"{generated_prefix}_feature_{row_number}" for row_number in source_row_numbers
    ]


def first_list_token(
    value: object,
    *,
    field_name: str,
    row_position: int,
) -> str:
    tokens = split_multi_value(value)
    if not tokens:
        raise PhosPyInputError(
            f"{field_name} must contain non-empty values; row_position={row_position}"
        )
    return tokens[0]


def multi_value_count(value: object) -> int:
    return len(split_multi_value(value))


def split_multi_value(value: object) -> list[str]:
    if is_missing(value):
        return []
    if isinstance(value, str):
        stripped = value.strip()
        if stripped == "":
            return []
        return [
            token.strip()
            for token in _MULTI_VALUE_SPLIT_PATTERN.split(stripped)
            if token.strip()
        ]
    return [str(value).strip()]


def required_text(
    value: object,
    *,
    field_name: str,
    row_position: int,
) -> str:
    if is_missing(value):
        raise PhosPyInputError(
            f"{field_name} must not contain missing values; row_position={row_position}"
        )
    token = str(value).strip()
    if token == "":
        raise PhosPyInputError(
            f"{field_name} must contain non-empty values; row_position={row_position}"
        )
    return token


def optional_text(value: object) -> str | None:
    if is_missing(value):
        return None
    token = str(value).strip()
    if token == "":
        return None
    return token


def is_missing(value: object) -> bool:
    if value is None or value is pd.NA or value is pd.NaT:
        return True
    if isinstance(value, float):
        return math.isnan(value)
    if isinstance(value, np.floating):
        scalar_value: object = value
        return str(scalar_value).lower() == "nan"
    if isinstance(value, (np.datetime64, np.timedelta64)):
        temporal_value: object = value
        return str(temporal_value) == "NaT"
    return False


__all__ = [
    "build_row_ids",
    "build_unique_feature_ids",
    "find_column",
    "first_list_token",
    "is_missing",
    "multi_value_count",
    "normalise_column_label",
    "optional_text",
    "parse_flag",
    "raise_for_forbidden_flags",
    "require_non_empty_unique_columns",
    "required_text",
    "resolve_column",
    "resolve_flag_series",
    "resolve_intensity_columns",
    "resolve_required_column",
    "split_multi_value",
]
