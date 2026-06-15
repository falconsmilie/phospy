"""Offline readers for enrichment set collection files."""

from __future__ import annotations

import csv
import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO, cast

from phospy.errors.input import PhosPyInputError, UnsupportedInputFormatError
from phospy.errors.validation import WorkflowValidationError
from phospy.science.enrichment.collections import (
    EnrichmentCollectionKind,
    EnrichmentIdentifierKind,
    EnrichmentSet,
    EnrichmentSetCollection,
)

_CSV = "csv"
_TSV = "tsv"
_REQUIRED_TABLE_COLUMNS = frozenset({"set_id", "name", "identifier"})
_IDENTIFIER_KIND_COLUMN = "identifier_kind"


@dataclass(slots=True)
class _PendingSet:
    set_id: str
    name: str
    identifier_kind: EnrichmentIdentifierKind
    identifiers: list[str]
    source_name: str | None
    source_version: str | None
    description: str | None


def read_enrichment_sets_gmt(
    path: Path,
    *,
    identifier_kind: EnrichmentIdentifierKind,
    collection_kind: EnrichmentCollectionKind | None = None,
    source_name: str | None = None,
    source_version: str | None = None,
) -> EnrichmentSetCollection:
    """Read a local GMT-like set file.

    The expected shape is ``set_id<TAB>description<TAB>identifier...``. GMT does
    not carry identifier semantics, so ``identifier_kind`` is required.
    """

    normalized_path = Path(path)
    enrichment_sets: list[EnrichmentSet] = []
    try:
        with normalized_path.open("r", encoding="utf-8", newline="") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.rstrip("\r\n")
                if line.strip() == "":
                    continue
                cells = line.split("\t")
                if len(cells) < 2:
                    raise PhosPyInputError(
                        "failed to parse GMT enrichment set file "
                        f"'{normalized_path}': line {line_number} must contain "
                        "set_id and description columns"
                    )
                set_id = cells[0]
                description = cells[1].strip() or None
                identifiers = cells[2:]
                enrichment_sets.append(
                    EnrichmentSet(
                        set_id=set_id,
                        name=set_id,
                        identifiers=identifiers,
                        identifier_kind=identifier_kind,
                        source_name=source_name,
                        source_version=source_version,
                        description=description,
                    )
                )
    except FileNotFoundError as exc:
        raise PhosPyInputError(f"input file does not exist: {normalized_path}") from exc
    except PermissionError as exc:
        raise PhosPyInputError(
            f"permission denied while reading input file: {normalized_path}"
        ) from exc
    except (OSError, UnicodeDecodeError) as exc:
        raise PhosPyInputError(
            f"failed to read GMT enrichment set file '{normalized_path}': {exc}"
        ) from exc
    except WorkflowValidationError as exc:
        raise PhosPyInputError(
            f"invalid GMT enrichment set file '{normalized_path}': {exc}"
        ) from exc
    return _build_collection(
        normalized_path,
        enrichment_sets=enrichment_sets,
        identifier_kind=identifier_kind,
        collection_kind=collection_kind,
        source_name=source_name,
        source_version=source_version,
    )


def read_enrichment_sets_table(
    path: Path,
    *,
    identifier_kind: EnrichmentIdentifierKind | None = None,
    collection_kind: EnrichmentCollectionKind | None = None,
    source_name: str | None = None,
    source_version: str | None = None,
) -> EnrichmentSetCollection:
    """Read local CSV/TSV enrichment set rows.

    Required columns are ``set_id``, ``name``, and ``identifier``. If the table
    lacks an ``identifier_kind`` column, ``identifier_kind`` must be provided.
    """

    normalized_path = Path(path)
    table_format = _table_format_from_path(normalized_path)
    sep = "," if table_format == _CSV else "\t"
    return _read_delimited_enrichment_sets(
        normalized_path,
        sep=sep,
        identifier_kind=identifier_kind,
        collection_kind=collection_kind,
        source_name=source_name,
        source_version=source_version,
    )


def read_enrichment_sets_csv(
    path: Path,
    *,
    identifier_kind: EnrichmentIdentifierKind | None = None,
    collection_kind: EnrichmentCollectionKind | None = None,
    source_name: str | None = None,
    source_version: str | None = None,
) -> EnrichmentSetCollection:
    """Read local CSV enrichment set rows."""

    return _read_delimited_enrichment_sets(
        Path(path),
        sep=",",
        identifier_kind=identifier_kind,
        collection_kind=collection_kind,
        source_name=source_name,
        source_version=source_version,
    )


def read_enrichment_sets_tsv(
    path: Path,
    *,
    identifier_kind: EnrichmentIdentifierKind | None = None,
    collection_kind: EnrichmentCollectionKind | None = None,
    source_name: str | None = None,
    source_version: str | None = None,
) -> EnrichmentSetCollection:
    """Read local TSV enrichment set rows."""

    return _read_delimited_enrichment_sets(
        Path(path),
        sep="\t",
        identifier_kind=identifier_kind,
        collection_kind=collection_kind,
        source_name=source_name,
        source_version=source_version,
    )


def load_enrichment_sets_gmt(
    path: Path,
    *,
    identifier_kind: EnrichmentIdentifierKind,
    collection_kind: EnrichmentCollectionKind | None = None,
    source_name: str | None = None,
    source_version: str | None = None,
) -> EnrichmentSetCollection:
    """Alias for ``read_enrichment_sets_gmt``."""

    _warn_deprecated_load_alias(
        "load_enrichment_sets_gmt",
        "read_enrichment_sets_gmt",
    )
    return read_enrichment_sets_gmt(
        path,
        identifier_kind=identifier_kind,
        collection_kind=collection_kind,
        source_name=source_name,
        source_version=source_version,
    )


def load_enrichment_sets_table(
    path: Path,
    *,
    identifier_kind: EnrichmentIdentifierKind | None = None,
    collection_kind: EnrichmentCollectionKind | None = None,
    source_name: str | None = None,
    source_version: str | None = None,
) -> EnrichmentSetCollection:
    """Alias for ``read_enrichment_sets_table``."""

    _warn_deprecated_load_alias(
        "load_enrichment_sets_table",
        "read_enrichment_sets_table",
    )
    return read_enrichment_sets_table(
        path,
        identifier_kind=identifier_kind,
        collection_kind=collection_kind,
        source_name=source_name,
        source_version=source_version,
    )


def load_enrichment_sets_csv(
    path: Path,
    *,
    identifier_kind: EnrichmentIdentifierKind | None = None,
    collection_kind: EnrichmentCollectionKind | None = None,
    source_name: str | None = None,
    source_version: str | None = None,
) -> EnrichmentSetCollection:
    """Alias for ``read_enrichment_sets_csv``."""

    _warn_deprecated_load_alias(
        "load_enrichment_sets_csv",
        "read_enrichment_sets_csv",
    )
    return read_enrichment_sets_csv(
        path,
        identifier_kind=identifier_kind,
        collection_kind=collection_kind,
        source_name=source_name,
        source_version=source_version,
    )


def load_enrichment_sets_tsv(
    path: Path,
    *,
    identifier_kind: EnrichmentIdentifierKind | None = None,
    collection_kind: EnrichmentCollectionKind | None = None,
    source_name: str | None = None,
    source_version: str | None = None,
) -> EnrichmentSetCollection:
    """Alias for ``read_enrichment_sets_tsv``."""

    _warn_deprecated_load_alias(
        "load_enrichment_sets_tsv",
        "read_enrichment_sets_tsv",
    )
    return read_enrichment_sets_tsv(
        path,
        identifier_kind=identifier_kind,
        collection_kind=collection_kind,
        source_name=source_name,
        source_version=source_version,
    )


def _warn_deprecated_load_alias(alias_name: str, replacement_name: str) -> None:
    warnings.warn(
        (
            f"{alias_name} is deprecated and will be removed in a future release; "
            f"use {replacement_name} instead."
        ),
        DeprecationWarning,
        stacklevel=3,
    )


def _read_delimited_enrichment_sets(
    path: Path,
    *,
    sep: str,
    identifier_kind: EnrichmentIdentifierKind | None,
    collection_kind: EnrichmentCollectionKind | None,
    source_name: str | None,
    source_version: str | None,
) -> EnrichmentSetCollection:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            return _read_delimited_enrichment_sets_from_handle(
                path,
                handle,
                sep=sep,
                identifier_kind=identifier_kind,
                collection_kind=collection_kind,
                source_name=source_name,
                source_version=source_version,
            )
    except FileNotFoundError as exc:
        raise PhosPyInputError(f"input file does not exist: {path}") from exc
    except PermissionError as exc:
        raise PhosPyInputError(
            f"permission denied while reading input file: {path}"
        ) from exc
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        raise PhosPyInputError(
            f"failed to read enrichment set table '{path}': {exc}"
        ) from exc
    except WorkflowValidationError as exc:
        raise PhosPyInputError(f"invalid enrichment set table '{path}': {exc}") from exc


def _read_delimited_enrichment_sets_from_handle(
    path: Path,
    handle: TextIO,
    *,
    sep: str,
    identifier_kind: EnrichmentIdentifierKind | None,
    collection_kind: EnrichmentCollectionKind | None,
    source_name: str | None,
    source_version: str | None,
) -> EnrichmentSetCollection:
    reader = csv.DictReader(handle, delimiter=sep)
    if reader.fieldnames is None:
        raise PhosPyInputError(f"enrichment set table '{path}' is empty")
    fieldnames = {field.strip() for field in reader.fieldnames if field is not None}
    missing_columns = sorted(_REQUIRED_TABLE_COLUMNS - fieldnames)
    if missing_columns:
        joined = ", ".join(missing_columns)
        raise PhosPyInputError(
            f"enrichment set table '{path}' is missing required columns: {joined}"
        )
    has_identifier_kind_column = _IDENTIFIER_KIND_COLUMN in fieldnames
    if not has_identifier_kind_column and identifier_kind is None:
        raise PhosPyInputError(
            "enrichment set table identifier_kind must be provided when the file "
            "does not contain an identifier_kind column"
        )

    pending_by_id: dict[str, _PendingSet] = {}
    row_count = 0
    for row_number, row in enumerate(reader, start=2):
        row_count += 1
        set_id = _required_cell(row, "set_id", path=path, row_number=row_number)
        name = _required_cell(row, "name", path=path, row_number=row_number)
        identifier = _required_cell(
            row,
            "identifier",
            path=path,
            row_number=row_number,
        )
        row_identifier_kind = _resolve_row_identifier_kind(
            row,
            explicit_identifier_kind=identifier_kind,
            path=path,
            row_number=row_number,
        )
        row_source_name = _optional_cell(row, "source_name") or source_name
        row_source_version = _optional_cell(row, "source_version") or source_version
        row_description = _optional_cell(row, "description")
        pending = pending_by_id.get(set_id)
        if pending is None:
            pending_by_id[set_id] = _PendingSet(
                set_id=set_id,
                name=name,
                identifier_kind=row_identifier_kind,
                identifiers=[identifier],
                source_name=row_source_name,
                source_version=row_source_version,
                description=row_description,
            )
            continue
        _require_matching_set_field(
            pending.name,
            name,
            path=path,
            row_number=row_number,
            set_id=set_id,
            column_name="name",
        )
        _require_matching_set_field(
            pending.identifier_kind,
            row_identifier_kind,
            path=path,
            row_number=row_number,
            set_id=set_id,
            column_name="identifier_kind",
        )
        pending.source_name = _merge_optional_set_field(
            pending.source_name,
            row_source_name,
            path=path,
            row_number=row_number,
            set_id=set_id,
            column_name="source_name",
        )
        pending.source_version = _merge_optional_set_field(
            pending.source_version,
            row_source_version,
            path=path,
            row_number=row_number,
            set_id=set_id,
            column_name="source_version",
        )
        pending.description = _merge_optional_set_field(
            pending.description,
            row_description,
            path=path,
            row_number=row_number,
            set_id=set_id,
            column_name="description",
        )
        pending.identifiers.append(identifier)
    if row_count == 0:
        raise PhosPyInputError(f"enrichment set table '{path}' contains no rows")
    enrichment_sets = [
        EnrichmentSet(
            set_id=pending.set_id,
            name=pending.name,
            identifiers=pending.identifiers,
            identifier_kind=pending.identifier_kind,
            source_name=pending.source_name,
            source_version=pending.source_version,
            description=pending.description,
        )
        for pending in pending_by_id.values()
    ]
    return _build_collection(
        path,
        enrichment_sets=enrichment_sets,
        identifier_kind=identifier_kind,
        collection_kind=collection_kind,
        source_name=source_name
        or _shared_optional(
            [enrichment_set.source_name for enrichment_set in enrichment_sets]
        ),
        source_version=source_version
        or _shared_optional(
            [enrichment_set.source_version for enrichment_set in enrichment_sets]
        ),
    )


def _build_collection(
    path: Path,
    *,
    enrichment_sets: Sequence[EnrichmentSet],
    identifier_kind: EnrichmentIdentifierKind | None,
    collection_kind: EnrichmentCollectionKind | None,
    source_name: str | None,
    source_version: str | None,
) -> EnrichmentSetCollection:
    try:
        return EnrichmentSetCollection(
            sets=tuple(enrichment_sets),
            identifier_kind=identifier_kind,
            collection_kind=collection_kind,
            source_name=source_name,
            source_version=source_version,
        )
    except WorkflowValidationError as exc:
        raise PhosPyInputError(
            f"invalid enrichment set collection in '{path}': {exc}"
        ) from exc


def _resolve_row_identifier_kind(
    row: dict[str, str],
    *,
    explicit_identifier_kind: EnrichmentIdentifierKind | None,
    path: Path,
    row_number: int,
) -> EnrichmentIdentifierKind:
    row_identifier_kind = _optional_cell(row, _IDENTIFIER_KIND_COLUMN)
    if row_identifier_kind is None:
        if explicit_identifier_kind is None:
            raise PhosPyInputError(
                f"enrichment set table '{path}' row {row_number} must include "
                "identifier_kind or the reader must be called with identifier_kind"
            )
        return explicit_identifier_kind
    if (
        explicit_identifier_kind is not None
        and row_identifier_kind != explicit_identifier_kind
    ):
        raise PhosPyInputError(
            f"enrichment set table '{path}' row {row_number} identifier_kind "
            f"does not match explicit identifier_kind; "
            f"observed={row_identifier_kind!r}, expected={explicit_identifier_kind!r}"
        )
    return cast(EnrichmentIdentifierKind, row_identifier_kind)


def _required_cell(
    row: dict[str, str],
    column_name: str,
    *,
    path: Path,
    row_number: int,
) -> str:
    value = row.get(column_name)
    if value is None or value.strip() == "":
        raise PhosPyInputError(
            f"enrichment set table '{path}' row {row_number} column "
            f"{column_name!r} must be non-empty"
        )
    return value.strip()


def _optional_cell(row: dict[str, str], column_name: str) -> str | None:
    value = row.get(column_name)
    if value is None:
        return None
    normalized = value.strip()
    if normalized == "":
        return None
    return normalized


def _require_matching_set_field(
    observed: object,
    expected: object,
    *,
    path: Path,
    row_number: int,
    set_id: str,
    column_name: str,
) -> None:
    if observed == expected:
        return
    raise PhosPyInputError(
        f"enrichment set table '{path}' row {row_number} has conflicting "
        f"{column_name} for set_id={set_id!r}; observed={expected!r}, "
        f"expected={observed!r}"
    )


def _merge_optional_set_field(
    observed: str | None,
    new_value: str | None,
    *,
    path: Path,
    row_number: int,
    set_id: str,
    column_name: str,
) -> str | None:
    if observed is None:
        return new_value
    if new_value is None or observed == new_value:
        return observed
    raise PhosPyInputError(
        f"enrichment set table '{path}' row {row_number} has conflicting "
        f"{column_name} for set_id={set_id!r}; observed={new_value!r}, "
        f"expected={observed!r}"
    )


def _shared_optional(values: Sequence[str | None]) -> str | None:
    non_empty = [value for value in values if value is not None]
    if not non_empty:
        return None
    first = non_empty[0]
    if all(value == first for value in non_empty):
        return first
    return None


def _table_format_from_path(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return _CSV
    if suffix in {".tsv", ".txt"}:
        return _TSV
    raise UnsupportedInputFormatError(
        "unsupported enrichment set table format for "
        f"'{path}'. supported formats: csv (.csv), tsv (.tsv), txt (.txt)"
    )


__all__ = [
    "load_enrichment_sets_csv",
    "load_enrichment_sets_gmt",
    "load_enrichment_sets_table",
    "load_enrichment_sets_tsv",
    "read_enrichment_sets_csv",
    "read_enrichment_sets_gmt",
    "read_enrichment_sets_table",
    "read_enrichment_sets_tsv",
]
