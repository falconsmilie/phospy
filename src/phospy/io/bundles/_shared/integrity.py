"""Bundle file-integrity records and verification helpers."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from pathlib import Path
from string import hexdigits

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.paths import resolve_bundle_relative_path
from phospy.io.bundles._shared.primitives import (
    require_int,
    require_mapping,
    require_str,
)
from phospy.provenance.hashing import DEFAULT_FILE_HASH_ALGORITHM, hash_file_bytes

_FILE_ENTRY_ALLOWED_FIELDS = frozenset({"path", "sha256", "byte_size", "logical_type"})
_TABLE_ENTRY_ALLOWED_FIELDS = frozenset(
    {"path", "sha256", "byte_size", "logical_type", "shape"}
)
_TABLE_SHAPE_ALLOWED_FIELDS = frozenset({"rows", "columns"})


def build_table_file_entry(
    *,
    bundle_root: Path,
    path: Path,
    table: pd.DataFrame,
    logical_type: str,
) -> dict[str, object]:
    """Build the manifest entry for a serialized table file."""

    entry = _build_file_entry(
        bundle_root=bundle_root,
        path=path,
        logical_type=logical_type,
    )
    entry["shape"] = {
        "rows": int(table.shape[0]),
        "columns": int(table.shape[1]),
    }
    return entry


def build_json_file_entry(
    *,
    bundle_root: Path,
    path: Path,
    logical_type: str,
) -> dict[str, object]:
    """Build the manifest entry for a serialized JSON sidecar."""

    return _build_file_entry(
        bundle_root=bundle_root,
        path=path,
        logical_type=logical_type,
    )


def require_file_entry(
    value: object,
    *,
    field_name: str,
    expected_logical_type: str | None = None,
) -> Mapping[str, object]:
    """Validate and return a manifest file entry."""

    entry = require_mapping(value, field_name=field_name)
    _reject_unsupported_fields(
        entry,
        field_name=field_name,
        allowed_fields=_FILE_ENTRY_ALLOWED_FIELDS,
    )
    _require_common_file_entry_fields(
        entry,
        field_name=field_name,
        expected_logical_type=expected_logical_type,
    )
    return entry


def require_table_entry(
    value: object,
    *,
    field_name: str,
    expected_logical_type: str | None = None,
) -> Mapping[str, object]:
    """Validate and return a manifest table-file entry."""

    entry = require_mapping(value, field_name=field_name)
    _reject_unsupported_fields(
        entry,
        field_name=field_name,
        allowed_fields=_TABLE_ENTRY_ALLOWED_FIELDS,
    )
    _require_common_file_entry_fields(
        entry,
        field_name=field_name,
        expected_logical_type=expected_logical_type,
    )
    shape = require_mapping(entry.get("shape"), field_name=f"{field_name}.shape")
    _reject_unsupported_fields(
        shape,
        field_name=f"{field_name}.shape",
        allowed_fields=_TABLE_SHAPE_ALLOWED_FIELDS,
    )
    rows = require_int(shape.get("rows"), field_name=f"{field_name}.shape.rows")
    columns = require_int(
        shape.get("columns"),
        field_name=f"{field_name}.shape.columns",
    )
    if rows < 0:
        raise PhosPyInputError(f"{field_name}.shape.rows must be non-negative")
    if columns < 0:
        raise PhosPyInputError(f"{field_name}.shape.columns must be non-negative")
    return entry


def require_optional_table_entry(
    value: object,
    *,
    field_name: str,
    expected_logical_type: str | None = None,
) -> Mapping[str, object] | None:
    """Validate an optional table entry, preserving explicit null."""

    if value is None:
        return None
    return require_table_entry(
        value,
        field_name=field_name,
        expected_logical_type=expected_logical_type,
    )


def file_entry_path(
    entry: Mapping[str, object],
    *,
    bundle_root: Path,
    field_name: str,
) -> Path:
    """Resolve a validated manifest file entry to an on-disk path."""

    return resolve_bundle_relative_path(
        bundle_root,
        require_str(entry.get("path"), field_name=f"{field_name}.path"),
        field_name=f"{field_name}.path",
    )


def table_entry_path(
    entry: Mapping[str, object],
    *,
    bundle_root: Path,
    field_name: str,
) -> Path:
    """Resolve a validated manifest table entry to an on-disk path."""

    return file_entry_path(entry, bundle_root=bundle_root, field_name=field_name)


def validate_table_entry_shape(
    table: pd.DataFrame,
    entry: Mapping[str, object],
    *,
    field_name: str,
) -> None:
    """Verify a loaded table matches the shape declared in the manifest."""

    shape = require_mapping(entry.get("shape"), field_name=f"{field_name}.shape")
    expected_rows = require_int(
        shape.get("rows"),
        field_name=f"{field_name}.shape.rows",
    )
    expected_columns = require_int(
        shape.get("columns"),
        field_name=f"{field_name}.shape.columns",
    )
    observed_shape = (int(table.shape[0]), int(table.shape[1]))
    expected_shape = (expected_rows, expected_columns)
    if observed_shape != expected_shape:
        raise PhosPyInputError(
            "bundle table shape mismatch: "
            f"{field_name}; expected rows={expected_rows}, "
            f"columns={expected_columns}; observed rows={observed_shape[0]}, "
            f"columns={observed_shape[1]}"
        )


def verify_bundle_integrity(
    *,
    bundle_root: Path,
    manifest_payload: Mapping[str, object],
    manifest_filename: str,
) -> None:
    """Verify manifest-declared file sizes/digests and reject extra files."""

    root = Path(bundle_root)
    records = list(
        _iter_declared_file_entries(
            manifest_payload,
            field_name="bundle manifest",
        )
    )
    declared_by_path: dict[str, Mapping[str, object]] = {}
    for field_name, entry in records:
        relative_path = require_str(
            entry.get("path"),
            field_name=f"{field_name}.path",
        )
        resolved = resolve_bundle_relative_path(
            root,
            relative_path,
            field_name=f"{field_name}.path",
        )
        normalized_relative_path = _relative_bundle_path(root, resolved)
        if normalized_relative_path in declared_by_path:
            raise PhosPyInputError(
                "bundle integrity check failed: duplicate manifest file entry "
                f"path={normalized_relative_path}"
            )
        declared_by_path[normalized_relative_path] = entry
        _verify_declared_file(
            path=resolved,
            entry=entry,
            field_name=field_name,
            relative_path=normalized_relative_path,
        )

    observed_files = {
        _relative_bundle_path(root, path)
        for path in root.rglob("*")
        if path.is_file() and _relative_bundle_path(root, path) != manifest_filename
    }
    declared_files = set(declared_by_path)
    extra_files = sorted(observed_files - declared_files)
    if extra_files:
        extra = ", ".join(extra_files)
        raise PhosPyInputError(
            "bundle integrity check failed: bundle contains undeclared file(s): "
            f"{extra}"
        )


def _build_file_entry(
    *,
    bundle_root: Path,
    path: Path,
    logical_type: str,
) -> dict[str, object]:
    relative_path = _relative_bundle_path(bundle_root, path)
    try:
        byte_size = int(path.stat().st_size)
        digest = hash_file_bytes(path, algorithm=DEFAULT_FILE_HASH_ALGORITHM)
    except OSError as exc:
        raise PhosPyInputError(f"failed to hash bundle file '{path}': {exc}") from exc
    return {
        "path": relative_path,
        "sha256": digest,
        "byte_size": byte_size,
        "logical_type": logical_type,
    }


def _iter_declared_file_entries(
    value: object,
    *,
    field_name: str,
) -> Iterator[tuple[str, Mapping[str, object]]]:
    if isinstance(value, Mapping):
        if {"path", "sha256", "byte_size"}.issubset(set(value.keys())):
            yield field_name, value
            return
        for key, item in value.items():
            yield from _iter_declared_file_entries(
                item,
                field_name=f"{field_name}.{str(key)}",
            )
        return
    if isinstance(value, list):
        for index, item in enumerate(value):
            yield from _iter_declared_file_entries(
                item,
                field_name=f"{field_name}[{index}]",
            )


def _require_common_file_entry_fields(
    entry: Mapping[str, object],
    *,
    field_name: str,
    expected_logical_type: str | None,
) -> None:
    require_str(entry.get("path"), field_name=f"{field_name}.path")
    digest = require_str(entry.get("sha256"), field_name=f"{field_name}.sha256")
    _require_sha256(digest, field_name=f"{field_name}.sha256")
    byte_size = require_int(
        entry.get("byte_size"),
        field_name=f"{field_name}.byte_size",
    )
    if byte_size < 0:
        raise PhosPyInputError(f"{field_name}.byte_size must be non-negative")
    logical_type = require_str(
        entry.get("logical_type"),
        field_name=f"{field_name}.logical_type",
    )
    if expected_logical_type is not None and logical_type != expected_logical_type:
        raise PhosPyInputError(
            f"{field_name}.logical_type must be '{expected_logical_type}'"
        )


def _verify_declared_file(
    *,
    path: Path,
    entry: Mapping[str, object],
    field_name: str,
    relative_path: str,
) -> None:
    if not path.is_file():
        raise PhosPyInputError(
            "bundle integrity check failed: declared file is missing: "
            f"path={relative_path}"
        )
    expected_size = require_int(
        entry.get("byte_size"),
        field_name=f"{field_name}.byte_size",
    )
    try:
        actual_size = int(path.stat().st_size)
    except OSError as exc:
        raise PhosPyInputError(
            "bundle integrity check failed: could not stat declared file "
            f"path={relative_path}: {exc}"
        ) from exc
    if actual_size != expected_size:
        raise PhosPyInputError(
            "bundle integrity check failed: declared file size mismatch: "
            f"path={relative_path}; expected byte_size={expected_size}; "
            f"actual byte_size={actual_size}"
        )
    expected_digest = require_str(
        entry.get("sha256"),
        field_name=f"{field_name}.sha256",
    )
    try:
        actual_digest = hash_file_bytes(path)
    except OSError as exc:
        raise PhosPyInputError(
            "bundle integrity check failed: could not hash declared file "
            f"path={relative_path}: {exc}"
        ) from exc
    if actual_digest != expected_digest:
        raise PhosPyInputError(
            "bundle integrity check failed: declared file digest mismatch: "
            f"path={relative_path}; expected sha256={expected_digest}; "
            f"actual sha256={actual_digest}"
        )


def _require_sha256(value: str, *, field_name: str) -> None:
    if len(value) != 64 or any(character not in hexdigits for character in value):
        raise PhosPyInputError(f"{field_name} must be a 64-character SHA-256 digest")


def _reject_unsupported_fields(
    payload: Mapping[str, object],
    *,
    field_name: str,
    allowed_fields: frozenset[str],
) -> None:
    unknown_fields = sorted(
        str(key) for key in payload.keys() if str(key) not in allowed_fields
    )
    if unknown_fields:
        unknown = ", ".join(unknown_fields)
        raise PhosPyInputError(f"{field_name} contains unsupported field(s): {unknown}")


def _relative_bundle_path(bundle_root: Path, path: Path) -> str:
    root = Path(bundle_root).resolve()
    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(root).as_posix()
    except ValueError as exc:
        raise PhosPyInputError(
            f"bundle file path points outside the bundle root: {path}"
        ) from exc
