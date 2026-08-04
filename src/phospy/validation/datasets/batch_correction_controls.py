"""Control-site provenance validation for applied batch correction."""
# pyright: reportUnknownMemberType=false

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import cast

from phospy.errors.input import PhosPyInputError
from phospy.provenance.models import BATCH_CORRECTION_SELECTED_SITE_KEY_ROW_SENTINELS
from phospy.validation.datasets._batch_correction_helpers import (
    NOT_PROVIDED_VALUES,
    duplicates_in_order,
    format_labels,
    format_positions,
    has_non_missing_text,
    is_missing_value,
    is_not_provided,
    reject_not_provided_required_mapping,
    require_non_empty_mapping,
)

_SELECTED_SITE_KEY_ROW_SENTINELS = (
    BATCH_CORRECTION_SELECTED_SITE_KEY_ROW_SENTINELS | NOT_PROVIDED_VALUES
)
_STRICT_CONTROL_SOURCE_TYPE_MARKERS = frozenset({"packaged", "reference", "external"})
_CALLER_CONTROL_SOURCE_AUDIT_FIELDS = (
    "organism",
    "identifier_namespace",
    "source_version",
    "license",
    "redistribution",
)
_STRICT_CONTROL_SOURCE_REQUIRED_FIELDS = (
    "organism",
    "identifier_namespace",
    "source_name",
    "source_version",
    "license",
    "redistribution",
    "selection_method",
)


def normalize_applied_selected_site_key_rows(rows: Sequence[object]) -> tuple[str, ...]:
    """Normalize and validate applied selected control row identifiers."""

    if not rows:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls/control provenance must include selected_site_key_rows"
        )

    normalized_rows: list[str] = []
    missing_rows: list[int] = []
    blank_rows: list[int] = []
    sentinel_rows: list[int] = []
    for position, row in enumerate(tuple(rows)):
        if is_missing_value(row):
            missing_rows.append(position)
            continue
        normalized = str(row).strip()
        if normalized == "":
            blank_rows.append(position)
            continue
        if _is_selected_site_key_row_sentinel(normalized):
            sentinel_rows.append(position)
            continue
        normalized_rows.append(normalized)

    if missing_rows:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls/control provenance selected_site_key_rows contains missing "
            f"site_key rows at positions {format_positions(missing_rows)}"
        )
    if blank_rows:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls/control provenance selected_site_key_rows contains blank "
            f"site_key rows at positions {format_positions(blank_rows)}"
        )
    if sentinel_rows:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls/control provenance selected_site_key_rows contains sentinel "
            f"site_key rows at positions {format_positions(sentinel_rows)}"
        )

    return tuple(normalized_rows)


def _is_selected_site_key_row_sentinel(value: object) -> bool:
    return str(value).strip().lower() in _SELECTED_SITE_KEY_ROW_SENTINELS


def require_unique_selected_control_site_rows(rows: Sequence[str]) -> None:
    duplicates = duplicates_in_order(rows)
    if duplicates:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance "
            "selected_site_key_rows contains duplicate selected control row "
            f"identifiers: {format_labels(duplicates)}"
        )


def require_control_site_source_metadata(
    source: Mapping[str, object],
    *,
    selected_site_key_rows: Sequence[str],
) -> None:
    require_non_empty_mapping(
        source,
        field_name="BatchCorrectionProvenance.control_site_source",
    )
    reject_not_provided_required_mapping(
        source,
        field_name="BatchCorrectionProvenance.control_site_source",
    )
    source_type = _source_type(source)
    if source_type is None or is_not_provided(source_type):
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance selected "
            "controls/control provenance must include control source metadata"
        )

    if _has_strict_control_source_type(source):
        missing = tuple(
            field_name
            for field_name in _STRICT_CONTROL_SOURCE_REQUIRED_FIELDS
            if not has_non_missing_text(source.get(field_name))
        )
        if missing:
            raise PhosPyInputError(
                "corrected_preprocessing_output BatchCorrectionProvenance "
                "packaged/reference/external control-source metadata is "
                f"incomplete; missing {format_labels(missing)}"
            )
        return

    missing_without_reason = tuple(
        field_name
        for field_name in _CALLER_CONTROL_SOURCE_AUDIT_FIELDS
        if not has_non_missing_text(source.get(field_name))
        and not _has_metadata_missing_reason(
            source,
            field_name,
            selected_site_key_rows=selected_site_key_rows,
        )
    )
    if missing_without_reason:
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance control "
            "source is missing "
            f"{format_labels(missing_without_reason)} without explicit rationale"
        )

    has_source_name = has_non_missing_text(source.get("source_name"))
    has_source_version = has_non_missing_text(source.get("source_version"))
    has_unavailable_reason = has_non_missing_text(
        source.get("source_version_unavailable_reason")
    )
    has_missing_reason = _has_metadata_missing_reason(
        source,
        "source_version",
        selected_site_key_rows=selected_site_key_rows,
    )
    if has_source_name and not (
        has_source_version or has_unavailable_reason or has_missing_reason
    ):
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance control "
            "source declares source_name without source_version or explicit "
            "source_version_unavailable_reason"
        )
    if source_type == "caller_supplied" and not (
        has_source_version or has_unavailable_reason or has_missing_reason
    ):
        raise PhosPyInputError(
            "corrected_preprocessing_output BatchCorrectionProvenance "
            "caller_supplied control source must record source_version or "
            "source_version_unavailable_reason"
        )


def _source_type(source: Mapping[str, object]) -> str | None:
    _require_consistent_source_aliases(source)
    for key in ("source_type", "source"):
        value = source.get(key)
        if has_non_missing_text(value):
            return str(value).strip().lower()
    return None


def _require_consistent_source_aliases(source: Mapping[str, object]) -> None:
    source_type = source.get("source_type")
    source_alias = source.get("source")
    if not has_non_missing_text(source_type) or not has_non_missing_text(source_alias):
        return
    normalized_source_type = str(source_type).strip().lower()
    normalized_source_alias = str(source_alias).strip().lower()
    if normalized_source_type == normalized_source_alias:
        return
    raise PhosPyInputError(
        "corrected_preprocessing_output BatchCorrectionProvenance control "
        "source has conflicting source declarations: "
        f"source_type={normalized_source_type!r}, source={normalized_source_alias!r}"
    )


def _has_strict_control_source_type(source: Mapping[str, object]) -> bool:
    for key in (
        "source_type",
        "source",
        "control_site_set_source_type",
        "source_name",
    ):
        value = source.get(key)
        if has_non_missing_text(value) and _is_strict_control_source_type(
            str(value).strip().lower()
        ):
            return True
    return False


def _is_strict_control_source_type(source_type: str | None) -> bool:
    if source_type is None:
        return False
    tokens = frozenset(source_type.replace("-", "_").split("_"))
    return bool(tokens & _STRICT_CONTROL_SOURCE_TYPE_MARKERS)


def _has_metadata_missing_reason(
    source: Mapping[str, object],
    field_name: str,
    *,
    selected_site_key_rows: Sequence[str],
) -> bool:
    reasons = source.get("metadata_missing_reason")
    if isinstance(reasons, Mapping) and has_non_missing_text(
        cast(Mapping[str, object], reasons).get(field_name)
    ):
        return True
    if has_non_missing_text(source.get(f"{field_name}_missing_reason")):
        return True
    if field_name == "source_version" and has_non_missing_text(
        source.get("source_version_unavailable_reason")
    ):
        return True
    return _has_metadata_missing_reason_by_site_key(
        source,
        field_name,
        selected_site_key_rows=selected_site_key_rows,
    )


def _has_metadata_missing_reason_by_site_key(
    source: Mapping[str, object],
    field_name: str,
    *,
    selected_site_key_rows: Sequence[str],
) -> bool:
    reasons_by_site_key = source.get("metadata_missing_reason_by_site_key")
    if not isinstance(reasons_by_site_key, Mapping):
        return False
    selected = tuple(str(site_key) for site_key in selected_site_key_rows)
    if not selected:
        return False
    by_site_key = cast(Mapping[str, object], reasons_by_site_key)
    for site_key in selected:
        site_reasons = by_site_key.get(site_key)
        if not isinstance(site_reasons, Mapping):
            return False
        if not has_non_missing_text(
            cast(Mapping[str, object], site_reasons).get(field_name)
        ):
            return False
    return True


__all__ = ["normalize_applied_selected_site_key_rows"]
