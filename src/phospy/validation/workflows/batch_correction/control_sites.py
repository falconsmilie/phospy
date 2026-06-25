"""Control-site eligibility validation for batch-correction workflows."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from enum import Enum

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.control_sites import (
    CONTROL_SITE_REASON_CONTROL_ANNOTATION_NOT_IN_SITE_ROWS,
    CONTROL_SITE_REASON_DUPLICATE_CONTROL_ANNOTATION,
    CONTROL_SITE_REASON_MISSING_CONTROL_ANNOTATION,
    ControlSiteEligibility,
    ControlSiteMapping,
    ControlSiteSet,
    ControlSiteStatus,
)

_DEFAULT_IDENTIFIER_NAMESPACE = "site_key"
_CALLER_SUPPLIED_SOURCE_TYPES = frozenset({"caller_supplied", "local", "local_file"})
_STRICT_SOURCE_TYPE_MARKERS = frozenset({"packaged", "reference", "external"})
_CALLER_AUDIT_FIELDS = (
    "organism",
    "identifier_namespace",
    "source_version",
    "license",
    "redistribution",
)
_STRICT_SOURCE_REQUIRED_FIELDS = (
    "organism",
    "identifier_namespace",
    "source_name",
    "source_version",
    "license",
    "redistribution",
    "selection_method",
)
_CONTROL_IDENTITY_FIELDS = (
    "source_type",
    "organism",
    "identifier_namespace",
    "source_name",
    "source_version",
    "license",
    "redistribution",
)


@dataclass(frozen=True, slots=True)
class _ControlSiteDatasetContext:
    site_keys: tuple[str, ...]
    site_metadata: pd.DataFrame | None
    dataset_organism: str | None


class ControlSiteMappingContractValidator:
    """Validate mapped control annotations without selecting replacement controls."""

    def run(self, *, mapping: ControlSiteMapping) -> None:
        if not mapping.row_eligibility and not mapping.unmapped_annotations:
            raise PhosPyInputError(
                "control-site validation requires at least one control mapping; "
                "got an empty control-site annotation set"
            )

        missing = tuple(
            row.site_key
            for row in mapping.row_eligibility
            if CONTROL_SITE_REASON_MISSING_CONTROL_ANNOTATION in row.reasons
            and row.site_key is not None
        )
        if missing:
            raise PhosPyInputError(
                "control-site validation failed: missing control mappings for "
                f"dataset site_key values {_format_labels(missing)}"
            )

        duplicates = tuple(
            row.site_key
            for row in mapping.row_eligibility
            if CONTROL_SITE_REASON_DUPLICATE_CONTROL_ANNOTATION in row.reasons
            and row.site_key is not None
        )
        if duplicates:
            raise PhosPyInputError(
                "control-site validation failed: duplicate control mappings for "
                f"site_key values {_format_labels(duplicates)}"
            )

        absent = tuple(
            row.site_key
            for row in mapping.unmapped_annotations
            if CONTROL_SITE_REASON_CONTROL_ANNOTATION_NOT_IN_SITE_ROWS in row.reasons
            and row.site_key is not None
        )
        if absent:
            raise PhosPyInputError(
                "control-site validation failed: control site_key values are absent "
                f"from the dataset: {_format_labels(absent)}"
            )

        invalid = tuple(row for row in mapping.all_eligibility if not row.valid)
        if invalid:
            details = tuple(
                f"{row.site_key or '<missing site_key>'}({', '.join(row.reasons)})"
                for row in invalid
            )
            raise PhosPyInputError(
                "control-site validation failed: invalid control mappings "
                f"{_format_labels(details)}"
            )

        ambiguous = _ambiguous_control_labels(mapping)
        if ambiguous:
            raise PhosPyInputError(
                "control-site validation failed: ambiguous control labels "
                f"{_format_labels(ambiguous)}. A control label must not be used for "
                "both eligible controls and non-control/excluded annotations."
            )


class ControlSiteMetadataCompatibilityValidator:
    """Validate control-source metadata against dataset metadata."""

    def run(
        self,
        *,
        mapping: ControlSiteMapping,
        site_metadata: pd.DataFrame | None,
        dataset_organism: str | None,
        expected_identifier_namespace: str | None,
        control_site_source_type: str | None = None,
    ) -> None:
        site_metadata_by_key = _site_metadata_by_key(site_metadata)
        controls = tuple(row for row in mapping.row_eligibility if row.is_control)
        _reject_ambiguous_control_metadata(controls)
        for row in mapping.row_eligibility:
            if not row.is_control:
                continue
            _validate_control_audit_metadata(
                row,
                control_site_source_type=control_site_source_type,
            )
            expected_organism = _expected_organism_for_row(
                row,
                site_metadata_by_key=site_metadata_by_key,
                dataset_organism=dataset_organism,
            )
            if (
                row.organism is not None
                and expected_organism is not None
                and _normalize_metadata_label(row.organism)
                != _normalize_metadata_label(expected_organism)
            ):
                raise PhosPyInputError(
                    "control-site validation failed: incompatible organism metadata "
                    f"for site_key {row.site_key!r}; control source declares "
                    f"{row.organism!r}, dataset expects {expected_organism!r}"
                )
            if (
                row.identifier_namespace is not None
                and expected_identifier_namespace is not None
                and _normalize_metadata_label(row.identifier_namespace)
                != _normalize_metadata_label(expected_identifier_namespace)
            ):
                raise PhosPyInputError(
                    "control-site validation failed: incompatible identifier "
                    f"namespace metadata for site_key {row.site_key!r}; control "
                    f"source declares {row.identifier_namespace!r}, expected "
                    f"{expected_identifier_namespace!r} because controls are mapped "
                    "by dataset site_key values"
                )


class ControlSiteMethodEligibilityValidator:
    """Validate method-level control count and weighting/grouping support."""

    def run(
        self,
        *,
        mapping: ControlSiteMapping,
        method: str,
        min_eligible_controls: int,
        n_unwanted_factors: int | None,
        supports_weights: bool,
        supports_groups: bool,
        supports_weighted_groups: bool,
    ) -> None:
        if isinstance(min_eligible_controls, bool) or not isinstance(
            min_eligible_controls, int
        ):
            raise PhosPyInputError(
                "control-site validation min_eligible_controls must be an int"
            )
        required_controls = int(min_eligible_controls)
        if required_controls < 1:
            raise PhosPyInputError(
                "control-site validation min_eligible_controls must be greater than "
                "or equal to 1"
            )
        if n_unwanted_factors is not None:
            if isinstance(n_unwanted_factors, bool) or not isinstance(
                n_unwanted_factors, int
            ):
                raise PhosPyInputError(
                    "control-site validation n_unwanted_factors must be an int when "
                    "provided"
                )
            if n_unwanted_factors < 1:
                raise PhosPyInputError(
                    "control-site validation n_unwanted_factors must be greater than "
                    "or equal to 1 when provided"
                )
            required_controls = max(required_controls, int(n_unwanted_factors) + 1)

        controls = tuple(row for row in mapping.row_eligibility if row.is_control)
        control_count = len(controls)
        if control_count < required_controls:
            factor_clause = (
                ""
                if n_unwanted_factors is None
                else f" and n_unwanted_factors={n_unwanted_factors}"
            )
            raise PhosPyInputError(
                "control-site validation failed: too few eligible controls for "
                f"method {method!r}{factor_clause}; required at least "
                f"{required_controls}, observed {control_count}"
            )

        weighted = tuple(row for row in controls if row.weight is not None)
        grouped = tuple(row for row in controls if row.group is not None)
        invalid_weights = tuple(
            row for row in weighted if _has_invalid_control_weight(row)
        )
        if invalid_weights:
            site_keys = tuple(
                row.site_key or "<missing site_key>" for row in invalid_weights
            )
            raise PhosPyInputError(
                "control-site validation failed: control weights must be positive "
                f"finite values for eligible controls; invalid site_key values "
                f"{_format_labels(site_keys)}"
            )
        if weighted and not supports_weights:
            raise PhosPyInputError(
                "control-site validation failed: unsupported control weighting for "
                f"method {method!r}; remove control weights or use a method that "
                "declares weight support"
            )
        if grouped and not supports_groups:
            raise PhosPyInputError(
                "control-site validation failed: unsupported control grouping for "
                f"method {method!r}; remove control groups or use a method that "
                "declares group support"
            )
        if weighted and len(weighted) != control_count:
            raise PhosPyInputError(
                "control-site validation failed: partial control weighting is "
                "unsupported; provide weights for every eligible control or none"
            )
        if grouped and len(grouped) != control_count:
            raise PhosPyInputError(
                "control-site validation failed: partial control grouping is "
                "unsupported; provide groups for every eligible control or none"
            )
        if weighted and grouped and not supports_weighted_groups:
            raise PhosPyInputError(
                "control-site validation failed: unsupported control weighting/"
                f"grouping combination for method {method!r}; weighted grouped "
                "controls are not declared as supported"
            )


class ControlSiteEligibilityValidator:
    """Validate caller-supplied control-site eligibility for batch correction.

    The validator inspects dataset metadata and control-site annotations only. It
    never mutates intensity matrices, drops controls, selects fallback controls,
    or estimates unwanted factors.
    """

    def __init__(
        self,
        *,
        mapping_validator: ControlSiteMappingContractValidator | None = None,
        metadata_validator: ControlSiteMetadataCompatibilityValidator | None = None,
        method_validator: ControlSiteMethodEligibilityValidator | None = None,
    ) -> None:
        self._mapping_validator = (
            mapping_validator or ControlSiteMappingContractValidator()
        )
        self._metadata_validator = (
            metadata_validator or ControlSiteMetadataCompatibilityValidator()
        )
        self._method_validator = (
            method_validator or ControlSiteMethodEligibilityValidator()
        )

    def run(
        self,
        *,
        control_set: ControlSiteSet | None,
        method: str,
        min_eligible_controls: int,
        site_keys: Iterable[object] | None = None,
        phospho: pd.DataFrame | None = None,
        site_metadata: pd.DataFrame | None = None,
        dataset_organism: object | None = None,
        n_unwanted_factors: int | None = None,
        expected_identifier_namespace: str | None = _DEFAULT_IDENTIFIER_NAMESPACE,
        control_site_source_type: str | None = "caller_supplied",
        supports_weights: bool = False,
        supports_groups: bool = False,
        supports_weighted_groups: bool = False,
    ) -> ControlSiteMapping:
        """Return a validated control-site mapping or raise ``PhosPyInputError``."""

        if control_set is None:
            raise PhosPyInputError(
                "control-site validation requires caller-supplied control mappings"
            )
        context = _resolve_dataset_context(
            site_keys=site_keys,
            phospho=phospho,
            site_metadata=site_metadata,
            dataset_organism=dataset_organism,
        )
        mapping = control_set.map_to_site_keys(context.site_keys)
        self._mapping_validator.run(mapping=mapping)
        self._metadata_validator.run(
            mapping=mapping,
            site_metadata=context.site_metadata,
            dataset_organism=context.dataset_organism,
            expected_identifier_namespace=expected_identifier_namespace,
            control_site_source_type=_normalize_metadata_label(
                control_site_source_type
            ),
        )
        self._method_validator.run(
            mapping=mapping,
            method=str(method),
            min_eligible_controls=min_eligible_controls,
            n_unwanted_factors=n_unwanted_factors,
            supports_weights=supports_weights,
            supports_groups=supports_groups,
            supports_weighted_groups=supports_weighted_groups,
        )
        return mapping


def _resolve_dataset_context(
    *,
    site_keys: Iterable[object] | None,
    phospho: pd.DataFrame | None,
    site_metadata: pd.DataFrame | None,
    dataset_organism: object | None,
) -> _ControlSiteDatasetContext:
    if site_keys is not None:
        resolved_site_keys = _normalize_site_key_axis(site_keys)
    elif phospho is not None:
        resolved_site_keys = _normalize_site_key_axis(phospho.index.tolist())
    elif site_metadata is not None:
        resolved_site_keys = _normalize_site_key_axis(site_metadata.index.tolist())
    else:
        raise PhosPyInputError(
            "control-site validation requires dataset site_key values via "
            "site_keys, phospho, or site_metadata"
        )
    resolved_metadata = (
        None if site_metadata is None else site_metadata.copy(deep=False)
    )
    return _ControlSiteDatasetContext(
        site_keys=resolved_site_keys,
        site_metadata=resolved_metadata,
        dataset_organism=_normalize_metadata_label(dataset_organism),
    )


def _normalize_site_key_axis(site_keys: Iterable[object]) -> tuple[str, ...]:
    resolved = tuple(str(site_key).strip() for site_key in site_keys)
    if not resolved:
        raise PhosPyInputError(
            "control-site validation requires at least one dataset site_key value"
        )
    missing_positions = tuple(
        position for position, site_key in enumerate(resolved) if site_key == ""
    )
    if missing_positions:
        raise PhosPyInputError(
            "control-site validation dataset site_key values must not be blank; "
            f"blank positions {_format_positions(missing_positions)}"
        )
    duplicates = _duplicates(resolved)
    if duplicates:
        raise PhosPyInputError(
            "control-site validation dataset site_key values must be unique; "
            f"duplicates {_format_labels(duplicates)}"
        )
    return resolved


def _ambiguous_control_labels(mapping: ControlSiteMapping) -> tuple[str, ...]:
    statuses_by_label: dict[str, set[ControlSiteStatus]] = {}
    for row in mapping.all_eligibility:
        if row.label is None:
            continue
        label = row.label.strip()
        if not label:
            continue
        statuses_by_label.setdefault(label, set()).add(row.control_status)
    ambiguous: list[str] = []
    for label, statuses in statuses_by_label.items():
        if ControlSiteStatus.CONTROL in statuses and len(statuses) > 1:
            ambiguous.append(label)
    return tuple(ambiguous)


def _has_invalid_control_weight(row: ControlSiteEligibility) -> bool:
    weight = row.weight
    return weight is not None and (not math.isfinite(weight) or weight <= 0.0)


def _reject_ambiguous_control_metadata(
    controls: Sequence[ControlSiteEligibility],
) -> None:
    for field_name in _CONTROL_IDENTITY_FIELDS:
        values = tuple(
            value
            for row in controls
            if (value := _normalize_metadata_label(getattr(row, field_name, None)))
            is not None
        )
        unique = tuple(dict.fromkeys(values))
        if len(unique) > 1:
            raise PhosPyInputError(
                "control-site validation failed: ambiguous control metadata for "
                f"{field_name!r}; accepted controls declare {_format_labels(unique)}"
            )


def _validate_control_audit_metadata(
    row: ControlSiteEligibility,
    *,
    control_site_source_type: str | None,
) -> None:
    source_type = _resolve_control_source_type(
        row,
        control_site_source_type=control_site_source_type,
    )
    if _is_strict_control_source_type(source_type):
        missing = tuple(
            field_name
            for field_name in _STRICT_SOURCE_REQUIRED_FIELDS
            if not _has_strict_metadata_field(row, field_name)
        )
        if missing:
            raise PhosPyInputError(
                "control-site validation failed: packaged/reference/external "
                "control-source metadata is incomplete for site_key "
                f"{row.site_key!r}; missing "
                f"{_format_labels(missing)}"
            )
        return

    if _has_external_or_formal_source_name(row) and not _has_metadata_field(
        row,
        "source_version",
    ):
        raise PhosPyInputError(
            "control-site validation failed: missing control metadata for "
            f"site_key {row.site_key!r}; external/formal source_name "
            "requires source_version"
        )

    missing_without_reason = tuple(
        field_name
        for field_name in _CALLER_AUDIT_FIELDS
        if not _has_metadata_field(row, field_name)
        and not _has_missing_metadata_reason(row, field_name)
    )
    has_source_identity = _has_metadata_field(
        row, "source_type"
    ) or _has_metadata_field(
        row,
        "source_name",
    )
    if not has_source_identity:
        missing_without_reason = ("source_type/source_name", *missing_without_reason)
    if missing_without_reason:
        raise PhosPyInputError(
            "control-site validation failed: missing control metadata for site_key "
            f"{row.site_key!r}; provide {_format_labels(missing_without_reason)} "
            "or metadata_missing_reason"
        )


def _resolve_control_source_type(
    row: ControlSiteEligibility,
    *,
    control_site_source_type: str | None,
) -> str:
    row_source_type = _normalize_metadata_label(getattr(row, "source_type", None))
    caller_source_type = _normalize_metadata_label(control_site_source_type)
    source_name = _normalize_metadata_label(getattr(row, "source_name", None))
    for value in (row_source_type, caller_source_type, source_name):
        if value is not None and _is_strict_control_source_type(value):
            return value
    for value in (
        row_source_type,
        caller_source_type,
        source_name,
    ):
        if value is not None:
            return value
    return "caller_supplied"


def _has_metadata_field(row: ControlSiteEligibility, field_name: str) -> bool:
    return _normalize_metadata_label(getattr(row, field_name, None)) is not None


def _has_strict_metadata_field(
    row: ControlSiteEligibility,
    field_name: str,
) -> bool:
    value = _normalize_metadata_label(getattr(row, field_name, None))
    if value is None:
        return False
    if field_name in {"source_name", "selection_method"}:
        return value not in _CALLER_SUPPLIED_SOURCE_TYPES
    return True


def _is_strict_control_source_type(source_type: str | None) -> bool:
    if source_type is None:
        return False
    tokens = frozenset(source_type.replace("-", "_").split("_"))
    return bool(tokens & _STRICT_SOURCE_TYPE_MARKERS)


def _has_missing_metadata_reason(
    row: ControlSiteEligibility,
    field_name: str,
) -> bool:
    reasons = getattr(row, "metadata_missing_reason", {})
    if not isinstance(reasons, Mapping):
        return False
    reason = reasons.get(field_name)
    return _normalize_metadata_label(reason) is not None


def _has_external_or_formal_source_name(row: ControlSiteEligibility) -> bool:
    source_type = _normalize_metadata_label(getattr(row, "source_type", None))
    source_name = _normalize_metadata_label(getattr(row, "source_name", None))
    if source_type is not None and source_type not in _CALLER_SUPPLIED_SOURCE_TYPES:
        return True
    if source_name is None:
        return False
    return source_name not in _CALLER_SUPPLIED_SOURCE_TYPES


def _site_metadata_by_key(
    site_metadata: pd.DataFrame | None,
) -> dict[str, dict[str, str]]:
    if site_metadata is None:
        return {}
    metadata_by_key: dict[str, dict[str, str]] = {}
    for position, index_value in enumerate(site_metadata.index.tolist()):
        site_key = str(index_value).strip()
        if site_key == "":
            continue
        row: dict[str, str] = {}
        for column in ("organism", "identifier_namespace"):
            if column not in site_metadata.columns:
                continue
            value = _normalize_metadata_label(site_metadata.iloc[position][column])
            if value is not None:
                row[column] = value
        metadata_by_key[site_key] = row
    return metadata_by_key


def _expected_organism_for_row(
    row: ControlSiteEligibility,
    *,
    site_metadata_by_key: dict[str, dict[str, str]],
    dataset_organism: str | None,
) -> str | None:
    if row.site_key is not None:
        site_row = site_metadata_by_key.get(row.site_key)
        if site_row is not None and "organism" in site_row:
            return site_row["organism"]
    return dataset_organism


def _normalize_metadata_label(value: object | None) -> str | None:
    if value is None:
        return None
    if isinstance(value, Enum):
        value = value.value
    if bool(pd.Series((value,), dtype="object").isna().iat[0]):
        return None
    text = str(value).strip()
    return text.lower() if text else None


def _duplicates(values: Sequence[str]) -> tuple[str, ...]:
    seen: set[str] = set()
    duplicates: list[str] = []
    for value in values:
        if value in seen and value not in duplicates:
            duplicates.append(value)
        seen.add(value)
    return tuple(duplicates)


def _format_positions(positions: Sequence[int]) -> str:
    preview = ", ".join(str(position) for position in tuple(positions)[:8])
    suffix = "" if len(positions) <= 8 else ", ..."
    return f"[{preview}{suffix}]"


def _format_labels(labels: Sequence[object]) -> str:
    preview = ", ".join(repr(value) for value in tuple(labels)[:5])
    suffix = "" if len(labels) <= 5 else " ..."
    return f"{preview}{suffix}"


__all__ = [
    "ControlSiteEligibilityValidator",
    "ControlSiteMappingContractValidator",
    "ControlSiteMetadataCompatibilityValidator",
    "ControlSiteMethodEligibilityValidator",
]
