"""Explicit control-site representation and site-key mapping models."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import TypeAlias, cast

from phospy.policies import PolicyEnum
from phospy.science.configs.preprocessing.control_sites import (
    CONTROL_SITE_SELECTION_METHOD_CALLER_SUPPLIED,
    CONTROL_SITE_SOURCE_CALLER_SUPPLIED,
    CONTROL_SITE_STATUS_CONTROL,
    CONTROL_SITE_STATUS_EXCLUDED,
    CONTROL_SITE_STATUS_INVALID,
    CONTROL_SITE_STATUS_NON_CONTROL,
    CONTROL_SITE_STATUS_UNKNOWN,
)

CONTROL_SITE_REASON_DUPLICATE_CONTROL_ANNOTATION = "duplicate_control_annotation"
CONTROL_SITE_REASON_INVALID_CONTROL_STATUS = "invalid_control_status"
CONTROL_SITE_REASON_INVALID_WEIGHT = "invalid_weight"
CONTROL_SITE_REASON_MISSING_CONTROL_ANNOTATION = "missing_control_annotation"
CONTROL_SITE_REASON_MISSING_SITE_KEY = "missing_site_key"
CONTROL_SITE_REASON_CONTROL_ANNOTATION_NOT_IN_SITE_ROWS = (
    "control_annotation_not_in_site_rows"
)

ControlSiteReason: TypeAlias = str


class ControlSiteStatus(PolicyEnum):
    """Stable per-site control annotation status."""

    CONTROL = CONTROL_SITE_STATUS_CONTROL
    NON_CONTROL = CONTROL_SITE_STATUS_NON_CONTROL
    EXCLUDED = CONTROL_SITE_STATUS_EXCLUDED
    UNKNOWN = CONTROL_SITE_STATUS_UNKNOWN
    INVALID = CONTROL_SITE_STATUS_INVALID


@dataclass(frozen=True, slots=True)
class ControlSiteSourceMetadata:
    """Source metadata attached to caller-supplied control-site annotations."""

    source_type: str | None = CONTROL_SITE_SOURCE_CALLER_SUPPLIED
    organism: str | None = None
    identifier_namespace: str | None = None
    source_name: str | None = CONTROL_SITE_SOURCE_CALLER_SUPPLIED
    source_version: str | None = None
    license: str | None = None
    redistribution: str | None = None
    selection_method: str | None = CONTROL_SITE_SELECTION_METHOD_CALLER_SUPPLIED
    metadata_missing_reason: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for field_name in (
            "source_type",
            "organism",
            "identifier_namespace",
            "source_name",
            "source_version",
            "license",
            "redistribution",
            "selection_method",
        ):
            object.__setattr__(
                self,
                field_name,
                _optional_non_empty_string(getattr(self, field_name)),
            )
        object.__setattr__(
            self,
            "metadata_missing_reason",
            _normalize_missing_reason_mapping(self.metadata_missing_reason),
        )

    def to_payload(self) -> dict[str, object]:
        """Return JSON-compatible metadata."""

        return {
            "source_type": self.source_type,
            "organism": self.organism,
            "identifier_namespace": self.identifier_namespace,
            "source_name": self.source_name,
            "source_version": self.source_version,
            "license": self.license,
            "redistribution": self.redistribution,
            "selection_method": self.selection_method,
            "metadata_missing_reason": dict(self.metadata_missing_reason),
        }


@dataclass(frozen=True, slots=True)
class ControlSiteAnnotation:
    """One caller-supplied control-site annotation.

    Local scalar issues are stored as ``structural_reasons`` so later validation
    can decide whether to reject, repair, or report them.
    """

    site_key: object | None
    control_status: object = ControlSiteStatus.CONTROL
    label: object | None = None
    weight: object | None = None
    group: object | None = None
    source_type: object | None = None
    organism: object | None = None
    identifier_namespace: object | None = None
    source_name: object | None = None
    source_version: object | None = None
    license: object | None = None
    redistribution: object | None = None
    selection_method: object | None = None
    metadata_missing_reason: Mapping[str, str] = field(default_factory=dict)
    exclusion_reason: object | None = None
    structural_reasons: Sequence[str] = ()

    def __post_init__(self) -> None:
        reasons = list(_normalize_reason_tuple(self.structural_reasons))
        site_key = _optional_non_empty_string(self.site_key)
        if site_key is None:
            reasons.append(CONTROL_SITE_REASON_MISSING_SITE_KEY)
        status, status_reasons = _coerce_control_site_status(self.control_status)
        reasons.extend(status_reasons)
        weight, weight_reasons = _coerce_optional_weight(self.weight)
        reasons.extend(weight_reasons)

        object.__setattr__(self, "site_key", site_key)
        object.__setattr__(self, "control_status", status)
        object.__setattr__(self, "label", _optional_non_empty_string(self.label))
        object.__setattr__(self, "weight", weight)
        object.__setattr__(self, "group", _optional_non_empty_string(self.group))
        object.__setattr__(
            self,
            "source_type",
            _optional_non_empty_string(self.source_type),
        )
        object.__setattr__(
            self,
            "organism",
            _optional_non_empty_string(self.organism),
        )
        object.__setattr__(
            self,
            "identifier_namespace",
            _optional_non_empty_string(self.identifier_namespace),
        )
        object.__setattr__(
            self,
            "source_name",
            _optional_non_empty_string(self.source_name),
        )
        object.__setattr__(
            self,
            "source_version",
            _optional_non_empty_string(self.source_version),
        )
        object.__setattr__(self, "license", _optional_non_empty_string(self.license))
        object.__setattr__(
            self,
            "redistribution",
            _optional_non_empty_string(self.redistribution),
        )
        object.__setattr__(
            self,
            "selection_method",
            _optional_non_empty_string(self.selection_method),
        )
        object.__setattr__(
            self,
            "metadata_missing_reason",
            _normalize_missing_reason_mapping(self.metadata_missing_reason),
        )
        object.__setattr__(
            self,
            "exclusion_reason",
            _optional_non_empty_string(self.exclusion_reason),
        )
        object.__setattr__(self, "structural_reasons", _dedupe_reasons(reasons))

    @property
    def structurally_valid(self) -> bool:
        return (
            not self.structural_reasons
            and _annotation_control_status(self) is not ControlSiteStatus.INVALID
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible annotation payload."""

        return {
            "site_key": self.site_key,
            "control_status": _annotation_control_status(self).value,
            "label": self.label,
            "weight": self.weight,
            "group": self.group,
            "source_type": self.source_type,
            "organism": self.organism,
            "identifier_namespace": self.identifier_namespace,
            "source_name": self.source_name,
            "source_version": self.source_version,
            "license": self.license,
            "redistribution": self.redistribution,
            "selection_method": self.selection_method,
            "metadata_missing_reason": dict(self.metadata_missing_reason),
            "exclusion_reason": self.exclusion_reason,
            "structural_reasons": list(self.structural_reasons),
            "structurally_valid": self.structurally_valid,
        }


@dataclass(frozen=True, slots=True)
class ControlSiteEligibility:
    """Mapped control-site state for one site-key row or one unmapped annotation."""

    site_key: str | None
    control_status: ControlSiteStatus
    valid: bool
    reasons: tuple[ControlSiteReason, ...] = ()
    label: str | None = None
    weight: float | None = None
    group: str | None = None
    source_type: str | None = None
    organism: str | None = None
    identifier_namespace: str | None = None
    source_name: str | None = None
    source_version: str | None = None
    license: str | None = None
    redistribution: str | None = None
    selection_method: str | None = None
    metadata_missing_reason: Mapping[str, str] = field(default_factory=dict)
    exclusion_reason: str | None = None
    row_position: int | None = None
    annotation_indices: tuple[int, ...] = ()

    @property
    def is_control(self) -> bool:
        return self.valid and self.control_status is ControlSiteStatus.CONTROL

    @property
    def is_weighted_control(self) -> bool:
        return self.is_control and self.weight is not None

    @property
    def annotation_count(self) -> int:
        return len(self.annotation_indices)

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible mapped eligibility payload."""

        return {
            "site_key": self.site_key,
            "control_status": self.control_status.value,
            "valid": self.valid,
            "reasons": list(self.reasons),
            "label": self.label,
            "weight": self.weight,
            "group": self.group,
            "source_type": self.source_type,
            "organism": self.organism,
            "identifier_namespace": self.identifier_namespace,
            "source_name": self.source_name,
            "source_version": self.source_version,
            "license": self.license,
            "redistribution": self.redistribution,
            "selection_method": self.selection_method,
            "metadata_missing_reason": dict(self.metadata_missing_reason),
            "exclusion_reason": self.exclusion_reason,
            "row_position": self.row_position,
            "annotation_indices": list(self.annotation_indices),
            "annotation_count": self.annotation_count,
        }


@dataclass(frozen=True, slots=True)
class ControlSiteMapping:
    """Control-site annotations mapped to a target site-key axis."""

    row_eligibility: tuple[ControlSiteEligibility, ...]
    unmapped_annotations: tuple[ControlSiteEligibility, ...] = ()

    @property
    def all_eligibility(self) -> tuple[ControlSiteEligibility, ...]:
        return self.row_eligibility + self.unmapped_annotations

    @property
    def control_status_by_site_key(self) -> dict[str, ControlSiteStatus]:
        return {
            row.site_key: row.control_status
            for row in self.row_eligibility
            if row.site_key is not None
        }

    @property
    def control_weight_by_site_key(self) -> dict[str, float]:
        return {
            row.site_key: row.weight
            for row in self.row_eligibility
            if row.site_key is not None and row.weight is not None
        }

    @property
    def grouped_control_site_keys(self) -> dict[str, tuple[str, ...]]:
        grouped: dict[str, list[str]] = defaultdict(list)
        for row in self.row_eligibility:
            if row.is_control and row.group is not None and row.site_key is not None:
                grouped[row.group].append(row.site_key)
        return {group: tuple(site_keys) for group, site_keys in grouped.items()}

    @property
    def invalid(self) -> tuple[ControlSiteEligibility, ...]:
        return tuple(row for row in self.all_eligibility if not row.valid)

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible mapping payload."""

        return {
            "row_eligibility": [row.to_payload() for row in self.row_eligibility],
            "unmapped_annotations": [
                row.to_payload() for row in self.unmapped_annotations
            ],
            "invalid": [row.to_payload() for row in self.invalid],
        }


@dataclass(frozen=True, slots=True)
class ControlSiteSet:
    """Caller-supplied control-site annotations and source metadata."""

    annotations: Sequence[ControlSiteAnnotation] = ()
    source_metadata: ControlSiteSourceMetadata = field(
        default_factory=ControlSiteSourceMetadata
    )
    unannotated_site_status: object = ControlSiteStatus.UNKNOWN

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "annotations",
            _normalize_annotations(self.annotations),
        )
        source_metadata = self.source_metadata
        if not isinstance(source_metadata, ControlSiteSourceMetadata):
            source_metadata = ControlSiteSourceMetadata()
        object.__setattr__(self, "source_metadata", source_metadata)
        status, _ = _coerce_control_site_status(self.unannotated_site_status)
        if status is ControlSiteStatus.INVALID:
            status = ControlSiteStatus.UNKNOWN
        object.__setattr__(self, "unannotated_site_status", status)

    @classmethod
    def from_site_keys(
        cls,
        site_keys: Iterable[object],
        *,
        source_metadata: ControlSiteSourceMetadata | None = None,
        label: str | None = "caller_supplied_control",
        group: str | None = None,
    ) -> ControlSiteSet:
        """Build a caller-supplied control set from control site keys only."""

        return cls(
            annotations=tuple(
                ControlSiteAnnotation(
                    site_key=site_key,
                    control_status=ControlSiteStatus.CONTROL,
                    label=label,
                    group=group,
                )
                for site_key in site_keys
            ),
            source_metadata=source_metadata or ControlSiteSourceMetadata(),
            unannotated_site_status=ControlSiteStatus.NON_CONTROL,
        )

    @classmethod
    def from_binary_controls(
        cls,
        controls: Mapping[object, object],
        *,
        source_metadata: ControlSiteSourceMetadata | None = None,
    ) -> ControlSiteSet:
        """Build a control set from ``site_key -> bool/status`` annotations."""

        return cls(
            annotations=tuple(
                ControlSiteAnnotation(site_key=site_key, control_status=status)
                for site_key, status in controls.items()
            ),
            source_metadata=source_metadata or ControlSiteSourceMetadata(),
            unannotated_site_status=ControlSiteStatus.UNKNOWN,
        )

    @classmethod
    def from_weighted_controls(
        cls,
        controls: Mapping[object, object],
        *,
        source_metadata: ControlSiteSourceMetadata | None = None,
        group: str | None = None,
    ) -> ControlSiteSet:
        """Build a caller-supplied weighted control set."""

        return cls(
            annotations=tuple(
                ControlSiteAnnotation(
                    site_key=site_key,
                    control_status=ControlSiteStatus.CONTROL,
                    weight=weight,
                    group=group,
                )
                for site_key, weight in controls.items()
            ),
            source_metadata=source_metadata or ControlSiteSourceMetadata(),
            unannotated_site_status=ControlSiteStatus.NON_CONTROL,
        )

    def map_to_site_keys(self, site_keys: Iterable[object]) -> ControlSiteMapping:
        """Map annotations onto an existing site-key axis without selecting sites."""

        target_site_keys = tuple(_optional_non_empty_string(site) for site in site_keys)
        annotations_by_site_key: dict[str, list[tuple[int, ControlSiteAnnotation]]] = (
            defaultdict(list)
        )
        for index, annotation in enumerate(self.annotations):
            annotation_site_key = _annotation_site_key(annotation)
            if annotation_site_key is None:
                continue
            annotations_by_site_key[annotation_site_key].append((index, annotation))

        row_eligibility = tuple(
            self._eligibility_for_site_row(
                site_key=site_key,
                row_position=row_position,
                annotations=annotations_by_site_key.get(site_key or "", []),
            )
            for row_position, site_key in enumerate(target_site_keys)
        )
        target_site_key_set = {
            site_key for site_key in target_site_keys if site_key is not None
        }
        unmapped_annotations = tuple(
            self._unmapped_annotation_eligibility(
                annotation=annotation,
                annotation_index=annotation_index,
            )
            for annotation_index, annotation in enumerate(self.annotations)
            if _annotation_site_key(annotation) is None
            or _annotation_site_key(annotation) not in target_site_key_set
        )
        return ControlSiteMapping(
            row_eligibility=row_eligibility,
            unmapped_annotations=unmapped_annotations,
        )

    def _eligibility_for_site_row(
        self,
        *,
        site_key: str | None,
        row_position: int,
        annotations: Sequence[tuple[int, ControlSiteAnnotation]],
    ) -> ControlSiteEligibility:
        if site_key is None:
            return ControlSiteEligibility(
                site_key=None,
                control_status=ControlSiteStatus.INVALID,
                valid=False,
                reasons=(CONTROL_SITE_REASON_MISSING_SITE_KEY,),
                row_position=row_position,
            )
        if not annotations:
            missing_status = self.unannotated_site_status
            missing_status = cast(ControlSiteStatus, missing_status)
            valid = missing_status is not ControlSiteStatus.UNKNOWN
            reasons: tuple[str, ...] = (
                () if valid else (CONTROL_SITE_REASON_MISSING_CONTROL_ANNOTATION,)
            )
            metadata = self._metadata_payload()
            return ControlSiteEligibility(
                site_key=site_key,
                control_status=missing_status,
                valid=valid,
                reasons=reasons,
                row_position=row_position,
                source_type=metadata.source_type,
                organism=metadata.organism,
                identifier_namespace=metadata.identifier_namespace,
                source_name=metadata.source_name,
                source_version=metadata.source_version,
                license=metadata.license,
                redistribution=metadata.redistribution,
                selection_method=metadata.selection_method,
                metadata_missing_reason=metadata.metadata_missing_reason,
            )
        if len(annotations) > 1:
            reasons = _dedupe_reasons(
                [
                    CONTROL_SITE_REASON_DUPLICATE_CONTROL_ANNOTATION,
                    *(
                        reason
                        for _, annotation in annotations
                        for reason in annotation.structural_reasons
                    ),
                ]
            )
            metadata = self._metadata_payload()
            return ControlSiteEligibility(
                site_key=site_key,
                control_status=ControlSiteStatus.INVALID,
                valid=False,
                reasons=reasons,
                row_position=row_position,
                annotation_indices=tuple(index for index, _ in annotations),
                source_type=metadata.source_type,
                organism=metadata.organism,
                identifier_namespace=metadata.identifier_namespace,
                source_name=metadata.source_name,
                source_version=metadata.source_version,
                license=metadata.license,
                redistribution=metadata.redistribution,
                selection_method=metadata.selection_method,
                metadata_missing_reason=metadata.metadata_missing_reason,
            )
        annotation_index, annotation = annotations[0]
        return self._eligibility_from_annotation(
            annotation=annotation,
            annotation_index=annotation_index,
            row_position=row_position,
        )

    def _eligibility_from_annotation(
        self,
        *,
        annotation: ControlSiteAnnotation,
        annotation_index: int,
        row_position: int | None,
    ) -> ControlSiteEligibility:
        valid = annotation.structurally_valid
        metadata = self._metadata_payload(annotation)
        return ControlSiteEligibility(
            site_key=_annotation_site_key(annotation),
            control_status=_annotation_control_status(annotation),
            valid=valid,
            reasons=_annotation_structural_reasons(annotation),
            label=_annotation_label(annotation),
            weight=_annotation_weight(annotation),
            group=_annotation_group(annotation),
            exclusion_reason=_annotation_exclusion_reason(annotation),
            row_position=row_position,
            annotation_indices=(annotation_index,),
            source_type=metadata.source_type,
            organism=metadata.organism,
            identifier_namespace=metadata.identifier_namespace,
            source_name=metadata.source_name,
            source_version=metadata.source_version,
            license=metadata.license,
            redistribution=metadata.redistribution,
            selection_method=metadata.selection_method,
            metadata_missing_reason=metadata.metadata_missing_reason,
        )

    def _unmapped_annotation_eligibility(
        self,
        *,
        annotation: ControlSiteAnnotation,
        annotation_index: int,
    ) -> ControlSiteEligibility:
        reasons = list(annotation.structural_reasons)
        annotation_site_key = _annotation_site_key(annotation)
        if annotation_site_key is not None:
            reasons.append(CONTROL_SITE_REASON_CONTROL_ANNOTATION_NOT_IN_SITE_ROWS)
        metadata = self._metadata_payload(annotation)
        return ControlSiteEligibility(
            site_key=annotation_site_key,
            control_status=_annotation_control_status(annotation),
            valid=False,
            reasons=_dedupe_reasons(reasons),
            label=_annotation_label(annotation),
            weight=_annotation_weight(annotation),
            group=_annotation_group(annotation),
            exclusion_reason=_annotation_exclusion_reason(annotation),
            row_position=None,
            annotation_indices=(annotation_index,),
            source_type=metadata.source_type,
            organism=metadata.organism,
            identifier_namespace=metadata.identifier_namespace,
            source_name=metadata.source_name,
            source_version=metadata.source_version,
            license=metadata.license,
            redistribution=metadata.redistribution,
            selection_method=metadata.selection_method,
            metadata_missing_reason=metadata.metadata_missing_reason,
        )

    def _metadata_payload(
        self,
        annotation: ControlSiteAnnotation | None = None,
    ) -> ControlSiteSourceMetadata:
        metadata = self.source_metadata
        return ControlSiteSourceMetadata(
            source_type=_first_non_none(
                annotation,
                "source_type",
                metadata.source_type,
            ),
            organism=_first_non_none(annotation, "organism", metadata.organism),
            identifier_namespace=_first_non_none(
                annotation,
                "identifier_namespace",
                metadata.identifier_namespace,
            ),
            source_name=_first_non_none(
                annotation,
                "source_name",
                metadata.source_name,
            ),
            source_version=_first_non_none(
                annotation,
                "source_version",
                metadata.source_version,
            ),
            license=_first_non_none(annotation, "license", metadata.license),
            redistribution=_first_non_none(
                annotation,
                "redistribution",
                metadata.redistribution,
            ),
            selection_method=_first_non_none(
                annotation,
                "selection_method",
                metadata.selection_method,
            ),
            metadata_missing_reason=_merged_missing_reason_mapping(
                metadata.metadata_missing_reason,
                () if annotation is None else annotation.metadata_missing_reason,
            ),
        )


def _normalize_annotations(
    annotations: Sequence[ControlSiteAnnotation],
) -> tuple[ControlSiteAnnotation, ...]:
    return tuple(
        annotation
        if isinstance(annotation, ControlSiteAnnotation)
        else ControlSiteAnnotation(
            site_key=None,
            control_status=ControlSiteStatus.INVALID,
            structural_reasons=("invalid_control_site_annotation",),
        )
        for annotation in annotations
    )


def _coerce_control_site_status(
    value: object,
) -> tuple[ControlSiteStatus, tuple[str, ...]]:
    if isinstance(value, ControlSiteStatus):
        return value, ()
    if isinstance(value, bool):
        return (
            ControlSiteStatus.CONTROL if value else ControlSiteStatus.NON_CONTROL,
            (),
        )
    if value is None:
        return ControlSiteStatus.UNKNOWN, ()
    normalized = str(value).strip().lower()
    aliases = {
        "true": ControlSiteStatus.CONTROL.value,
        "false": ControlSiteStatus.NON_CONTROL.value,
        "stable": ControlSiteStatus.CONTROL.value,
        "unstable": ControlSiteStatus.NON_CONTROL.value,
        "non-control": ControlSiteStatus.NON_CONTROL.value,
        "not_control": ControlSiteStatus.NON_CONTROL.value,
    }
    normalized = aliases.get(normalized, normalized)
    try:
        return ControlSiteStatus(normalized), ()
    except ValueError:
        return (
            ControlSiteStatus.INVALID,
            (CONTROL_SITE_REASON_INVALID_CONTROL_STATUS,),
        )


def _coerce_optional_weight(
    value: object | None,
) -> tuple[float | None, tuple[str, ...]]:
    if value is None:
        return None, ()
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None, (CONTROL_SITE_REASON_INVALID_WEIGHT,)
    weight = float(value)
    if not math.isfinite(weight) or weight < 0:
        return None, (CONTROL_SITE_REASON_INVALID_WEIGHT,)
    return weight, ()


def _optional_non_empty_string(value: object | None) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if text == "":
        return None
    return text


def _normalize_reason_tuple(reasons: Sequence[str]) -> tuple[str, ...]:
    if isinstance(reasons, str):
        return (reasons.strip(),) if reasons.strip() else ()
    return _dedupe_reasons(
        str(reason).strip() for reason in tuple(reasons) if str(reason).strip()
    )


def _normalize_missing_reason_mapping(
    reasons: Mapping[str, str] | Sequence[tuple[str, str]] | None,
) -> dict[str, str]:
    if reasons is None:
        return {}
    items = reasons.items() if isinstance(reasons, Mapping) else tuple(reasons)
    normalized: dict[str, str] = {}
    for key, reason in items:
        field_name = _optional_non_empty_string(key)
        reason_text = _optional_non_empty_string(reason)
        if field_name is None or reason_text is None:
            continue
        normalized[field_name] = reason_text
    return normalized


def _merged_missing_reason_mapping(
    base: Mapping[str, str],
    override: Mapping[str, str] | Sequence[tuple[str, str]],
) -> dict[str, str]:
    merged = dict(base)
    merged.update(_normalize_missing_reason_mapping(override))
    return merged


def _dedupe_reasons(reasons: Iterable[str]) -> tuple[str, ...]:
    seen: set[str] = set()
    deduped: list[str] = []
    for reason in reasons:
        if reason in seen:
            continue
        seen.add(reason)
        deduped.append(reason)
    return tuple(deduped)


def _first_non_none(
    annotation: ControlSiteAnnotation | None,
    field_name: str,
    fallback: str | None,
) -> str | None:
    if annotation is None:
        return fallback
    value = getattr(annotation, field_name)
    return fallback if value is None else value


def _annotation_site_key(annotation: ControlSiteAnnotation) -> str | None:
    return cast(str | None, annotation.site_key)


def _annotation_control_status(annotation: ControlSiteAnnotation) -> ControlSiteStatus:
    return cast(ControlSiteStatus, annotation.control_status)


def _annotation_label(annotation: ControlSiteAnnotation) -> str | None:
    return cast(str | None, annotation.label)


def _annotation_weight(annotation: ControlSiteAnnotation) -> float | None:
    return cast(float | None, annotation.weight)


def _annotation_group(annotation: ControlSiteAnnotation) -> str | None:
    return cast(str | None, annotation.group)


def _annotation_exclusion_reason(annotation: ControlSiteAnnotation) -> str | None:
    return cast(str | None, annotation.exclusion_reason)


def _annotation_structural_reasons(
    annotation: ControlSiteAnnotation,
) -> tuple[ControlSiteReason, ...]:
    return cast(tuple[ControlSiteReason, ...], annotation.structural_reasons)


__all__ = [
    "CONTROL_SITE_REASON_CONTROL_ANNOTATION_NOT_IN_SITE_ROWS",
    "CONTROL_SITE_REASON_DUPLICATE_CONTROL_ANNOTATION",
    "CONTROL_SITE_REASON_INVALID_CONTROL_STATUS",
    "CONTROL_SITE_REASON_INVALID_WEIGHT",
    "CONTROL_SITE_REASON_MISSING_CONTROL_ANNOTATION",
    "CONTROL_SITE_REASON_MISSING_SITE_KEY",
    "ControlSiteAnnotation",
    "ControlSiteEligibility",
    "ControlSiteMapping",
    "ControlSiteReason",
    "ControlSiteSet",
    "ControlSiteSourceMetadata",
    "ControlSiteStatus",
]
