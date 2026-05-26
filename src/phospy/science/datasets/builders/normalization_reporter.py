"""Diagnostics and label-normalisation helpers for dataset convention normalisation."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import UnsupportedInputFormatError
from phospy.science.sites.identifiers import (
    SiteIdentifierNormalisationRecord,
    SiteIdentifierNormalisationReport,
    build_site_identifier_normalisation_report,
    canonicalize_site_index,
)


@dataclass(frozen=True, slots=True)
class IndexLabelNormalizationPolicy:
    allow_non_string_labels: bool = True
    detect_duplicate_labels_after_normalisation: bool = True


@dataclass(frozen=True, slots=True)
class _NormalizedIndexLabels:
    index: pd.Index
    raw_labels: tuple[str, ...]


SITE_IDENTIFIER_INDEX_POLICY = IndexLabelNormalizationPolicy(
    allow_non_string_labels=True,
    detect_duplicate_labels_after_normalisation=True,
)
SAMPLE_LABEL_INDEX_POLICY = IndexLabelNormalizationPolicy(
    allow_non_string_labels=True,
    detect_duplicate_labels_after_normalisation=True,
)


class DatasetConventionNormalisationReporter:
    """Own index-label normalisation and site-identifier diagnostics payloads."""

    def normalize_index_labels(
        self,
        index: pd.Index,
        *,
        field_name: str,
        policy: IndexLabelNormalizationPolicy,
    ) -> pd.Index:
        normalized = self._validate_and_normalize_index_labels(
            index,
            field_name=field_name,
            policy=policy,
        )
        return normalized.index

    def normalize_supported_site_index_if_present(
        self,
        index: pd.Index,
        *,
        field_name: str,
        site_identifier_records: list[SiteIdentifierNormalisationRecord],
    ) -> pd.Index:
        normalized = self._validate_and_normalize_index_labels(
            index,
            field_name=field_name,
            policy=SITE_IDENTIFIER_INDEX_POLICY,
        )
        if not _has_site_like_tokens(normalized.index):
            return normalized.index
        canonical = canonicalize_site_index(
            normalized.index,
            field_name=field_name,
            error_type=UnsupportedInputFormatError,
            require_unique=False,
        )
        site_identifier_records.extend(
            _site_identifier_normalisation_changes(
                raw_labels=normalized.raw_labels,
                normalized_labels=canonical,
                field_name=field_name,
            )
        )
        return canonical

    def canonicalize_site_index_with_label_validation(
        self,
        index: pd.Index,
        *,
        field_name: str,
        site_identifier_records: list[SiteIdentifierNormalisationRecord],
        index_name: str | None = None,
    ) -> pd.Index:
        normalized = self._validate_and_normalize_index_labels(
            index,
            field_name=field_name,
            policy=SITE_IDENTIFIER_INDEX_POLICY,
        )
        canonical = canonicalize_site_index(
            normalized.index,
            field_name=field_name,
            error_type=UnsupportedInputFormatError,
            require_unique=False,
            index_name=index_name,
        )
        site_identifier_records.extend(
            _site_identifier_normalisation_changes(
                raw_labels=normalized.raw_labels,
                normalized_labels=canonical,
                field_name=field_name,
            )
        )
        return canonical

    @staticmethod
    def build_site_identifier_report(
        records: list[SiteIdentifierNormalisationRecord],
    ) -> SiteIdentifierNormalisationReport | None:
        return build_site_identifier_normalisation_report(records)

    def _validate_and_normalize_index_labels(
        self,
        index: pd.Index,
        *,
        field_name: str,
        policy: IndexLabelNormalizationPolicy,
    ) -> _NormalizedIndexLabels:
        raw_objects = index.tolist()
        if not raw_objects:
            return _NormalizedIndexLabels(index=index.copy(), raw_labels=())

        raw_labels: list[str] = []
        normalized_labels: list[str] = []
        missing_positions: list[int] = []
        blank_positions: list[int] = []
        non_string_positions: list[int] = []

        for position, value in enumerate(raw_objects):
            if _is_missing_label(value):
                missing_positions.append(position)
                continue
            if not isinstance(value, str) and not policy.allow_non_string_labels:
                non_string_positions.append(position)
                continue
            raw_label = str(value)
            normalized_label = raw_label.strip()
            if normalized_label == "":
                blank_positions.append(position)
                continue
            raw_labels.append(raw_label)
            normalized_labels.append(normalized_label)

        if missing_positions:
            raise UnsupportedInputFormatError(
                f"{field_name} must not contain missing labels; found missing labels "
                f"at positions: {_position_preview(missing_positions)}"
            )
        if non_string_positions:
            raise UnsupportedInputFormatError(
                f"{field_name} must contain string labels; found non-string labels at "
                f"positions: {_position_preview(non_string_positions)}"
            )
        if blank_positions:
            raise UnsupportedInputFormatError(
                f"{field_name} must contain non-blank labels; found blank labels at "
                f"positions: {_position_preview(blank_positions)}"
            )

        normalized_index = pd.Index(normalized_labels, name=index.name)
        if policy.detect_duplicate_labels_after_normalisation:
            _raise_if_duplicate_labels_introduced_by_normalisation(
                raw_labels=raw_labels,
                normalized_labels=normalized_labels,
                field_name=field_name,
            )
        return _NormalizedIndexLabels(
            index=normalized_index,
            raw_labels=tuple(raw_labels),
        )


def _raise_if_duplicate_labels_introduced_by_normalisation(
    *,
    raw_labels: list[str],
    normalized_labels: list[str],
    field_name: str,
) -> None:
    normalized_index = pd.Index(normalized_labels)
    if normalized_index.is_unique:
        return
    duplicate_labels = list(
        dict.fromkeys(normalized_index[normalized_index.duplicated()])
    )
    introduced_by_normalisation: list[str] = []
    for label in duplicate_labels:
        raw_variants = {
            raw_label
            for raw_label, normalized_label in zip(
                raw_labels, normalized_labels, strict=False
            )
            if normalized_label == label
        }
        if raw_variants != {label}:
            introduced_by_normalisation.append(label)
    if not introduced_by_normalisation:
        return
    preview = ", ".join(repr(label) for label in introduced_by_normalisation[:5])
    suffix = "" if len(introduced_by_normalisation) <= 5 else " ..."
    raise UnsupportedInputFormatError(
        f"{field_name} contains duplicate labels introduced by normalization: "
        f"{preview}{suffix}. Provide unique labels after trimming whitespace."
    )


def _site_identifier_normalisation_changes(
    *,
    raw_labels: tuple[str, ...],
    normalized_labels: pd.Index,
    field_name: str,
) -> tuple[SiteIdentifierNormalisationRecord, ...]:
    records: list[SiteIdentifierNormalisationRecord] = []
    for position, (raw_label, normalized_label) in enumerate(
        zip(raw_labels, normalized_labels.tolist(), strict=False)
    ):
        if raw_label == normalized_label:
            continue
        records.append(
            SiteIdentifierNormalisationRecord(
                field_name=field_name,
                row_position=position,
                original_value=raw_label,
                normalised_value=str(normalized_label),
            )
        )
    return tuple(records)


def _position_preview(positions: list[int]) -> str:
    preview = ", ".join(str(position) for position in positions[:8])
    suffix = "" if len(positions) <= 8 else ", ..."
    return f"[{preview}{suffix}]"


def _is_missing_label(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _has_site_like_tokens(index: pd.Index) -> bool:
    for value in index.tolist():
        if ";" in str(value):
            return True
    return False
