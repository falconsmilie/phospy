"""Resolve batch-correction sample metadata below the dataset boundary."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import PhosPyInputError


@dataclass(frozen=True, slots=True)
class ResolvedBatchCorrectionMetadata:
    """Batch and condition labels aligned to phospho matrix sample order."""

    batch_by_sample: Mapping[str, str]
    condition_by_sample: Mapping[str, str]
    sample_order: tuple[str, ...]

    @property
    def batch_labels(self) -> tuple[str, ...]:
        return tuple(self.batch_by_sample[sample] for sample in self.sample_order)

    @property
    def condition_labels(self) -> tuple[str, ...]:
        return tuple(self.condition_by_sample[sample] for sample in self.sample_order)


class BatchCorrectionMetadataResolver:
    """Resolve explicit batch-correction labels from sample metadata."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        batch_column: str,
        condition_column: str,
    ) -> ResolvedBatchCorrectionMetadata:
        if sample_metadata is None:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.batch_correction "
                "requires sample_metadata input data"
            )

        sample_order = _normalize_label_index(
            phospho.columns,
            field_name="dataset build request phospho.columns",
        )
        metadata_index = _normalize_label_index(
            sample_metadata.index,
            field_name="dataset build request sample_metadata.index",
        )
        _require_unique_labels(
            sample_order,
            field_name="dataset build request phospho.columns",
        )
        _require_unique_labels(
            metadata_index,
            field_name="dataset build request sample_metadata.index",
        )

        batch_column = str(batch_column).strip()
        condition_column = str(condition_column).strip()
        _require_column(
            sample_metadata,
            column=batch_column,
            field_name=(
                "dataset build request preprocessing_config.batch_correction."
                "batch_column"
            ),
        )
        _require_column(
            sample_metadata,
            column=condition_column,
            field_name=(
                "dataset build request preprocessing_config.batch_correction."
                "condition_column"
            ),
        )
        _require_no_missing_or_extra_metadata_rows(
            metadata_index=metadata_index,
            sample_order=sample_order,
        )

        aligned = sample_metadata.copy(deep=False)
        aligned.index = metadata_index
        aligned = aligned.reindex(sample_order)

        return ResolvedBatchCorrectionMetadata(
            batch_by_sample=_resolve_labels_by_sample(
                aligned,
                column=batch_column,
                sample_order=sample_order,
                label_kind="batch",
            ),
            condition_by_sample=_resolve_labels_by_sample(
                aligned,
                column=condition_column,
                sample_order=sample_order,
                label_kind="condition",
            ),
            sample_order=tuple(sample_order.tolist()),
        )


def _normalize_label_index(index: pd.Index, *, field_name: str) -> pd.Index:
    normalized: list[str] = []
    missing_positions: list[int] = []
    blank_positions: list[int] = []
    for position, value in enumerate(index.tolist()):
        if _is_missing_value(value):
            missing_positions.append(position)
            continue
        label = str(value).strip()
        if label == "":
            blank_positions.append(position)
            continue
        normalized.append(label)

    if missing_positions:
        raise PhosPyInputError(
            f"{field_name} must not contain missing labels for batch-correction "
            f"metadata resolution; found positions: "
            f"{_format_positions(missing_positions)}"
        )
    if blank_positions:
        raise PhosPyInputError(
            f"{field_name} must not contain blank labels for batch-correction "
            f"metadata resolution; found positions: {_format_positions(blank_positions)}"
        )
    return pd.Index(normalized, name=index.name)


def _require_unique_labels(index: pd.Index, *, field_name: str) -> None:
    if index.is_unique:
        return
    duplicates = list(dict.fromkeys(index[index.duplicated()].astype(str).tolist()))
    preview = ", ".join(repr(value) for value in duplicates[:5])
    suffix = "" if len(duplicates) <= 5 else " ..."
    raise PhosPyInputError(
        f"{field_name} contains duplicate sample labels for batch-correction "
        f"metadata resolution: {preview}{suffix}"
    )


def _require_column(
    sample_metadata: pd.DataFrame,
    *,
    column: str,
    field_name: str,
) -> None:
    matches = [
        candidate for candidate in sample_metadata.columns if candidate == column
    ]
    if len(matches) == 1:
        return
    if len(matches) > 1:
        raise PhosPyInputError(
            f"{field_name} resolves to duplicate sample_metadata column {column!r}"
        )
    raise PhosPyInputError(
        f"{field_name} references missing sample_metadata column {column!r}"
    )


def _require_no_missing_or_extra_metadata_rows(
    *,
    metadata_index: pd.Index,
    sample_order: pd.Index,
) -> None:
    metadata_samples = set(metadata_index.tolist())
    matrix_samples = set(sample_order.tolist())
    missing_samples = [
        sample for sample in sample_order.tolist() if sample not in metadata_samples
    ]
    extra_samples = [
        sample for sample in metadata_index.tolist() if sample not in matrix_samples
    ]
    if not missing_samples and not extra_samples:
        return

    details: list[str] = []
    if missing_samples:
        details.append(
            "missing rows for batch-correction samples: "
            f"{_format_label_preview(missing_samples)}"
        )
    if extra_samples:
        details.append(
            "rows not present in phospho columns for batch-correction metadata "
            f"resolution: {_format_label_preview(extra_samples)}"
        )
    raise PhosPyInputError(
        "dataset build request sample_metadata has incompatible sample rows for "
        f"batch-correction metadata resolution; {'; '.join(details)}"
    )


def _resolve_labels_by_sample(
    frame: pd.DataFrame,
    *,
    column: str,
    sample_order: pd.Index,
    label_kind: str,
) -> dict[str, str]:
    labels: dict[str, str] = {}
    missing_samples: list[str] = []
    blank_samples: list[str] = []
    for sample, value in zip(
        sample_order.tolist(),
        frame.loc[:, column].tolist(),
        strict=True,
    ):
        if _is_missing_value(value):
            missing_samples.append(str(sample))
            continue
        label = str(value).strip()
        if label == "":
            blank_samples.append(str(sample))
            continue
        labels[str(sample)] = label

    if missing_samples:
        preview = ", ".join(repr(value) for value in missing_samples[:5])
        suffix = "" if len(missing_samples) <= 5 else " ..."
        raise PhosPyInputError(
            f"dataset build request sample_metadata column {column!r} contains "
            f"missing {label_kind} labels for samples: {preview}{suffix}"
        )
    if blank_samples:
        preview = ", ".join(repr(value) for value in blank_samples[:5])
        suffix = "" if len(blank_samples) <= 5 else " ..."
        raise PhosPyInputError(
            f"dataset build request sample_metadata column {column!r} contains "
            f"blank {label_kind} labels for samples: {preview}{suffix}"
        )
    return labels


def _is_missing_value(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _format_positions(positions: list[int]) -> str:
    preview = ", ".join(str(position) for position in positions[:8])
    suffix = "" if len(positions) <= 8 else ", ..."
    return f"[{preview}{suffix}]"


def _format_label_preview(labels: list[str]) -> str:
    preview = ", ".join(repr(value) for value in labels[:5])
    suffix = "" if len(labels) <= 5 else " ..."
    return f"{preview}{suffix}"


def levels_in_sample_order(
    labels_by_sample: Mapping[str, str],
    *,
    sample_order: tuple[str, ...],
) -> tuple[str, ...]:
    """Return distinct labels in resolved sample order."""

    levels: list[str] = []
    seen: set[str] = set()
    for sample in sample_order:
        level = labels_by_sample[sample]
        if level in seen:
            continue
        seen.add(level)
        levels.append(level)
    return tuple(levels)


__all__ = [
    "BatchCorrectionMetadataResolver",
    "ResolvedBatchCorrectionMetadata",
    "levels_in_sample_order",
]
