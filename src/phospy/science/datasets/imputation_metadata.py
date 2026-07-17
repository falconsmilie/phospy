"""Imputation observation-mask metadata for analysis-ready datasets."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt
import pandas as pd

from phospy.errors.validation import DatasetValidationError
from phospy.frames.ownership import borrow_dataframe, export_dataframe, own_dataframe
from phospy.frames.validation import (
    require_dataframe,
    require_exact_index_match,
)

IMPUTATION_FEATURE_METADATA_COLUMNS = (
    "imputed_cell_count",
    "observed_cell_count",
    "imputed_fraction",
)
IMPUTATION_OBSERVATION_SUMMARY_COLUMNS = (
    "feature_id",
    "observed_cell_count",
    "imputed_cell_count",
    "total_analysed_cell_count",
    "imputed_fraction",
)


@dataclass(frozen=True, slots=True, init=False)
class ImputationObservationMetadata:
    """Dataset-owned originally-observed vs imputed-cell metadata.

    The internal mask is aligned to `dataset.phospho`: rows are phosphosite
    features, columns are samples, and True means the value was originally
    observed rather than imputed. Public accessors return defensive snapshots.
    """

    _observed_mask: pd.DataFrame = field(init=False, repr=False)
    _feature_summary: pd.DataFrame = field(init=False, repr=False)

    def __init__(
        self,
        *,
        observed_mask: pd.DataFrame,
        phospho_index: pd.Index,
        sample_index: pd.Index,
        _assume_owned: bool = False,
    ) -> None:
        mask = own_dataframe(
            observed_mask,
            field_name="dataset.imputation_observation_mask",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        require_dataframe(
            mask,
            field_name="dataset.imputation_observation_mask",
            allow_empty=False,
            error_type=DatasetValidationError,
        )
        require_exact_index_match(
            left=mask.index,
            right=phospho_index,
            left_name="dataset.imputation_observation_mask.index",
            right_name="dataset.phospho.index",
            error_type=DatasetValidationError,
        )
        require_exact_index_match(
            left=mask.columns,
            right=sample_index,
            left_name="dataset.imputation_observation_mask.columns",
            right_name="dataset.phospho.columns",
            error_type=DatasetValidationError,
        )
        _require_boolean_observation_mask(mask)

        observed = pd.DataFrame(
            mask.to_numpy(dtype=bool),
            index=mask.index,
            columns=mask.columns,
        )
        observed.index.name = mask.index.name
        observed.columns.name = mask.columns.name
        feature_summary = _build_imputation_feature_summary(observed)
        object.__setattr__(self, "_observed_mask", observed)
        object.__setattr__(self, "_feature_summary", feature_summary)

    @property
    def feature_summary(self) -> pd.DataFrame:
        return export_dataframe(self._feature_summary)

    @property
    def observed_mask(self) -> pd.DataFrame:
        return export_dataframe(self._observed_mask)

    def feature_summary_dataframe(self) -> pd.DataFrame:
        """Return per-feature imputation counts isolated from this metadata."""

        return export_dataframe(self._feature_summary)

    def observed_mask_dataframe(self) -> pd.DataFrame:
        """Return a defensive observed-cell mask snapshot."""

        return export_dataframe(self._observed_mask)

    def feature_observation_summary_dataframe(
        self,
        *,
        feature_ids: Sequence[object],
        sample_ids: Sequence[object],
    ) -> pd.DataFrame:
        """Return feature-level observation counts for a requested subset."""

        requested_feature_ids = _requested_label_list(
            feature_ids,
            field_name="dataset.imputation_observation_summary.feature_ids",
        )
        requested_sample_ids = _requested_label_list(
            sample_ids,
            field_name="dataset.imputation_observation_summary.sample_ids",
        )
        _require_requested_labels_present(
            requested_feature_ids,
            available_labels=self._observed_mask.index,
            field_name="dataset.imputation_observation_summary.feature_ids",
            available_field_name="dataset.imputation_observation_mask.index",
        )
        _require_requested_labels_present(
            requested_sample_ids,
            available_labels=self._observed_mask.columns,
            field_name="dataset.imputation_observation_summary.sample_ids",
            available_field_name="dataset.imputation_observation_mask.columns",
        )
        observed_values = self._observed_mask.to_numpy(dtype=bool)
        feature_positions = _label_positions(
            requested_feature_ids,
            available_labels=self._observed_mask.index,
        )
        sample_positions = _label_positions(
            requested_sample_ids,
            available_labels=self._observed_mask.columns,
        )
        observed_subset = pd.DataFrame(
            observed_values[np.ix_(feature_positions, sample_positions)],
            index=pd.Index(
                requested_feature_ids,
                name=self._observed_mask.index.name,
            ),
            columns=pd.Index(
                requested_sample_ids,
                name=self._observed_mask.columns.name,
            ),
        )
        summary = _build_imputation_observation_summary(observed_subset)
        return export_dataframe(summary)

    def aggregated_observed_mask_dataframe(
        self,
        *,
        sample_groups: Sequence[tuple[object, Sequence[object]]],
    ) -> pd.DataFrame:
        """Return an observed-cell mask collapsed to requested sample groups."""

        if not sample_groups:
            raise DatasetValidationError(
                "dataset.imputation_observation_mask sample_groups must contain "
                "at least one sample group"
            )
        observed_values = self._observed_mask.to_numpy(dtype=bool)
        aggregated_column_values: list[np.ndarray] = []
        output_labels: list[object] = []
        for output_label, input_sample_ids in sample_groups:
            requested_sample_ids = _requested_label_list(
                input_sample_ids,
                field_name=(
                    "dataset.imputation_observation_mask sample_groups."
                    f"{output_label!r}.sample_ids"
                ),
            )
            _require_requested_labels_present(
                requested_sample_ids,
                available_labels=self._observed_mask.columns,
                field_name=(
                    "dataset.imputation_observation_mask sample_groups."
                    f"{output_label!r}.sample_ids"
                ),
                available_field_name="dataset.imputation_observation_mask.columns",
            )
            sample_positions = _label_positions(
                requested_sample_ids,
                available_labels=self._observed_mask.columns,
            )
            collapsed = np.all(observed_values[:, sample_positions], axis=1)
            aggregated_column_values.append(collapsed.astype(bool))
            output_labels.append(output_label)
        aggregated_values = np.column_stack(aggregated_column_values).astype(bool)
        aggregated = pd.DataFrame(
            aggregated_values,
            index=pd.Index(
                self._observed_mask.index.tolist(),
                name=self._observed_mask.index.name,
            ),
            columns=pd.Index(output_labels, name=self._observed_mask.columns.name),
        )
        return export_dataframe(aggregated)

    def _borrow_observed_mask_frame(self) -> pd.DataFrame:
        """Package-private defensive mask snapshot for internal read paths."""

        return borrow_dataframe(self._observed_mask)


def _build_imputation_observation_metadata_or_none(
    *,
    imputation_observation_mask: pd.DataFrame | None,
    phospho_index: pd.Index,
    sample_index: pd.Index,
) -> ImputationObservationMetadata | None:
    if imputation_observation_mask is None:
        return None
    return ImputationObservationMetadata(
        observed_mask=imputation_observation_mask,
        phospho_index=phospho_index,
        sample_index=sample_index,
        _assume_owned=True,
    )


def _require_boolean_observation_mask(mask: pd.DataFrame) -> None:
    values = mask.to_numpy(dtype="object", copy=False)
    missing_values: npt.NDArray[np.bool_] = np.asarray(
        pd.isna(values),
        dtype=bool,
    )
    if bool(missing_values.any()):
        raise DatasetValidationError(
            "dataset.imputation_observation_mask must not contain missing values"
        )
    if all(pd.api.types.is_bool_dtype(dtype) for dtype in mask.dtypes):
        return

    boolean_result = np.frompyfunc(
        lambda value: isinstance(value, (bool, np.bool_)),
        1,
        1,
    )(values)
    boolean_cells: npt.NDArray[np.bool_] = np.asarray(boolean_result, dtype=bool)
    invalid_locations = np.argwhere(~boolean_cells)
    if invalid_locations.size == 0:
        return
    row_index, column_index = invalid_locations[0]
    raise DatasetValidationError(
        "dataset.imputation_observation_mask must contain only boolean "
        "values; "
        f"invalid_cell=({mask.index[int(row_index)]!r}, "
        f"{mask.columns[int(column_index)]!r})"
    )


def _build_imputation_feature_summary(observed_mask: pd.DataFrame) -> pd.DataFrame:
    sample_count = int(observed_mask.shape[1])
    observed_values = observed_mask.to_numpy(dtype=bool)
    observed_counts = observed_values.sum(axis=1).astype(np.int64)
    imputed_counts = (sample_count - observed_counts).astype(np.int64)
    imputed_fraction = imputed_counts.astype(float) / float(sample_count)
    summary = pd.DataFrame(
        {
            "imputed_cell_count": imputed_counts,
            "observed_cell_count": observed_counts,
            "imputed_fraction": imputed_fraction,
        },
        index=observed_mask.index,
    )
    summary.index.name = observed_mask.index.name
    return summary


def _build_imputation_observation_summary(observed_mask: pd.DataFrame) -> pd.DataFrame:
    sample_count = int(observed_mask.shape[1])
    observed_values = observed_mask.to_numpy(dtype=bool)
    observed_counts = observed_values.sum(axis=1).astype(np.int64)
    imputed_counts = (sample_count - observed_counts).astype(np.int64)
    imputed_fraction = imputed_counts.astype(float) / float(sample_count)
    feature_ids = observed_mask.index.tolist()
    summary = pd.DataFrame(
        {
            "feature_id": feature_ids,
            "observed_cell_count": observed_counts,
            "imputed_cell_count": imputed_counts,
            "total_analysed_cell_count": np.full(
                shape=observed_counts.shape,
                fill_value=sample_count,
                dtype=np.int64,
            ),
            "imputed_fraction": imputed_fraction,
        },
        index=pd.Index(feature_ids, name=observed_mask.index.name),
        columns=list(IMPUTATION_OBSERVATION_SUMMARY_COLUMNS),
    )
    return summary


def _requested_label_list(
    labels: Sequence[object],
    *,
    field_name: str,
) -> list[object]:
    requested: list[object]
    if isinstance(labels, str | bytes):
        requested = [labels]
    else:
        try:
            requested = list(labels)
        except TypeError as exc:
            raise DatasetValidationError(
                f"{field_name} must be a sequence of labels"
            ) from exc
    if not requested:
        raise DatasetValidationError(f"{field_name} must contain at least one label")
    return requested


def _label_positions(
    requested_labels: Sequence[object],
    *,
    available_labels: pd.Index,
) -> list[int]:
    positions_by_label = {
        label: int(position) for position, label in enumerate(available_labels.tolist())
    }
    return [positions_by_label[label] for label in requested_labels]


def _require_requested_labels_present(
    requested_labels: Sequence[object],
    *,
    available_labels: pd.Index,
    field_name: str,
    available_field_name: str,
) -> None:
    available = set(available_labels.tolist())
    missing: list[object] = []
    for label in requested_labels:
        if label in available:
            continue
        missing.append(label)
    if missing:
        raise DatasetValidationError(
            f"{field_name} contains labels absent from {available_field_name}: "
            + _format_label_preview(missing)
        )


def _format_label_preview(labels: Sequence[object]) -> str:
    preview = [repr(label) for label in labels[:5]]
    if len(labels) > 5:
        preview.append("...")
    return ", ".join(preview)


build_imputation_observation_metadata_or_none = (
    _build_imputation_observation_metadata_or_none
)
require_boolean_observation_mask = _require_boolean_observation_mask


__all__ = [
    "IMPUTATION_FEATURE_METADATA_COLUMNS",
    "IMPUTATION_OBSERVATION_SUMMARY_COLUMNS",
    "ImputationObservationMetadata",
    "build_imputation_observation_metadata_or_none",
    "require_boolean_observation_mask",
]
