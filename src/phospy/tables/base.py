"""Lightweight internal schema wrappers for scientific tables."""

from __future__ import annotations

from dataclasses import InitVar, dataclass
from typing import ClassVar

import pandas as pd

from phospy._frame_ownership import export_dataframe, own_dataframe
from phospy.errors.validation import PhosPyValidationError
from phospy.validation.common.dataframes import (
    require_no_duplicate_labels,
    require_string_index,
)

ValidationErrorType = type[PhosPyValidationError]


@dataclass(frozen=True, slots=True)
class TableSchema:
    """Base wrapper for one owned, validated DataFrame contract."""

    frame: pd.DataFrame
    _assume_owned: InitVar[bool] = False

    _field_name: ClassVar[str] = "table.frame"
    _error_type: ClassVar[ValidationErrorType] = PhosPyValidationError

    def __post_init__(self, _assume_owned: bool) -> None:
        frame = self._own_frame(_assume_owned)
        validated = self._validate_frame(frame)
        object.__setattr__(self, "frame", validated)

    def _own_frame(self, assume_owned: bool) -> pd.DataFrame:
        return own_dataframe(
            self.frame,
            field_name=self._field_name,
            error_type=self._error_type,
            assume_owned=assume_owned,
        )

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        return frame

    def to_pandas(self) -> pd.DataFrame:
        """Return a table snapshot; mutating it does not mutate this object."""

        return export_dataframe(self.frame)

    @classmethod
    def _from_owned(cls, *, frame: pd.DataFrame) -> TableSchema:
        return cls(frame=frame, _assume_owned=True)


def require_canonical_label_index(
    index: pd.Index,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> pd.Index:
    """Require canonical non-empty stripped string labels for one index."""

    require_string_index(
        index,
        field_name=field_name,
        error_type=error_type,
    )
    require_no_duplicate_labels(
        index,
        field_name=field_name,
        error_type=error_type,
    )
    values = index.tolist()
    stripped_values = [value.strip() for value in values]
    if any(value == "" for value in stripped_values):
        raise error_type(f"{field_name} must contain non-empty string labels")
    collisions: dict[str, set[str]] = {}
    for raw_value, stripped_value in zip(values, stripped_values, strict=False):
        collisions.setdefault(stripped_value, set()).add(raw_value)
    colliding = [
        value for value, raw_values in collisions.items() if len(raw_values) > 1
    ]
    if colliding:
        preview = ", ".join(repr(value) for value in colliding[:5])
        suffix = "" if len(colliding) <= 5 else " ..."
        raise error_type(
            f"{field_name} contains colliding labels when stripped: {preview}{suffix}"
        )
    if any(
        raw_value != stripped_value
        for raw_value, stripped_value in zip(values, stripped_values, strict=False)
    ):
        raise error_type(f"{field_name} must contain canonical non-empty string labels")
    return index
