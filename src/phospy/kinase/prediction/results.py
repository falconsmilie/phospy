"""Kinase prediction-stage result model."""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from phospy.errors.validation import PhosPyValidationError
from phospy.frames.comparison import dataframe_equals, optional_dataframe_equals
from phospy.frames.ownership import (
    _borrow_dataframe,
    _borrow_optional_dataframe,
    export_dataframe,
    export_optional_dataframe,
    own_optional_dataframe,
)


@dataclass(frozen=True, slots=True, init=False, eq=False)
class KinasePredictionResult:
    """Prediction-stage outputs.

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for explicit prediction-content comparison.
    """

    __hash__ = object.__hash__

    _pred_mat: pd.DataFrame = field(init=False, repr=False)
    _substrate_list: pd.DataFrame | None = field(init=False, repr=False)

    def __init__(
        self,
        pred_mat: pd.DataFrame,
        substrate_list: pd.DataFrame | None = None,
    ) -> None:
        from phospy.kinase.prediction.tables import KinasePredictionMatrix

        pred_mat = KinasePredictionMatrix(
            frame=pred_mat,
            field_name="prediction_result.pred_mat",
            _assume_owned=False,
        ).frame
        substrate_list = own_optional_dataframe(
            substrate_list,
            field_name="prediction_result.substrate_list",
            error_type=PhosPyValidationError,
            assume_owned=False,
        )
        object.__setattr__(self, "_pred_mat", pred_mat)
        object.__setattr__(self, "_substrate_list", substrate_list)

    @property
    def pred_mat(self) -> pd.DataFrame:
        return export_dataframe(self._pred_mat)

    @property
    def substrate_list(self) -> pd.DataFrame | None:
        return _export_public_substrate_list(self._substrate_list)

    def _borrow_pred_mat_frame(self) -> pd.DataFrame:
        """Package-private prediction-matrix snapshot for internal views."""

        return _borrow_dataframe(self._pred_mat)

    def _borrow_substrate_list_frame(self) -> pd.DataFrame | None:
        """Package-private substrate-list snapshot for internal views."""

        return _borrow_optional_dataframe(self._substrate_list)

    @classmethod
    def _from_owned(
        cls,
        *,
        pred_mat: pd.DataFrame,
        substrate_list: pd.DataFrame | None = None,
    ) -> KinasePredictionResult:
        _require_frame_type(pred_mat, field_name="prediction_result.pred_mat")
        _require_optional_frame_type(
            substrate_list,
            field_name="prediction_result.substrate_list",
        )
        result = object.__new__(cls)
        object.__setattr__(result, "_pred_mat", pred_mat)
        object.__setattr__(result, "_substrate_list", substrate_list)
        return result

    @classmethod
    def from_trusted_owned(
        cls,
        *,
        pred_mat: pd.DataFrame,
        substrate_list: pd.DataFrame | None = None,
    ) -> KinasePredictionResult:
        """Construct from already-owned tables at trusted workflow boundaries."""

        return cls._from_owned(
            pred_mat=pred_mat,
            substrate_list=substrate_list,
        )

    def to_dataframe(self) -> pd.DataFrame:
        """Return a `pred_mat` snapshot isolated from this result."""

        return export_dataframe(self._pred_mat)

    def substrate_list_dataframe(self) -> pd.DataFrame | None:
        """Return an optional substrate-list snapshot isolated from this result."""

        return _export_public_substrate_list(self._substrate_list)

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another prediction result has the same content."""

        if not isinstance(other, KinasePredictionResult):
            return False
        return dataframe_equals(self._pred_mat, other._pred_mat) and (
            optional_dataframe_equals(self._substrate_list, other._substrate_list)
        )


def _export_public_substrate_list(table: pd.DataFrame | None) -> pd.DataFrame | None:
    exported = export_optional_dataframe(table)
    if exported is None:
        return None
    legacy_columns = ["kinase", "substrate_site", "score", "rank"]
    if not all(column in exported.columns for column in legacy_columns):
        return exported
    if {"site_key", "display_id"}.issubset(set(exported.columns)):
        return exported
    return exported.loc[:, legacy_columns]


def _require_frame_type(value: object, *, field_name: str) -> None:
    if not isinstance(value, pd.DataFrame):
        raise PhosPyValidationError(f"{field_name} must be a pandas DataFrame")


def _require_optional_frame_type(value: object | None, *, field_name: str) -> None:
    if value is not None and not isinstance(value, pd.DataFrame):
        raise PhosPyValidationError(
            f"{field_name} must be a pandas DataFrame when provided"
        )


def _substrate_list_uses_encoded_site_keys(table: pd.DataFrame) -> bool:
    if "site_key" not in table.columns:
        return False
    values = table.loc[:, "site_key"].astype(str).tolist()
    if not values:
        return False
    return all(value.startswith("phospy:v1|") for value in values)


__all__ = ["KinasePredictionResult"]
