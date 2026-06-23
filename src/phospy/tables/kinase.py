"""Kinase/prediction scientific table wrappers."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import pandas as pd

from phospy.errors.validation import PhosPyValidationError
from phospy.frames.ownership import own_dataframe
from phospy.science.sites.validation import require_site_key_index
from phospy.tables.base import TableSchema, require_canonical_label_index
from phospy.validation.common.dataframes import (
    require_dataframe,
    require_finite_numeric_dataframe,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
)
from phospy.validation.common.numeric_frames import require_numeric_unit_interval

KINASE_SUBSTRATE_CONTRIBUTION_STATUS_INCLUDED = "included"
KINASE_SUBSTRATE_CONTRIBUTION_STATUS_EXCLUDED = "excluded"
KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_BELOW_MIN_SUBSTRATES = (
    "kinase_below_min_substrates"
)
KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_MISSING_SCORE_VALUE = "missing_score_value"
KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_NO_SCORE_COLUMN = "kinase_score_unavailable"
KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_NOT_IN_PROFILE_SUPPORT = (
    "substrate_not_in_profile_support"
)
KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_NOT_QUANTIFIED = "substrate_not_quantified"
KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS: tuple[str, ...] = (
    "kinase",
    "substrate_site",
    "substrate_identifier",
    "value_used_in_scoring",
    "score_component",
    "score_source",
    "reference_source_name",
    "reference_source_version",
    "reference_bundle_id",
    "reference_identifier_namespace",
    "status",
    "exclusion_reason",
    "ambiguous",
)
_KINASE_SUBSTRATE_CONTRIBUTION_STATUSES = frozenset(
    {
        KINASE_SUBSTRATE_CONTRIBUTION_STATUS_INCLUDED,
        KINASE_SUBSTRATE_CONTRIBUTION_STATUS_EXCLUDED,
    }
)
_CONTRIBUTION_REQUIRED_TEXT_COLUMNS = (
    "kinase",
    "substrate_site",
    "score_component",
)
_CONTRIBUTION_OPTIONAL_TEXT_COLUMNS = (
    "substrate_identifier",
    "score_source",
    "reference_source_name",
    "reference_source_version",
    "reference_bundle_id",
    "reference_identifier_namespace",
    "exclusion_reason",
)


@dataclass(frozen=True, slots=True)
class KinaseScoreMatrix(TableSchema):
    """Schema wrapper for kinase scoring matrices."""

    field_name: str = field(
        default="scoring_result.profile_scores",
        repr=False,
        compare=False,
    )
    allow_missing: bool = field(default=True, repr=False, compare=False)
    enforce_unit_interval: bool = field(default=False, repr=False, compare=False)
    require_site_index: bool = field(default=True, repr=False, compare=False)

    _field_name = "scoring_result.profile_scores"
    _error_type = PhosPyValidationError

    def _own_frame(self, assume_owned: bool) -> pd.DataFrame:
        return own_dataframe(
            self.frame,
            field_name=self.field_name,
            error_type=self._error_type,
            assume_owned=assume_owned,
        )

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        _validate_kinase_score_like_matrix(
            frame=frame,
            field_name=self.field_name,
            error_type=self._error_type,
            allow_missing=self.allow_missing,
            enforce_unit_interval=self.enforce_unit_interval,
            require_site_index=self.require_site_index,
        )
        return frame


@dataclass(frozen=True, slots=True)
class KinasePredictionMatrix(TableSchema):
    """Schema wrapper for ``prediction_result.pred_mat``."""

    field_name: str = field(
        default="prediction_result.pred_mat",
        repr=False,
        compare=False,
    )
    allow_missing: bool = field(default=True, repr=False, compare=False)
    enforce_unit_interval: bool = field(default=True, repr=False, compare=False)

    _field_name = "prediction_result.pred_mat"
    _error_type = PhosPyValidationError

    def _own_frame(self, assume_owned: bool) -> pd.DataFrame:
        return own_dataframe(
            self.frame,
            field_name=self.field_name,
            error_type=self._error_type,
            assume_owned=assume_owned,
        )

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        _validate_kinase_score_like_matrix(
            frame=frame,
            field_name=self.field_name,
            error_type=self._error_type,
            allow_missing=self.allow_missing,
            enforce_unit_interval=self.enforce_unit_interval,
            require_site_index=True,
        )
        return frame


@dataclass(frozen=True, slots=True)
class KinaseSubstrateContributionTable(TableSchema):
    """Schema wrapper for optional substrate contribution rows."""

    field_name: str = field(
        default="kinase_result.substrate_contributions",
        repr=False,
        compare=False,
    )

    _field_name = "kinase_result.substrate_contributions"
    _error_type = PhosPyValidationError

    def _own_frame(self, assume_owned: bool) -> pd.DataFrame:
        return own_dataframe(
            self.frame,
            field_name=self.field_name,
            error_type=self._error_type,
            assume_owned=assume_owned,
        )

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        require_dataframe(
            frame,
            field_name=self.field_name,
            allow_empty=True,
            error_type=self._error_type,
        )
        require_unique_columns(
            frame,
            field_name=self.field_name,
            error_type=self._error_type,
        )
        observed_columns = tuple(str(column) for column in frame.columns.tolist())
        if observed_columns != KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS:
            raise self._error_type(
                f"{self.field_name} columns must exactly match "
                f"{list(KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS)}"
            )
        if frame.empty:
            return frame
        _require_contribution_required_text(frame, field_name=self.field_name)
        _require_contribution_optional_text(frame, field_name=self.field_name)
        _require_contribution_numeric_or_missing(frame, field_name=self.field_name)
        _require_contribution_status(frame, field_name=self.field_name)
        _require_contribution_ambiguous_flags(frame, field_name=self.field_name)
        return frame


def _validate_kinase_score_like_matrix(
    *,
    frame: pd.DataFrame,
    field_name: str,
    error_type: type[PhosPyValidationError],
    allow_missing: bool,
    enforce_unit_interval: bool,
    require_site_index: bool,
) -> None:
    require_dataframe(
        frame,
        field_name=field_name,
        allow_empty=False,
        error_type=error_type,
    )
    require_unique_columns(
        frame,
        field_name=field_name,
        error_type=error_type,
    )
    require_numeric_dataframe(
        frame,
        field_name=field_name,
        error_type=error_type,
    )
    require_finite_numeric_dataframe(
        frame,
        field_name=field_name,
        error_type=error_type,
        allow_missing=allow_missing,
    )
    require_unique_index(
        frame,
        field_name=field_name,
        error_type=error_type,
    )
    if require_site_index:
        _require_site_index_identity(
            frame.index,
            field_name=f"{field_name}.index",
            error_type=error_type,
        )
    else:
        require_canonical_label_index(
            frame.index,
            field_name=f"{field_name}.index",
            error_type=error_type,
        )
    require_canonical_label_index(
        frame.columns,
        field_name=f"{field_name}.columns",
        error_type=error_type,
    )
    if enforce_unit_interval:
        require_numeric_unit_interval(
            frame,
            field_name=field_name,
            error_type=error_type,
        )


def _require_contribution_required_text(
    frame: pd.DataFrame,
    *,
    field_name: str,
) -> None:
    for column_name in _CONTRIBUTION_REQUIRED_TEXT_COLUMNS:
        invalid = [
            index
            for index, value in frame.loc[:, column_name].items()
            if not isinstance(value, str) or value.strip() == ""
        ]
        if invalid:
            raise PhosPyValidationError(
                f"{field_name}.{column_name} must contain non-empty strings"
            )


def _require_contribution_optional_text(
    frame: pd.DataFrame,
    *,
    field_name: str,
) -> None:
    for column_name in _CONTRIBUTION_OPTIONAL_TEXT_COLUMNS:
        invalid = [
            index
            for index, value in frame.loc[:, column_name].items()
            if not _is_missing_value(value)
            and (not isinstance(value, str) or value.strip() == "")
        ]
        if invalid:
            raise PhosPyValidationError(
                f"{field_name}.{column_name} must contain non-empty strings or "
                "missing values"
            )


def _require_contribution_numeric_or_missing(
    frame: pd.DataFrame,
    *,
    field_name: str,
) -> None:
    invalid = []
    for index, value in frame.loc[:, "value_used_in_scoring"].items():
        if _is_missing_value(value):
            continue
        try:
            numeric_value = float(value)
        except (TypeError, ValueError):
            invalid.append(index)
            continue
        if not math.isfinite(numeric_value):
            invalid.append(index)
    if invalid:
        raise PhosPyValidationError(
            f"{field_name}.value_used_in_scoring must contain finite numeric "
            "values or missing values"
        )


def _require_contribution_status(
    frame: pd.DataFrame,
    *,
    field_name: str,
) -> None:
    invalid_statuses = [
        value
        for value in frame.loc[:, "status"].tolist()
        if not isinstance(value, str)
        or value not in _KINASE_SUBSTRATE_CONTRIBUTION_STATUSES
    ]
    if invalid_statuses:
        allowed = ", ".join(sorted(_KINASE_SUBSTRATE_CONTRIBUTION_STATUSES))
        raise PhosPyValidationError(f"{field_name}.status must be one of: {allowed}")
    for row in frame.loc[:, ["status", "exclusion_reason"]].itertuples(index=False):
        status = str(row.status)
        exclusion_reason = row.exclusion_reason
        if status == KINASE_SUBSTRATE_CONTRIBUTION_STATUS_INCLUDED:
            if not _is_missing_value(exclusion_reason):
                raise PhosPyValidationError(
                    f"{field_name}.exclusion_reason must be missing for included rows"
                )
            continue
        if _is_missing_value(exclusion_reason):
            raise PhosPyValidationError(
                f"{field_name}.exclusion_reason must be set for excluded rows"
            )


def _require_contribution_ambiguous_flags(
    frame: pd.DataFrame,
    *,
    field_name: str,
) -> None:
    invalid = [
        value
        for value in frame.loc[:, "ambiguous"].tolist()
        if not _is_bool_value(value)
    ]
    if invalid:
        raise PhosPyValidationError(f"{field_name}.ambiguous must contain bool values")


def _is_bool_value(value: object) -> bool:
    if isinstance(value, bool):
        return True
    return type(value).__name__ == "bool_"


def _is_missing_value(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _require_site_index_identity(
    index: pd.Index,
    *,
    field_name: str,
    error_type: type[PhosPyValidationError],
) -> None:
    require_site_key_index(
        index,
        field_name=field_name,
        error_type=error_type,
    )
