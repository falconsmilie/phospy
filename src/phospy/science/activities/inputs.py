"""Validated activity-stage input models."""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import pandas as pd

from phospy.errors.validation import PhosPyValidationError
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.frames.comparison import dataframe_equals
from phospy.science.activities.membership import ActivityMembershipSelection
from phospy.science.activities.semantics import (
    ActivityInputMatrix,
    ActivityInputSemantics,
    ActivityProfileMetadata,
)
from phospy.science.tables.activity import ActivityMatrix
from phospy.science.tables.kinase import KinasePredictionMatrix


@dataclass(frozen=True, slots=True)
class PredMatOverlapSummary:
    """Resolved overlap diagnostics between prediction and phospho matrices."""

    overlap_count: int
    pred_mat_rows: int
    phospho_rows: int


@dataclass(frozen=True, slots=True, eq=False)
class KinaseActivityInputs:
    """Trusted activity-like score inputs resolved by workflow validation.

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for explicit activity-input content comparison.
    """

    __hash__ = object.__hash__

    pred_mat: pd.DataFrame
    phospho_matrix: pd.DataFrame
    threshold: float
    min_substrates: int
    top_n_substrates: int
    overlap_summary: PredMatOverlapSummary
    activity_input: ActivityInputMatrix | None = None
    membership_selection: ActivityMembershipSelection | None = None

    def __post_init__(self) -> None:
        if self.activity_input is None:
            warnings.warn(
                (
                    "KinaseActivityInputs constructed without typed activity_input "
                    "semantics is deprecated; treating phospho_matrix columns as "
                    "sample/profile labels with sample-level abundance semantics."
                ),
                DeprecationWarning,
                stacklevel=2,
            )
            activity_input = ActivityInputMatrix.sample_level_abundance(
                self.phospho_matrix,
                field_name="dataset.phospho",
                _assume_owned=True,
            )
        elif isinstance(self.activity_input, ActivityInputMatrix):
            activity_input = self.activity_input
        else:
            raise WorkflowBoundaryError(
                "activity input activity_input must be ActivityInputMatrix or None"
            )
        try:
            pred_mat = KinasePredictionMatrix(
                frame=self.pred_mat,
                field_name="prediction_result.pred_mat",
                _assume_owned=True,
            ).frame
            phospho_matrix = ActivityMatrix(
                frame=activity_input.frame,
                field_name="dataset.phospho",
                _assume_owned=True,
            ).frame
        except PhosPyValidationError as exc:
            raise WorkflowBoundaryError(
                seam="kinase.activity.input_schema",
                next_action=(
                    "ensure prediction_result.pred_mat and the typed activity input "
                    "matrix satisfy activity-stage input table schema requirements"
                ),
                details={"schema_error": str(exc)},
                message_prefix="kinase workflow boundary validation failed",
            ) from exc
        if not isinstance(self.overlap_summary, PredMatOverlapSummary):
            raise WorkflowBoundaryError(
                "activity input overlap_summary must be PredMatOverlapSummary"
            )
        if self.membership_selection is None:
            membership_selection = ActivityMembershipSelection.missing(
                selected_kinase_universe=pred_mat.columns.astype(str).tolist(),
                selected_substrate_universe=pred_mat.index.astype(str).tolist(),
            )
        elif isinstance(self.membership_selection, ActivityMembershipSelection):
            membership_selection = self.membership_selection
        else:
            raise WorkflowBoundaryError(
                "activity input membership_selection must be "
                "ActivityMembershipSelection or None"
            )
        object.__setattr__(self, "pred_mat", pred_mat)
        object.__setattr__(self, "phospho_matrix", phospho_matrix)
        object.__setattr__(
            self,
            "activity_input",
            ActivityInputMatrix(
                frame=phospho_matrix,
                semantics=activity_input.semantics,
                profile_metadata=activity_input.profile_metadata,
                field_name="dataset.phospho",
                _assume_owned=True,
            ),
        )
        object.__setattr__(self, "membership_selection", membership_selection)

    @property
    def input_semantics(self) -> ActivityInputSemantics:
        if self.activity_input is None:
            raise RuntimeError("KinaseActivityInputs.activity_input was not resolved")
        return self.activity_input.semantics

    @property
    def profile_metadata(self) -> ActivityProfileMetadata:
        if self.activity_input is None:
            raise RuntimeError("KinaseActivityInputs.activity_input was not resolved")
        return self.activity_input.profile_metadata

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another validated input has the same content."""

        if not isinstance(other, KinaseActivityInputs):
            return False
        if self.activity_input is None or other.activity_input is None:
            same_activity_input = self.activity_input is other.activity_input
        else:
            same_activity_input = self.activity_input.scientifically_equals(
                other.activity_input
            )
        return (
            dataframe_equals(self.pred_mat, other.pred_mat)
            and dataframe_equals(self.phospho_matrix, other.phospho_matrix)
            and self.threshold == other.threshold
            and self.min_substrates == other.min_substrates
            and self.top_n_substrates == other.top_n_substrates
            and self.overlap_summary == other.overlap_summary
            and same_activity_input
            and self.membership_selection == other.membership_selection
        )


__all__ = [
    "KinaseActivityInputs",
    "PredMatOverlapSummary",
]
