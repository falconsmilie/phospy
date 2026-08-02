"""Activity-like score stage result models."""

from __future__ import annotations

import warnings
from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from phospy.errors.validation import PhosPyValidationError
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.frames.ownership import (
    export_dataframe,
    export_optional_dataframe,
    export_series,
)
from phospy.provenance.scientific_policy_models import ScientificPolicyRecord
from phospy.science.activities.semantics import (
    ActivityInputMatrix,
    ActivityInputSemantics,
    ActivityProfileAxis,
    ActivityProfileMetadata,
    ActivityQuantitativeSemantics,
)
from phospy.science.activities.threshold_membership import (
    ActivityThresholdMembershipDiagnostics,
)
from phospy.science.tables.activity import (
    ActivityCountMatrix,
    ActivityCountSeries,
    ActivityMatrix,
    ActivityStatisticsTable,
    ActivityTargetTable,
)
from phospy.science.tables.kinase import KinasePredictionMatrix


@dataclass(frozen=True, slots=True)
class ActivityMethodMetadata:
    """Stable scientific identity metadata for an activity-like scoring method."""

    activity_method_id: str
    activity_method_family: str
    activity_method_label: str
    is_ksea: bool
    is_phosr_kinase_activity_equivalent: bool

    def to_payload(self) -> dict[str, object]:
        """Return a scalar metadata snapshot, not an export/report payload."""

        return {
            "activity_method_id": self.activity_method_id,
            "activity_method_family": self.activity_method_family,
            "activity_method_label": self.activity_method_label,
            "is_ksea": bool(self.is_ksea),
            "is_phosr_kinase_activity_equivalent": bool(
                self.is_phosr_kinase_activity_equivalent
            ),
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> ActivityMethodMetadata:
        method_id = str(payload.get("activity_method_id", "")).strip()
        method_family = str(payload.get("activity_method_family", "")).strip()
        method_label = str(payload.get("activity_method_label", "")).strip()
        is_ksea = payload.get("is_ksea")
        is_phosr_equivalent = payload.get("is_phosr_kinase_activity_equivalent")
        if not method_id:
            raise ValueError("activity_method_id must be a non-empty string")
        if not method_family:
            raise ValueError("activity_method_family must be a non-empty string")
        if not method_label:
            raise ValueError("activity_method_label must be a non-empty string")
        if not isinstance(is_ksea, bool):
            raise ValueError("is_ksea must be a bool")
        if not isinstance(is_phosr_equivalent, bool):
            raise ValueError("is_phosr_kinase_activity_equivalent must be a bool")
        return cls(
            activity_method_id=method_id,
            activity_method_family=method_family,
            activity_method_label=method_label,
            is_ksea=is_ksea,
            is_phosr_kinase_activity_equivalent=is_phosr_equivalent,
        )


SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD = ActivityMethodMetadata(
    activity_method_id="simplified_weighted_substrate_activity_v1",
    activity_method_family="heuristic_weighted_substrate_score",
    activity_method_label="simplified weighted substrate activity-like score",
    is_ksea=False,
    is_phosr_kinase_activity_equivalent=False,
)

KSEA_ZSCORE_ACTIVITY_METHOD = ActivityMethodMetadata(
    activity_method_id="ksea_zscore_v1",
    activity_method_family="substrate_set_enrichment",
    activity_method_label="KSEA-style z-score kinase activity score",
    is_ksea=True,
    is_phosr_kinase_activity_equivalent=False,
)

SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_METHOD = ActivityMethodMetadata(
    activity_method_id="ssgsea_substrate_enrichment_activity_v1",
    activity_method_family="substrate_set_enrichment",
    activity_method_label="ssGSEA substrate enrichment activity-like score",
    is_ksea=False,
    is_phosr_kinase_activity_equivalent=False,
)


@dataclass(frozen=True, slots=True)
class ActivityMethodSummary:
    """Method-level score computability counters."""

    kinases_evaluated: int
    kinase_condition_pairs_evaluated: int
    kinase_condition_pairs_computed: int
    kinase_condition_pairs_insufficient_substrates: int
    kinase_condition_pairs_invalid_background_variance: int
    kinase_condition_pairs_no_finite_background_values: int
    kinase_condition_pairs_no_finite_substrate_values: int

    def to_payload(self) -> dict[str, int]:
        """Return computability counters as a plain defensive snapshot."""

        return {
            "kinases_evaluated": int(self.kinases_evaluated),
            "kinase_profile_pairs_evaluated": int(
                self.kinase_condition_pairs_evaluated
            ),
            "kinase_profile_pairs_computed": int(self.kinase_condition_pairs_computed),
            "kinase_profile_pairs_insufficient_substrates": int(
                self.kinase_condition_pairs_insufficient_substrates
            ),
            "kinase_profile_pairs_invalid_background_variance": int(
                self.kinase_condition_pairs_invalid_background_variance
            ),
            "kinase_profile_pairs_no_finite_background_values": int(
                self.kinase_condition_pairs_no_finite_background_values
            ),
            "kinase_profile_pairs_no_finite_substrate_values": int(
                self.kinase_condition_pairs_no_finite_substrate_values
            ),
            "kinase_condition_pairs_evaluated": int(
                self.kinase_condition_pairs_evaluated
            ),
            "kinase_condition_pairs_computed": int(
                self.kinase_condition_pairs_computed
            ),
            "kinase_condition_pairs_insufficient_substrates": int(
                self.kinase_condition_pairs_insufficient_substrates
            ),
            "kinase_condition_pairs_invalid_background_variance": int(
                self.kinase_condition_pairs_invalid_background_variance
            ),
            "kinase_condition_pairs_no_finite_background_values": int(
                self.kinase_condition_pairs_no_finite_background_values
            ),
            "kinase_condition_pairs_no_finite_substrate_values": int(
                self.kinase_condition_pairs_no_finite_substrate_values
            ),
        }

    @property
    def kinase_profile_pairs_evaluated(self) -> int:
        return self.kinase_condition_pairs_evaluated

    @property
    def kinase_profile_pairs_computed(self) -> int:
        return self.kinase_condition_pairs_computed

    @property
    def kinase_profile_pairs_insufficient_substrates(self) -> int:
        return self.kinase_condition_pairs_insufficient_substrates

    @property
    def kinase_profile_pairs_invalid_background_variance(self) -> int:
        return self.kinase_condition_pairs_invalid_background_variance

    @property
    def kinase_profile_pairs_no_finite_background_values(self) -> int:
        return self.kinase_condition_pairs_no_finite_background_values

    @property
    def kinase_profile_pairs_no_finite_substrate_values(self) -> int:
        return self.kinase_condition_pairs_no_finite_substrate_values

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> ActivityMethodSummary:
        return cls(
            kinases_evaluated=_coerce_payload_int(
                payload=payload,
                field_name="kinases_evaluated",
            ),
            kinase_condition_pairs_evaluated=_coerce_payload_int_with_fallback(
                payload=payload,
                field_name="kinase_profile_pairs_evaluated",
                fallback_field_name="kinase_condition_pairs_evaluated",
            ),
            kinase_condition_pairs_computed=_coerce_payload_int_with_fallback(
                payload=payload,
                field_name="kinase_profile_pairs_computed",
                fallback_field_name="kinase_condition_pairs_computed",
            ),
            kinase_condition_pairs_insufficient_substrates=(
                _coerce_payload_int_with_fallback(
                    payload=payload,
                    field_name="kinase_profile_pairs_insufficient_substrates",
                    fallback_field_name=(
                        "kinase_condition_pairs_insufficient_substrates"
                    ),
                )
            ),
            kinase_condition_pairs_invalid_background_variance=(
                _coerce_payload_int_with_fallback(
                    payload=payload,
                    field_name="kinase_profile_pairs_invalid_background_variance",
                    fallback_field_name=(
                        "kinase_condition_pairs_invalid_background_variance"
                    ),
                )
            ),
            kinase_condition_pairs_no_finite_background_values=(
                _coerce_payload_int_with_fallback(
                    payload=payload,
                    field_name="kinase_profile_pairs_no_finite_background_values",
                    fallback_field_name=(
                        "kinase_condition_pairs_no_finite_background_values"
                    ),
                )
            ),
            kinase_condition_pairs_no_finite_substrate_values=(
                _coerce_payload_int_with_fallback(
                    payload=payload,
                    field_name="kinase_profile_pairs_no_finite_substrate_values",
                    fallback_field_name=(
                        "kinase_condition_pairs_no_finite_substrate_values"
                    ),
                )
            ),
        )


def _coerce_payload_int_with_fallback(
    *,
    payload: Mapping[str, object],
    field_name: str,
    fallback_field_name: str,
) -> int:
    if field_name in payload:
        return _coerce_payload_int(payload=payload, field_name=field_name)
    return _coerce_payload_int(payload=payload, field_name=fallback_field_name)


@dataclass(frozen=True, slots=True)
class PredMatOverlapSummary:
    """Resolved overlap diagnostics between prediction and phospho matrices."""

    overlap_count: int
    pred_mat_rows: int
    phospho_rows: int


@dataclass(frozen=True, slots=True)
class KinaseActivityInputs:
    """Trusted activity-like score inputs resolved by workflow validation."""

    pred_mat: pd.DataFrame
    phospho_matrix: pd.DataFrame
    threshold: float
    min_substrates: int
    top_n_substrates: int
    overlap_summary: PredMatOverlapSummary
    activity_input: ActivityInputMatrix | None = None

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


@dataclass(frozen=True, slots=True, init=False)
class ActivityMethodDiagnostics:
    """Typed method diagnostics carried alongside activity result matrices."""

    method_summary: ActivityMethodSummary | None
    threshold_membership_diagnostics: ActivityThresholdMembershipDiagnostics | None
    _statistics_table: pd.DataFrame | None = field(init=False, repr=False)

    def __init__(
        self,
        *,
        method_summary: ActivityMethodSummary | None = None,
        threshold_membership_diagnostics: ActivityThresholdMembershipDiagnostics
        | None = None,
        statistics_table: pd.DataFrame | None = None,
        _assume_owned: bool = False,
    ) -> None:
        _validate_optional_method_summary(method_summary)
        _validate_optional_threshold_membership_diagnostics(
            threshold_membership_diagnostics
        )
        if statistics_table is not None:
            statistics_table = ActivityStatisticsTable(
                frame=statistics_table,
                _assume_owned=_assume_owned,
            ).frame
        object.__setattr__(self, "method_summary", method_summary)
        object.__setattr__(
            self,
            "threshold_membership_diagnostics",
            threshold_membership_diagnostics,
        )
        object.__setattr__(self, "_statistics_table", statistics_table)

    @property
    def statistics_table(self) -> pd.DataFrame | None:
        """Return an optional statistics-table snapshot, not a report export."""

        return export_optional_dataframe(self._statistics_table)


class WeightedSubstrateActivityDiagnostics(ActivityMethodDiagnostics):
    """Diagnostics for simplified weighted substrate activity-like scores."""

    __slots__ = ()


class KseaZScoreActivityDiagnostics(ActivityMethodDiagnostics):
    """Diagnostics for KSEA-style kinase activity score outputs."""

    __slots__ = ()


class SsgseaSubstrateEnrichmentActivityDiagnostics(ActivityMethodDiagnostics):
    """Diagnostics for ssGSEA-style substrate enrichment activity-like scores."""

    __slots__ = ()


def _validate_optional_method_summary(
    method_summary: ActivityMethodSummary | None,
) -> None:
    if method_summary is not None and not isinstance(
        method_summary,
        ActivityMethodSummary,
    ):
        raise WorkflowBoundaryError(
            "activity_result.method_summary must be ActivityMethodSummary or None"
        )


def _validate_optional_threshold_membership_diagnostics(
    threshold_membership_diagnostics: ActivityThresholdMembershipDiagnostics | None,
) -> None:
    if threshold_membership_diagnostics is not None and not isinstance(
        threshold_membership_diagnostics,
        ActivityThresholdMembershipDiagnostics,
    ):
        raise WorkflowBoundaryError(
            "activity_result.threshold_membership_diagnostics must be "
            "ActivityThresholdMembershipDiagnostics or None"
        )


def _coerce_payload_int(*, payload: Mapping[str, object], field_name: str) -> int:
    value = payload.get(field_name, 0)
    if isinstance(value, bool | int | float | str | bytes | bytearray):
        return int(value)
    raise ValueError(f"{field_name} must be int-compatible")


@dataclass(frozen=True, slots=True, init=False)
class KinaseActivityResult:
    """Activity-like score stage outputs.

    Outputs are deliberately table-first. ``activity_matrix`` and
    ``substrate_count_matrix`` are the stable method-neutral matrices. Legacy
    weighted/KSEA sidecars remain available for existing callers but are no
    longer required by the core result contract.

    ``activity_matrix`` contains exploratory kinase activity scores or
    activity-like substrate summaries. Scores depend on substrate coverage,
    reference evidence, threshold rules, and the selected method. Sparse or
    missing substrate support weakens interpretation, and causal kinase activity
    claims require external validation.

    - ``activity_matrix``: primary kinase activity score matrix for the selected method
    - ``activity_scores``: deprecated compatibility alias for ``activity_matrix``
    - ``weighted_activity``: deprecated compatibility alias for ``activity_matrix``
    - ``p_value_matrix``: optional activity p-value matrix
    - ``q_value_matrix``: optional multiple-testing-adjusted activity q-value matrix
    - ``confidence_interval_low``: optional lower confidence interval matrix
    - ``confidence_interval_high``: optional upper confidence interval matrix
    - ``substrate_count_matrix``: method-neutral count matrix for substrates used
      per kinase-condition score when defined by the method
    - ``thresholded_substrate_mean_activity``: sample-by-kinase mean phospho
      signal over predicted substrates above threshold for legacy weighted output
    - ``thresholded_substrate_counts``: compatibility sidecar count series
    - ``activity_substrate_counts``: method-neutral condition-specific count matrix
      legacy accessor retained for KSEA-style count matrices
    - ``target_counts``: thresholded predicted target counts per kinase
    - ``target_table``: thresholded kinase-target edge table
    - ``method_diagnostics``: typed method diagnostics, not an arbitrary mapping
    - ``policy_provenance``: scientific policy records attached to this result
    - ``threshold_membership_diagnostics``: threshold inclusion rule metadata used
      by thresholded substrate membership diagnostics
    - ``activity_method``: stable method identity metadata for these outputs

    Public DataFrame helpers are defensive in-memory snapshots. They do not
    write files, format reports, plot figures, or run additional science.
    """

    activity_method: ActivityMethodMetadata
    _activity_matrix: pd.DataFrame = field(init=False, repr=False)
    _p_value_matrix: pd.DataFrame | None = field(init=False, repr=False)
    _q_value_matrix: pd.DataFrame | None = field(init=False, repr=False)
    _confidence_interval_low: pd.DataFrame | None = field(init=False, repr=False)
    _confidence_interval_high: pd.DataFrame | None = field(init=False, repr=False)
    _substrate_count_matrix: pd.DataFrame = field(init=False, repr=False)
    _thresholded_substrate_mean_activity: pd.DataFrame = field(init=False, repr=False)
    _thresholded_substrate_counts: pd.Series = field(init=False, repr=False)
    _activity_substrate_counts: pd.DataFrame | None = field(init=False, repr=False)
    _target_counts: pd.Series = field(init=False, repr=False)
    _target_table: pd.DataFrame = field(init=False, repr=False)
    _statistics_table: pd.DataFrame | None = field(init=False, repr=False)
    method_diagnostics: ActivityMethodDiagnostics
    policy_provenance: tuple[ScientificPolicyRecord, ...]
    threshold_membership_diagnostics: ActivityThresholdMembershipDiagnostics | None
    method_summary: ActivityMethodSummary | None
    input_semantics: ActivityInputSemantics
    profile_metadata: ActivityProfileMetadata

    def __init__(
        self,
        weighted_activity: pd.DataFrame | None = None,
        thresholded_substrate_mean_activity: pd.DataFrame | None = None,
        thresholded_substrate_counts: pd.Series | None = None,
        target_counts: pd.Series | None = None,
        target_table: pd.DataFrame | None = None,
        threshold_membership_diagnostics: ActivityThresholdMembershipDiagnostics
        | None = None,
        activity_substrate_counts: pd.DataFrame | None = None,
        statistics_table: pd.DataFrame | None = None,
        method_summary: ActivityMethodSummary | None = None,
        activity_method: ActivityMethodMetadata = (
            SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD
        ),
        activity_matrix: pd.DataFrame | None = None,
        p_value_matrix: pd.DataFrame | None = None,
        q_value_matrix: pd.DataFrame | None = None,
        confidence_interval_low: pd.DataFrame | None = None,
        confidence_interval_high: pd.DataFrame | None = None,
        substrate_count_matrix: pd.DataFrame | None = None,
        method_diagnostics: ActivityMethodDiagnostics | None = None,
        policy_provenance: tuple[ScientificPolicyRecord, ...]
        | list[ScientificPolicyRecord]
        | ScientificPolicyRecord
        | None = None,
        input_semantics: ActivityInputSemantics | None = None,
        profile_metadata: ActivityProfileMetadata | None = None,
    ) -> None:
        self._init_activity_result(
            weighted_activity=weighted_activity,
            thresholded_substrate_mean_activity=thresholded_substrate_mean_activity,
            thresholded_substrate_counts=thresholded_substrate_counts,
            target_counts=target_counts,
            target_table=target_table,
            threshold_membership_diagnostics=threshold_membership_diagnostics,
            activity_substrate_counts=activity_substrate_counts,
            statistics_table=statistics_table,
            method_summary=method_summary,
            activity_method=activity_method,
            activity_matrix=activity_matrix,
            p_value_matrix=p_value_matrix,
            q_value_matrix=q_value_matrix,
            confidence_interval_low=confidence_interval_low,
            confidence_interval_high=confidence_interval_high,
            substrate_count_matrix=substrate_count_matrix,
            method_diagnostics=method_diagnostics,
            policy_provenance=policy_provenance,
            input_semantics=input_semantics,
            profile_metadata=profile_metadata,
            assume_owned=False,
        )

    def _init_activity_result(
        self,
        *,
        weighted_activity: pd.DataFrame | None = None,
        thresholded_substrate_mean_activity: pd.DataFrame | None = None,
        thresholded_substrate_counts: pd.Series | None = None,
        target_counts: pd.Series | None = None,
        target_table: pd.DataFrame | None = None,
        threshold_membership_diagnostics: ActivityThresholdMembershipDiagnostics
        | None = None,
        activity_substrate_counts: pd.DataFrame | None = None,
        statistics_table: pd.DataFrame | None = None,
        method_summary: ActivityMethodSummary | None = None,
        activity_method: ActivityMethodMetadata = (
            SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD
        ),
        activity_matrix: pd.DataFrame | None = None,
        p_value_matrix: pd.DataFrame | None = None,
        q_value_matrix: pd.DataFrame | None = None,
        confidence_interval_low: pd.DataFrame | None = None,
        confidence_interval_high: pd.DataFrame | None = None,
        substrate_count_matrix: pd.DataFrame | None = None,
        method_diagnostics: ActivityMethodDiagnostics | None = None,
        policy_provenance: tuple[ScientificPolicyRecord, ...]
        | list[ScientificPolicyRecord]
        | ScientificPolicyRecord
        | None = None,
        input_semantics: ActivityInputSemantics | None = None,
        profile_metadata: ActivityProfileMetadata | None = None,
        assume_owned: bool,
    ) -> None:
        if not isinstance(activity_method, ActivityMethodMetadata):
            raise WorkflowBoundaryError(
                "activity_result.activity_method must be ActivityMethodMetadata"
            )
        if activity_matrix is not None and weighted_activity is not None:
            raise WorkflowBoundaryError(
                "activity_result must receive either activity_matrix or "
                "weighted_activity, not both"
            )
        raw_activity_matrix = (
            activity_matrix if activity_matrix is not None else weighted_activity
        )
        if raw_activity_matrix is None:
            raise WorkflowBoundaryError(
                "activity_result.activity_matrix must be provided"
            )

        activity_matrix = ActivityMatrix(
            frame=raw_activity_matrix,
            field_name="activity_result.activity_matrix",
            _assume_owned=assume_owned,
        ).frame
        should_apply_profile_axis_name = (
            input_semantics is not None and profile_metadata is not None
        )
        input_semantics, profile_metadata = _resolve_activity_result_semantics(
            activity_matrix=activity_matrix,
            input_semantics=input_semantics,
            profile_metadata=profile_metadata,
        )
        if should_apply_profile_axis_name:
            activity_matrix = _apply_activity_profile_axis_name(
                activity_matrix,
                input_semantics=input_semantics,
            )
        p_value_matrix = _validate_optional_probability_matrix(
            p_value_matrix,
            field_name="activity_result.p_value_matrix",
            assume_owned=assume_owned,
        )
        if should_apply_profile_axis_name:
            p_value_matrix = _apply_optional_activity_profile_axis_name(
                p_value_matrix,
                input_semantics=input_semantics,
            )
        q_value_matrix = _validate_optional_probability_matrix(
            q_value_matrix,
            field_name="activity_result.q_value_matrix",
            assume_owned=assume_owned,
        )
        if should_apply_profile_axis_name:
            q_value_matrix = _apply_optional_activity_profile_axis_name(
                q_value_matrix,
                input_semantics=input_semantics,
            )
        confidence_interval_low = _validate_optional_activity_matrix(
            confidence_interval_low,
            field_name="activity_result.confidence_interval_low",
            assume_owned=assume_owned,
        )
        if should_apply_profile_axis_name:
            confidence_interval_low = _apply_optional_activity_profile_axis_name(
                confidence_interval_low,
                input_semantics=input_semantics,
            )
        confidence_interval_high = _validate_optional_activity_matrix(
            confidence_interval_high,
            field_name="activity_result.confidence_interval_high",
            assume_owned=assume_owned,
        )
        if should_apply_profile_axis_name:
            confidence_interval_high = _apply_optional_activity_profile_axis_name(
                confidence_interval_high,
                input_semantics=input_semantics,
            )
        if substrate_count_matrix is None:
            substrate_count_matrix = activity_substrate_counts
        if substrate_count_matrix is None:
            substrate_count_matrix = _empty_count_matrix()
        substrate_count_matrix = ActivityCountMatrix(
            frame=substrate_count_matrix,
            field_name="activity_result.substrate_count_matrix",
            _assume_owned=assume_owned,
        ).frame
        if should_apply_profile_axis_name:
            substrate_count_matrix = _apply_activity_profile_axis_name(
                substrate_count_matrix,
                input_semantics=input_semantics,
            )
        if thresholded_substrate_mean_activity is None:
            thresholded_substrate_mean_activity = _empty_activity_matrix()
        thresholded_substrate_mean_activity = ActivityMatrix(
            frame=thresholded_substrate_mean_activity,
            field_name="activity_result.thresholded_substrate_mean_activity",
            _assume_owned=assume_owned,
        ).frame
        if should_apply_profile_axis_name:
            thresholded_substrate_mean_activity = _apply_activity_profile_axis_name(
                thresholded_substrate_mean_activity,
                input_semantics=input_semantics,
            )
        if thresholded_substrate_counts is None:
            thresholded_substrate_counts = _empty_count_series("n_substrates")
        thresholded_substrate_counts = ActivityCountSeries(
            series=thresholded_substrate_counts,
            field_name="activity_result.thresholded_substrate_counts",
            _assume_owned=assume_owned,
        ).series
        if activity_substrate_counts is not None:
            activity_substrate_counts = ActivityCountMatrix(
                frame=activity_substrate_counts,
                field_name="activity_result.activity_substrate_counts",
                _assume_owned=assume_owned,
            ).frame
            if should_apply_profile_axis_name:
                activity_substrate_counts = _apply_activity_profile_axis_name(
                    activity_substrate_counts,
                    input_semantics=input_semantics,
                )
        if target_counts is None:
            target_counts = _empty_count_series("n_targets")
        target_counts = ActivityCountSeries(
            series=target_counts,
            field_name="activity_result.target_counts",
            _assume_owned=assume_owned,
        ).series
        if target_table is None:
            target_table = _empty_target_table()
        target_table = ActivityTargetTable(
            frame=target_table,
            _assume_owned=assume_owned,
        ).frame
        if method_diagnostics is not None:
            if not isinstance(method_diagnostics, ActivityMethodDiagnostics):
                raise WorkflowBoundaryError(
                    "activity_result.method_diagnostics must be "
                    "ActivityMethodDiagnostics or None"
                )
            if method_summary is None:
                method_summary = method_diagnostics.method_summary
            if threshold_membership_diagnostics is None:
                threshold_membership_diagnostics = (
                    method_diagnostics.threshold_membership_diagnostics
                )
            if statistics_table is None:
                statistics_table = method_diagnostics.statistics_table
        if statistics_table is not None:
            statistics_table = ActivityStatisticsTable(
                frame=statistics_table,
                _assume_owned=assume_owned,
            ).frame
        _validate_optional_method_summary(method_summary)
        _validate_optional_threshold_membership_diagnostics(
            threshold_membership_diagnostics
        )
        if method_diagnostics is None:
            method_diagnostics = _build_activity_method_diagnostics(
                activity_method=activity_method,
                method_summary=method_summary,
                threshold_membership_diagnostics=threshold_membership_diagnostics,
                statistics_table=statistics_table,
            )
        policy_provenance_tuple = _coerce_policy_provenance(policy_provenance)

        object.__setattr__(self, "_activity_matrix", activity_matrix)
        object.__setattr__(self, "_p_value_matrix", p_value_matrix)
        object.__setattr__(self, "_q_value_matrix", q_value_matrix)
        object.__setattr__(
            self,
            "_confidence_interval_low",
            confidence_interval_low,
        )
        object.__setattr__(
            self,
            "_confidence_interval_high",
            confidence_interval_high,
        )
        object.__setattr__(self, "_substrate_count_matrix", substrate_count_matrix)
        object.__setattr__(
            self,
            "_thresholded_substrate_mean_activity",
            thresholded_substrate_mean_activity,
        )
        object.__setattr__(
            self,
            "_thresholded_substrate_counts",
            thresholded_substrate_counts,
        )
        object.__setattr__(
            self,
            "_activity_substrate_counts",
            activity_substrate_counts,
        )
        object.__setattr__(self, "_target_counts", target_counts)
        object.__setattr__(self, "activity_method", activity_method)
        object.__setattr__(self, "_target_table", target_table)
        object.__setattr__(self, "_statistics_table", statistics_table)
        object.__setattr__(self, "method_diagnostics", method_diagnostics)
        object.__setattr__(self, "policy_provenance", policy_provenance_tuple)
        object.__setattr__(
            self,
            "threshold_membership_diagnostics",
            threshold_membership_diagnostics,
        )
        object.__setattr__(self, "method_summary", method_summary)
        object.__setattr__(self, "input_semantics", input_semantics)
        object.__setattr__(self, "profile_metadata", profile_metadata)

    @property
    def activity_matrix(self) -> pd.DataFrame:
        """Return the primary kinase activity score matrix for the selected method."""

        return export_dataframe(self._activity_matrix)

    @property
    def activity_scores(self) -> pd.DataFrame:
        """Deprecated compatibility alias for :attr:`activity_matrix`."""

        warnings.warn(
            (
                "KinaseActivityResult.activity_scores is deprecated and will be "
                "removed in a future release; use "
                "KinaseActivityResult.activity_matrix instead."
            ),
            DeprecationWarning,
            stacklevel=2,
        )
        return self.activity_matrix

    @property
    def weighted_activity(self) -> pd.DataFrame:
        """Deprecated compatibility alias for :attr:`activity_matrix`."""

        warnings.warn(
            (
                "KinaseActivityResult.weighted_activity is deprecated and will be "
                "removed in a future release; use "
                "KinaseActivityResult.activity_matrix instead."
            ),
            DeprecationWarning,
            stacklevel=2,
        )
        return self.activity_matrix

    @property
    def p_value_matrix(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._p_value_matrix)

    @property
    def q_value_matrix(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._q_value_matrix)

    @property
    def confidence_interval_low(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._confidence_interval_low)

    @property
    def confidence_interval_high(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._confidence_interval_high)

    @property
    def substrate_count_matrix(self) -> pd.DataFrame:
        return export_dataframe(self._substrate_count_matrix)

    @property
    def thresholded_substrate_mean_activity(self) -> pd.DataFrame:
        return export_dataframe(self._thresholded_substrate_mean_activity)

    @property
    def thresholded_substrate_counts(self) -> pd.Series:
        return export_series(self._thresholded_substrate_counts)

    @property
    def activity_substrate_counts(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._activity_substrate_counts)

    @property
    def target_counts(self) -> pd.Series:
        return export_series(self._target_counts)

    @property
    def target_table(self) -> pd.DataFrame:
        return _export_public_target_table(self._target_table)

    @property
    def statistics_table(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._statistics_table)

    @property
    def count_field_semantics(self) -> dict[str, str]:
        """Return count-field meaning text as a fresh mapping snapshot."""

        if self.activity_method.is_ksea:
            return {
                "substrate_count_matrix": (
                    "condition-specific finite substrate count used for each "
                    "kinase-condition kinase activity score"
                ),
                "activity_substrate_counts": (
                    "condition-specific finite substrate count used for each "
                    "kinase-condition kinase activity score; legacy KSEA accessor"
                ),
                "thresholded_substrate_counts": (
                    "global post-threshold evidence membership count before "
                    "condition-specific finite-value filtering"
                ),
                "target_counts": (
                    "global post-threshold predicted target membership count"
                ),
            }
        if self.activity_method.activity_method_family == "substrate_set_enrichment":
            return {
                "substrate_count_matrix": (
                    "condition-specific finite substrate count used for each "
                    "kinase-condition substrate-set enrichment score"
                ),
                "activity_substrate_counts": (
                    "condition-specific finite substrate count used for each "
                    "kinase-condition substrate-set enrichment score"
                ),
                "thresholded_substrate_counts": (
                    "global kinase-substrate membership count before "
                    "condition-specific finite-value filtering"
                ),
                "target_counts": "global kinase-substrate membership count",
            }
        return {
            "substrate_count_matrix": (
                "condition-specific finite substrate count used for each primary "
                "kinase-condition activity-like score when supplied by the method"
            ),
            "thresholded_substrate_counts": (
                "global thresholded substrate membership count per kinase"
            ),
            "target_counts": "global thresholded predicted target count per kinase",
        }

    @classmethod
    def _from_owned(
        cls,
        *,
        weighted_activity: pd.DataFrame | None = None,
        thresholded_substrate_mean_activity: pd.DataFrame | None = None,
        thresholded_substrate_counts: pd.Series | None = None,
        target_counts: pd.Series | None = None,
        target_table: pd.DataFrame | None = None,
        threshold_membership_diagnostics: ActivityThresholdMembershipDiagnostics
        | None = None,
        activity_substrate_counts: pd.DataFrame | None = None,
        statistics_table: pd.DataFrame | None = None,
        method_summary: ActivityMethodSummary | None = None,
        activity_method: ActivityMethodMetadata = (
            SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD
        ),
        activity_matrix: pd.DataFrame | None = None,
        p_value_matrix: pd.DataFrame | None = None,
        q_value_matrix: pd.DataFrame | None = None,
        confidence_interval_low: pd.DataFrame | None = None,
        confidence_interval_high: pd.DataFrame | None = None,
        substrate_count_matrix: pd.DataFrame | None = None,
        method_diagnostics: ActivityMethodDiagnostics | None = None,
        policy_provenance: tuple[ScientificPolicyRecord, ...]
        | list[ScientificPolicyRecord]
        | ScientificPolicyRecord
        | None = None,
        input_semantics: ActivityInputSemantics | None = None,
        profile_metadata: ActivityProfileMetadata | None = None,
    ) -> KinaseActivityResult:
        result = object.__new__(cls)
        KinaseActivityResult._init_activity_result(
            result,
            weighted_activity=weighted_activity,
            thresholded_substrate_mean_activity=thresholded_substrate_mean_activity,
            thresholded_substrate_counts=thresholded_substrate_counts,
            activity_substrate_counts=activity_substrate_counts,
            target_counts=target_counts,
            target_table=target_table,
            threshold_membership_diagnostics=threshold_membership_diagnostics,
            statistics_table=statistics_table,
            method_summary=method_summary,
            activity_method=activity_method,
            activity_matrix=activity_matrix,
            p_value_matrix=p_value_matrix,
            q_value_matrix=q_value_matrix,
            confidence_interval_low=confidence_interval_low,
            confidence_interval_high=confidence_interval_high,
            substrate_count_matrix=substrate_count_matrix,
            method_diagnostics=method_diagnostics,
            policy_provenance=policy_provenance,
            input_semantics=input_semantics,
            profile_metadata=profile_metadata,
            assume_owned=True,
        )
        return result

    def to_dataframe(self) -> pd.DataFrame:
        """Return a primary kinase activity score snapshot isolated from this result."""

        return self.activity_matrix

    def p_value_matrix_dataframe(self) -> pd.DataFrame | None:
        """Return a p-value matrix snapshot when available, not an export."""

        return export_optional_dataframe(self._p_value_matrix)

    def q_value_matrix_dataframe(self) -> pd.DataFrame | None:
        """Return a q-value matrix snapshot when available, not an export."""

        return export_optional_dataframe(self._q_value_matrix)

    def confidence_interval_low_dataframe(self) -> pd.DataFrame | None:
        """Return a lower confidence interval snapshot when available."""

        return export_optional_dataframe(self._confidence_interval_low)

    def confidence_interval_high_dataframe(self) -> pd.DataFrame | None:
        """Return an upper confidence interval snapshot when available."""

        return export_optional_dataframe(self._confidence_interval_high)

    def substrate_count_matrix_dataframe(self) -> pd.DataFrame:
        """Return a substrate-count matrix snapshot isolated from this result."""

        return export_dataframe(self._substrate_count_matrix)

    def thresholded_substrate_mean_activity_dataframe(self) -> pd.DataFrame:
        """Return a thresholded-mean snapshot isolated from this result."""

        return export_dataframe(self._thresholded_substrate_mean_activity)

    def target_table_dataframe(self) -> pd.DataFrame:
        """Return a target-table snapshot isolated from this result."""

        return _export_public_target_table(self._target_table)

    def statistics_table_dataframe(self) -> pd.DataFrame | None:
        """Return a statistics-table snapshot, not a formatted report."""

        return export_optional_dataframe(self._statistics_table)


def _resolve_activity_result_semantics(
    *,
    activity_matrix: pd.DataFrame,
    input_semantics: ActivityInputSemantics | None,
    profile_metadata: ActivityProfileMetadata | None,
) -> tuple[ActivityInputSemantics, ActivityProfileMetadata]:
    profile_ids = tuple(str(column) for column in activity_matrix.columns)
    if input_semantics is None and profile_metadata is None:
        warnings.warn(
            (
                "KinaseActivityResult constructed without explicit activity input "
                "semantics is deprecated; treating activity columns as "
                "sample/profile labels with sample-level abundance semantics."
            ),
            DeprecationWarning,
            stacklevel=3,
        )
        input_semantics = ActivityInputSemantics(
            profile_axis=ActivityProfileAxis.SAMPLE,
            quantitative_semantics=(
                ActivityQuantitativeSemantics.SAMPLE_LEVEL_ABUNDANCE
            ),
        )
        profile_metadata = ActivityProfileMetadata(
            axis=ActivityProfileAxis.SAMPLE,
            profile_ids=profile_ids,
            sample_ids=profile_ids,
        )
    elif input_semantics is None or profile_metadata is None:
        raise WorkflowBoundaryError(
            "activity_result.input_semantics and activity_result.profile_metadata "
            "must be provided together"
        )
    if not isinstance(input_semantics, ActivityInputSemantics):
        raise WorkflowBoundaryError(
            "activity_result.input_semantics must be ActivityInputSemantics"
        )
    if not isinstance(profile_metadata, ActivityProfileMetadata):
        raise WorkflowBoundaryError(
            "activity_result.profile_metadata must be ActivityProfileMetadata"
        )
    if profile_metadata.axis is not input_semantics.profile_axis:
        raise WorkflowBoundaryError(
            "activity_result.profile_metadata.axis must match "
            "activity_result.input_semantics.profile_axis"
        )
    observed_profile_ids = tuple(str(value) for value in profile_metadata.profile_ids)
    if observed_profile_ids != profile_ids:
        raise WorkflowBoundaryError(
            "activity_result.profile_metadata.profile_ids must exactly match "
            "activity_result.activity_matrix columns; "
            f"expected={profile_ids!r}, got={observed_profile_ids!r}"
        )
    return input_semantics, profile_metadata


def _activity_profile_axis_name(
    input_semantics: ActivityInputSemantics,
) -> str:
    if input_semantics.has_real_condition_contract:
        return "condition"
    return "profile_id"


def _apply_activity_profile_axis_name(
    frame: pd.DataFrame,
    *,
    input_semantics: ActivityInputSemantics,
) -> pd.DataFrame:
    if input_semantics.profile_axis is ActivityProfileAxis.SAMPLE:
        return frame
    renamed = frame.copy(deep=False)
    renamed.columns = renamed.columns.copy()
    renamed.columns.name = _activity_profile_axis_name(input_semantics)
    return renamed


def _apply_optional_activity_profile_axis_name(
    frame: pd.DataFrame | None,
    *,
    input_semantics: ActivityInputSemantics,
) -> pd.DataFrame | None:
    if frame is None:
        return None
    return _apply_activity_profile_axis_name(
        frame,
        input_semantics=input_semantics,
    )


def _validate_optional_activity_matrix(
    matrix: pd.DataFrame | None,
    *,
    field_name: str,
    assume_owned: bool,
) -> pd.DataFrame | None:
    if matrix is None:
        return None
    return ActivityMatrix(
        frame=matrix,
        field_name=field_name,
        _assume_owned=assume_owned,
    ).frame


def _validate_optional_probability_matrix(
    matrix: pd.DataFrame | None,
    *,
    field_name: str,
    assume_owned: bool,
) -> pd.DataFrame | None:
    matrix = _validate_optional_activity_matrix(
        matrix,
        field_name=field_name,
        assume_owned=assume_owned,
    )
    if matrix is None:
        return None
    values = matrix.to_numpy(dtype="float64", copy=False)
    finite_mask = np.isfinite(values)
    finite_values = values[finite_mask]
    if ((finite_values < 0.0) | (finite_values > 1.0)).any():
        raise PhosPyValidationError(
            f"{field_name} must be between 0.0 and 1.0 when present"
        )
    return matrix


def _build_activity_method_diagnostics(
    *,
    activity_method: ActivityMethodMetadata,
    method_summary: ActivityMethodSummary | None,
    threshold_membership_diagnostics: ActivityThresholdMembershipDiagnostics | None,
    statistics_table: pd.DataFrame | None,
) -> ActivityMethodDiagnostics:
    diagnostics_cls: type[ActivityMethodDiagnostics]
    if activity_method.is_ksea:
        diagnostics_cls = KseaZScoreActivityDiagnostics
    elif (
        activity_method.activity_method_id
        == SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_METHOD.activity_method_id
    ):
        diagnostics_cls = SsgseaSubstrateEnrichmentActivityDiagnostics
    elif (
        activity_method.activity_method_id
        == SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD.activity_method_id
    ):
        diagnostics_cls = WeightedSubstrateActivityDiagnostics
    else:
        diagnostics_cls = ActivityMethodDiagnostics
    return diagnostics_cls(
        method_summary=method_summary,
        threshold_membership_diagnostics=threshold_membership_diagnostics,
        statistics_table=statistics_table,
        _assume_owned=True,
    )


def _coerce_policy_provenance(
    policy_provenance: tuple[ScientificPolicyRecord, ...]
    | list[ScientificPolicyRecord]
    | ScientificPolicyRecord
    | None,
) -> tuple[ScientificPolicyRecord, ...]:
    if policy_provenance is None:
        return ()
    if isinstance(policy_provenance, ScientificPolicyRecord):
        return (policy_provenance,)
    if not isinstance(policy_provenance, tuple | list):
        raise WorkflowBoundaryError(
            "activity_result.policy_provenance must contain ScientificPolicyRecord "
            "objects"
        )
    for record in policy_provenance:
        if not isinstance(record, ScientificPolicyRecord):
            raise WorkflowBoundaryError(
                "activity_result.policy_provenance must contain "
                "ScientificPolicyRecord objects"
            )
    return tuple(policy_provenance)


def _empty_activity_matrix() -> pd.DataFrame:
    return pd.DataFrame(dtype=float)


def _empty_count_matrix() -> pd.DataFrame:
    return pd.DataFrame(dtype="int64")


def _empty_count_series(name: str) -> pd.Series:
    series = pd.Series(dtype="int64", name=name)
    series.index.name = "kinase"
    return series


def _empty_target_table() -> pd.DataFrame:
    return pd.DataFrame(columns=["site_id", "kinase", "score"])


def _export_public_target_table(table: pd.DataFrame) -> pd.DataFrame:
    exported = export_dataframe(table)
    if {"site_key", "display_id"}.issubset(exported.columns):
        return exported
    legacy_columns = ["site_id", "kinase", "score"]
    if not all(column in exported.columns for column in legacy_columns):
        return exported
    return exported.loc[:, legacy_columns]
