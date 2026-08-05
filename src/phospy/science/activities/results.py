"""Activity result-table model."""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import pandas as pd

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.frames.ownership import (
    export_dataframe,
    export_optional_dataframe,
    export_series,
)
from phospy.provenance.scientific_policy_models import ScientificPolicyRecord
from phospy.science.activities.diagnostics import (
    ActivityMethodDiagnostics,
    _build_activity_method_diagnostics,
    _validate_optional_method_summary,
    _validate_optional_threshold_membership_diagnostics,
)
from phospy.science.activities.membership import ActivityMembershipSelection
from phospy.science.activities.method_models import (
    SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD,
    ActivityMethodMetadata,
    ActivityMethodSummary,
)
from phospy.science.activities.result_validation import (
    _apply_activity_profile_axis_name,
    _apply_optional_activity_profile_axis_name,
    _empty_activity_matrix,
    _empty_count_matrix,
    _empty_count_series,
    _empty_target_table,
    _export_public_target_table,
    _resolve_activity_result_semantics,
    _validate_activity_statistics_profile_contract,
    _validate_optional_activity_matrix,
    _validate_optional_probability_matrix,
)
from phospy.science.activities.semantics import (
    ActivityInputSemantics,
    ActivityProfileMetadata,
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
      per kinase-profile score when defined by the method
    - ``thresholded_substrate_mean_activity``: sample-by-kinase mean phospho
      signal over predicted substrates above threshold for legacy weighted output
    - ``thresholded_substrate_counts``: compatibility sidecar count series
    - ``activity_substrate_counts``: method-neutral profile-specific count matrix
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
    membership_selection: ActivityMembershipSelection | None

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
        membership_selection: ActivityMembershipSelection | None = None,
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
            membership_selection=membership_selection,
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
        membership_selection: ActivityMembershipSelection | None = None,
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
            _validate_activity_statistics_profile_contract(
                statistics_table=statistics_table,
                input_semantics=input_semantics,
                profile_metadata=profile_metadata,
            )
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
        if membership_selection is not None and not isinstance(
            membership_selection,
            ActivityMembershipSelection,
        ):
            raise WorkflowBoundaryError(
                "activity_result.membership_selection must be "
                "ActivityMembershipSelection or None"
            )

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
        object.__setattr__(self, "membership_selection", membership_selection)

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

        profile_specific = (
            "condition-specific"
            if self.input_semantics.has_real_condition_contract
            else "profile-specific"
        )
        kinase_profile_pair = (
            "kinase-condition"
            if self.input_semantics.has_real_condition_contract
            else "kinase-profile"
        )
        if self.activity_method.is_ksea:
            return {
                "substrate_count_matrix": (
                    f"{profile_specific} finite substrate count used for each "
                    f"{kinase_profile_pair} kinase activity score"
                ),
                "activity_substrate_counts": (
                    f"{profile_specific} finite substrate count used for each "
                    f"{kinase_profile_pair} kinase activity score; legacy KSEA "
                    "accessor"
                ),
                "thresholded_substrate_counts": (
                    "global post-threshold evidence membership count before "
                    f"{profile_specific} finite-value filtering"
                ),
                "target_counts": (
                    "global post-threshold predicted target membership count"
                ),
            }
        if self.activity_method.activity_method_family == "substrate_set_enrichment":
            return {
                "substrate_count_matrix": (
                    f"{profile_specific} finite substrate count used for each "
                    f"{kinase_profile_pair} substrate-set enrichment score"
                ),
                "activity_substrate_counts": (
                    f"{profile_specific} finite substrate count used for each "
                    f"{kinase_profile_pair} substrate-set enrichment score"
                ),
                "thresholded_substrate_counts": (
                    "global kinase-substrate membership count before "
                    f"{profile_specific} finite-value filtering"
                ),
                "target_counts": "global kinase-substrate membership count",
            }
        return {
            "substrate_count_matrix": (
                f"{profile_specific} finite substrate count used for each primary "
                f"{kinase_profile_pair} activity-like score when supplied by the "
                "method"
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
        membership_selection: ActivityMembershipSelection | None = None,
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
            membership_selection=membership_selection,
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

    def legacy_condition_statistics_table_dataframe(self) -> pd.DataFrame | None:
        """Return an old condition-shaped statistics-table compatibility snapshot."""

        warnings.warn(
            (
                "KinaseActivityResult.legacy_condition_statistics_table_dataframe() "
                "is deprecated; use statistics_table_dataframe() or statistics_table "
                "and the canonical profile_id column. This adapter adds "
                "condition=profile_id only for compatibility and does not establish "
                "a biological condition contract."
            ),
            DeprecationWarning,
            stacklevel=2,
        )
        table = export_optional_dataframe(self._statistics_table)
        if table is None:
            return None
        legacy_table = table.copy(deep=True)
        legacy_table.loc[:, "condition"] = legacy_table.loc[:, "profile_id"]
        return legacy_table


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


__all__ = [
    "KinaseActivityResult",
]
