"""Signalome domain models."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Literal

import pandas as pd

from phospy.errors.validation import PhosPyValidationError
from phospy.frames.comparison import dataframe_equals, optional_dataframe_equals
from phospy.frames.ownership import (
    export_dataframe,
    export_optional_dataframe,
    own_dataframe,
    own_optional_dataframe,
)
from phospy.provenance.models import TableFingerprint
from phospy.science.configs import (
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
    SignalomeScorePreconditioningPolicy,
)

SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS = "correlation_thresholds"
SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT = "explicit_module_count"
SignalomeModuleSelectionStrategy = Literal[
    "correlation_thresholds",
    "explicit_module_count",
]
SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_STABLE = "stable"
SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_UNSTABLE = "unstable"
SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_NOT_COMPUTABLE = "not_computable"
SignalomeModuleSelectionStabilityStatus = Literal[
    "stable",
    "unstable",
    "not_computable",
]
SIGNALOME_MODULE_SELECTION_STABILITY_METHOD = (
    "seeded_score_perturbation_and_threshold_grid"
)
SIGNALOME_MODULE_SELECTION_STABILITY_VERSION = "1"
SIGNALOME_CORRELATION_STATUS_FINITE = "finite"
SIGNALOME_CORRELATION_STATUS_CONSTANT_PROFILE = "constant_profile"
SIGNALOME_CORRELATION_STATUS_INSUFFICIENT_OBSERVATIONS = "insufficient_observations"
SIGNALOME_CORRELATION_STATUS_MISSING_VALUES = "missing_values"
SIGNALOME_CORRELATION_STATUS_NON_FINITE_VALUES = "non_finite_values"
SIGNALOME_CORRELATION_STATUS_UNDEFINED = "undefined"
SignalomeCorrelationStatus = Literal[
    "finite",
    "constant_profile",
    "insufficient_observations",
    "missing_values",
    "non_finite_values",
    "undefined",
]


@dataclass(frozen=True, slots=True)
class SignalomeClusterCandidateScore:
    """One candidate module-count score summary."""

    min_median_correlation: float
    mean_median_correlation: float


@dataclass(frozen=True, slots=True)
class SignalomeModuleSelectionThresholdSensitivityRecord:
    """Selected count for one threshold-grid sensitivity point."""

    primary_threshold: float
    fallback_threshold: float
    selected_module_count: int
    threshold_used: float | None


@dataclass(frozen=True, slots=True)
class SignalomeModuleSelectionThresholdSensitivity:
    """Structured threshold-sensitivity summary for automatic module selection."""

    method: str
    version: str
    records: tuple[SignalomeModuleSelectionThresholdSensitivityRecord, ...]
    selected_count_frequency: dict[int, int]
    disagrees_with_selected_count: bool

    def __post_init__(self) -> None:
        records = tuple(self.records)
        frequencies = {
            int(key): int(value) for key, value in self.selected_count_frequency.items()
        }
        object.__setattr__(self, "records", records)
        object.__setattr__(self, "selected_count_frequency", frequencies)


@dataclass(frozen=True, slots=True)
class SignalomeModuleSelectionAssignmentSimilaritySummary:
    """Summary of partition agreement across seeded perturbations."""

    metric: str
    evaluated_perturbations: int
    minimum: float | None
    median: float | None
    mean: float | None
    maximum: float | None


@dataclass(frozen=True, slots=True)
class SignalomeModuleSelectionStabilityReport:
    """Descriptive stability/sensitivity report for module-count selection."""

    evaluation_method: str
    evaluation_version: str
    seed_policy: str
    random_seed: int | None
    perturbation_count: int
    selected_count_frequency: dict[int, int]
    assignment_similarity_metric: str
    assignment_similarity: SignalomeModuleSelectionAssignmentSimilaritySummary
    threshold_sensitivity: SignalomeModuleSelectionThresholdSensitivity
    status: SignalomeModuleSelectionStabilityStatus
    limitations: tuple[str, ...]
    not_computable_reason: str | None = None
    base_selected_module_count: int = 1
    input_site_count: int = 0
    input_dimension_count: int = 0

    def __post_init__(self) -> None:
        frequencies = {
            int(key): int(value) for key, value in self.selected_count_frequency.items()
        }
        limitations = tuple(str(value) for value in self.limitations)
        object.__setattr__(self, "selected_count_frequency", frequencies)
        object.__setattr__(self, "limitations", limitations)


def default_signalome_module_selection_stability_report() -> (
    SignalomeModuleSelectionStabilityReport
):
    """Return a stable placeholder module-selection stability report."""

    threshold_sensitivity = SignalomeModuleSelectionThresholdSensitivity(
        method="not_captured",
        version="0",
        records=(),
        selected_count_frequency={},
        disagrees_with_selected_count=False,
    )
    assignment_similarity = SignalomeModuleSelectionAssignmentSimilaritySummary(
        metric="not_captured",
        evaluated_perturbations=0,
        minimum=None,
        median=None,
        mean=None,
        maximum=None,
    )
    return SignalomeModuleSelectionStabilityReport(
        evaluation_method="not_captured",
        evaluation_version="0",
        seed_policy="not_captured",
        random_seed=None,
        perturbation_count=0,
        selected_count_frequency={},
        assignment_similarity_metric="not_captured",
        assignment_similarity=assignment_similarity,
        threshold_sensitivity=threshold_sensitivity,
        status=SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_NOT_COMPUTABLE,
        limitations=(
            "Module-selection stability diagnostics were not captured for this result.",
        ),
        not_computable_reason="module selection stability diagnostics were not captured",
        base_selected_module_count=1,
        input_site_count=0,
        input_dimension_count=0,
    )


@dataclass(frozen=True, slots=True)
class SignalomeModuleSelectionDiagnostics:
    """Structured module-selection diagnostics."""

    strategy: SignalomeModuleSelectionStrategy
    selected_module_count: int
    requested_module_count: int | None
    threshold_used: float | None
    max_clusters_evaluated: int
    candidate_scores: dict[int, SignalomeClusterCandidateScore]
    reason: str
    zero_variance_profile_count: int = 0
    near_constant_profile_count: int = 0
    excluded_from_correlation_count: int = 0
    stability_report: SignalomeModuleSelectionStabilityReport = field(
        default_factory=default_signalome_module_selection_stability_report
    )

    @property
    def used_automatic_selection(self) -> bool:
        return self.requested_module_count is None


@dataclass(frozen=True, slots=True)
class SignalomeClusteringPreparationDiagnostics:
    """Structured diagnostics for signalome clustering matrix preparation."""

    preparation_policy_id: str
    input_dimension_count: int
    retained_dimension_count: int
    retained_dimension_labels: tuple[str, ...]
    dropped_fully_missing_dimension_count: int
    dropped_fully_missing_dimension_labels: tuple[str, ...]
    dropped_fully_missing_dimension_preview: tuple[str, ...]
    dropped_fully_missing_value_count: int
    non_finite_input_value_count: int
    missing_after_non_finite_normalization_count: int
    imputed_value_count: int
    imputed_value_counts_by_dimension: Mapping[str, int]
    prepared_matrix_fingerprint: TableFingerprint | None = None

    def __post_init__(self) -> None:
        retained = tuple(str(value) for value in self.retained_dimension_labels)
        dropped = tuple(
            str(value) for value in self.dropped_fully_missing_dimension_labels
        )
        preview = tuple(
            str(value) for value in self.dropped_fully_missing_dimension_preview
        )
        imputation_counts = {
            str(key): int(value)
            for key, value in self.imputed_value_counts_by_dimension.items()
        }
        object.__setattr__(self, "retained_dimension_labels", retained)
        object.__setattr__(self, "dropped_fully_missing_dimension_labels", dropped)
        object.__setattr__(self, "dropped_fully_missing_dimension_preview", preview)
        object.__setattr__(
            self,
            "imputed_value_counts_by_dimension",
            MappingProxyType(imputation_counts),
        )


@dataclass(frozen=True, slots=True)
class SignalomeScorePreconditioningDiagnostics:
    """Structured score preconditioning diagnostics for signalome inputs."""

    input_row_count: int
    dropped_all_missing_row_count: int
    retained_row_count: int
    policy: SignalomeScorePreconditioningPolicy = (
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
    )


@dataclass(frozen=True, slots=True)
class SignalomeAlignmentInputDiagnostics:
    """Structured provided/retained/dropped counts for one aligned input lane."""

    provided_count: int
    retained_count: int
    dropped_count: int
    dropped_reasons: dict[str, int] = field(default_factory=dict)


@dataclass(frozen=True, slots=True, init=False)
class SignalomeAlignmentDiagnostics:
    """Structured alignment diagnostics across signalome scientific inputs."""

    dataset_sites: SignalomeAlignmentInputDiagnostics
    prediction_score_sites: SignalomeAlignmentInputDiagnostics
    downstream_score_sites: SignalomeAlignmentInputDiagnostics
    kinases: SignalomeAlignmentInputDiagnostics
    protein_group_ids: SignalomeAlignmentInputDiagnostics

    def __init__(
        self,
        *,
        dataset_sites: SignalomeAlignmentInputDiagnostics,
        prediction_score_sites: SignalomeAlignmentInputDiagnostics,
        downstream_score_sites: SignalomeAlignmentInputDiagnostics,
        kinases: SignalomeAlignmentInputDiagnostics,
        protein_group_ids: SignalomeAlignmentInputDiagnostics | None = None,
        protein_identifiers: SignalomeAlignmentInputDiagnostics | None = None,
    ) -> None:
        """Create diagnostics, accepting legacy protein_identifiers alias."""

        if protein_group_ids is None:
            if protein_identifiers is None:
                raise TypeError(
                    "SignalomeAlignmentDiagnostics requires protein_group_ids"
                )
            protein_group_ids = protein_identifiers
        elif (
            protein_identifiers is not None and protein_identifiers != protein_group_ids
        ):
            raise ValueError(
                "SignalomeAlignmentDiagnostics received conflicting "
                "protein_group_ids and legacy protein_identifiers diagnostics"
            )
        object.__setattr__(self, "dataset_sites", dataset_sites)
        object.__setattr__(self, "prediction_score_sites", prediction_score_sites)
        object.__setattr__(self, "downstream_score_sites", downstream_score_sites)
        object.__setattr__(self, "kinases", kinases)
        object.__setattr__(self, "protein_group_ids", protein_group_ids)

    @property
    def protein_identifiers(self) -> SignalomeAlignmentInputDiagnostics:
        """Legacy alias for Signalome protein grouping diagnostics."""

        return self.protein_group_ids


@dataclass(frozen=True, slots=True)
class SignalomeNetworkCorrelationDiagnostics:
    """Structured diagnostics for kinase-network correlation eligibility."""

    total_candidate_correlations: int
    finite_correlations: int
    undefined_correlations: int
    constant_profile_correlations: int
    insufficient_observation_correlations: int
    missing_value_correlations: int
    non_finite_value_correlations: int
    edges_created: int
    edges_skipped_non_finite_correlation: int
    edges_skipped_below_threshold: int = 0
    edges_skipped_insufficient_paired_observations: int = 0
    edges_skipped_constant_profile: int = 0
    edges_skipped_missing_score: int = 0
    edges_skipped_non_finite_score: int = 0
    edges_skipped_undefined_correlation: int = 0


def default_signalome_module_selection_diagnostics() -> (
    SignalomeModuleSelectionDiagnostics
):
    """Return a stable placeholder diagnostics payload."""

    return SignalomeModuleSelectionDiagnostics(
        strategy=SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
        selected_module_count=1,
        requested_module_count=None,
        threshold_used=None,
        max_clusters_evaluated=1,
        candidate_scores={},
        reason="module selection diagnostics were not captured",
    )


def default_signalome_clustering_preparation_diagnostics() -> (
    SignalomeClusteringPreparationDiagnostics
):
    """Return a stable placeholder clustering-preparation diagnostics payload."""

    return SignalomeClusteringPreparationDiagnostics(
        preparation_policy_id="not_captured",
        input_dimension_count=0,
        retained_dimension_count=0,
        retained_dimension_labels=(),
        dropped_fully_missing_dimension_count=0,
        dropped_fully_missing_dimension_labels=(),
        dropped_fully_missing_dimension_preview=(),
        dropped_fully_missing_value_count=0,
        non_finite_input_value_count=0,
        missing_after_non_finite_normalization_count=0,
        imputed_value_count=0,
        imputed_value_counts_by_dimension={},
        prepared_matrix_fingerprint=None,
    )


def default_signalome_score_preconditioning_diagnostics() -> (
    SignalomeScorePreconditioningDiagnostics
):
    """Return a stable placeholder score-preconditioning diagnostics payload."""

    return SignalomeScorePreconditioningDiagnostics(
        input_row_count=0,
        dropped_all_missing_row_count=0,
        retained_row_count=0,
        policy=SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    )


def default_signalome_alignment_diagnostics() -> SignalomeAlignmentDiagnostics:
    """Return stable placeholder alignment diagnostics payload."""

    def _empty() -> SignalomeAlignmentInputDiagnostics:
        return SignalomeAlignmentInputDiagnostics(
            provided_count=0,
            retained_count=0,
            dropped_count=0,
            dropped_reasons={},
        )

    return SignalomeAlignmentDiagnostics(
        dataset_sites=_empty(),
        prediction_score_sites=_empty(),
        downstream_score_sites=_empty(),
        kinases=_empty(),
        protein_group_ids=_empty(),
    )


def default_signalome_network_correlation_diagnostics() -> (
    SignalomeNetworkCorrelationDiagnostics
):
    """Return a stable placeholder network-correlation diagnostics payload."""

    return SignalomeNetworkCorrelationDiagnostics(
        total_candidate_correlations=0,
        finite_correlations=0,
        undefined_correlations=0,
        constant_profile_correlations=0,
        insufficient_observation_correlations=0,
        missing_value_correlations=0,
        non_finite_value_correlations=0,
        edges_created=0,
        edges_skipped_non_finite_correlation=0,
        edges_skipped_below_threshold=0,
        edges_skipped_insufficient_paired_observations=0,
        edges_skipped_constant_profile=0,
        edges_skipped_missing_score=0,
        edges_skipped_non_finite_score=0,
        edges_skipped_undefined_correlation=0,
    )


@dataclass(frozen=True, slots=True, init=False, eq=False)
class SignalomeAssignments:
    """Score-derived signalome module assignment table.

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for explicit table-content comparison.
    """

    __hash__ = object.__hash__

    _table: pd.DataFrame = field(init=False, repr=False)

    def __init__(self, table: pd.DataFrame, _assume_owned: bool = False) -> None:
        from phospy.science.tables.signalome import SignalomeAssignmentsTable

        table = own_dataframe(
            table,
            field_name="signalome_result.module_assignments.table",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        object.__setattr__(
            self,
            "_table",
            SignalomeAssignmentsTable(frame=table, _assume_owned=True).frame,
        )

    @classmethod
    def _from_owned(cls, *, table: pd.DataFrame) -> SignalomeAssignments:
        return cls(table=table, _assume_owned=True)

    @classmethod
    def from_trusted_owned(cls, *, table: pd.DataFrame) -> SignalomeAssignments:
        """Construct from an already-owned assignment table."""

        return cls(table=table, _assume_owned=True)

    def to_pandas(self) -> pd.DataFrame:
        """Return an assignments snapshot isolated from this object."""

        return export_dataframe(self._table)

    @property
    def table(self) -> pd.DataFrame:
        return export_dataframe(self._table)

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another assignment container has the same table."""

        if not isinstance(other, SignalomeAssignments):
            return False
        return dataframe_equals(self._table, other._table)


@dataclass(frozen=True, slots=True, init=False, eq=False)
class SignalomeModules:
    """Candidate kinase-supported module summary table.

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for explicit table-content comparison.
    """

    __hash__ = object.__hash__

    _table: pd.DataFrame = field(init=False, repr=False)

    def __init__(self, table: pd.DataFrame, _assume_owned: bool = False) -> None:
        from phospy.science.tables.signalome import SignalomeModulesTable

        table = own_dataframe(
            table,
            field_name="signalome_result.signalome_modules.table",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        object.__setattr__(
            self,
            "_table",
            SignalomeModulesTable(frame=table, _assume_owned=True).frame,
        )

    @classmethod
    def _from_owned(cls, *, table: pd.DataFrame) -> SignalomeModules:
        return cls(table=table, _assume_owned=True)

    @classmethod
    def from_trusted_owned(cls, *, table: pd.DataFrame) -> SignalomeModules:
        """Construct from an already-owned module table."""

        return cls(table=table, _assume_owned=True)

    def to_pandas(self) -> pd.DataFrame:
        """Return a modules snapshot isolated from this object."""

        return export_dataframe(self._table)

    @property
    def table(self) -> pd.DataFrame:
        return export_dataframe(self._table)

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another module container has the same table."""

        if not isinstance(other, SignalomeModules):
            return False
        return dataframe_equals(self._table, other._table)


@dataclass(frozen=True, slots=True, init=False, eq=False)
class KinaseNetwork:
    """Kinase score-profile association tables derived from signalome analysis.

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for explicit network-content comparison.
    """

    __hash__ = object.__hash__

    correlation_diagnostics: SignalomeNetworkCorrelationDiagnostics = field(
        default_factory=default_signalome_network_correlation_diagnostics
    )
    _edges: pd.DataFrame = field(init=False, repr=False)
    _nodes: pd.DataFrame | None = field(init=False, repr=False)
    _candidate_correlations: pd.DataFrame | None = field(init=False, repr=False)

    def __init__(
        self,
        edges: pd.DataFrame,
        nodes: pd.DataFrame | None = None,
        candidate_correlations: pd.DataFrame | None = None,
        correlation_diagnostics: SignalomeNetworkCorrelationDiagnostics | None = None,
        _assume_owned: bool = False,
    ) -> None:
        from phospy.science.tables.signalome import (
            KinaseNetworkCandidateCorrelationsTable,
            KinaseNetworkEdgesTable,
            KinaseNetworkNodesTable,
        )

        object.__setattr__(
            self,
            "correlation_diagnostics",
            (
                default_signalome_network_correlation_diagnostics()
                if correlation_diagnostics is None
                else correlation_diagnostics
            ),
        )

        edges = own_dataframe(
            edges,
            field_name="signalome_result.kinase_network.edges",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        nodes = own_optional_dataframe(
            nodes,
            field_name="signalome_result.kinase_network.nodes",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        candidate_correlations = own_optional_dataframe(
            candidate_correlations,
            field_name="signalome_result.kinase_network.candidate_correlations",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        if not isinstance(
            self.correlation_diagnostics,
            SignalomeNetworkCorrelationDiagnostics,
        ):
            raise PhosPyValidationError(
                "signalome_result.kinase_network.correlation_diagnostics must be "
                "SignalomeNetworkCorrelationDiagnostics"
            )
        edges = KinaseNetworkEdgesTable(frame=edges, _assume_owned=True).frame
        if nodes is not None:
            nodes = KinaseNetworkNodesTable(frame=nodes, _assume_owned=True).frame
        if candidate_correlations is not None:
            candidate_correlations = KinaseNetworkCandidateCorrelationsTable(
                frame=candidate_correlations,
                _assume_owned=True,
            ).frame
        object.__setattr__(self, "_edges", edges)
        object.__setattr__(self, "_nodes", nodes)
        object.__setattr__(self, "_candidate_correlations", candidate_correlations)

    @property
    def edges(self) -> pd.DataFrame:
        return export_dataframe(self._edges)

    @property
    def nodes(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._nodes)

    @property
    def candidate_correlations(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._candidate_correlations)

    @classmethod
    def _from_owned(
        cls,
        *,
        edges: pd.DataFrame,
        nodes: pd.DataFrame | None = None,
        candidate_correlations: pd.DataFrame | None = None,
        correlation_diagnostics: SignalomeNetworkCorrelationDiagnostics | None = None,
    ) -> KinaseNetwork:
        return cls(
            edges=edges,
            nodes=nodes,
            candidate_correlations=candidate_correlations,
            correlation_diagnostics=correlation_diagnostics,
            _assume_owned=True,
        )

    @classmethod
    def from_trusted_owned(
        cls,
        *,
        edges: pd.DataFrame,
        nodes: pd.DataFrame | None = None,
        candidate_correlations: pd.DataFrame | None = None,
        correlation_diagnostics: SignalomeNetworkCorrelationDiagnostics | None = None,
    ) -> KinaseNetwork:
        """Construct from already-owned network tables."""

        return cls(
            edges=edges,
            nodes=nodes,
            candidate_correlations=candidate_correlations,
            correlation_diagnostics=correlation_diagnostics,
            _assume_owned=True,
        )

    def to_pandas(self) -> pd.DataFrame:
        """Return an edges snapshot isolated from this network."""

        return export_dataframe(self._edges)

    def nodes_dataframe(self) -> pd.DataFrame | None:
        """Return an optional nodes snapshot isolated from this network."""

        return export_optional_dataframe(self._nodes)

    def candidate_correlations_dataframe(self) -> pd.DataFrame | None:
        """Return optional candidate correlations isolated from this network."""

        return export_optional_dataframe(self._candidate_correlations)

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another network has the same scientific content."""

        if not isinstance(other, KinaseNetwork):
            return False
        return (
            self.correlation_diagnostics == other.correlation_diagnostics
            and dataframe_equals(self._edges, other._edges)
            and optional_dataframe_equals(self._nodes, other._nodes)
            and optional_dataframe_equals(
                self._candidate_correlations,
                other._candidate_correlations,
            )
        )


__all__ = [
    "KinaseNetwork",
    "SignalomeAlignmentDiagnostics",
    "SignalomeAlignmentInputDiagnostics",
    "SignalomeClusteringPreparationDiagnostics",
    "SIGNALOME_CORRELATION_STATUS_CONSTANT_PROFILE",
    "SIGNALOME_CORRELATION_STATUS_FINITE",
    "SIGNALOME_CORRELATION_STATUS_INSUFFICIENT_OBSERVATIONS",
    "SIGNALOME_CORRELATION_STATUS_MISSING_VALUES",
    "SIGNALOME_CORRELATION_STATUS_NON_FINITE_VALUES",
    "SIGNALOME_CORRELATION_STATUS_UNDEFINED",
    "SignalomeAssignments",
    "SignalomeModuleSelectionAssignmentSimilaritySummary",
    "SignalomeClusterCandidateScore",
    "SignalomeCorrelationStatus",
    "SignalomeModuleSelectionDiagnostics",
    "SignalomeModuleSelectionStabilityReport",
    "SignalomeModuleSelectionStabilityStatus",
    "SignalomeModuleSelectionStrategy",
    "SignalomeModuleSelectionThresholdSensitivity",
    "SignalomeModuleSelectionThresholdSensitivityRecord",
    "SignalomeNetworkCorrelationDiagnostics",
    "SignalomeScorePreconditioningDiagnostics",
    "SignalomeScorePreconditioningPolicy",
    "SignalomeModules",
    "SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT",
    "SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP",
    "default_signalome_alignment_diagnostics",
    "default_signalome_clustering_preparation_diagnostics",
    "default_signalome_score_preconditioning_diagnostics",
    "SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS",
    "SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT",
    "SIGNALOME_MODULE_SELECTION_STABILITY_METHOD",
    "SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_NOT_COMPUTABLE",
    "SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_STABLE",
    "SIGNALOME_MODULE_SELECTION_STABILITY_STATUS_UNSTABLE",
    "SIGNALOME_MODULE_SELECTION_STABILITY_VERSION",
    "default_signalome_network_correlation_diagnostics",
    "default_signalome_module_selection_diagnostics",
    "default_signalome_module_selection_stability_report",
]
