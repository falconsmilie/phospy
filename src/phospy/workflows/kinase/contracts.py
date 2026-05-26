"""Kinase workflow internal contracts."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Protocol

import pandas as pd

from phospy.contracts.configs import (
    KinaseActivityMethod,
    KinaseActivityPValueMethod,
    KinaseAdaptivePolicy,
    KinasePredictionMode,
    KinaseProfileMissingValueStrategy,
)
from phospy.contracts.requests import KinaseWorkflowRequest
from phospy.contracts.results import KinaseWorkflowResult
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.prediction.policies import (
    DEFAULT_PREDICTION_SAMPLING_POLICY,
    PredictionSamplingPolicy,
)
from phospy.science.references.models import ReferenceBundle
from phospy.science.sites.validation import (
    require_canonical_site_index,
    require_canonical_site_series,
)
from phospy.tables.datasets import PhosphoIntensityMatrix
from phospy.validation.common.dataframes import (
    require_canonical_string_column,
    require_columns,
    require_dataframe,
    require_non_empty_string_column,
    require_unique_index,
    require_unique_row_pairs,
)

if TYPE_CHECKING:
    from phospy.science.references.resolution import ReferenceResolverContract


@dataclass(frozen=True, slots=True)
class ResolvedKinaseWorkflowRequest:
    """Interpreter output for kinase workflow execution."""

    dataset: AnalysisReadyPhosphoDataset
    references: ReferenceBundle
    kinase_substrate_map: pd.DataFrame
    site_sequences: pd.DataFrame
    scoring_site_index: pd.Index
    activity_phospho_matrix: pd.DataFrame
    execution_config: ResolvedKinaseExecutionConfig
    site_sequence_merge_diagnostics: dict[str, object] = field(default_factory=dict)
    _kinase_substrate_reference: pd.DataFrame = field(
        init=False,
        repr=False,
        compare=False,
    )
    _site_sequence_reference: pd.DataFrame = field(
        init=False,
        repr=False,
        compare=False,
    )
    _activity_phospho_table: PhosphoIntensityMatrix = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        kinase_substrate_map = self._validate_kinase_substrate_map(
            self.kinase_substrate_map
        )
        site_sequences = self._validate_site_sequences(self.site_sequences)
        activity_phospho_table = self._validate_activity_phospho_table(
            self.activity_phospho_matrix
        )
        activity_phospho_matrix = activity_phospho_table.frame
        if not isinstance(self.scoring_site_index, pd.Index):
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.scoring_site_index_type; "
                "scoring_site_index must be a pandas Index; "
                "next_action=ensure interpreter passes a pandas Index for "
                "resolved scoring-site alignment"
            )
        if not activity_phospho_matrix.index.equals(self.scoring_site_index):
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.activity_site_alignment; "
                "activity_phospho_matrix.index must exactly match scoring_site_index; "
                "next_action=ensure interpreted activity phospho rows are aligned "
                "to the resolved scoring-site index"
            )
        if not kinase_substrate_map.equals(self.references.kinase_substrate_map):
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.kinase_substrate_reference_alignment; "
                "kinase_substrate_map must match "
                "references.kinase_substrate_map exactly; "
                "next_action=consume the normalized reference table from "
                "ReferenceBundle without workflow-local identifier cleanup"
            )
        if not self._has_compatible_site_sequence_reference_alignment(
            merged_site_sequences=site_sequences,
            reference_site_sequences=self.references.site_sequences,
            merge_diagnostics=self.site_sequence_merge_diagnostics,
        ):
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.site_sequence_reference_alignment; "
                "site_sequences must preserve every references.site_sequences row "
                "exactly unless conflict_policy='prefer_dataset' and each changed "
                "row is declared in site_sequence_merge_diagnostics.conflict_diagnostics; "
                "next_action=ensure execution-time site-sequence merge diagnostics "
                "match any reference-row overrides"
            )
        object.__setattr__(self, "kinase_substrate_map", kinase_substrate_map)
        object.__setattr__(self, "site_sequences", site_sequences)
        object.__setattr__(self, "activity_phospho_matrix", activity_phospho_matrix)
        object.__setattr__(self, "_kinase_substrate_reference", kinase_substrate_map)
        object.__setattr__(self, "_site_sequence_reference", site_sequences)
        object.__setattr__(self, "_activity_phospho_table", activity_phospho_table)

    @staticmethod
    def _validate_kinase_substrate_map(value: object) -> pd.DataFrame:
        try:
            frame = require_dataframe(
                value,
                field_name="kinase_request.kinase_substrate_map",
                allow_empty=False,
                error_type=WorkflowBoundaryError,
            )
            require_columns(
                frame,
                field_name="kinase_request.kinase_substrate_map",
                required_columns=("kinase", "substrate_site"),
                error_type=WorkflowBoundaryError,
            )
            require_non_empty_string_column(
                frame,
                field_name="kinase_request.kinase_substrate_map",
                column_name="kinase",
                error_type=WorkflowBoundaryError,
            )
            require_canonical_string_column(
                frame,
                field_name="kinase_request.kinase_substrate_map",
                column_name="kinase",
                error_type=WorkflowBoundaryError,
            )
            substrate_site = frame.loc[:, "substrate_site"]
            require_canonical_site_series(
                substrate_site,
                field_name="kinase_request.kinase_substrate_map.substrate_site",
                error_type=WorkflowBoundaryError,
            )
            require_unique_row_pairs(
                frame,
                field_name="kinase_request.kinase_substrate_map",
                column_names=("kinase", "substrate_site"),
                error_type=WorkflowBoundaryError,
            )
            return frame
        except WorkflowBoundaryError as exc:
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.kinase_substrate_map_schema; "
                f"{exc}; "
                "next_action=provide a prepared, already-normalized "
                "kinase_substrate_map from ReferenceBundle"
            ) from exc

    @staticmethod
    def _validate_site_sequences(value: object) -> pd.DataFrame:
        try:
            frame = require_dataframe(
                value,
                field_name="kinase_request.site_sequences",
                allow_empty=False,
                error_type=WorkflowBoundaryError,
            )
            require_columns(
                frame,
                field_name="kinase_request.site_sequences",
                required_columns=("site_sequence",),
                error_type=WorkflowBoundaryError,
            )
            require_non_empty_string_column(
                frame,
                field_name="kinase_request.site_sequences",
                column_name="site_sequence",
                error_type=WorkflowBoundaryError,
            )
            require_canonical_string_column(
                frame,
                field_name="kinase_request.site_sequences",
                column_name="site_sequence",
                error_type=WorkflowBoundaryError,
            )
            require_unique_index(
                frame,
                field_name="kinase_request.site_sequences",
                error_type=WorkflowBoundaryError,
            )
            require_canonical_site_index(
                frame.index,
                field_name="kinase_request.site_sequences.index",
                error_type=WorkflowBoundaryError,
            )
            return frame
        except WorkflowBoundaryError as exc:
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.site_sequence_schema; "
                f"{exc}; "
                "next_action=provide prepared, already-normalized "
                "site_sequences from ReferenceBundle or interpreter merge output"
            ) from exc

    @staticmethod
    def _validate_activity_phospho_table(value: object) -> PhosphoIntensityMatrix:
        frame = require_dataframe(
            value,
            field_name="kinase_request.activity_phospho_matrix",
            allow_empty=False,
            error_type=WorkflowBoundaryError,
        )
        return PhosphoIntensityMatrix(
            frame=frame,
            allow_missing=True,
            _assume_owned=True,
        )

    @property
    def kinase_substrate_reference(self) -> pd.DataFrame:
        return self._kinase_substrate_reference

    @property
    def site_sequence_reference(self) -> pd.DataFrame:
        return self._site_sequence_reference

    @property
    def activity_phospho_table(self) -> PhosphoIntensityMatrix:
        return self._activity_phospho_table

    @staticmethod
    def _has_compatible_site_sequence_reference_alignment(
        *,
        merged_site_sequences: pd.DataFrame,
        reference_site_sequences: pd.DataFrame,
        merge_diagnostics: Mapping[str, object] | None = None,
    ) -> bool:
        if not {"site_sequence"}.issubset(set(merged_site_sequences.columns)):
            return False
        if not {"site_sequence"}.issubset(set(reference_site_sequences.columns)):
            return False
        if merged_site_sequences.index.has_duplicates:
            return False
        conflict_policy = None
        prefer_dataset_conflicts: dict[str, Mapping[str, object]] = {}
        if isinstance(merge_diagnostics, Mapping):
            policy_value = merge_diagnostics.get("conflict_policy")
            if isinstance(policy_value, str):
                conflict_policy = policy_value
            conflict_rows = merge_diagnostics.get("conflict_diagnostics")
            if isinstance(conflict_rows, list):
                for row in conflict_rows:
                    if not isinstance(row, Mapping):
                        continue
                    site_value = row.get("site_id")
                    if not isinstance(site_value, str):
                        continue
                    prefer_dataset_conflicts[site_value] = row
        for site_id, sequence in reference_site_sequences.loc[
            :, "site_sequence"
        ].items():
            if site_id not in merged_site_sequences.index:
                return False
            merged_value = merged_site_sequences.at[site_id, "site_sequence"]
            merged_sequence = str(merged_value)
            reference_sequence = str(sequence)
            if merged_sequence == reference_sequence:
                continue
            if conflict_policy != "prefer_dataset":
                return False
            conflict_row = prefer_dataset_conflicts.get(str(site_id))
            if conflict_row is None:
                return False
            selected_source = conflict_row.get("selected_sequence_source")
            if str(selected_source) != "dataset":
                return False
            selected_sequence = conflict_row.get("selected_sequence")
            if str(selected_sequence) != merged_sequence:
                return False
            conflict_reference = conflict_row.get("reference_sequence")
            if str(conflict_reference) != reference_sequence:
                return False
        return True


@dataclass(frozen=True, slots=True)
class ResolvedKinaseActivityExecutionConfig:
    """Execution-ready kinase activity-stage config."""

    method: KinaseActivityMethod
    threshold: float
    min_substrates: int
    top_n_substrates: int
    ksea_min_substrates: int
    ksea_evidence_threshold: float
    ksea_p_value_method: KinaseActivityPValueMethod
    ksea_adjust_p_values: bool


@dataclass(frozen=True, slots=True)
class ResolvedKinaseExecutionConfig:
    """Execution-ready kinase workflow config resolved by the interpreter."""

    scoring_min_substrates: int
    include_diagnostic_scoring_tables: bool
    profile_missing_value_strategy: KinaseProfileMissingValueStrategy
    prediction_top_k: int
    prediction_deterministic_max_selected_kinases: int
    prediction_adaptive_ensemble_runs: int
    prediction_mode: KinasePredictionMode
    prediction_adaptive_policy: KinaseAdaptivePolicy
    prediction_n_iterations: int
    prediction_random_state: int | None
    prediction_sampling_policy: PredictionSamplingPolicy = (
        DEFAULT_PREDICTION_SAMPLING_POLICY
    )
    activity: ResolvedKinaseActivityExecutionConfig | None = None


class KinaseWorkflowValidatorContract(Protocol):
    """Internal contract for kinase workflow request validation."""

    def run(self, request: KinaseWorkflowRequest) -> KinaseWorkflowRequest: ...


class KinaseWorkflowInterpreterContract(Protocol):
    """Internal contract for kinase workflow request interpretation."""

    _reference_resolver: ReferenceResolverContract

    def run(self, request: KinaseWorkflowRequest) -> ResolvedKinaseWorkflowRequest: ...


class KinaseWorkflowExecutorContract(Protocol):
    """Internal contract for kinase workflow execution."""

    def run(self, request: ResolvedKinaseWorkflowRequest) -> KinaseWorkflowResult: ...


__all__ = [
    "KinaseWorkflowExecutorContract",
    "KinaseWorkflowInterpreterContract",
    "KinaseWorkflowValidatorContract",
    "ResolvedKinaseActivityExecutionConfig",
    "ResolvedKinaseExecutionConfig",
    "ResolvedKinaseWorkflowRequest",
]
