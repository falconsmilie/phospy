"""Kinase workflow internal contracts."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Protocol

import pandas as pd

from phospy.contracts.configs import (
    KINASE_RELIABILITY_PROFILE_EXPLORATORY,
    KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET,
    KinaseActivityMethod,
    KinaseActivityPValueMethod,
    KinaseActivitySsgseaRankingDirection,
    KinaseAdaptivePolicy,
    KinaseAttritionPolicy,
    KinasePredictionMode,
    KinaseProfileMissingValueStrategy,
    KinaseReliabilityProfile,
    KinaseScoringMode,
    LocalisationRequirement,
    ProfileSelfInclusionPolicy,
    ReferenceContextCompatibilityPolicy,
    normalize_kinase_scoring_mode,
)
from phospy.contracts.requests import KinaseWorkflowRequest
from phospy.contracts.results import KinaseWorkflowResult
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.provenance.models import RowAttritionRecord
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.prediction.policies import (
    DEFAULT_PREDICTION_SAMPLING_POLICY,
    PredictionSamplingPolicy,
)
from phospy.science.quantitative_method_contracts import (
    ResolvedMethodQuantitativeInputContract,
)
from phospy.science.references.models import ReferenceBundle
from phospy.science.sites.validation import (
    require_site_key_index,
    require_site_key_series,
)
from phospy.validation.common.config_values import coerce_policy_enum
from phospy.validation.common.dataframes import (
    require_canonical_string_column,
    require_columns,
    require_dataframe,
    require_finite_numeric_dataframe,
    require_non_empty_string_column,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
    require_unique_row_pairs,
)
from phospy.validation.datasets.protein_scoped_site_identity import (
    enforce_display_id_column,
    enforce_site_key_column_matches_index,
)
from phospy.validation.datasets.site_metadata import (
    enforce_site_sequence_context_contract,
)
from phospy.workflows.kinase.attrition_metrics import (
    KinaseAttritionMetrics,
    KinaseAttritionPolicyViolation,
    build_kinase_attrition_metrics,
)
from phospy.workflows.kinase.scoring_mode_contracts import (
    kinase_scoring_mode_input_contract,
)
from phospy.workflows.kinase.sequence_contracts import (
    dataset_sequence_source_label,
    kinase_sequence_context_contract,
)
from phospy.workflows.kinase.site_sequence_policy import (
    resolve_site_sequence_conflict_policy,
)

if TYPE_CHECKING:
    from phospy.science.references.kinase_library import KinaseLibraryResource
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
    kinase_library_resource: KinaseLibraryResource | None = None
    attrition_metrics: KinaseAttritionMetrics | None = None
    attrition_policy_violations: tuple[KinaseAttritionPolicyViolation, ...] = ()
    row_attrition_records: tuple[RowAttritionRecord, ...] = ()
    site_identity_map: pd.DataFrame | None = None
    site_sequence_merge_diagnostics: dict[str, object] = field(default_factory=dict)
    reference_resolution_details: dict[str, object] = field(default_factory=dict)
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
    _activity_phospho_table: pd.DataFrame = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        kinase_substrate_map = self._validate_kinase_substrate_map(
            self.kinase_substrate_map,
            scoring_mode=self.execution_config.scoring_mode,
        )
        site_sequences = self._validate_site_sequences(self.site_sequences)
        activity_phospho_matrix = self._validate_activity_phospho_table(
            self.activity_phospho_matrix
        )
        if not isinstance(self.scoring_site_index, pd.Index):
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.scoring_site_index_type; "
                "scoring_site_index must be a pandas Index; "
                "next_action=ensure interpreter passes a pandas Index for "
                "resolved scoring-site alignment"
            )
        _require_site_index_identity(
            self.scoring_site_index,
            field_name="kinase_request.scoring_site_index",
            error_type=WorkflowBoundaryError,
        )
        site_identity_map = self._validate_site_identity_map(
            self.site_identity_map,
            scoring_site_index=self.scoring_site_index,
        )
        if not activity_phospho_matrix.index.equals(self.scoring_site_index):
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.activity_site_alignment; "
                "activity_phospho_matrix.index must exactly match scoring_site_index; "
                "next_action=ensure interpreted activity phospho rows are aligned "
                "to the resolved scoring-site index"
            )
        if not self.scoring_site_index.equals(site_sequences.index):
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.scoring_site_sequence_alignment; "
                "scoring_site_index must exactly match site_sequences.index; "
                "next_action=ensure sequence support is projected onto the same "
                "site_key rows selected for scoring"
            )
        if not self.scoring_site_index.equals(site_identity_map.index):
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.scoring_site_identity_alignment; "
                "scoring_site_index must exactly match site_identity_map.index; "
                "next_action=ensure scoring site selection preserves "
                "site_key/display_id identity mapping"
            )
        self._validate_sequence_context_contract(
            site_sequences=site_sequences,
            site_identity_map=site_identity_map,
            scoring_mode=self.execution_config.scoring_mode,
            kinase_library_resource=self.kinase_library_resource,
            references=self.references,
            dataset=self.dataset,
            merge_diagnostics=self.site_sequence_merge_diagnostics,
        )
        if not set(kinase_substrate_map.loc[:, "substrate_site"].astype(str)).issubset(
            set(self.scoring_site_index.astype(str))
        ):
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.kinase_substrate_reference_alignment; "
                "kinase_substrate_map.substrate_site must be projected to scoring "
                "site_key rows; next_action=map reference display IDs onto "
                "dataset site_key before scoring"
            )
        attrition_metrics = self.attrition_metrics
        if attrition_metrics is None:
            attrition_metrics = build_kinase_attrition_metrics(
                dataset_site_index=self.dataset.phospho.index,
                reference_site_index=kinase_substrate_map.loc[:, "substrate_site"],
                sequence_supported_site_index=self.scoring_site_index,
            )
        elif not isinstance(attrition_metrics, KinaseAttritionMetrics):
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.attrition_metrics_type; "
                "attrition_metrics must be KinaseAttritionMetrics; "
                "next_action=ensure interpreter attaches structured attrition "
                "metrics before executor scoring"
            )
        attrition_policy_violations = tuple(self.attrition_policy_violations)
        for violation in attrition_policy_violations:
            if isinstance(violation, KinaseAttritionPolicyViolation):
                continue
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.attrition_policy_violations_type; "
                "attrition_policy_violations must contain "
                "KinaseAttritionPolicyViolation values; "
                "next_action=ensure interpreter attaches structured attrition "
                "policy violations before executor scoring"
            )
        row_attrition_records = tuple(self.row_attrition_records)
        for record in row_attrition_records:
            if isinstance(record, RowAttritionRecord):
                continue
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.row_attrition_records_type; "
                "row_attrition_records must contain RowAttritionRecord values; "
                "next_action=ensure interpreter attaches structured causal "
                "row-attrition records before executor scoring"
            )
        if not self._has_compatible_site_sequence_reference_alignment(
            merged_site_sequences=site_sequences,
            site_identity_map=site_identity_map,
            reference_site_sequences=self.references.site_sequences,
            merge_diagnostics=self.site_sequence_merge_diagnostics,
        ):
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.site_sequence_reference_alignment; "
                "site_sequences must project reference display rows onto scoring "
                "site_key rows; sequence overrides are only allowed when "
                "conflict_policy='prefer_dataset' and declared in "
                "site_sequence_merge_diagnostics.conflict_diagnostics; "
                "next_action=ensure execution-time site-sequence projection "
                "diagnostics match reference-row overrides"
            )
        object.__setattr__(self, "kinase_substrate_map", kinase_substrate_map)
        object.__setattr__(self, "site_sequences", site_sequences)
        object.__setattr__(self, "site_identity_map", site_identity_map)
        object.__setattr__(self, "activity_phospho_matrix", activity_phospho_matrix)
        object.__setattr__(self, "attrition_metrics", attrition_metrics)
        object.__setattr__(
            self,
            "attrition_policy_violations",
            attrition_policy_violations,
        )
        object.__setattr__(self, "row_attrition_records", row_attrition_records)
        object.__setattr__(self, "_kinase_substrate_reference", kinase_substrate_map)
        object.__setattr__(self, "_site_sequence_reference", site_sequences)
        object.__setattr__(self, "_activity_phospho_table", activity_phospho_matrix)

    @staticmethod
    def _validate_kinase_substrate_map(
        value: object,
        *,
        scoring_mode: str,
    ) -> pd.DataFrame:
        mode_contract = kinase_scoring_mode_input_contract(scoring_mode)
        allow_empty = not mode_contract.requires_substrate_reference_overlap
        try:
            frame = require_dataframe(
                value,
                field_name="kinase_request.kinase_substrate_map",
                allow_empty=allow_empty,
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
            _require_site_series_identity(
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
                "next_action=provide a projected kinase_substrate_map "
                "indexed by dataset site_key rows"
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
                required_columns=("site_sequence", "display_id"),
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
            require_non_empty_string_column(
                frame,
                field_name="kinase_request.site_sequences",
                column_name="display_id",
                error_type=WorkflowBoundaryError,
            )
            require_unique_index(
                frame,
                field_name="kinase_request.site_sequences",
                error_type=WorkflowBoundaryError,
            )
            _require_site_index_identity(
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
                "next_action=provide site_sequences projected to site_key with "
                "display_id support metadata"
            ) from exc

    @staticmethod
    def _validate_activity_phospho_table(value: object) -> pd.DataFrame:
        frame = require_dataframe(
            value,
            field_name="kinase_request.activity_phospho_matrix",
            allow_empty=False,
            error_type=WorkflowBoundaryError,
        )
        require_unique_columns(
            frame,
            field_name="kinase_request.activity_phospho_matrix",
            error_type=WorkflowBoundaryError,
        )
        require_numeric_dataframe(
            frame,
            field_name="kinase_request.activity_phospho_matrix",
            error_type=WorkflowBoundaryError,
        )
        require_finite_numeric_dataframe(
            frame,
            field_name="kinase_request.activity_phospho_matrix",
            error_type=WorkflowBoundaryError,
            allow_missing=True,
        )
        require_unique_index(
            frame,
            field_name="kinase_request.activity_phospho_matrix",
            error_type=WorkflowBoundaryError,
        )
        _require_site_index_identity(
            frame.index,
            field_name="kinase_request.activity_phospho_matrix.index",
            error_type=WorkflowBoundaryError,
        )
        return frame

    @staticmethod
    def _validate_site_identity_map(
        value: object | None, *, scoring_site_index: pd.Index
    ) -> pd.DataFrame:
        try:
            if value is None:
                raise WorkflowBoundaryError(
                    "kinase_request.site_identity_map is required"
                )
            frame = require_dataframe(
                value,
                field_name="kinase_request.site_identity_map",
                allow_empty=False,
                error_type=WorkflowBoundaryError,
            )
            require_columns(
                frame,
                field_name="kinase_request.site_identity_map",
                required_columns=("site_key", "display_id"),
                error_type=WorkflowBoundaryError,
            )
            require_non_empty_string_column(
                frame,
                field_name="kinase_request.site_identity_map",
                column_name="site_key",
                error_type=WorkflowBoundaryError,
            )
            require_non_empty_string_column(
                frame,
                field_name="kinase_request.site_identity_map",
                column_name="display_id",
                error_type=WorkflowBoundaryError,
            )
            require_unique_index(
                frame,
                field_name="kinase_request.site_identity_map",
                error_type=WorkflowBoundaryError,
            )
            _require_site_index_identity(
                frame.index,
                field_name="kinase_request.site_identity_map.index",
                error_type=WorkflowBoundaryError,
            )
            display_ids = enforce_display_id_column(
                site_metadata=frame,
                field_name="kinase_request.site_identity_map",
                error_type=WorkflowBoundaryError,
                column_name="display_id",
            )
            site_keys = enforce_site_key_column_matches_index(
                site_metadata=frame,
                field_name="kinase_request.site_identity_map",
                error_type=WorkflowBoundaryError,
                site_key_column="site_key",
            )
            if list(site_keys.astype(str)) != list(scoring_site_index.astype(str)):
                raise WorkflowBoundaryError(
                    "kinase_request.site_identity_map.site_key must exactly match "
                    "kinase_request.scoring_site_index"
                )
            _ = display_ids
            return frame
        except WorkflowBoundaryError as exc:
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.contracts.site_identity_map_schema; "
                f"{exc}; "
                "next_action=provide one row per scoring site_key with mapped "
                "display_id values"
            ) from exc

    @property
    def kinase_substrate_reference(self) -> pd.DataFrame:
        return self._kinase_substrate_reference

    @property
    def site_sequence_reference(self) -> pd.DataFrame:
        return self._site_sequence_reference

    @property
    def activity_phospho_table(self) -> pd.DataFrame:
        return self._activity_phospho_table

    @staticmethod
    def _has_compatible_site_sequence_reference_alignment(
        *,
        merged_site_sequences: pd.DataFrame,
        site_identity_map: pd.DataFrame,
        reference_site_sequences: pd.DataFrame,
        merge_diagnostics: Mapping[str, object] | None = None,
    ) -> bool:
        if not {"site_sequence", "display_id"}.issubset(
            set(merged_site_sequences.columns)
        ):
            return False
        if not {"site_key", "display_id"}.issubset(set(site_identity_map.columns)):
            return False
        if not merged_site_sequences.index.equals(site_identity_map.index):
            return False
        if not {"site_sequence"}.issubset(set(reference_site_sequences.columns)):
            return False
        if merged_site_sequences.index.has_duplicates:
            return False
        conflict_policy = None
        prefer_dataset_conflicts: dict[str, Mapping[str, object]] = {}
        if isinstance(merge_diagnostics, Mapping):
            policy_value = merge_diagnostics.get("conflict_policy")
            try:
                conflict_policy = resolve_site_sequence_conflict_policy(
                    policy_value,
                    field_name="site_sequence_merge_diagnostics.conflict_policy",
                    error_type=WorkflowBoundaryError,
                )
            except WorkflowBoundaryError:
                conflict_policy = None
            conflict_rows = merge_diagnostics.get("conflict_diagnostics")
            if isinstance(conflict_rows, list):
                for row in conflict_rows:
                    if not isinstance(row, Mapping):
                        continue
                    site_value = row.get("site_key")
                    if not isinstance(site_value, str):
                        site_value = row.get("site_id")
                    if not isinstance(site_value, str):
                        continue
                    prefer_dataset_conflicts[site_value] = row
        for site_key, row in merged_site_sequences.loc[
            :, ["site_sequence", "display_id"]
        ].iterrows():
            display_id = str(row["display_id"])
            if display_id not in reference_site_sequences.index:
                continue
            merged_sequence = str(row["site_sequence"])
            reference_sequence = str(
                reference_site_sequences.at[display_id, "site_sequence"]
            )
            if merged_sequence == reference_sequence:
                continue
            if conflict_policy != KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET:
                return False
            conflict_row = prefer_dataset_conflicts.get(str(site_key))
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

    @staticmethod
    def _validate_sequence_context_contract(
        *,
        site_sequences: pd.DataFrame,
        site_identity_map: pd.DataFrame,
        scoring_mode: str,
        kinase_library_resource: KinaseLibraryResource | None,
        references: ReferenceBundle,
        dataset: AnalysisReadyPhosphoDataset,
        merge_diagnostics: Mapping[str, object] | None,
    ) -> None:
        contract = kinase_sequence_context_contract(
            scoring_mode=scoring_mode,
            kinase_library_resource=kinase_library_resource,
        )
        if contract is None:
            return
        context_frame = pd.DataFrame(
            {
                "site_sequence": site_sequences.loc[:, "site_sequence"].tolist(),
                "display_id": site_sequences.loc[:, "display_id"].astype(str).tolist(),
                "site": [
                    _site_token_from_display_id(display_id)
                    for display_id in site_sequences.loc[:, "display_id"]
                    .astype(str)
                    .tolist()
                ],
            },
            index=site_sequences.index.copy(),
        )
        source_by_site = _resolved_sequence_source_by_site(
            site_sequences=site_sequences,
            site_identity_map=site_identity_map,
            references=references,
            dataset=dataset,
            merge_diagnostics=merge_diagnostics,
        )
        enforce_site_sequence_context_contract(
            site_metadata=context_frame,
            field_name="kinase workflow interpreted site_sequences",
            workflow_name="kinase workflow interpreted request",
            scoring_mode=scoring_mode,
            contract=contract,
            error_type=WorkflowBoundaryError,
            sequence_source_by_site=source_by_site,
            allow_unknown_site_residue=False,
        )


def _require_site_index_identity(
    index: pd.Index,
    *,
    field_name: str,
    error_type: type[WorkflowBoundaryError],
) -> None:
    require_site_key_index(
        index,
        field_name=field_name,
        error_type=error_type,
    )


def _require_site_series_identity(
    series: pd.Series,
    *,
    field_name: str,
    error_type: type[WorkflowBoundaryError],
) -> None:
    require_site_key_series(
        series,
        field_name=field_name,
        error_type=error_type,
    )


def _site_token_from_display_id(display_id: object) -> str:
    parts = str(display_id).split(";")
    if len(parts) >= 2:
        site_token = parts[1].strip()
        if site_token:
            return site_token
    return str(display_id).strip()


def _resolved_sequence_source_by_site(
    *,
    site_sequences: pd.DataFrame,
    site_identity_map: pd.DataFrame,
    references: ReferenceBundle,
    dataset: AnalysisReadyPhosphoDataset,
    merge_diagnostics: Mapping[str, object] | None,
) -> dict[str, str]:
    dataset_source = dataset_sequence_source_label(dataset)
    conflict_source_by_site = _selected_conflict_source_by_site(merge_diagnostics)
    reference_display_ids = set(references.site_sequences.index.astype(str).tolist())
    source_by_site: dict[str, str] = {}
    for site_id, row in site_sequences.loc[:, ["display_id"]].iterrows():
        site_key = str(site_id)
        selected_conflict_source = conflict_source_by_site.get(site_key)
        if selected_conflict_source == "dataset":
            source_by_site[site_key] = dataset_source or "unknown"
            continue
        if selected_conflict_source == "reference":
            source_by_site[site_key] = "reference"
            continue
        display_id = str(row["display_id"])
        if display_id in reference_display_ids:
            source_by_site[site_key] = "reference"
            continue
        if site_key in site_identity_map.index:
            identity_display_id = str(site_identity_map.at[site_key, "display_id"])
            if identity_display_id in reference_display_ids:
                source_by_site[site_key] = "reference"
                continue
        source_by_site[site_key] = dataset_source or "unknown"
    return source_by_site


def _selected_conflict_source_by_site(
    merge_diagnostics: Mapping[str, object] | None,
) -> dict[str, str]:
    if not isinstance(merge_diagnostics, Mapping):
        return {}
    rows = merge_diagnostics.get("conflict_diagnostics")
    if not isinstance(rows, list):
        return {}
    selected_by_site: dict[str, str] = {}
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        site_key = row.get("site_key")
        if not isinstance(site_key, str):
            site_key = row.get("site_id")
        if not isinstance(site_key, str) or site_key.strip() == "":
            continue
        selected_source = row.get("selected_sequence_source")
        if not isinstance(selected_source, str):
            continue
        selected_source_text = selected_source.strip().lower()
        if selected_source_text not in {"dataset", "reference"}:
            continue
        selected_by_site[site_key] = selected_source_text
    return selected_by_site


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
    ssgsea_min_substrates: int = 5
    ssgsea_ranking_direction: KinaseActivitySsgseaRankingDirection = "descending"
    ssgsea_permutations: int = 0
    ssgsea_random_seed: int | None = 0
    ssgsea_adjust_p_values: bool = True
    method_input_contract: ResolvedMethodQuantitativeInputContract | None = None


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
    profile_self_inclusion_policy: ProfileSelfInclusionPolicy = (
        ProfileSelfInclusionPolicy.ALLOW
    )
    attrition_policy: KinaseAttritionPolicy = field(
        default_factory=KinaseAttritionPolicy
    )
    prediction_sampling_policy: PredictionSamplingPolicy = (
        DEFAULT_PREDICTION_SAMPLING_POLICY
    )
    activity: ResolvedKinaseActivityExecutionConfig | None = None
    scoring_mode: KinaseScoringMode = KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED
    include_substrate_contributions: bool = False
    requested_reliability_profile: KinaseReliabilityProfile | None = None
    effective_reliability_profile: KinaseReliabilityProfile = (
        KINASE_RELIABILITY_PROFILE_EXPLORATORY
    )
    localisation_requirement: LocalisationRequirement = field(
        default_factory=LocalisationRequirement
    )
    reference_context_compatibility_policy: ReferenceContextCompatibilityPolicy = (
        ReferenceContextCompatibilityPolicy.REQUIRE_KNOWN_MATCH
    )
    scoring_method_input_contract: ResolvedMethodQuantitativeInputContract | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "scoring_mode",
            normalize_kinase_scoring_mode(self.scoring_mode),
        )
        object.__setattr__(
            self,
            "profile_self_inclusion_policy",
            coerce_policy_enum(
                ProfileSelfInclusionPolicy,
                self.profile_self_inclusion_policy,
                field_name="kinase.execution_config.profile_self_inclusion_policy",
                error_type=WorkflowBoundaryError,
            ),
        )
        object.__setattr__(
            self,
            "reference_context_compatibility_policy",
            coerce_policy_enum(
                ReferenceContextCompatibilityPolicy,
                self.reference_context_compatibility_policy,
                field_name=(
                    "kinase.execution_config.reference_context_compatibility_policy"
                ),
                error_type=WorkflowBoundaryError,
            ),
        )
        if self.requested_reliability_profile is not None:
            object.__setattr__(
                self,
                "requested_reliability_profile",
                coerce_policy_enum(
                    KinaseReliabilityProfile,
                    self.requested_reliability_profile,
                    field_name=(
                        "kinase.execution_config.requested_reliability_profile"
                    ),
                    error_type=WorkflowBoundaryError,
                ),
            )
        object.__setattr__(
            self,
            "effective_reliability_profile",
            coerce_policy_enum(
                KinaseReliabilityProfile,
                self.effective_reliability_profile,
                field_name="kinase.execution_config.effective_reliability_profile",
                error_type=WorkflowBoundaryError,
            ),
        )


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
