"""Internal interpreter for kinase workflow requests."""

from __future__ import annotations

from typing import NoReturn

import pandas as pd

from phospy.contracts.configs import ReferenceContextCompatibilityPolicy
from phospy.contracts.requests import KinaseWorkflowRequest
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.prediction.policies import resolve_prediction_sampling_policy
from phospy.science.references.kinase_library import KinaseLibraryResource
from phospy.science.references.models import ReferenceBundle, ReferencePreset
from phospy.science.references.resolution import (
    BundledReferenceProvider,
    ReferenceResolver,
    ReferenceResolverContract,
)
from phospy.validation.identity_contracts import (
    validate_reference_context_compatibility,
)
from phospy.validation.references.compatibility import ReferenceCompatibilityValidator
from phospy.workflows._pandas_typing import (
    dataframe_column,
    dataframe_copy,
    dataframe_loc,
    dataframe_reindex,
    dataframe_reset_index,
    index_as_strings,
    series_as_strings,
)
from phospy.workflows._row_attrition import make_row_attrition_record
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseActivityExecutionConfig,
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.reference_projection import KinaseReferenceProjector
from phospy.workflows.kinase.resolved_validator import (
    ResolvedKinaseEligibilityValidator,
    ResolvedKinaseInputs,
)
from phospy.workflows.kinase.scoring_mode_contracts import (
    kinase_scoring_mode_input_contract,
)
from phospy.workflows.kinase.site_sequence_policy import (
    resolve_site_sequence_conflict_policy,
)
from phospy.workflows.kinase.site_sequence_support import (
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
    KinaseSiteSequenceSupportBuilder,
)


class KinaseWorkflowInterpreter:
    """Resolve workflow request defaults and references for execution."""

    def __init__(
        self,
        *,
        reference_resolver: ReferenceResolverContract | None = None,
        reference_projector: KinaseReferenceProjector | None = None,
        site_sequence_support_builder: KinaseSiteSequenceSupportBuilder | None = None,
        resolved_validator: ResolvedKinaseEligibilityValidator | None = None,
    ) -> None:
        self._reference_resolver = reference_resolver or ReferenceResolver(
            provider=BundledReferenceProvider(),
            compatibility_validator=ReferenceCompatibilityValidator(),
        )
        self._reference_projector = reference_projector or KinaseReferenceProjector()
        self._site_sequence_support_builder = (
            site_sequence_support_builder or KinaseSiteSequenceSupportBuilder()
        )
        self._resolved_validator = (
            resolved_validator or ResolvedKinaseEligibilityValidator()
        )

    def run(self, request: KinaseWorkflowRequest) -> ResolvedKinaseWorkflowRequest:
        dataset_view = DatasetInternalView(request.dataset)
        dataset_phospho = dataset_view.phospho
        dataset_site_metadata = dataset_view.site_metadata
        site_identity_map = self._build_site_identity_map(
            dataset=dataset_phospho,
            site_metadata=dataset_site_metadata,
        )
        site_sequence_conflict_policy = resolve_site_sequence_conflict_policy(
            request.site_sequence_conflict_policy,
            field_name="site_sequence_conflict_policy",
            error_type=WorkflowBoundaryError,
        )
        references = self._reference_resolver.run(
            request.references,
            dataset_organism=request.dataset.organism,
        )
        validate_reference_context_compatibility(
            request.dataset.reference_context,
            None
            if references.provenance is None
            else references.provenance.reference_context,
            operation="kinase workflow resolved dataset/reference bundle",
            allow_unknown=(
                request.scoring_config.reference_context_compatibility_policy
                is ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
            error_type=WorkflowBoundaryError,
        )
        kinase_library_resource = self._resolve_kinase_library_resource(
            request=request,
        )
        reference_site_count = len(
            frozenset(
                series_as_strings(
                    dataframe_column(
                        references.kinase_substrate_map,
                        "substrate_site",
                    )
                )
            )
        )
        projection_result = self._reference_projector.run(
            reference_kinase_substrate_map=references.kinase_substrate_map,
            site_identity_map=site_identity_map,
            ambiguity_policy=request.reference_display_ambiguity_policy,
        )
        kinase_substrate_map = projection_result.kinase_substrate_map
        display_reference_matching = (
            projection_result.display_reference_matching_payload()
        )
        merge_result = self._site_sequence_support_builder.run(
            dataset=dataset_phospho,
            site_metadata=dataset_site_metadata,
            reference_site_sequences=references.site_sequences,
            conflict_policy=site_sequence_conflict_policy,
        )
        if (
            merge_result.dataset_reference_conflict_count > 0
            and site_sequence_conflict_policy
            == KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR
        ):
            self._raise_boundary_error(
                seam="kinase.interpreter.site_sequence_conflict",
                next_action=(
                    "fix dataset site_sequence values for conflicting sites or use "
                    "site_sequence_conflict_policy='prefer_reference' or "
                    "'prefer_dataset' on KinaseWorkflowRequest"
                ),
                conflict_policy=site_sequence_conflict_policy,
                dataset_reference_conflict_count=int(
                    merge_result.dataset_reference_conflict_count
                ),
                conflict_diagnostics=[
                    item.to_payload() for item in merge_result.conflicts
                ],
            )
        site_sequences = merge_result.site_sequences
        scoring_site_index = self._resolve_scoring_site_index(
            dataset=dataset_phospho,
            site_sequences=site_sequences,
        )
        row_attrition_records = tuple(
            record
            for record in (
                make_row_attrition_record(
                    workflow="kinase",
                    stage="kinase_sequence_context",
                    reason="sites_missing_valid_centered_sequence",
                    input_site_ids=dataset_phospho.index,
                    output_site_ids=scoring_site_index,
                ),
            )
            if record is not None
        )
        activity_phospho_matrix = dataframe_loc(
            dataset_phospho,
            rows=scoring_site_index,
        )
        execution_config = self._resolve_execution_config(request)
        resolved_inputs = ResolvedKinaseInputs(
            dataset=request.dataset,
            dataset_phospho=dataset_phospho,
            references=references,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            site_identity_map=site_identity_map,
            scoring_site_index=scoring_site_index,
            activity_phospho_matrix=activity_phospho_matrix,
            execution_config=execution_config,
            reference_site_count=reference_site_count,
            kinase_library_resource=kinase_library_resource,
        )
        self._resolved_validator.run(resolved_inputs)
        scoring_site_keys = set(index_as_strings(scoring_site_index))
        site_sequences = dataframe_reindex(site_sequences, scoring_site_index)
        site_identity_map = dataframe_reindex(site_identity_map, scoring_site_index)
        substrate_sites = series_as_strings(
            dataframe_column(kinase_substrate_map, "substrate_site")
        )
        substrate_site_mask = [site in scoring_site_keys for site in substrate_sites]
        kinase_substrate_map = dataframe_copy(
            dataframe_reset_index(
                dataframe_loc(kinase_substrate_map, rows=substrate_site_mask),
                drop=True,
            ),
            deep=True,
        )
        return ResolvedKinaseWorkflowRequest(
            dataset=request.dataset,
            references=references,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            site_identity_map=site_identity_map,
            scoring_site_index=scoring_site_index,
            activity_phospho_matrix=activity_phospho_matrix,
            execution_config=execution_config,
            kinase_library_resource=kinase_library_resource,
            attrition_metrics=resolved_inputs.attrition_metrics,
            attrition_policy_violations=(resolved_inputs.attrition_policy_violations),
            row_attrition_records=row_attrition_records,
            site_sequence_merge_diagnostics={
                **merge_result.diagnostics_payload(),
                "reference_sequence_count": int(references.site_sequences.shape[0]),
                "execution_sequence_count": int(site_sequences.shape[0]),
                "reference_substrate_map_count": int(
                    references.kinase_substrate_map.shape[0]
                ),
                "execution_substrate_map_count": int(kinase_substrate_map.shape[0]),
                "display_reference_matching": display_reference_matching,
            },
            reference_resolution_details=self._build_reference_resolution_details(
                request=request,
                references=references,
            ),
        )

    @staticmethod
    def _build_site_identity_map(
        *,
        dataset: pd.DataFrame,
        site_metadata: pd.DataFrame,
    ) -> pd.DataFrame:
        missing = [
            column
            for column in ("site_key", "display_id")
            if column not in site_metadata.columns
        ]
        if missing:
            joined = ", ".join(missing)
            raise WorkflowBoundaryError(
                seam="kinase.interpreter.site_identity_map",
                next_action=(
                    "ensure the analysis-ready dataset metadata includes site_key "
                    "and display_id"
                ),
                details={"missing_columns": missing},
                message_prefix=(
                    "kinase workflow boundary validation failed: missing required "
                    f"site identity columns: {joined}"
                ),
            )
        aligned_metadata = dataframe_reindex(site_metadata, dataset.index)
        site_key_values = series_as_strings(
            dataframe_column(aligned_metadata, "site_key"),
            fill_missing="",
            strip=True,
        )
        display_id_values = series_as_strings(
            dataframe_column(aligned_metadata, "display_id"),
            fill_missing="",
            strip=True,
        )
        site_keys = index_as_strings(dataset.index)
        if site_key_values != site_keys:
            raise WorkflowBoundaryError(
                seam="kinase.interpreter.site_identity_map",
                next_action=(
                    "ensure dataset.site_metadata.site_key exactly matches "
                    "dataset.phospho.index"
                ),
                details={"dataset_site_count": int(dataset.index.size)},
                message_prefix="kinase workflow boundary validation failed",
            )
        empty_display_id_count = sum(
            1 for display_id in display_id_values if display_id == ""
        )
        if empty_display_id_count:
            raise WorkflowBoundaryError(
                seam="kinase.interpreter.site_identity_map",
                next_action=(
                    "ensure dataset.site_metadata.display_id is present for every "
                    "scoring site"
                ),
                details={"empty_display_id_count": empty_display_id_count},
                message_prefix="kinase workflow boundary validation failed",
            )
        return pd.DataFrame(
            {"site_key": site_keys, "display_id": display_id_values},
            index=pd.Index(site_keys, name="site_key"),
        )

    @staticmethod
    def _resolve_execution_config(
        request: KinaseWorkflowRequest,
    ) -> ResolvedKinaseExecutionConfig:
        prediction_sampling_policy = resolve_prediction_sampling_policy(
            request.prediction_config.adaptive_policy
        )
        activity = (
            None
            if request.activity_config is None or not request.activity_config.enabled
            else ResolvedKinaseActivityExecutionConfig(
                method=request.activity_config.method,
                threshold=float(request.activity_config.threshold),
                min_substrates=int(request.activity_config.min_substrates),
                top_n_substrates=int(request.activity_config.top_n_substrates),
                ksea_min_substrates=int(request.activity_config.ksea_min_substrates),
                ksea_evidence_threshold=float(
                    request.activity_config.threshold
                    if request.activity_config.ksea_evidence_threshold is None
                    else request.activity_config.ksea_evidence_threshold
                ),
                ksea_p_value_method=request.activity_config.ksea_p_value_method,
                ksea_adjust_p_values=bool(request.activity_config.ksea_adjust_p_values),
                ssgsea_min_substrates=int(
                    request.activity_config.ssgsea_min_substrates
                ),
                ssgsea_ranking_direction=(
                    request.activity_config.ssgsea_ranking_direction
                ),
                ssgsea_permutations=int(request.activity_config.ssgsea_permutations),
                ssgsea_random_seed=(
                    None
                    if request.activity_config.ssgsea_random_seed is None
                    else int(request.activity_config.ssgsea_random_seed)
                ),
                ssgsea_adjust_p_values=bool(
                    request.activity_config.ssgsea_adjust_p_values
                ),
            )
        )
        return ResolvedKinaseExecutionConfig(
            scoring_min_substrates=int(request.scoring_config.min_substrates),
            scoring_mode=request.scoring_config.scoring_mode,
            include_diagnostic_scoring_tables=bool(
                request.scoring_config.include_diagnostic_scoring_tables
            ),
            include_substrate_contributions=bool(
                request.scoring_config.include_substrate_contributions
            ),
            profile_missing_value_strategy=request.scoring_config.profile_missing_value_strategy,
            profile_self_inclusion_policy=(
                request.scoring_config.profile_self_inclusion_policy
            ),
            prediction_top_k=int(request.prediction_config.top_k),
            prediction_deterministic_max_selected_kinases=int(
                request.prediction_config.deterministic_max_selected_kinases
            ),
            prediction_adaptive_ensemble_runs=int(
                request.prediction_config.adaptive_ensemble_runs
            ),
            prediction_mode=request.prediction_config.mode,
            prediction_adaptive_policy=request.prediction_config.adaptive_policy,
            prediction_sampling_policy=prediction_sampling_policy,
            prediction_n_iterations=int(request.prediction_config.n_iterations),
            prediction_random_state=(
                None
                if request.prediction_config.random_state is None
                else int(request.prediction_config.random_state)
            ),
            attrition_policy=request.scoring_config.attrition_policy,
            activity=activity,
            requested_reliability_profile=(
                request.scoring_config.requested_reliability_profile
            ),
            effective_reliability_profile=(
                request.scoring_config.effective_reliability_profile
            ),
            localisation_requirement=request.scoring_config.localisation_requirement,
            reference_context_compatibility_policy=(
                request.scoring_config.reference_context_compatibility_policy
            ),
        )

    @staticmethod
    def _build_reference_resolution_details(
        *,
        request: KinaseWorkflowRequest,
        references: ReferenceBundle,
    ) -> dict[str, object]:
        reference_input = request.references
        reference_input_kind = "reference_bundle"
        reference_input_value = "explicit_bundle"
        if isinstance(reference_input, ReferencePreset):
            reference_input_kind = "reference_preset"
            reference_input_value = reference_input.value
        provenance = references.provenance
        return {
            "reference_input_kind": reference_input_kind,
            "reference_input_value": reference_input_value,
            "dataset_organism": (
                None
                if request.dataset.organism is None
                else request.dataset.organism.value
            ),
            "resolved_reference_organism": references.organism.value,
            "resolved_reference_source_type": (
                "unknown" if provenance is None else provenance.source_type
            ),
            "resolved_reference_bundle_id": (
                None if provenance is None else provenance.bundle_id
            ),
            "resolved_reference_source_name": (
                None if provenance is None else provenance.source_name
            ),
            "resolved_reference_source_version": (
                None if provenance is None else provenance.source_version
            ),
        }

    def _resolve_kinase_library_resource(
        self,
        *,
        request: KinaseWorkflowRequest,
    ) -> KinaseLibraryResource | None:
        scoring_mode = request.scoring_config.scoring_mode
        mode_contract = kinase_scoring_mode_input_contract(scoring_mode)
        if not mode_contract.requires_kinase_library_resource:
            return None
        resource = request.kinase_library_resource
        if not isinstance(resource, KinaseLibraryResource):
            self._raise_boundary_error(
                seam="kinase.interpreter.kinase_library_resource_type",
                next_action=(
                    "provide a KinaseLibraryResource via "
                    "KinaseWorkflowRequest.kinase_library_resource"
                ),
                scoring_mode=scoring_mode,
                observed_type=type(resource).__name__,
            )
        dataset_organism = request.dataset.organism
        if dataset_organism is not None and dataset_organism.value not in {
            str(value).lower() for value in resource.organisms
        }:
            self._raise_boundary_error(
                seam="kinase.interpreter.kinase_library_resource_organism",
                next_action=(
                    "provide a Kinase Library resource whose organisms include "
                    "the analysis-ready dataset organism"
                ),
                scoring_mode=scoring_mode,
                dataset_organism=dataset_organism.value,
                resource_organisms=tuple(resource.organisms),
            )
        if not resource.sequence_window.central_residue_required:
            self._raise_boundary_error(
                seam="kinase.interpreter.kinase_library_resource_window",
                next_action=(
                    "provide a Kinase Library resource with a central "
                    "phospho-residue window"
                ),
                scoring_mode=scoring_mode,
                sequence_window=resource.sequence_window.to_payload(),
            )
        return resource

    @staticmethod
    def _resolve_scoring_site_index(
        *,
        dataset: pd.DataFrame,
        site_sequences: pd.DataFrame,
    ) -> pd.Index:
        sequence_sites = set(site_sequences.index)
        scoring_sites = [
            site_id for site_id in dataset.index if site_id in sequence_sites
        ]
        return pd.Index(scoring_sites, name=dataset.index.name)

    @staticmethod
    def _raise_boundary_error(
        *,
        seam: str,
        next_action: str,
        **details: object,
    ) -> NoReturn:
        raise WorkflowBoundaryError(
            seam=seam,
            next_action=next_action,
            details=details,
            message_prefix="kinase workflow boundary validation failed",
        )
