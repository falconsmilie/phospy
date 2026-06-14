"""Internal interpreter for kinase workflow requests."""

from __future__ import annotations

from typing import NoReturn, TypedDict, cast

import pandas as pd

from phospy.contracts.configs import (
    KINASE_SCORING_MODES_REQUIRING_KINASE_LIBRARY,
)
from phospy.contracts.requests import KinaseWorkflowRequest
from phospy.errors.workflows import PhosPyWorkflowError, WorkflowBoundaryError
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.prediction.motif_scoring import (
    KINASE_LIBRARY_RESIDUE_CLASS_SER_THR,
    KINASE_LIBRARY_RESIDUE_CLASS_TYR,
)
from phospy.science.prediction.policies import resolve_prediction_sampling_policy
from phospy.science.references.kinase_library import KinaseLibraryResource
from phospy.science.references.resolution import (
    BundledReferenceProvider,
    ReferenceResolver,
    ReferenceResolverContract,
)
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseActivityExecutionConfig,
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.reference_projection import KinaseReferenceProjector
from phospy.workflows.kinase.site_sequence_support import (
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
    KinaseSiteSequenceSupportBuilder,
)


class _OverlapSummary(TypedDict):
    dataset_sites: int
    reference_sites: int
    overlap_sites: int
    reference_kinases: int
    kinases_with_overlap: int
    max_quantified_sites_per_kinase: int
    per_kinase_quantified: pd.Series


class KinaseWorkflowInterpreter:
    """Resolve workflow request defaults and references for execution."""

    _KINASE_COLUMN = "kinase"
    _SUBSTRATE_COLUMN = "substrate_site"

    def __init__(
        self,
        *,
        reference_resolver: ReferenceResolverContract | None = None,
        reference_projector: KinaseReferenceProjector | None = None,
        site_sequence_support_builder: KinaseSiteSequenceSupportBuilder | None = None,
    ) -> None:
        self._reference_resolver = reference_resolver or ReferenceResolver(
            provider=BundledReferenceProvider()
        )
        self._reference_projector = reference_projector or KinaseReferenceProjector()
        self._site_sequence_support_builder = (
            site_sequence_support_builder or KinaseSiteSequenceSupportBuilder()
        )

    def run(self, request: KinaseWorkflowRequest) -> ResolvedKinaseWorkflowRequest:
        dataset_view = DatasetInternalView(request.dataset)
        dataset_phospho = dataset_view.phospho
        dataset_site_metadata = dataset_view.site_metadata
        site_identity_map = self._build_site_identity_map(
            dataset=dataset_phospho,
            site_metadata=dataset_site_metadata,
        )
        site_sequence_conflict_policy = request.site_sequence_conflict_policy
        references = self._reference_resolver.run(
            request.references,
            dataset_organism=request.dataset.organism,
        )
        kinase_library_resource = self._resolve_kinase_library_resource(
            request=request,
        )
        reference_site_count = int(
            references.kinase_substrate_map.loc[:, "substrate_site"]
            .astype(str)
            .nunique()
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
        overlap_counts = self._summarize_overlap(
            dataset=dataset_phospho,
            kinase_substrate_map=kinase_substrate_map,
        )
        self._validate_reference_coverage(
            overlap_counts=overlap_counts,
            request=request,
            reference_site_count=reference_site_count,
        )
        self._validate_eligible_kinases(
            overlap_counts=overlap_counts,
            request=request,
        )
        scoring_site_index = self._resolve_scoring_site_index(
            dataset=dataset_phospho,
            site_sequences=site_sequences,
        )
        self._validate_scoring_site_support(
            scoring_site_index=scoring_site_index,
            dataset=dataset_phospho,
            site_sequences=site_sequences,
        )
        scoring_site_keys = set(scoring_site_index.astype(str))
        site_sequences = site_sequences.reindex(scoring_site_index)
        site_identity_map = site_identity_map.reindex(scoring_site_index)
        self._validate_kinase_library_resource_usability(
            resource=kinase_library_resource,
            site_sequences=site_sequences,
            request=request,
        )
        kinase_substrate_map = (
            kinase_substrate_map.loc[
                kinase_substrate_map.loc[:, "substrate_site"]
                .astype(str)
                .isin(scoring_site_keys),
                :,
            ]
            .reset_index(drop=True)
            .copy(deep=True)
        )
        activity_phospho_matrix = dataset_phospho.loc[scoring_site_index, :]
        execution_config = self._resolve_execution_config(request)
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
        aligned_metadata = site_metadata.reindex(dataset.index)
        site_key_values = (
            aligned_metadata.loc[:, "site_key"].fillna("").astype(str).str.strip()
        )
        display_id_values = (
            aligned_metadata.loc[:, "display_id"].fillna("").astype(str).str.strip()
        )
        site_keys = dataset.index.astype(str).tolist()
        if site_key_values.tolist() != site_keys:
            raise WorkflowBoundaryError(
                seam="kinase.interpreter.site_identity_map",
                next_action=(
                    "ensure dataset.site_metadata.site_key exactly matches "
                    "dataset.phospho.index"
                ),
                details={"dataset_site_count": int(dataset.index.size)},
                message_prefix="kinase workflow boundary validation failed",
            )
        if (display_id_values == "").any():
            raise WorkflowBoundaryError(
                seam="kinase.interpreter.site_identity_map",
                next_action=(
                    "ensure dataset.site_metadata.display_id is present for every "
                    "scoring site"
                ),
                details={
                    "empty_display_id_count": int((display_id_values == "").sum())
                },
                message_prefix="kinase workflow boundary validation failed",
            )
        return pd.DataFrame(
            {"site_key": site_keys, "display_id": display_id_values.tolist()},
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
            profile_missing_value_strategy=request.scoring_config.profile_missing_value_strategy,
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
            activity=activity,
        )

    def _resolve_kinase_library_resource(
        self,
        *,
        request: KinaseWorkflowRequest,
    ) -> KinaseLibraryResource | None:
        scoring_mode = request.scoring_config.scoring_mode
        if scoring_mode not in KINASE_SCORING_MODES_REQUIRING_KINASE_LIBRARY:
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

    def _validate_kinase_library_resource_usability(
        self,
        *,
        resource: KinaseLibraryResource | None,
        site_sequences: pd.DataFrame,
        request: KinaseWorkflowRequest,
    ) -> None:
        if resource is None:
            return
        resource_residue_classes = {
            matrix.residue_class.value for matrix in resource.matrices
        }
        site_residue_classes = {
            residue_class
            for residue_class in (
                _central_residue_class(
                    value,
                    upstream=int(resource.sequence_window.upstream_residues),
                    downstream=int(resource.sequence_window.downstream_residues),
                )
                for value in site_sequences.loc[:, "site_sequence"].tolist()
            )
            if residue_class is not None
        }
        if site_residue_classes.intersection(resource_residue_classes):
            return
        self._raise_boundary_error(
            seam="kinase.interpreter.kinase_library_resource_usability",
            next_action=(
                "provide a Kinase Library resource with residue_class lanes that "
                "match at least one resolved scoring-site sequence"
            ),
            scoring_mode=request.scoring_config.scoring_mode,
            resource_residue_classes=tuple(sorted(resource_residue_classes)),
            scoring_site_residue_classes=tuple(sorted(site_residue_classes)),
            scoring_site_count=int(site_sequences.shape[0]),
        )

    @classmethod
    def _summarize_overlap(
        cls,
        *,
        dataset: pd.DataFrame,
        kinase_substrate_map: pd.DataFrame,
    ) -> _OverlapSummary:
        dataset_sites = set(dataset.index.tolist())
        reference_sites = set(
            kinase_substrate_map.loc[:, cls._SUBSTRATE_COLUMN].tolist()
        )
        overlapping_sites = dataset_sites.intersection(reference_sites)
        overlapping_map = kinase_substrate_map[
            kinase_substrate_map.loc[:, cls._SUBSTRATE_COLUMN].isin(overlapping_sites)
        ]
        per_kinase_quantified = cast(
            pd.Series,
            overlapping_map.groupby(cls._KINASE_COLUMN, sort=False)[
                cls._SUBSTRATE_COLUMN
            ]
            .nunique()
            .astype("int64"),
        )
        max_quantified_sites_per_kinase = (
            0
            if per_kinase_quantified.empty
            else int(per_kinase_quantified.to_numpy(dtype="int64", copy=False).max())
        )
        return {
            "dataset_sites": len(dataset_sites),
            "reference_sites": len(reference_sites),
            "overlap_sites": len(overlapping_sites),
            "reference_kinases": int(
                kinase_substrate_map.loc[:, cls._KINASE_COLUMN].nunique()
            ),
            "kinases_with_overlap": int(per_kinase_quantified.size),
            "max_quantified_sites_per_kinase": max_quantified_sites_per_kinase,
            "per_kinase_quantified": per_kinase_quantified,
        }

    def _validate_reference_coverage(
        self,
        *,
        overlap_counts: _OverlapSummary,
        request: KinaseWorkflowRequest,
        reference_site_count: int,
    ) -> None:
        overlap_sites = int(overlap_counts["overlap_sites"])
        if overlap_sites > 0:
            return
        self._raise_boundary_error(
            seam="kinase.interpreter.reference_coverage",
            next_action=(
                "use references whose display IDs match dataset display_id metadata "
                "or verify dataset site_key/display_id identity mapping"
            ),
            dataset_sites=overlap_counts["dataset_sites"],
            reference_sites=reference_site_count,
            overlap_sites=overlap_sites,
            scoring_config_min_substrates=request.scoring_config.min_substrates,
        )

    def _validate_eligible_kinases(
        self,
        *,
        overlap_counts: _OverlapSummary,
        request: KinaseWorkflowRequest,
    ) -> None:
        per_kinase_quantified = overlap_counts["per_kinase_quantified"]
        if not isinstance(per_kinase_quantified, pd.Series):
            raise PhosPyWorkflowError(
                "kinase workflow interpreter expected "
                "overlap_counts['per_kinase_quantified'] to be a pandas Series, "
                f"got {type(per_kinase_quantified).__name__}"
            )
        eligible_mask = per_kinase_quantified >= request.scoring_config.min_substrates
        eligible_kinases = per_kinase_quantified[eligible_mask]
        if not eligible_kinases.empty:
            return
        self._raise_boundary_error(
            seam="kinase.interpreter.eligible_kinases",
            next_action=(
                "lower scoring_config.min_substrates or provide references with "
                "deeper overlap for the current dataset "
                "(scientific floor: min_substrates >= 2)"
            ),
            reference_kinases=overlap_counts["reference_kinases"],
            kinases_with_overlap=overlap_counts["kinases_with_overlap"],
            eligible_kinases=int(eligible_kinases.size),
            max_quantified_sites_per_kinase=overlap_counts[
                "max_quantified_sites_per_kinase"
            ],
            scoring_config_min_substrates=request.scoring_config.min_substrates,
            prediction_config_deterministic_max_selected_kinases=(
                request.prediction_config.deterministic_max_selected_kinases
            ),
            prediction_config_adaptive_ensemble_runs=(
                request.prediction_config.adaptive_ensemble_runs
            ),
            prediction_config_mode=request.prediction_config.mode,
        )

    @staticmethod
    def _resolve_scoring_site_index(
        *,
        dataset: pd.DataFrame,
        site_sequences: pd.DataFrame,
    ) -> pd.Index:
        sequence_sites = set(site_sequences.index.tolist())
        scoring_sites = [
            site_id for site_id in dataset.index if site_id in sequence_sites
        ]
        return pd.Index(scoring_sites, name=dataset.index.name)

    def _validate_scoring_site_support(
        self,
        *,
        scoring_site_index: pd.Index,
        dataset: pd.DataFrame,
        site_sequences: pd.DataFrame,
    ) -> None:
        if not scoring_site_index.empty:
            return
        self._raise_boundary_error(
            seam="kinase.interpreter.sequence_support",
            next_action=(
                "ensure references.site_sequences display IDs can be projected to "
                "dataset site_key rows"
            ),
            dataset_sites=int(dataset.index.size),
            reference_sequence_sites=int(site_sequences.index.size),
            sequence_supported_sites=0,
        )

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


def _central_residue_class(
    value: object,
    *,
    upstream: int,
    downstream: int,
) -> str | None:
    if not isinstance(value, str):
        return None
    sequence = value.strip().upper()
    if not sequence:
        return None
    window_size = int(upstream) + 1 + int(downstream)
    if len(sequence) == window_size:
        centre_index = int(upstream)
    elif len(sequence) >= window_size and len(sequence) % 2 == 1:
        centre_index = len(sequence) // 2
    else:
        return None
    if centre_index < 0 or centre_index >= len(sequence):
        return None
    residue = sequence[centre_index]
    if residue in {"S", "T"}:
        return KINASE_LIBRARY_RESIDUE_CLASS_SER_THR
    if residue == "Y":
        return KINASE_LIBRARY_RESIDUE_CLASS_TYR
    return None
