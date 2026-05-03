"""Internal interpreter for kinase workflow requests."""

from __future__ import annotations

from typing import NoReturn, TypedDict, cast

import pandas as pd

from phospy.api.requests import KinaseWorkflowRequest
from phospy.errors.workflows import PhosPyWorkflowError, WorkflowBoundaryError
from phospy.prediction.policies import resolve_prediction_sampling_policy
from phospy.references.resolution import (
    BundledReferenceProvider,
    ReferenceResolver,
    ReferenceResolverContract,
)
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseActivityExecutionConfig,
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)


class _OverlapSummary(TypedDict):
    dataset_sites: int
    reference_sites: int
    overlap_sites: int
    reference_kinases: int
    kinases_with_overlap: int
    max_quantified_sites_per_kinase: int
    per_kinase_quantified: pd.Series


class _SiteSequenceMergeSummary(TypedDict):
    site_sequences: pd.DataFrame
    dataset_sequences_added: int
    dataset_reference_conflict_count: int
    dataset_sequences_missing: int
    dataset_sequences_available: int


class KinaseWorkflowInterpreter:
    """Resolve workflow request defaults and references for execution."""

    _KINASE_COLUMN = "kinase"
    _SUBSTRATE_COLUMN = "substrate_site"
    _SITE_SEQUENCE_COLUMN = "site_sequence"

    def __init__(
        self, *, reference_resolver: ReferenceResolverContract | None = None
    ) -> None:
        self._reference_resolver = reference_resolver or ReferenceResolver(
            provider=BundledReferenceProvider()
        )

    def run(self, request: KinaseWorkflowRequest) -> ResolvedKinaseWorkflowRequest:
        references = self._reference_resolver.run(
            request.references,
            dataset_organism=request.dataset.organism,
        )
        # ReferenceBundle construction is the sole identifier-normalisation
        # boundary. Downstream workflow stages consume these tables as-is.
        kinase_substrate_map = references.kinase_substrate_map
        merge_summary = self._build_execution_site_sequences(
            dataset=request.dataset.phospho,
            site_metadata=request.dataset.site_metadata,
            reference_site_sequences=references.site_sequences,
        )
        site_sequences = merge_summary["site_sequences"]
        overlap_counts = self._summarize_overlap(
            dataset=request.dataset.phospho,
            kinase_substrate_map=kinase_substrate_map,
        )
        self._validate_reference_coverage(
            overlap_counts=overlap_counts,
            request=request,
        )
        self._validate_eligible_kinases(
            overlap_counts=overlap_counts,
            request=request,
        )
        scoring_site_index = self._resolve_scoring_site_index(
            dataset=request.dataset.phospho,
            site_sequences=site_sequences,
        )
        self._validate_scoring_site_support(
            scoring_site_index=scoring_site_index,
            dataset=request.dataset.phospho,
            site_sequences=site_sequences,
        )
        activity_phospho_matrix = request.dataset.phospho.loc[scoring_site_index, :]
        execution_config = self._resolve_execution_config(request)
        return ResolvedKinaseWorkflowRequest(
            dataset=request.dataset,
            references=references,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            scoring_site_index=scoring_site_index,
            activity_phospho_matrix=activity_phospho_matrix,
            execution_config=execution_config,
            site_sequence_merge_diagnostics={
                "dataset_sequences_added": int(
                    merge_summary["dataset_sequences_added"]
                ),
                "dataset_reference_conflict_count": int(
                    merge_summary["dataset_reference_conflict_count"]
                ),
                "dataset_sequences_missing": int(
                    merge_summary["dataset_sequences_missing"]
                ),
                "dataset_sequences_available": int(
                    merge_summary["dataset_sequences_available"]
                ),
                "reference_sequence_count": int(references.site_sequences.shape[0]),
                "execution_sequence_count": int(site_sequences.shape[0]),
            },
        )

    @classmethod
    def _build_execution_site_sequences(
        cls,
        *,
        dataset: pd.DataFrame,
        site_metadata: pd.DataFrame,
        reference_site_sequences: pd.DataFrame,
    ) -> _SiteSequenceMergeSummary:
        merged = reference_site_sequences.copy(deep=True)
        dataset_sequences_available = 0
        dataset_sequences_missing = 0
        dataset_sequences_added = 0
        dataset_reference_conflict_count = 0
        if "site_sequence" not in site_metadata.columns:
            return {
                "site_sequences": merged,
                "dataset_sequences_added": 0,
                "dataset_reference_conflict_count": 0,
                "dataset_sequences_missing": int(dataset.shape[0]),
                "dataset_sequences_available": 0,
            }
        dataset_sequence_series = (
            site_metadata.reindex(dataset.index)
            .loc[:, "site_sequence"]
            .astype("string")
            .str.strip()
        )
        for site_id in dataset.index.tolist():
            site_key = str(site_id)
            sequence_value = dataset_sequence_series.loc[site_id]
            has_sequence = bool(pd.notna(sequence_value)) and str(sequence_value) != ""
            if not has_sequence:
                dataset_sequences_missing += 1
                continue
            dataset_sequences_available += 1
            sequence_text = str(sequence_value)
            if site_key in merged.index:
                reference_value = str(merged.at[site_key, cls._SITE_SEQUENCE_COLUMN])
                if reference_value != sequence_text:
                    dataset_reference_conflict_count += 1
                continue
            merged.loc[site_key, cls._SITE_SEQUENCE_COLUMN] = sequence_text
            dataset_sequences_added += 1
        if merged.index.name != reference_site_sequences.index.name:
            merged.index.name = reference_site_sequences.index.name
        return {
            "site_sequences": merged,
            "dataset_sequences_added": dataset_sequences_added,
            "dataset_reference_conflict_count": dataset_reference_conflict_count,
            "dataset_sequences_missing": dataset_sequences_missing,
            "dataset_sequences_available": dataset_sequences_available,
        }

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
            )
        )
        return ResolvedKinaseExecutionConfig(
            scoring_min_substrates=int(request.scoring_config.min_substrates),
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
    ) -> None:
        overlap_sites = int(overlap_counts["overlap_sites"])
        if overlap_sites > 0:
            return
        self._raise_boundary_error(
            seam="kinase.interpreter.reference_coverage",
            next_action=(
                "use references that contain dataset phosphosite IDs or verify site "
                "identifier formatting in dataset.phospho.index"
            ),
            dataset_sites=overlap_counts["dataset_sites"],
            reference_sites=overlap_counts["reference_sites"],
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
        eligible_kinases = cast(pd.Series, per_kinase_quantified[eligible_mask])
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
                "ensure references.site_sequences contains sequence entries for "
                "dataset phosphosite IDs"
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
