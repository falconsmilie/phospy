"""Scoring orchestration for kinase workflow execution."""

from __future__ import annotations

from collections.abc import Callable

import pandas as pd

from phospy.contracts.configs import (
    KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
    KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF,
)
from phospy.errors.workflows import WorkflowStageError
from phospy.science.prediction.models import KinaseScoringResult
from phospy.science.prediction.motif_scoring import (
    DEFAULT_MOTIF_FLANK_SIZE,
    SEQUENCE_SEMANTICS_CENTRED_SEQUENCE,
    MotifScoringResult,
    build_motif_library,
    get_motif_library_validation,
    score_phosphosite_motifs,
)
from phospy.science.prediction.scoring import (
    SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY,
    DownstreamScoreSelectionPolicy,
    build_kinase_score_source_diagnostics,
    build_kinase_score_source_summary,
    fuse_profile_and_motif_scores_by_rank_weight,
    resolve_downstream_score_matrix,
)
from phospy.science.scoring.policy_models import DownstreamScoreSource
from phospy.workflows.kinase.component_models import KinaseScoringRunResult
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.kinase_library_scoring import (
    KINASE_LIBRARY_WORKFLOW_SCORE_SCALE,
    KinaseLibraryWorkflowScorer,
)
from phospy.workflows.kinase.science import (
    KinaseProfileBuild,
    build_kinase_profiles,
    score_profile_correlations,
)


class KinaseScoringRunner:
    """Run profile/motif scoring and resolve the downstream score lane."""

    def __init__(
        self,
        *,
        build_profiles: Callable[..., KinaseProfileBuild] = build_kinase_profiles,
        score_profiles: Callable[..., pd.DataFrame] = score_profile_correlations,
        build_motif_library_fn: Callable[
            ..., tuple[dict[str, pd.DataFrame], pd.Series]
        ] = build_motif_library,
        get_motif_library_validation_fn: Callable[..., object | None] = (
            get_motif_library_validation
        ),
        score_motifs: Callable[..., MotifScoringResult] = score_phosphosite_motifs,
        fuse_scores: Callable[..., tuple[pd.DataFrame, pd.DataFrame | None]] = (
            fuse_profile_and_motif_scores_by_rank_weight
        ),
        select_downstream: Callable[
            ...,
            tuple[pd.DataFrame, DownstreamScoreSource]
            | tuple[
                pd.DataFrame,
                DownstreamScoreSource,
                DownstreamScoreSelectionPolicy,
            ],
        ] = (resolve_downstream_score_matrix),
        kinase_library_scorer: KinaseLibraryWorkflowScorer | None = None,
    ) -> None:
        self._build_profiles = build_profiles
        self._score_profiles = score_profiles
        self._build_motif_library = build_motif_library_fn
        self._get_motif_library_validation = get_motif_library_validation_fn
        self._score_motifs = score_motifs
        self._fuse_scores = fuse_scores
        self._select_downstream = select_downstream
        self._kinase_library_scorer = (
            kinase_library_scorer or KinaseLibraryWorkflowScorer()
        )

    def run(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        config: ResolvedKinaseExecutionConfig,
    ) -> KinaseScoringRunResult:
        include_diagnostic_tables = config.include_diagnostic_scoring_tables
        scoring_min_substrates = int(config.scoring_min_substrates)
        scoring_phospho = request.activity_phospho_matrix
        profile_build = self._build_profiles(
            phospho=scoring_phospho,
            kinase_substrate_map=request.kinase_substrate_map,
            min_substrates=scoring_min_substrates,
            allow_single_substrate_profiles=False,
            profile_missing_value_strategy=config.profile_missing_value_strategy,
        )
        if profile_build.profile_matrix.empty:
            raise WorkflowStageError(
                "kinase workflow internal invariant failed at seam="
                "kinase.executor.scoring_profiles; interpreter should reject "
                "requests with zero eligible kinases before scoring"
            )
        profile_scores = self._score_profiles(
            phospho=scoring_phospho,
            profile_matrix=profile_build.profile_matrix,
        )
        sequence_series = request.site_sequences.loc[:, "site_sequence"]
        site_identity_series = request.site_sequences.loc[:, "display_id"]
        if config.scoring_mode == KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF:
            return self._run_kinase_library_motif_mode(
                request=request,
                config=config,
                profile_scores=profile_scores,
                profile_build=profile_build,
                sequence_series=sequence_series,
                site_identity_series=site_identity_series,
            )
        if config.scoring_mode == KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF:
            return self._run_combined_profile_motif_mode(
                request=request,
                config=config,
                profile_scores=profile_scores,
                profile_build=profile_build,
                sequence_series=sequence_series,
                site_identity_series=site_identity_series,
            )
        eligible_kinases = set(profile_scores.columns.astype(str))
        # `request.kinase_substrate_map` comes from ReferenceBundle and is
        # already normalized at reference-ingestion time.
        motif_kinase_substrate_map = request.kinase_substrate_map.loc[
            request.kinase_substrate_map.loc[:, "kinase"]
            .astype(str)
            .isin(eligible_kinases)
        ]
        motif_result = self._resolve_motif_scores(
            scoring_min_substrates=scoring_min_substrates,
            scoring_phospho=scoring_phospho,
            motif_kinase_substrate_map=motif_kinase_substrate_map,
            sequence_series=sequence_series,
            site_identity_series=site_identity_series,
        )
        try:
            (
                rank_weighted_fusion_scores,
                score_fusion_weights,
            ) = self._fuse_scores(
                motif_scores=motif_result.motif_scores,
                profile_scores=profile_scores,
                motif_sizes=motif_result.motif_sizes,
                profile_sizes=profile_build.substrate_counts.astype(float),
                allow_profile_only_fallback=True,
                emit_weights=include_diagnostic_tables,
            )
        except ValueError as exc:
            raise WorkflowStageError(
                "kinase workflow internal invariant failed at seam="
                "kinase.executor.rank_weighted_fusion_scoring; "
                f"{exc}"
            ) from exc
        diagnostic_motif_scores: pd.DataFrame | None = None
        diagnostic_score_source_matrix: pd.DataFrame | None = None
        if include_diagnostic_tables:
            score_source_matrix, score_source_summary = (
                build_kinase_score_source_diagnostics(
                    motif_scores=motif_result.motif_scores,
                    profile_scores=profile_scores,
                    rank_weighted_fusion_scores=rank_weighted_fusion_scores,
                )
            )
            diagnostic_motif_scores = motif_result.motif_scores
            if diagnostic_motif_scores.empty:
                diagnostic_motif_scores = pd.DataFrame(
                    index=motif_result.motif_scores.index.copy(),
                    columns=profile_scores.columns.copy(),
                    dtype=float,
                )
            diagnostic_score_source_matrix = score_source_matrix
        else:
            score_source_summary = build_kinase_score_source_summary(
                motif_scores=motif_result.motif_scores,
                profile_scores=profile_scores,
                rank_weighted_fusion_scores=rank_weighted_fusion_scores,
            )
        scoring_result = KinaseScoringResult._from_owned(
            profile_scores=profile_scores,
            motif_scores=diagnostic_motif_scores,
            rank_weighted_fusion_scores=rank_weighted_fusion_scores,
            score_fusion_weights=score_fusion_weights,
            score_source_matrix=diagnostic_score_source_matrix,
            score_source_summary=score_source_summary,
            motif_sequence_validation=motif_result.sequence_validation,
            motif_library_validation=motif_result.library_validation,
            scoring_mode=config.scoring_mode,
            score_source=DownstreamScoreSource.RANK_WEIGHTED_FUSION_SCORES,
            score_scale="relative_support_score_unit_interval",
        )
        downstream_selection = self._select_downstream(
            profile_scores=profile_scores,
            rank_weighted_fusion_scores=rank_weighted_fusion_scores,
        )
        if len(downstream_selection) == 3:
            (
                downstream_score_matrix,
                downstream_score_source,
                downstream_score_selection_policy,
            ) = downstream_selection
        else:
            downstream_score_matrix, downstream_score_source = downstream_selection
            downstream_score_selection_policy = (
                SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY
            )
        return KinaseScoringRunResult(
            scoring_result=scoring_result,
            downstream_score_matrix=downstream_score_matrix,
            downstream_score_source=downstream_score_source,
            quantified_substrates=profile_build.quantified_substrates,
            downstream_score_selection_policy=downstream_score_selection_policy,
        )

    def _run_kinase_library_motif_mode(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        config: ResolvedKinaseExecutionConfig,
        profile_scores: pd.DataFrame,
        profile_build: KinaseProfileBuild,
        sequence_series: pd.Series,
        site_identity_series: pd.Series,
    ) -> KinaseScoringRunResult:
        library_resource = request.kinase_library_resource
        if library_resource is None:
            raise WorkflowStageError(
                "kinase workflow internal invariant failed at seam="
                "kinase.executor.kinase_library_resource; interpreter should "
                "resolve a KinaseLibraryResource for Kinase Library scoring modes"
            )
        library_result = self._kinase_library_scorer.run(
            resource=library_resource,
            site_sequences=sequence_series.loc[request.scoring_site_index],
            site_identities=site_identity_series.loc[request.scoring_site_index],
            site_index=tuple(str(site_id) for site_id in request.scoring_site_index),
        )
        scoring_result = KinaseScoringResult._from_owned(
            profile_scores=profile_scores,
            kinase_library_motif_scores=library_result.scores,
            kinase_library_site_diagnostics=library_result.site_diagnostics,
            kinase_library_kinase_diagnostics=library_result.kinase_diagnostics,
            scoring_mode=config.scoring_mode,
            score_source=DownstreamScoreSource.KINASE_LIBRARY_MOTIF_SCORES,
            score_scale=KINASE_LIBRARY_WORKFLOW_SCORE_SCALE,
            score_scale_metadata=library_result.score_scale_metadata,
        )
        return KinaseScoringRunResult(
            scoring_result=scoring_result,
            downstream_score_matrix=library_result.scores,
            downstream_score_source=DownstreamScoreSource.KINASE_LIBRARY_MOTIF_SCORES,
            quantified_substrates=profile_build.quantified_substrates,
            downstream_score_selection_policy=None,
        )

    def _run_combined_profile_motif_mode(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        config: ResolvedKinaseExecutionConfig,
        profile_scores: pd.DataFrame,
        profile_build: KinaseProfileBuild,
        sequence_series: pd.Series,
        site_identity_series: pd.Series,
    ) -> KinaseScoringRunResult:
        library_resource = request.kinase_library_resource
        if library_resource is None:
            raise WorkflowStageError(
                "kinase workflow internal invariant failed at seam="
                "kinase.executor.combined_profile_motif_resource; interpreter "
                "should resolve a KinaseLibraryResource for combined scoring"
            )
        library_result = self._kinase_library_scorer.run(
            resource=library_resource,
            site_sequences=sequence_series.loc[request.scoring_site_index],
            site_identities=site_identity_series.loc[request.scoring_site_index],
            site_index=tuple(str(site_id) for site_id in request.scoring_site_index),
        )
        motif_sizes = library_result.scores.notna().sum(axis=0).astype(float)
        motif_sizes.index.name = "kinase"
        try:
            combined_scores, score_fusion_weights = self._fuse_scores(
                motif_scores=library_result.scores,
                profile_scores=profile_scores,
                motif_sizes=motif_sizes,
                profile_sizes=profile_build.substrate_counts.astype(float),
                allow_profile_only_fallback=True,
                emit_weights=bool(config.include_diagnostic_scoring_tables),
            )
        except ValueError as exc:
            raise WorkflowStageError(
                "kinase workflow internal invariant failed at seam="
                "kinase.executor.combined_profile_motif_scoring; "
                f"{exc}"
            ) from exc
        scoring_result = KinaseScoringResult._from_owned(
            profile_scores=profile_scores,
            kinase_library_motif_scores=library_result.scores,
            combined_profile_motif_scores=combined_scores,
            score_fusion_weights=score_fusion_weights,
            kinase_library_site_diagnostics=library_result.site_diagnostics,
            kinase_library_kinase_diagnostics=library_result.kinase_diagnostics,
            scoring_mode=config.scoring_mode,
            score_source=DownstreamScoreSource.COMBINED_PROFILE_MOTIF_SCORES,
            score_scale="combined_profile_kinase_library_motif_unit_interval",
            score_scale_metadata={
                **library_result.score_scale_metadata,
                "authoritative_matrix": "combined_profile_motif_scores",
                "profile_score_scale": "relative_support_score_unit_interval",
                "fusion_policy": "rank_weighted_motif_profile_fusion_v1",
            },
        )
        return KinaseScoringRunResult(
            scoring_result=scoring_result,
            downstream_score_matrix=combined_scores,
            downstream_score_source=DownstreamScoreSource.COMBINED_PROFILE_MOTIF_SCORES,
            quantified_substrates=profile_build.quantified_substrates,
            downstream_score_selection_policy=None,
        )

    def _resolve_motif_scores(
        self,
        *,
        scoring_min_substrates: int,
        scoring_phospho: pd.DataFrame,
        motif_kinase_substrate_map: pd.DataFrame,
        sequence_series: pd.Series,
        site_identity_series: pd.Series,
    ) -> MotifScoringResult:
        motif_frequency_matrices, motif_sizes = self._build_motif_library(
            kinase_substrate_map=motif_kinase_substrate_map,
            site_sequences=sequence_series,
            site_identities=site_identity_series,
            flank_size=DEFAULT_MOTIF_FLANK_SIZE,
        )
        motif_library_validation = self._get_motif_library_validation(motif_sizes)
        return self._score_motifs(
            site_sequences=sequence_series.loc[scoring_phospho.index],
            site_identities=site_identity_series.loc[scoring_phospho.index],
            motif_frequency_matrices=motif_frequency_matrices,
            motif_sizes=motif_sizes,
            site_index=tuple(str(site_id) for site_id in scoring_phospho.index),
            min_motif_size=scoring_min_substrates,
            flank_size=DEFAULT_MOTIF_FLANK_SIZE,
            sequence_semantics=SEQUENCE_SEMANTICS_CENTRED_SEQUENCE,
            library_validation=motif_library_validation,
        )


__all__ = ["KinaseScoringRunner"]
