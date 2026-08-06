"""ssGSEA-style substrate-set enrichment activity-like score method."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.frames.validation import (
    require_canonical_string_column,
    require_columns,
    require_dataframe,
    require_non_empty_string_column,
    require_unique_columns,
    require_unique_index,
    require_unique_row_pairs,
)
from phospy.science.activities.method_contracts import (
    ssgsea_substrate_enrichment_activity_input_contract,
)
from phospy.science.activities.models import (
    SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_METHOD,
    ActivityMethodSummary,
    KinaseActivityResult,
    SsgseaSubstrateEnrichmentActivityDiagnostics,
)
from phospy.science.activities.scientific_policies import (
    SSGSEA_PERMUTATION_RNG_SEED_POLICY_VERSION,
    SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_POLICY_VERSION,
    SSGSEA_TIE_POLICY,
    build_ssgsea_substrate_enrichment_activity_policy,
)
from phospy.science.activities.semantics import (
    ActivityInputMatrix,
    ActivityInputSemantics,
    ActivityProfileAxis,
    ActivityQuantitativeSemantics,
    normalize_activity_input_matrix,
)
from phospy.science.activities.statistics import benjamini_hochberg_q_values

SSGSEA_RANKING_DIRECTION_DESCENDING = "descending"
SSGSEA_RANKING_DIRECTION_ASCENDING = "ascending"
SSGSEA_RANKING_DIRECTIONS = frozenset(
    {
        SSGSEA_RANKING_DIRECTION_DESCENDING,
        SSGSEA_RANKING_DIRECTION_ASCENDING,
    }
)
SSGSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG = "benjamini_hochberg"

SSGSEA_STATUS_COMPUTED = "computed"
SSGSEA_STATUS_INSUFFICIENT_SUBSTRATES = "insufficient_substrates"
SSGSEA_STATUS_NO_FINITE_BACKGROUND_VALUES = "no_finite_background_values"
SSGSEA_STATUS_INSUFFICIENT_BACKGROUND_SITES = "insufficient_background_sites"
SSGSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES = "no_finite_substrate_values"
SSGSEA_SIGNIFICANCE_STATUS_AVAILABLE = "permutation_significance_available"
SSGSEA_SIGNIFICANCE_STATUS_P_VALUE_AVAILABLE_Q_VALUE_DISABLED = (
    "permutation_p_value_available_q_value_unavailable_adjustment_disabled"
)
SSGSEA_SIGNIFICANCE_STATUS_UNAVAILABLE_NO_PERMUTATIONS = (
    "significance_unavailable_no_permutations"
)
SSGSEA_SIGNIFICANCE_STATUS_UNAVAILABLE_NOT_COMPUTABLE = (
    "significance_unavailable_not_computable"
)

_KINASE_COLUMN = "kinase"
_SUBSTRATE_COLUMN = "substrate_site"
_SSGSEA_PERMUTATION_STREAM_NAME = "substrate_label_permutation"
_SSGSEA_PERMUTATION_SEED_DIGEST_SIZE_BYTES = 16
_SSGSEA_PERMUTATION_RNG_SEED_HASH_POLICY_TOKEN = "stable_by_method_condition_kinase"
_SSGSEA_PERMUTATION_EXTREME_ATOL = 1e-12


@dataclass(frozen=True, slots=True)
class _RankedTieBlocks:
    site_labels: np.ndarray
    block_starts: np.ndarray
    block_sizes: np.ndarray

    @property
    def n_background(self) -> int:
        return int(self.site_labels.size)

    @property
    def n_blocks(self) -> int:
        return int(self.block_sizes.size)

    @property
    def has_ties(self) -> bool:
        return bool((self.block_sizes > 1).any())

    @property
    def n_tie_blocks(self) -> int:
        return int((self.block_sizes > 1).sum())

    @property
    def n_tied_sites(self) -> int:
        tie_mask = self.block_sizes > 1
        if not bool(tie_mask.any()):
            return 0
        return int(self.block_sizes[tie_mask].sum())

    @property
    def max_tie_block_size(self) -> int:
        if not self.has_ties:
            return 0
        return int(self.block_sizes.max())


@dataclass(frozen=True, slots=True)
class _SsgseaNullScoreCacheKey:
    n_background: int
    n_substrates: int
    tie_block_sizes: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class _SsgseaNullScoreEngine:
    """Reusable scoring constants for equivalent ssGSEA permutation nulls."""

    n_background: int
    n_substrates: int
    has_ties: bool
    rank_weights: np.ndarray | None
    rank_weight_multiplier: float
    rank_weight_constant: float
    block_index_by_position: np.ndarray | None
    block_hit_coefficients: np.ndarray | None
    block_constant: float

    @classmethod
    def from_ranked_blocks(
        cls,
        *,
        ranked_blocks: _RankedTieBlocks,
        n_substrates: int,
    ) -> _SsgseaNullScoreEngine:
        n_background = int(ranked_blocks.n_background)
        n_substrates = int(n_substrates)
        n_misses = n_background - n_substrates
        if n_background == 0 or n_substrates == 0 or n_misses == 0:
            raise ValueError(
                "ssgsea null score engine requires at least one background, "
                "substrate, and non-substrate site"
            )

        multiplier = (1.0 / float(n_substrates)) + (1.0 / float(n_misses))
        if not ranked_blocks.has_ties:
            rank_weights = (n_background - np.arange(n_background, dtype=float)).astype(
                float, copy=False
            )
            total_rank_weight = float(n_background * (n_background + 1) / 2)
            return cls(
                n_background=n_background,
                n_substrates=n_substrates,
                has_ties=False,
                rank_weights=rank_weights,
                rank_weight_multiplier=float(multiplier),
                rank_weight_constant=float(total_rank_weight / float(n_misses)),
                block_index_by_position=None,
                block_hit_coefficients=None,
                block_constant=0.0,
            )

        block_sizes = ranked_blocks.block_sizes.astype(float, copy=False)
        suffix_site_counts = (
            n_background - np.cumsum(ranked_blocks.block_sizes)
        ).astype(float, copy=False)
        block_area_coefficients = suffix_site_counts + ((block_sizes + 1.0) / 2.0)
        block_index_by_position = np.empty(n_background, dtype=np.int64)
        for block_index, block_start, block_size in zip(
            range(ranked_blocks.n_blocks),
            ranked_blocks.block_starts.tolist(),
            ranked_blocks.block_sizes.tolist(),
            strict=True,
        ):
            block_index_by_position[
                int(block_start) : int(block_start) + int(block_size)
            ] = int(block_index)

        return cls(
            n_background=n_background,
            n_substrates=n_substrates,
            has_ties=True,
            rank_weights=None,
            rank_weight_multiplier=0.0,
            rank_weight_constant=0.0,
            block_index_by_position=block_index_by_position,
            block_hit_coefficients=(multiplier * block_area_coefficients).astype(
                float,
                copy=False,
            ),
            block_constant=float(
                np.dot(block_sizes / float(n_misses), block_area_coefficients)
            ),
        )

    def score_selected_positions(self, hit_positions: np.ndarray) -> float:
        positions = np.asarray(hit_positions, dtype=np.int64)
        if int(positions.size) != int(self.n_substrates):
            return np.nan
        if not self.has_ties:
            if self.rank_weights is None:
                raise RuntimeError("ssgsea untied null score engine is incomplete")
            hit_weight_sum = float(self.rank_weights[positions].sum())
            area = hit_weight_sum * float(self.rank_weight_multiplier) - float(
                self.rank_weight_constant
            )
            return float(area / float(self.n_background))

        if self.block_index_by_position is None or self.block_hit_coefficients is None:
            raise RuntimeError("ssgsea tied null score engine is incomplete")
        block_indices = self.block_index_by_position[positions]
        block_hit_counts = np.bincount(
            block_indices,
            minlength=int(self.block_hit_coefficients.size),
        )
        area = float(np.dot(block_hit_counts, self.block_hit_coefficients)) - float(
            self.block_constant
        )
        return float(area / float(self.n_background))


class _SsgseaNullScoreCache:
    """Cache null scoring engines for mathematically equivalent cases.

    For untied profiles the permutation null depends only on the finite
    background size and selected substrate count. When tied rank blocks are
    present, the tie-block size sequence is part of the null definition.
    """

    def __init__(self) -> None:
        self._engines: dict[
            _SsgseaNullScoreCacheKey,
            _SsgseaNullScoreEngine,
        ] = {}

    @property
    def engine_count(self) -> int:
        return len(self._engines)

    def get(
        self,
        *,
        ranked_blocks: _RankedTieBlocks,
        n_substrates: int,
    ) -> _SsgseaNullScoreEngine:
        key = _null_score_cache_key(
            ranked_blocks=ranked_blocks,
            n_substrates=int(n_substrates),
        )
        engine = self._engines.get(key)
        if engine is None:
            engine = _SsgseaNullScoreEngine.from_ranked_blocks(
                ranked_blocks=ranked_blocks,
                n_substrates=int(n_substrates),
            )
            self._engines[key] = engine
        return engine


@dataclass(frozen=True, slots=True)
class SsgseaSubstrateEnrichmentActivityMethod:
    """PhosPy ssGSEA-style kinase substrate enrichment score."""

    min_substrates: int
    ranking_direction: str = SSGSEA_RANKING_DIRECTION_DESCENDING
    permutation_count: int = 0
    random_seed: int | None = 0
    adjust_p_values: bool = True
    q_value_method: str = SSGSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG

    def __post_init__(self) -> None:
        if isinstance(self.min_substrates, bool) or not isinstance(
            self.min_substrates,
            int,
        ):
            raise ValueError("ssgsea min_substrates must be an int")
        if int(self.min_substrates) < 1:
            raise ValueError("ssgsea min_substrates must be greater than or equal to 1")
        if self.ranking_direction not in SSGSEA_RANKING_DIRECTIONS:
            allowed = ", ".join(sorted(SSGSEA_RANKING_DIRECTIONS))
            raise ValueError(f"ssgsea ranking_direction must be one of: {allowed}")
        if isinstance(self.permutation_count, bool) or not isinstance(
            self.permutation_count,
            int,
        ):
            raise ValueError("ssgsea permutation_count must be an int")
        if int(self.permutation_count) < 0:
            raise ValueError(
                "ssgsea permutation_count must be greater than or equal to 0"
            )
        if self.random_seed is not None and (
            isinstance(self.random_seed, bool) or not isinstance(self.random_seed, int)
        ):
            raise ValueError("ssgsea random_seed must be an int or None")
        if self.random_seed is not None and int(self.random_seed) < 0:
            raise ValueError("ssgsea random_seed must be greater than or equal to 0")
        if int(self.permutation_count) > 0 and self.random_seed is None:
            raise ValueError(
                "ssgsea random_seed must be set when permutation_count is positive"
            )
        if not isinstance(self.adjust_p_values, bool):
            raise ValueError("ssgsea adjust_p_values must be a bool")
        if self.q_value_method != SSGSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG:
            raise ValueError("ssgsea q_value_method must be 'benjamini_hochberg'")

    def run(
        self,
        *,
        activity_input: ActivityInputMatrix | None = None,
        effect_matrix: pd.DataFrame | None = None,
        kinase_substrate_membership: pd.DataFrame,
    ) -> KinaseActivityResult:
        if activity_input is not None and effect_matrix is not None:
            raise WorkflowBoundaryError(
                "ssgsea activity requires either activity_input or legacy "
                "effect_matrix, not both"
            )
        if activity_input is None:
            if effect_matrix is None:
                raise WorkflowBoundaryError(
                    "ssgsea activity requires ActivityInputMatrix with explicit "
                    "contrast/effect semantics"
                )
            activity_input = normalize_activity_input_matrix(
                effect_matrix,
                field_name="activity_inputs.effect_matrix",
                legacy_dataframe_semantics=ActivityInputSemantics(
                    profile_axis=ActivityProfileAxis.EFFECT,
                    quantitative_semantics=(
                        ActivityQuantitativeSemantics.STANDARDISED_EFFECT
                    ),
                ),
                legacy_dataframe_warning=(
                    "Passing a raw DataFrame as effect_matrix is deprecated; "
                    "provide ActivityInputMatrix.contrast_log_fold_change(...) or "
                    "ActivityInputMatrix.standardised_effect(...) so activity "
                    "input semantics are explicit."
                ),
            )
        from phospy.science.quantitative_method_contracts import (
            resolve_activity_input_contract,
        )

        resolve_activity_input_contract(
            activity_input=activity_input,
            contract=ssgsea_substrate_enrichment_activity_input_contract(),
            context=(
                "ssGSEA substrate enrichment activity input requires explicit "
                "contrast/effect input"
            ),
        )
        effects = activity_input.frame
        membership = _validate_membership_table(kinase_substrate_membership)
        site_labels = np.asarray(effects.index.astype(str).tolist(), dtype=object)
        site_universe = set(str(value) for value in site_labels.tolist())
        aligned_membership = membership.loc[
            membership.loc[:, _SUBSTRATE_COLUMN].astype(str).isin(site_universe),
            :,
        ].copy(deep=True)
        kinases = _ordered_unique_strings(membership.loc[:, _KINASE_COLUMN])
        kinase_index = pd.Index(kinases, name="kinase")
        profile_index = pd.Index(
            effects.columns.astype(str).tolist(),
            name=effects.columns.name,
        )
        membership_by_kinase = _build_membership_lookup(
            kinases=kinases,
            membership=aligned_membership,
        )

        activity_scores = pd.DataFrame(
            np.nan,
            index=kinase_index,
            columns=profile_index,
            dtype=float,
        )
        substrate_count_table = pd.DataFrame(
            0,
            index=kinase_index,
            columns=profile_index,
            dtype=int,
        )
        p_value_matrix = (
            pd.DataFrame(
                np.nan,
                index=kinase_index,
                columns=profile_index,
                dtype=float,
            )
            if int(self.permutation_count) > 0
            else None
        )
        q_value_matrix = (
            pd.DataFrame(
                np.nan,
                index=kinase_index,
                columns=profile_index,
                dtype=float,
            )
            if int(self.permutation_count) > 0 and bool(self.adjust_p_values)
            else None
        )

        random_seed: int | None = None
        if int(self.permutation_count) > 0:
            if self.random_seed is None:
                raise ValueError(
                    "ssgsea random_seed must be set when permutation_count is positive"
                )
            random_seed = int(self.random_seed)
        effect_values = effects.to_numpy(dtype=float, copy=False)
        rows: list[dict[str, object]] = []
        counts = {
            SSGSEA_STATUS_COMPUTED: 0,
            SSGSEA_STATUS_INSUFFICIENT_SUBSTRATES: 0,
            SSGSEA_STATUS_NO_FINITE_BACKGROUND_VALUES: 0,
            SSGSEA_STATUS_INSUFFICIENT_BACKGROUND_SITES: 0,
            SSGSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES: 0,
        }
        null_score_cache = _SsgseaNullScoreCache()

        for profile_position, profile_id in enumerate(profile_index):
            profile_values = effect_values[:, profile_position]
            finite_positions = np.flatnonzero(np.isfinite(profile_values))
            ranked_blocks = _rank_site_blocks(
                site_labels=site_labels,
                values=profile_values,
                finite_positions=finite_positions,
                ranking_direction=str(self.ranking_direction),
            )
            n_background = ranked_blocks.n_background
            ranked_position_by_site = {
                str(site_id): int(position)
                for position, site_id in enumerate(ranked_blocks.site_labels.tolist())
            }

            for kinase_position, kinase_name in enumerate(kinase_index):
                substrate_sites = membership_by_kinase[str(kinase_name)]
                hit_positions = np.fromiter(
                    (
                        ranked_position_by_site[site_id]
                        for site_id in substrate_sites
                        if site_id in ranked_position_by_site
                    ),
                    dtype=np.int64,
                )
                hit_mask = np.zeros(n_background, dtype=bool)
                hit_mask[hit_positions] = True
                n_substrates = int(hit_mask.sum())
                tie_diagnostics = _build_kinase_tie_diagnostics(
                    ranked_blocks=ranked_blocks,
                    hit_mask=hit_mask,
                )
                substrate_count_table.iat[kinase_position, profile_position] = (
                    n_substrates
                )

                status, reason = _resolve_status(
                    n_background=n_background,
                    n_substrates=n_substrates,
                    min_substrates=int(self.min_substrates),
                )
                enrichment_score = np.nan
                p_value = np.nan
                if status == SSGSEA_STATUS_COMPUTED:
                    null_score_engine = null_score_cache.get(
                        ranked_blocks=ranked_blocks,
                        n_substrates=n_substrates,
                    )
                    enrichment_score = null_score_engine.score_selected_positions(
                        hit_positions,
                    )
                    if not np.isfinite(enrichment_score):
                        status = SSGSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES
                        reason = "enrichment score is not finite"
                    else:
                        activity_scores.iat[
                            kinase_position,
                            profile_position,
                        ] = float(enrichment_score)
                        if random_seed is not None and p_value_matrix is not None:
                            permutation_rng = _make_ssgsea_permutation_rng(
                                random_seed=int(random_seed),
                                profile_id=str(profile_id),
                                kinase_name=str(kinase_name),
                            )
                            p_value = _permutation_p_value(
                                observed_score=float(enrichment_score),
                                null_score_engine=null_score_engine,
                                permutation_count=int(self.permutation_count),
                                rng=permutation_rng,
                            )
                            p_value_matrix.iat[
                                kinase_position,
                                profile_position,
                            ] = float(p_value)

                counts[status] += 1
                significance_status = _resolve_significance_status(
                    computability_status=status,
                    permutation_count=int(self.permutation_count),
                    adjust_p_values=bool(self.adjust_p_values),
                )
                rows.append(
                    {
                        "kinase": str(kinase_name),
                        "profile_id": str(profile_id),
                        "z_score": np.nan,
                        "enrichment_score": (
                            float(enrichment_score)
                            if np.isfinite(enrichment_score)
                            else np.nan
                        ),
                        "p_value": float(p_value) if np.isfinite(p_value) else np.nan,
                        "q_value": np.nan,
                        "significance_status": significance_status,
                        "n_substrates": int(n_substrates),
                        "n_background_sites": int(n_background),
                        "evidence_threshold": np.nan,
                        "evidence_threshold_operator": "",
                        "evidence_threshold_description": (
                            "not thresholded; explicit kinase-substrate membership"
                        ),
                        "min_substrates": int(self.min_substrates),
                        "ranking_direction": str(self.ranking_direction),
                        "tie_policy": SSGSEA_TIE_POLICY,
                        "n_tie_blocks": ranked_blocks.n_tie_blocks,
                        "n_tied_sites": ranked_blocks.n_tied_sites,
                        "max_tie_block_size": ranked_blocks.max_tie_block_size,
                        "substrate_only_tie_blocks": tie_diagnostics[
                            "substrate_only_tie_blocks"
                        ],
                        "non_substrate_only_tie_blocks": tie_diagnostics[
                            "non_substrate_only_tie_blocks"
                        ],
                        "mixed_substrate_tie_blocks": tie_diagnostics[
                            "mixed_substrate_tie_blocks"
                        ],
                        "permutation_count": int(self.permutation_count),
                        "random_seed": (
                            np.nan
                            if self.random_seed is None
                            else int(self.random_seed)
                        ),
                        "computability_status": status,
                        "reason": reason,
                    }
                )

        statistics_table = pd.DataFrame.from_records(
            rows,
            columns=[
                "kinase",
                "profile_id",
                "z_score",
                "enrichment_score",
                "p_value",
                "q_value",
                "significance_status",
                "n_substrates",
                "n_background_sites",
                "evidence_threshold",
                "evidence_threshold_operator",
                "evidence_threshold_description",
                "min_substrates",
                "ranking_direction",
                "tie_policy",
                "n_tie_blocks",
                "n_tied_sites",
                "max_tie_block_size",
                "substrate_only_tie_blocks",
                "non_substrate_only_tie_blocks",
                "mixed_substrate_tie_blocks",
                "permutation_count",
                "random_seed",
                "computability_status",
                "reason",
            ],
        )
        if p_value_matrix is not None and q_value_matrix is not None:
            _apply_q_values(
                statistics_table=statistics_table,
                q_value_matrix=q_value_matrix,
                profile_index=profile_index,
            )

        target_counts = _build_target_counts(
            aligned_membership=aligned_membership,
            kinase_index=kinase_index,
        )
        thresholded_substrate_counts = target_counts.rename("n_substrates")
        target_table = _build_target_table(aligned_membership=aligned_membership)
        summary = ActivityMethodSummary(
            kinases_evaluated=int(len(kinase_index)),
            kinase_profile_pairs_evaluated=int(len(kinase_index) * len(profile_index)),
            kinase_profile_pairs_computed=counts[SSGSEA_STATUS_COMPUTED],
            kinase_profile_pairs_insufficient_substrates=counts[
                SSGSEA_STATUS_INSUFFICIENT_SUBSTRATES
            ],
            kinase_profile_pairs_invalid_background_variance=(
                counts[SSGSEA_STATUS_NO_FINITE_BACKGROUND_VALUES]
                + counts[SSGSEA_STATUS_INSUFFICIENT_BACKGROUND_SITES]
            ),
            kinase_profile_pairs_no_finite_background_values=counts[
                SSGSEA_STATUS_NO_FINITE_BACKGROUND_VALUES
            ],
            kinase_profile_pairs_no_finite_substrate_values=counts[
                SSGSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES
            ],
        )
        diagnostics = SsgseaSubstrateEnrichmentActivityDiagnostics(
            method_summary=summary,
            threshold_membership_diagnostics=None,
            statistics_table=statistics_table,
        )
        policy = build_ssgsea_substrate_enrichment_activity_policy(
            min_substrates=int(self.min_substrates),
            ranking_direction=str(self.ranking_direction),
            permutation_count=int(self.permutation_count),
            random_seed=self.random_seed,
            adjust_p_values=bool(self.adjust_p_values),
            q_value_method=(
                str(self.q_value_method)
                if int(self.permutation_count) > 0 and bool(self.adjust_p_values)
                else None
            ),
        )

        return KinaseActivityResult._from_owned(
            weighted_activity=activity_scores,
            p_value_matrix=p_value_matrix,
            q_value_matrix=q_value_matrix,
            thresholded_substrate_counts=thresholded_substrate_counts,
            activity_substrate_counts=substrate_count_table,
            substrate_count_matrix=substrate_count_table,
            target_counts=target_counts,
            target_table=target_table,
            statistics_table=statistics_table,
            method_summary=summary,
            method_diagnostics=diagnostics,
            policy_provenance=(policy,),
            activity_method=SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_METHOD,
            input_semantics=activity_input.semantics,
            profile_metadata=activity_input.profile_metadata,
        )


def _validate_membership_table(value: pd.DataFrame) -> pd.DataFrame:
    frame = require_dataframe(
        value,
        field_name="activity_inputs.kinase_substrate_membership",
        allow_empty=False,
        error_type=WorkflowBoundaryError,
    )
    require_columns(
        frame,
        field_name="activity_inputs.kinase_substrate_membership",
        required_columns=(_KINASE_COLUMN, _SUBSTRATE_COLUMN),
        error_type=WorkflowBoundaryError,
    )
    require_non_empty_string_column(
        frame,
        field_name="activity_inputs.kinase_substrate_membership",
        column_name=_KINASE_COLUMN,
        error_type=WorkflowBoundaryError,
    )
    require_canonical_string_column(
        frame,
        field_name="activity_inputs.kinase_substrate_membership",
        column_name=_KINASE_COLUMN,
        error_type=WorkflowBoundaryError,
    )
    require_non_empty_string_column(
        frame,
        field_name="activity_inputs.kinase_substrate_membership",
        column_name=_SUBSTRATE_COLUMN,
        error_type=WorkflowBoundaryError,
    )
    require_unique_row_pairs(
        frame,
        field_name="activity_inputs.kinase_substrate_membership",
        column_names=(_KINASE_COLUMN, _SUBSTRATE_COLUMN),
        error_type=WorkflowBoundaryError,
    )
    require_unique_index(
        frame.reset_index(drop=True),
        field_name="activity_inputs.kinase_substrate_membership",
        error_type=WorkflowBoundaryError,
    )
    require_unique_columns(
        frame,
        field_name="activity_inputs.kinase_substrate_membership",
        error_type=WorkflowBoundaryError,
    )
    return frame.copy(deep=True)


def _ordered_unique_strings(values: pd.Series) -> list[str]:
    return list(dict.fromkeys(str(value) for value in values.tolist()))


def _build_membership_lookup(
    *,
    kinases: list[str],
    membership: pd.DataFrame,
) -> dict[str, set[str]]:
    result = {kinase: set() for kinase in kinases}
    if membership.empty:
        return result
    for kinase, substrate_site in membership.loc[
        :,
        [_KINASE_COLUMN, _SUBSTRATE_COLUMN],
    ].itertuples(index=False):
        result.setdefault(str(kinase), set()).add(str(substrate_site))
    return result


def _rank_site_blocks(
    *,
    site_labels: np.ndarray,
    values: np.ndarray,
    finite_positions: np.ndarray,
    ranking_direction: str,
) -> _RankedTieBlocks:
    if finite_positions.size == 0:
        return _RankedTieBlocks(
            site_labels=np.asarray([], dtype=object),
            block_starts=np.asarray([], dtype=np.int64),
            block_sizes=np.asarray([], dtype=np.int64),
        )
    finite_values = values[finite_positions]
    # Sorting establishes only the order of distinct value blocks. The order of
    # rows inside an equal-value block is intentionally ignored by
    # _score_from_ranked_hit_mask.
    if ranking_direction == SSGSEA_RANKING_DIRECTION_ASCENDING:
        order = np.argsort(finite_values, kind="mergesort")
    else:
        order = np.argsort(-finite_values, kind="mergesort")
    ranked_positions = finite_positions[order]
    ranked_values = values[ranked_positions]
    block_starts = np.concatenate(
        (
            np.asarray([0], dtype=np.int64),
            np.flatnonzero(ranked_values[1:] != ranked_values[:-1]).astype(np.int64)
            + 1,
        )
    )
    block_ends = np.concatenate(
        (
            block_starts[1:],
            np.asarray([ranked_values.size], dtype=np.int64),
        )
    )
    return _RankedTieBlocks(
        site_labels=site_labels[ranked_positions],
        block_starts=block_starts,
        block_sizes=(block_ends - block_starts).astype(np.int64),
    )


def _resolve_status(
    *,
    n_background: int,
    n_substrates: int,
    min_substrates: int,
) -> tuple[str, str]:
    if n_background == 0:
        return (
            SSGSEA_STATUS_NO_FINITE_BACKGROUND_VALUES,
            "background requires at least one finite phosphosite effect value",
        )
    if n_substrates < int(min_substrates):
        return (
            SSGSEA_STATUS_INSUFFICIENT_SUBSTRATES,
            f"n_substrates={n_substrates} is below min_substrates={min_substrates}",
        )
    if n_background - n_substrates < 1:
        return (
            SSGSEA_STATUS_INSUFFICIENT_BACKGROUND_SITES,
            "enrichment requires at least one finite non-substrate background site",
        )
    return SSGSEA_STATUS_COMPUTED, ""


def _resolve_significance_status(
    *,
    computability_status: str,
    permutation_count: int,
    adjust_p_values: bool,
) -> str:
    if int(permutation_count) <= 0:
        return SSGSEA_SIGNIFICANCE_STATUS_UNAVAILABLE_NO_PERMUTATIONS
    if computability_status != SSGSEA_STATUS_COMPUTED:
        return SSGSEA_SIGNIFICANCE_STATUS_UNAVAILABLE_NOT_COMPUTABLE
    if not bool(adjust_p_values):
        return SSGSEA_SIGNIFICANCE_STATUS_P_VALUE_AVAILABLE_Q_VALUE_DISABLED
    return SSGSEA_SIGNIFICANCE_STATUS_AVAILABLE


def _score_from_hit_mask(hit_mask: np.ndarray) -> float:
    n_background = int(hit_mask.size)
    n_substrates = int(hit_mask.sum())
    n_misses = n_background - n_substrates
    if n_background == 0 or n_substrates == 0 or n_misses == 0:
        return np.nan
    hit_increment = 1.0 / float(n_substrates)
    miss_increment = 1.0 / float(n_misses)
    running = np.cumsum(
        np.where(hit_mask, hit_increment, -miss_increment).astype(float),
    )
    return float(running.sum() / float(n_background))


def _score_from_ranked_hit_mask(
    *,
    hit_mask: np.ndarray,
    ranked_blocks: _RankedTieBlocks,
) -> float:
    if not ranked_blocks.has_ties:
        return _score_from_hit_mask(hit_mask)
    return _score_from_block_hit_counts(
        block_hit_counts=_block_hit_counts(
            ranked_blocks=ranked_blocks,
            hit_mask=hit_mask,
        ),
        block_sizes=ranked_blocks.block_sizes,
    )


def _score_from_block_hit_counts(
    *,
    block_hit_counts: np.ndarray,
    block_sizes: np.ndarray,
) -> float:
    n_background = int(block_sizes.sum())
    n_substrates = int(block_hit_counts.sum())
    n_misses = n_background - n_substrates
    if n_background == 0 or n_substrates == 0 or n_misses == 0:
        return np.nan

    hit_increment = 1.0 / float(n_substrates)
    miss_increment = 1.0 / float(n_misses)
    running = 0.0
    area = 0.0
    for block_hits_raw, block_size_raw in zip(
        block_hit_counts.tolist(),
        block_sizes.tolist(),
        strict=True,
    ):
        block_size = int(block_size_raw)
        block_hits = int(block_hits_raw)
        block_misses = block_size - block_hits
        block_delta = (
            float(block_hits) * hit_increment - float(block_misses) * miss_increment
        )
        area += (
            float(block_size) * running
            + ((float(block_size) + 1.0) / 2.0) * block_delta
        )
        running += block_delta
    return float(area / float(n_background))


def _block_hit_counts(
    *,
    ranked_blocks: _RankedTieBlocks,
    hit_mask: np.ndarray,
) -> np.ndarray:
    if ranked_blocks.n_blocks == 0:
        return np.asarray([], dtype=np.int64)
    return np.add.reduceat(
        hit_mask.astype(np.int64, copy=False),
        ranked_blocks.block_starts,
    ).astype(np.int64, copy=False)


def _null_score_cache_key(
    *,
    ranked_blocks: _RankedTieBlocks,
    n_substrates: int,
) -> _SsgseaNullScoreCacheKey:
    return _SsgseaNullScoreCacheKey(
        n_background=int(ranked_blocks.n_background),
        n_substrates=int(n_substrates),
        tie_block_sizes=(
            tuple(int(value) for value in ranked_blocks.block_sizes.tolist())
            if ranked_blocks.has_ties
            else ()
        ),
    )


def _build_kinase_tie_diagnostics(
    *,
    ranked_blocks: _RankedTieBlocks,
    hit_mask: np.ndarray,
) -> dict[str, int]:
    if not ranked_blocks.has_ties:
        return {
            "substrate_only_tie_blocks": 0,
            "non_substrate_only_tie_blocks": 0,
            "mixed_substrate_tie_blocks": 0,
        }
    block_hit_counts = _block_hit_counts(
        ranked_blocks=ranked_blocks,
        hit_mask=hit_mask,
    )
    tied = ranked_blocks.block_sizes > 1
    tied_block_sizes = ranked_blocks.block_sizes[tied]
    tied_hit_counts = block_hit_counts[tied]
    substrate_only = tied_hit_counts == tied_block_sizes
    non_substrate_only = tied_hit_counts == 0
    mixed = (tied_hit_counts > 0) & (tied_hit_counts < tied_block_sizes)
    return {
        "substrate_only_tie_blocks": int(substrate_only.sum()),
        "non_substrate_only_tie_blocks": int(non_substrate_only.sum()),
        "mixed_substrate_tie_blocks": int(mixed.sum()),
    }


def _permutation_p_value(
    *,
    observed_score: float,
    null_score_engine: _SsgseaNullScoreEngine,
    permutation_count: int,
    rng: np.random.Generator,
) -> float:
    extreme_count = 0
    observed_abs = abs(float(observed_score))
    n_background = int(null_score_engine.n_background)
    n_substrates = int(null_score_engine.n_substrates)
    for _ in range(int(permutation_count)):
        hit_positions = rng.choice(
            int(n_background),
            size=int(n_substrates),
            replace=False,
        )
        score = null_score_engine.score_selected_positions(hit_positions)
        score_abs = abs(float(score))
        if score_abs > observed_abs or np.isclose(
            score_abs,
            observed_abs,
            rtol=0.0,
            atol=_SSGSEA_PERMUTATION_EXTREME_ATOL,
        ):
            extreme_count += 1
    return float((extreme_count + 1) / (int(permutation_count) + 1))


def _make_ssgsea_permutation_rng(
    *,
    random_seed: int,
    profile_id: str,
    kinase_name: str,
) -> np.random.Generator:
    return np.random.default_rng(
        _derive_ssgsea_permutation_seed(
            random_seed=int(random_seed),
            profile_id=str(profile_id),
            kinase_name=str(kinase_name),
        )
    )


def _derive_ssgsea_permutation_seed(
    *,
    random_seed: int,
    profile_id: str | None = None,
    condition_name: str | None = None,
    kinase_name: str,
) -> int:
    if profile_id is not None and condition_name is not None:
        if str(profile_id) != str(condition_name):
            raise ValueError("profile_id conflicts with legacy condition_name")
    resolved_profile_id = (
        str(profile_id)
        if profile_id is not None
        else ""
        if condition_name is None
        else str(condition_name)
    )
    seed_material = {
        "kinase": str(kinase_name),
        "method_id": SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_METHOD.activity_method_id,
        "method_version": SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_POLICY_VERSION,
        "profile_id": resolved_profile_id,
        "random_seed": int(random_seed),
        "seed_policy": _SSGSEA_PERMUTATION_RNG_SEED_HASH_POLICY_TOKEN,
        "seed_policy_version": SSGSEA_PERMUTATION_RNG_SEED_POLICY_VERSION,
        "stream": _SSGSEA_PERMUTATION_STREAM_NAME,
    }
    encoded = json.dumps(
        seed_material,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    digest = hashlib.blake2b(
        encoded,
        digest_size=_SSGSEA_PERMUTATION_SEED_DIGEST_SIZE_BYTES,
    ).digest()
    return int.from_bytes(digest, "little")


def _apply_q_values(
    *,
    statistics_table: pd.DataFrame,
    q_value_matrix: pd.DataFrame,
    profile_index: pd.Index,
) -> None:
    for profile_id in profile_index:
        profile_mask = statistics_table.loc[:, "profile_id"].astype(str) == str(
            profile_id
        )
        computed_mask = (
            statistics_table.loc[:, "computability_status"] == SSGSEA_STATUS_COMPUTED
        )
        selected = profile_mask & computed_mask
        if not bool(selected.any()):
            continue
        profile_p_values = statistics_table.loc[selected, "p_value"].astype(float)
        q_values = benjamini_hochberg_q_values(profile_p_values)
        statistics_table.loc[selected, "q_value"] = q_values.to_numpy(
            dtype=float,
            copy=False,
        )
        profile_rows = statistics_table.loc[selected, "kinase"].astype(str)
        for kinase_name, q_value in zip(
            profile_rows.tolist(),
            q_values.to_numpy(dtype=float, copy=False).tolist(),
            strict=True,
        ):
            q_value_matrix.at[str(kinase_name), str(profile_id)] = float(q_value)


def _build_target_counts(
    *,
    aligned_membership: pd.DataFrame,
    kinase_index: pd.Index,
) -> pd.Series:
    if aligned_membership.empty:
        counts = pd.Series(0, index=kinase_index.copy(), name="n_targets")
    else:
        counts = (
            aligned_membership.groupby(_KINASE_COLUMN, sort=False)[_SUBSTRATE_COLUMN]
            .nunique()
            .reindex(kinase_index, fill_value=0)
            .astype("int64")
            .rename("n_targets")
        )
    counts.index.name = "kinase"
    return counts.sort_values(ascending=False)


def _build_target_table(*, aligned_membership: pd.DataFrame) -> pd.DataFrame:
    if aligned_membership.empty:
        return pd.DataFrame(columns=["site_id", "kinase", "score"])
    target_table = aligned_membership.loc[
        :,
        [_SUBSTRATE_COLUMN, _KINASE_COLUMN],
    ].copy(deep=True)
    target_table.columns = ["site_id", "kinase"]
    target_table.loc[:, "score"] = 1.0
    return target_table.sort_values(["kinase", "site_id"], kind="mergesort")


__all__ = [
    "SSGSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG",
    "SSGSEA_RANKING_DIRECTION_ASCENDING",
    "SSGSEA_RANKING_DIRECTION_DESCENDING",
    "SSGSEA_RANKING_DIRECTIONS",
    "SSGSEA_SIGNIFICANCE_STATUS_AVAILABLE",
    "SSGSEA_SIGNIFICANCE_STATUS_P_VALUE_AVAILABLE_Q_VALUE_DISABLED",
    "SSGSEA_SIGNIFICANCE_STATUS_UNAVAILABLE_NO_PERMUTATIONS",
    "SSGSEA_SIGNIFICANCE_STATUS_UNAVAILABLE_NOT_COMPUTABLE",
    "SSGSEA_STATUS_COMPUTED",
    "SSGSEA_STATUS_INSUFFICIENT_BACKGROUND_SITES",
    "SSGSEA_STATUS_INSUFFICIENT_SUBSTRATES",
    "SSGSEA_STATUS_NO_FINITE_BACKGROUND_VALUES",
    "SSGSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES",
    "SsgseaSubstrateEnrichmentActivityMethod",
]
