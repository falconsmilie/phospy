"""ssGSEA-style substrate-set enrichment activity-like score method."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
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
# Historical v1 salt retained solely to preserve the existing deterministic
# ssGSEA permutation stream. The active public policy name is
# ``stable_by_method_profile_kinase`` and is emitted from
# ``scientific_policies.py``.
_SSGSEA_PERMUTATION_RNG_SEED_V1_COMPATIBILITY_TOKEN = (
    "stable_by_method_condition_kinase"
)
_SSGSEA_PERMUTATION_EXTREME_ATOL = 1e-12
_BoolArray = npt.NDArray[np.bool_]
_FloatArray = npt.NDArray[np.float64]
_IntArray = npt.NDArray[np.int64]
_ObjectArray = npt.NDArray[np.object_]


def _int_tuple(values: _IntArray) -> tuple[int, ...]:
    return tuple(int(values[index]) for index in range(int(values.size)))


def _object_tuple(values: _ObjectArray) -> tuple[object, ...]:
    return tuple(values[index] for index in range(int(values.size)))


@dataclass(frozen=True, slots=True)
class _RankedTieBlocks:
    site_labels: _ObjectArray
    block_starts: _IntArray
    block_sizes: _IntArray

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
    rank_weights: _FloatArray | None
    rank_weight_multiplier: float
    rank_weight_constant: float
    block_index_by_position: _IntArray | None
    block_hit_coefficients: _FloatArray | None
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
            rank_weights: _FloatArray = (
                n_background - np.arange(n_background, dtype=float)
            ).astype(
                float,
                copy=False,
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

        block_sizes: _FloatArray = ranked_blocks.block_sizes.astype(float, copy=False)
        suffix_site_counts: _FloatArray = (
            n_background - np.cumsum(ranked_blocks.block_sizes)
        ).astype(float, copy=False)
        block_area_coefficients: _FloatArray = suffix_site_counts + (
            (block_sizes + 1.0) / 2.0
        )
        block_index_by_position: _IntArray = np.empty(n_background, dtype=np.int64)
        block_starts = _int_tuple(ranked_blocks.block_starts)
        integer_block_sizes = _int_tuple(ranked_blocks.block_sizes)
        for block_index, block_start, block_size in zip(
            range(ranked_blocks.n_blocks),
            block_starts,
            integer_block_sizes,
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

    def score_selected_positions(self, hit_positions: _IntArray) -> float:
        positions: _IntArray = np.asarray(hit_positions, dtype=np.int64)
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
        block_indices: _IntArray = self.block_index_by_position[positions].astype(
            np.int64,
            copy=False,
        )
        block_hit_counts: _IntArray = np.bincount(
            block_indices,
            minlength=int(self.block_hit_coefficients.size),
        ).astype(np.int64, copy=False)
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
class _SsgseaPreparedInputs:
    activity_input: ActivityInputMatrix
    aligned_membership: pd.DataFrame
    kinase_index: pd.Index
    profile_index: pd.Index
    site_labels: _ObjectArray
    membership_by_kinase: dict[str, set[str]]
    effect_values: _FloatArray
    random_seed: int | None


@dataclass(frozen=True, slots=True)
class _SsgseaMutableScoringTables:
    activity_scores: pd.DataFrame
    substrate_count_table: pd.DataFrame
    p_value_matrix: pd.DataFrame | None
    q_value_matrix: pd.DataFrame | None
    rows: list[dict[str, object]]
    status_counts: dict[str, int]
    null_score_cache: _SsgseaNullScoreCache


@dataclass(frozen=True, slots=True)
class _SsgseaProfileRanking:
    position: int
    profile_id: str
    ranked_blocks: _RankedTieBlocks
    ranked_position_by_site: dict[str, int]


@dataclass(frozen=True, slots=True)
class _SsgseaPairScore:
    n_substrates: int
    enrichment_score: float
    p_value: float
    status: str
    reason: str
    tie_diagnostics: dict[str, int]


@dataclass(frozen=True, slots=True)
class _SsgseaScoredOutputs:
    activity_scores: pd.DataFrame
    substrate_count_table: pd.DataFrame
    p_value_matrix: pd.DataFrame | None
    q_value_matrix: pd.DataFrame | None
    statistics_table: pd.DataFrame
    status_counts: dict[str, int]


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
        min_substrates = _require_int(
            self.min_substrates,
            field_name="ssgsea min_substrates",
        )
        if min_substrates < 1:
            raise ValueError("ssgsea min_substrates must be greater than or equal to 1")
        if self.ranking_direction not in SSGSEA_RANKING_DIRECTIONS:
            allowed = ", ".join(sorted(SSGSEA_RANKING_DIRECTIONS))
            raise ValueError(f"ssgsea ranking_direction must be one of: {allowed}")
        permutation_count = _require_int(
            self.permutation_count,
            field_name="ssgsea permutation_count",
        )
        if permutation_count < 0:
            raise ValueError(
                "ssgsea permutation_count must be greater than or equal to 0"
            )
        random_seed = _require_optional_int(
            self.random_seed,
            field_name="ssgsea random_seed",
        )
        if random_seed is not None and random_seed < 0:
            raise ValueError("ssgsea random_seed must be greater than or equal to 0")
        if permutation_count > 0 and random_seed is None:
            raise ValueError(
                "ssgsea random_seed must be set when permutation_count is positive"
            )
        adjust_p_values = _require_bool(
            self.adjust_p_values,
            field_name="ssgsea adjust_p_values",
        )
        if self.q_value_method != SSGSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG:
            raise ValueError("ssgsea q_value_method must be 'benjamini_hochberg'")
        object.__setattr__(self, "min_substrates", min_substrates)
        object.__setattr__(self, "permutation_count", permutation_count)
        object.__setattr__(self, "random_seed", random_seed)
        object.__setattr__(self, "adjust_p_values", adjust_p_values)

    def run(
        self,
        *,
        activity_input: ActivityInputMatrix | None = None,
        effect_matrix: pd.DataFrame | None = None,
        kinase_substrate_membership: pd.DataFrame,
    ) -> KinaseActivityResult:
        prepared = _prepare_ssgsea_inputs(
            activity_input=activity_input,
            effect_matrix=effect_matrix,
            kinase_substrate_membership=kinase_substrate_membership,
            method=self,
        )
        scored = _score_ssgsea_profiles(method=self, prepared=prepared)
        corrected = _apply_ssgsea_multiple_testing(
            prepared=prepared,
            scored=scored,
        )
        return _assemble_ssgsea_result(
            method=self,
            prepared=prepared,
            scored=corrected,
        )


def _prepare_ssgsea_inputs(
    *,
    method: SsgseaSubstrateEnrichmentActivityMethod,
    activity_input: ActivityInputMatrix | None,
    effect_matrix: pd.DataFrame | None,
    kinase_substrate_membership: pd.DataFrame,
) -> _SsgseaPreparedInputs:
    resolved_activity_input = _resolve_ssgsea_activity_input(
        activity_input=activity_input,
        effect_matrix=effect_matrix,
    )
    from phospy.science.quantitative_method_contracts import (
        resolve_activity_input_contract,
    )

    resolve_activity_input_contract(
        activity_input=resolved_activity_input,
        contract=ssgsea_substrate_enrichment_activity_input_contract(),
        context=(
            "ssGSEA substrate enrichment activity input requires explicit "
            "contrast/effect input"
        ),
    )
    effects = resolved_activity_input.frame
    membership = _validate_membership_table(kinase_substrate_membership)
    site_labels: _ObjectArray = np.asarray(
        effects.index.astype(str).tolist(),
        dtype=object,
    )
    site_universe = set(str(value) for value in _object_tuple(site_labels))
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
    effect_values: _FloatArray = effects.to_numpy(dtype=float, copy=False)
    return _SsgseaPreparedInputs(
        activity_input=resolved_activity_input,
        aligned_membership=aligned_membership,
        kinase_index=kinase_index,
        profile_index=profile_index,
        site_labels=site_labels,
        membership_by_kinase=_build_membership_lookup(
            kinases=kinases,
            membership=aligned_membership,
        ),
        effect_values=effect_values,
        random_seed=_resolve_ssgsea_random_seed(method),
    )


def _resolve_ssgsea_activity_input(
    *,
    activity_input: ActivityInputMatrix | None,
    effect_matrix: pd.DataFrame | None,
) -> ActivityInputMatrix:
    if activity_input is not None and effect_matrix is not None:
        raise WorkflowBoundaryError(
            "ssgsea activity requires either activity_input or legacy "
            "effect_matrix, not both"
        )
    if activity_input is not None:
        return activity_input
    if effect_matrix is None:
        raise WorkflowBoundaryError(
            "ssgsea activity requires ActivityInputMatrix with explicit "
            "contrast/effect semantics"
        )
    return normalize_activity_input_matrix(
        effect_matrix,
        field_name="activity_inputs.effect_matrix",
        legacy_dataframe_semantics=ActivityInputSemantics(
            profile_axis=ActivityProfileAxis.EFFECT,
            quantitative_semantics=ActivityQuantitativeSemantics.STANDARDISED_EFFECT,
        ),
        legacy_dataframe_warning=(
            "Passing a raw DataFrame as effect_matrix is deprecated; "
            "provide ActivityInputMatrix.contrast_log_fold_change(...) or "
            "ActivityInputMatrix.standardised_effect(...) so activity "
            "input semantics are explicit."
        ),
    )


def _resolve_ssgsea_random_seed(
    method: SsgseaSubstrateEnrichmentActivityMethod,
) -> int | None:
    if int(method.permutation_count) <= 0:
        return None
    if method.random_seed is None:
        raise ValueError(
            "ssgsea random_seed must be set when permutation_count is positive"
        )
    return int(method.random_seed)


def _score_ssgsea_profiles(
    *,
    method: SsgseaSubstrateEnrichmentActivityMethod,
    prepared: _SsgseaPreparedInputs,
) -> _SsgseaScoredOutputs:
    tables = _initialise_ssgsea_scoring_tables(method=method, prepared=prepared)
    for profile_position, profile_id in enumerate(prepared.profile_index):
        ranking = _rank_ssgsea_profile(
            method=method,
            prepared=prepared,
            profile_position=profile_position,
            profile_id=str(profile_id),
        )
        for kinase_position, kinase_name in enumerate(prepared.kinase_index):
            score = _score_ssgsea_kinase_profile(
                method=method,
                prepared=prepared,
                ranking=ranking,
                tables=tables,
                kinase_name=str(kinase_name),
            )
            _record_ssgsea_pair_score(
                method=method,
                prepared=prepared,
                ranking=ranking,
                tables=tables,
                score=score,
                kinase_position=kinase_position,
                kinase_name=str(kinase_name),
            )
    return _SsgseaScoredOutputs(
        activity_scores=tables.activity_scores,
        substrate_count_table=tables.substrate_count_table,
        p_value_matrix=tables.p_value_matrix,
        q_value_matrix=tables.q_value_matrix,
        statistics_table=pd.DataFrame.from_records(
            tables.rows,
            columns=_ssgsea_statistics_columns(),
        ),
        status_counts=tables.status_counts,
    )


def _initialise_ssgsea_scoring_tables(
    *,
    method: SsgseaSubstrateEnrichmentActivityMethod,
    prepared: _SsgseaPreparedInputs,
) -> _SsgseaMutableScoringTables:
    p_value_matrix = (
        pd.DataFrame(
            np.nan,
            index=prepared.kinase_index,
            columns=prepared.profile_index,
            dtype=float,
        )
        if int(method.permutation_count) > 0
        else None
    )
    q_value_matrix = (
        pd.DataFrame(
            np.nan,
            index=prepared.kinase_index,
            columns=prepared.profile_index,
            dtype=float,
        )
        if int(method.permutation_count) > 0 and bool(method.adjust_p_values)
        else None
    )
    return _SsgseaMutableScoringTables(
        activity_scores=pd.DataFrame(
            np.nan,
            index=prepared.kinase_index,
            columns=prepared.profile_index,
            dtype=float,
        ),
        substrate_count_table=pd.DataFrame(
            0,
            index=prepared.kinase_index,
            columns=prepared.profile_index,
            dtype=int,
        ),
        p_value_matrix=p_value_matrix,
        q_value_matrix=q_value_matrix,
        rows=[],
        status_counts={
            SSGSEA_STATUS_COMPUTED: 0,
            SSGSEA_STATUS_INSUFFICIENT_SUBSTRATES: 0,
            SSGSEA_STATUS_NO_FINITE_BACKGROUND_VALUES: 0,
            SSGSEA_STATUS_INSUFFICIENT_BACKGROUND_SITES: 0,
            SSGSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES: 0,
        },
        null_score_cache=_SsgseaNullScoreCache(),
    )


def _rank_ssgsea_profile(
    *,
    method: SsgseaSubstrateEnrichmentActivityMethod,
    prepared: _SsgseaPreparedInputs,
    profile_position: int,
    profile_id: str,
) -> _SsgseaProfileRanking:
    profile_values: _FloatArray = prepared.effect_values[:, profile_position]
    finite_positions: _IntArray = np.flatnonzero(np.isfinite(profile_values)).astype(
        np.int64,
        copy=False,
    )
    ranked_blocks = _rank_site_blocks(
        site_labels=prepared.site_labels,
        values=profile_values,
        finite_positions=finite_positions,
        ranking_direction=str(method.ranking_direction),
    )
    return _SsgseaProfileRanking(
        position=profile_position,
        profile_id=profile_id,
        ranked_blocks=ranked_blocks,
        ranked_position_by_site={
            str(site_id): int(position)
            for position, site_id in enumerate(_object_tuple(ranked_blocks.site_labels))
        },
    )


def _score_ssgsea_kinase_profile(
    *,
    method: SsgseaSubstrateEnrichmentActivityMethod,
    prepared: _SsgseaPreparedInputs,
    ranking: _SsgseaProfileRanking,
    tables: _SsgseaMutableScoringTables,
    kinase_name: str,
) -> _SsgseaPairScore:
    hit_positions: _IntArray = np.fromiter(
        (
            ranking.ranked_position_by_site[site_id]
            for site_id in prepared.membership_by_kinase[kinase_name]
            if site_id in ranking.ranked_position_by_site
        ),
        dtype=np.int64,
    )
    hit_mask: _BoolArray = np.zeros(ranking.ranked_blocks.n_background, dtype=bool)
    hit_mask[hit_positions] = True
    n_substrates = int(hit_mask.sum())
    status, reason = _resolve_status(
        n_background=ranking.ranked_blocks.n_background,
        n_substrates=n_substrates,
        min_substrates=int(method.min_substrates),
    )
    enrichment_score = np.nan
    p_value = np.nan
    if status == SSGSEA_STATUS_COMPUTED:
        null_score_engine = tables.null_score_cache.get(
            ranked_blocks=ranking.ranked_blocks,
            n_substrates=n_substrates,
        )
        enrichment_score = null_score_engine.score_selected_positions(hit_positions)
        if not np.isfinite(enrichment_score):
            status = SSGSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES
            reason = "enrichment score is not finite"
        else:
            if prepared.random_seed is not None and tables.p_value_matrix is not None:
                p_value = _ssgsea_null_p_value(
                    method=method,
                    prepared=prepared,
                    ranking=ranking,
                    kinase_name=kinase_name,
                    observed_score=float(enrichment_score),
                    null_score_engine=null_score_engine,
                )
    return _SsgseaPairScore(
        n_substrates=n_substrates,
        enrichment_score=(
            float(enrichment_score) if np.isfinite(enrichment_score) else np.nan
        ),
        p_value=float(p_value) if np.isfinite(p_value) else np.nan,
        status=status,
        reason=reason,
        tie_diagnostics=_build_kinase_tie_diagnostics(
            ranked_blocks=ranking.ranked_blocks,
            hit_mask=hit_mask,
        ),
    )


def _ssgsea_null_p_value(
    *,
    method: SsgseaSubstrateEnrichmentActivityMethod,
    prepared: _SsgseaPreparedInputs,
    ranking: _SsgseaProfileRanking,
    kinase_name: str,
    observed_score: float,
    null_score_engine: _SsgseaNullScoreEngine,
) -> float:
    if prepared.random_seed is None:
        return np.nan
    permutation_rng = _make_ssgsea_permutation_rng(
        random_seed=int(prepared.random_seed),
        profile_id=ranking.profile_id,
        kinase_name=kinase_name,
    )
    return _permutation_p_value(
        observed_score=float(observed_score),
        null_score_engine=null_score_engine,
        permutation_count=int(method.permutation_count),
        rng=permutation_rng,
    )


def _record_ssgsea_pair_score(
    *,
    method: SsgseaSubstrateEnrichmentActivityMethod,
    prepared: _SsgseaPreparedInputs,
    ranking: _SsgseaProfileRanking,
    tables: _SsgseaMutableScoringTables,
    score: _SsgseaPairScore,
    kinase_position: int,
    kinase_name: str,
) -> None:
    tables.substrate_count_table.iat[kinase_position, ranking.position] = (
        score.n_substrates
    )
    if score.status == SSGSEA_STATUS_COMPUTED:
        tables.activity_scores.iat[kinase_position, ranking.position] = (
            score.enrichment_score
        )
        if tables.p_value_matrix is not None:
            tables.p_value_matrix.iat[kinase_position, ranking.position] = score.p_value
    tables.status_counts[score.status] += 1
    tables.rows.append(
        _ssgsea_statistics_row(
            method=method,
            ranking=ranking,
            score=score,
            kinase_name=kinase_name,
        )
    )


def _ssgsea_statistics_row(
    *,
    method: SsgseaSubstrateEnrichmentActivityMethod,
    ranking: _SsgseaProfileRanking,
    score: _SsgseaPairScore,
    kinase_name: str,
) -> dict[str, object]:
    significance_status = _resolve_significance_status(
        computability_status=score.status,
        permutation_count=int(method.permutation_count),
        adjust_p_values=bool(method.adjust_p_values),
    )
    return {
        "kinase": kinase_name,
        "profile_id": ranking.profile_id,
        "z_score": np.nan,
        "enrichment_score": score.enrichment_score,
        "p_value": score.p_value,
        "q_value": np.nan,
        "significance_status": significance_status,
        "n_substrates": int(score.n_substrates),
        "n_background_sites": int(ranking.ranked_blocks.n_background),
        "evidence_threshold": np.nan,
        "evidence_threshold_operator": "",
        "evidence_threshold_description": (
            "not thresholded; explicit kinase-substrate membership"
        ),
        "min_substrates": int(method.min_substrates),
        "ranking_direction": str(method.ranking_direction),
        "tie_policy": SSGSEA_TIE_POLICY,
        "n_tie_blocks": ranking.ranked_blocks.n_tie_blocks,
        "n_tied_sites": ranking.ranked_blocks.n_tied_sites,
        "max_tie_block_size": ranking.ranked_blocks.max_tie_block_size,
        "substrate_only_tie_blocks": score.tie_diagnostics["substrate_only_tie_blocks"],
        "non_substrate_only_tie_blocks": score.tie_diagnostics[
            "non_substrate_only_tie_blocks"
        ],
        "mixed_substrate_tie_blocks": score.tie_diagnostics[
            "mixed_substrate_tie_blocks"
        ],
        "permutation_count": int(method.permutation_count),
        "random_seed": (
            np.nan if method.random_seed is None else int(method.random_seed)
        ),
        "computability_status": score.status,
        "reason": score.reason,
    }


def _ssgsea_statistics_columns() -> list[str]:
    return [
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
    ]


def _apply_ssgsea_multiple_testing(
    *,
    prepared: _SsgseaPreparedInputs,
    scored: _SsgseaScoredOutputs,
) -> _SsgseaScoredOutputs:
    if scored.p_value_matrix is not None and scored.q_value_matrix is not None:
        _apply_q_values(
            statistics_table=scored.statistics_table,
            q_value_matrix=scored.q_value_matrix,
            profile_index=prepared.profile_index,
        )
    return scored


def _assemble_ssgsea_result(
    *,
    method: SsgseaSubstrateEnrichmentActivityMethod,
    prepared: _SsgseaPreparedInputs,
    scored: _SsgseaScoredOutputs,
) -> KinaseActivityResult:
    target_counts = _build_target_counts(
        aligned_membership=prepared.aligned_membership,
        kinase_index=prepared.kinase_index,
    )
    thresholded_substrate_counts = target_counts.rename("n_substrates")
    target_table = _build_target_table(aligned_membership=prepared.aligned_membership)
    summary = ActivityMethodSummary(
        kinases_evaluated=int(len(prepared.kinase_index)),
        kinase_profile_pairs_evaluated=int(
            len(prepared.kinase_index) * len(prepared.profile_index)
        ),
        kinase_profile_pairs_computed=scored.status_counts[SSGSEA_STATUS_COMPUTED],
        kinase_profile_pairs_insufficient_substrates=scored.status_counts[
            SSGSEA_STATUS_INSUFFICIENT_SUBSTRATES
        ],
        kinase_profile_pairs_invalid_background_variance=(
            scored.status_counts[SSGSEA_STATUS_NO_FINITE_BACKGROUND_VALUES]
            + scored.status_counts[SSGSEA_STATUS_INSUFFICIENT_BACKGROUND_SITES]
        ),
        kinase_profile_pairs_no_finite_background_values=scored.status_counts[
            SSGSEA_STATUS_NO_FINITE_BACKGROUND_VALUES
        ],
        kinase_profile_pairs_no_finite_substrate_values=scored.status_counts[
            SSGSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES
        ],
    )
    diagnostics = SsgseaSubstrateEnrichmentActivityDiagnostics(
        method_summary=summary,
        threshold_membership_diagnostics=None,
        statistics_table=scored.statistics_table,
    )
    policy = build_ssgsea_substrate_enrichment_activity_policy(
        min_substrates=int(method.min_substrates),
        ranking_direction=str(method.ranking_direction),
        permutation_count=int(method.permutation_count),
        random_seed=method.random_seed,
        adjust_p_values=bool(method.adjust_p_values),
        q_value_method=(
            str(method.q_value_method)
            if int(method.permutation_count) > 0 and bool(method.adjust_p_values)
            else None
        ),
    )

    return KinaseActivityResult.from_trusted_owned(
        weighted_activity=scored.activity_scores,
        p_value_matrix=scored.p_value_matrix,
        q_value_matrix=scored.q_value_matrix,
        thresholded_substrate_counts=thresholded_substrate_counts,
        activity_substrate_counts=scored.substrate_count_table,
        substrate_count_matrix=scored.substrate_count_table,
        target_counts=target_counts,
        target_table=target_table,
        statistics_table=scored.statistics_table,
        method_summary=summary,
        method_diagnostics=diagnostics,
        policy_provenance=(policy,),
        activity_method=SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_METHOD,
        input_semantics=prepared.activity_input.semantics,
        profile_metadata=prepared.activity_input.profile_metadata,
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


def _require_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{field_name} must be an int")
    return int(value)


def _require_optional_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{field_name} must be an int or None")
    return int(value)


def _require_bool(value: object, *, field_name: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{field_name} must be a bool")
    return value


def _ordered_unique_strings(values: pd.Series) -> list[str]:
    return list(dict.fromkeys(str(value) for value in values.tolist()))


def _build_membership_lookup(
    *,
    kinases: list[str],
    membership: pd.DataFrame,
) -> dict[str, set[str]]:
    result: dict[str, set[str]] = {kinase: set() for kinase in kinases}
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
    site_labels: _ObjectArray,
    values: _FloatArray,
    finite_positions: _IntArray,
    ranking_direction: str,
) -> _RankedTieBlocks:
    if finite_positions.size == 0:
        return _RankedTieBlocks(
            site_labels=np.asarray([], dtype=object),
            block_starts=np.asarray([], dtype=np.int64),
            block_sizes=np.asarray([], dtype=np.int64),
        )
    finite_values: _FloatArray = values[finite_positions]
    # Sorting establishes only the order of distinct value blocks. The order of
    # rows inside an equal-value block is intentionally ignored by block-count
    # scoring.
    if ranking_direction == SSGSEA_RANKING_DIRECTION_ASCENDING:
        order: _IntArray = np.argsort(finite_values, kind="mergesort").astype(
            np.int64,
            copy=False,
        )
    else:
        order = np.argsort(-finite_values, kind="mergesort").astype(
            np.int64,
            copy=False,
        )
    ranked_positions: _IntArray = finite_positions[order]
    ranked_values: _FloatArray = values[ranked_positions]
    block_starts: _IntArray = np.concatenate(
        (
            np.asarray([0], dtype=np.int64),
            np.flatnonzero(ranked_values[1:] != ranked_values[:-1]).astype(np.int64)
            + 1,
        )
    )
    block_ends: _IntArray = np.concatenate(
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


def _score_from_hit_mask(hit_mask: _BoolArray) -> float:
    n_background = int(hit_mask.size)
    n_substrates = int(hit_mask.sum())
    n_misses = n_background - n_substrates
    if n_background == 0 or n_substrates == 0 or n_misses == 0:
        return np.nan
    hit_increment = 1.0 / float(n_substrates)
    miss_increment = 1.0 / float(n_misses)
    steps: _FloatArray = np.where(
        hit_mask,
        hit_increment,
        -miss_increment,
    ).astype(float)
    running: _FloatArray = np.cumsum(steps)
    return float(running.sum() / float(n_background))


def _score_from_ranked_hit_mask(  # pyright: ignore[reportUnusedFunction] - imported by activity regression tests as a private numerical seam
    *,
    hit_mask: _BoolArray,
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
    block_hit_counts: _IntArray,
    block_sizes: _IntArray,
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
    for block_hits, block_size in zip(
        _int_tuple(block_hit_counts),
        _int_tuple(block_sizes),
        strict=True,
    ):
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
    hit_mask: _BoolArray,
) -> _IntArray:
    if ranked_blocks.n_blocks == 0:
        return np.asarray([], dtype=np.int64)
    hit_values: _IntArray = hit_mask.astype(np.int64, copy=False)
    counts: _IntArray = np.add.reduceat(
        hit_values,
        ranked_blocks.block_starts,
    ).astype(np.int64, copy=False)
    return counts


def _null_score_cache_key(
    *,
    ranked_blocks: _RankedTieBlocks,
    n_substrates: int,
) -> _SsgseaNullScoreCacheKey:
    return _SsgseaNullScoreCacheKey(
        n_background=int(ranked_blocks.n_background),
        n_substrates=int(n_substrates),
        tie_block_sizes=(
            _int_tuple(ranked_blocks.block_sizes) if ranked_blocks.has_ties else ()
        ),
    )


def _build_kinase_tie_diagnostics(
    *,
    ranked_blocks: _RankedTieBlocks,
    hit_mask: _BoolArray,
) -> dict[str, int]:
    if not ranked_blocks.has_ties:
        return {
            "substrate_only_tie_blocks": 0,
            "non_substrate_only_tie_blocks": 0,
            "mixed_substrate_tie_blocks": 0,
        }
    block_hit_counts: _IntArray = _block_hit_counts(
        ranked_blocks=ranked_blocks,
        hit_mask=hit_mask,
    )
    tied: _BoolArray = ranked_blocks.block_sizes > 1
    tied_block_sizes: _IntArray = ranked_blocks.block_sizes[tied]
    tied_hit_counts: _IntArray = block_hit_counts[tied]
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
        hit_positions: _IntArray = np.asarray(
            rng.choice(
                int(n_background),
                size=int(n_substrates),
                replace=False,
            ),
            dtype=np.int64,
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
    """Derive the deterministic child RNG seed for one profile/kinase stream.

    ``condition_name`` is retained only as a private compatibility alias for
    ``profile_id``. The active scientific identity is profile-based.
    """
    seed_material = _ssgsea_permutation_seed_material(
        random_seed=random_seed,
        profile_id=profile_id,
        condition_name=condition_name,
        kinase_name=kinase_name,
    )
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


def _ssgsea_permutation_seed_material(
    *,
    random_seed: int,
    profile_id: str | None = None,
    condition_name: str | None = None,
    kinase_name: str,
) -> dict[str, object]:
    """Return the v1-encoded seed material used for ssGSEA permutations."""
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
    seed_material: dict[str, object] = {
        "kinase": str(kinase_name),
        "method_id": SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_METHOD.activity_method_id,
        "method_version": SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_POLICY_VERSION,
        "profile_id": resolved_profile_id,
        "random_seed": int(random_seed),
        # Keep both the historical JSON key and value for byte-for-byte v1 stream
        # compatibility. Public provenance reports the current profile-based
        # policy, not this compatibility salt.
        "seed_policy": _SSGSEA_PERMUTATION_RNG_SEED_V1_COMPATIBILITY_TOKEN,
        "seed_policy_version": SSGSEA_PERMUTATION_RNG_SEED_POLICY_VERSION,
        "stream": _SSGSEA_PERMUTATION_STREAM_NAME,
    }
    return seed_material


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
