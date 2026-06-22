"""ssGSEA-style substrate-set enrichment activity-like score method."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.activities.models import (
    SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_METHOD,
    ActivityMethodSummary,
    KinaseActivityResult,
    SsgseaSubstrateEnrichmentActivityDiagnostics,
)
from phospy.science.activities.scientific_policies import (
    build_ssgsea_substrate_enrichment_activity_policy,
)
from phospy.science.activities.statistics import benjamini_hochberg_q_values
from phospy.tables.activity import ActivityMatrix
from phospy.validation.common.dataframes import (
    require_canonical_string_column,
    require_columns,
    require_dataframe,
    require_non_empty_string_column,
    require_unique_columns,
    require_unique_index,
    require_unique_row_pairs,
)

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

_KINASE_COLUMN = "kinase"
_SUBSTRATE_COLUMN = "substrate_site"


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
        effect_matrix: pd.DataFrame,
        kinase_substrate_membership: pd.DataFrame,
    ) -> KinaseActivityResult:
        effects = ActivityMatrix(
            frame=effect_matrix,
            field_name="activity_inputs.effect_matrix",
        ).frame
        membership = _validate_membership_table(kinase_substrate_membership)
        site_labels = np.asarray(effects.index.astype(str).tolist(), dtype=object)
        site_universe = set(str(value) for value in site_labels.tolist())
        aligned_membership = membership.loc[
            membership.loc[:, _SUBSTRATE_COLUMN].astype(str).isin(site_universe),
            :,
        ].copy(deep=True)
        kinases = _ordered_unique_strings(membership.loc[:, _KINASE_COLUMN])
        kinase_index = pd.Index(kinases, name="kinase")
        condition_index = pd.Index(
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
            columns=condition_index,
            dtype=float,
        )
        substrate_count_table = pd.DataFrame(
            0,
            index=kinase_index,
            columns=condition_index,
            dtype=int,
        )
        p_value_matrix = (
            pd.DataFrame(
                np.nan,
                index=kinase_index,
                columns=condition_index,
                dtype=float,
            )
            if int(self.permutation_count) > 0
            else None
        )
        q_value_matrix = (
            pd.DataFrame(
                np.nan,
                index=kinase_index,
                columns=condition_index,
                dtype=float,
            )
            if int(self.permutation_count) > 0 and bool(self.adjust_p_values)
            else None
        )

        rng: np.random.Generator | None = None
        if int(self.permutation_count) > 0:
            random_seed = self.random_seed
            if random_seed is None:
                raise ValueError(
                    "ssgsea random_seed must be set when permutation_count is positive"
                )
            rng = np.random.default_rng(int(random_seed))
        effect_values = effects.to_numpy(dtype=float, copy=False)
        rows: list[dict[str, object]] = []
        counts = {
            SSGSEA_STATUS_COMPUTED: 0,
            SSGSEA_STATUS_INSUFFICIENT_SUBSTRATES: 0,
            SSGSEA_STATUS_NO_FINITE_BACKGROUND_VALUES: 0,
            SSGSEA_STATUS_INSUFFICIENT_BACKGROUND_SITES: 0,
            SSGSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES: 0,
        }

        for condition_position, condition_name in enumerate(condition_index):
            condition_values = effect_values[:, condition_position]
            finite_positions = np.flatnonzero(np.isfinite(condition_values))
            ranked_site_labels = _rank_sites(
                site_labels=site_labels,
                values=condition_values,
                finite_positions=finite_positions,
                ranking_direction=str(self.ranking_direction),
            )
            n_background = int(ranked_site_labels.size)

            for kinase_position, kinase_name in enumerate(kinase_index):
                substrate_sites = membership_by_kinase[str(kinase_name)]
                hit_mask = np.fromiter(
                    (str(site_id) in substrate_sites for site_id in ranked_site_labels),
                    dtype=bool,
                    count=n_background,
                )
                n_substrates = int(hit_mask.sum())
                substrate_count_table.iat[kinase_position, condition_position] = (
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
                    enrichment_score = _score_from_hit_mask(hit_mask)
                    if not np.isfinite(enrichment_score):
                        status = SSGSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES
                        reason = "enrichment score is not finite"
                    else:
                        activity_scores.iat[
                            kinase_position,
                            condition_position,
                        ] = float(enrichment_score)
                        if rng is not None and p_value_matrix is not None:
                            p_value = _permutation_p_value(
                                observed_score=float(enrichment_score),
                                n_background=n_background,
                                n_substrates=n_substrates,
                                permutation_count=int(self.permutation_count),
                                rng=rng,
                            )
                            p_value_matrix.iat[
                                kinase_position,
                                condition_position,
                            ] = float(p_value)

                counts[status] += 1
                rows.append(
                    {
                        "kinase": str(kinase_name),
                        "condition": str(condition_name),
                        "z_score": np.nan,
                        "enrichment_score": (
                            float(enrichment_score)
                            if np.isfinite(enrichment_score)
                            else np.nan
                        ),
                        "p_value": float(p_value) if np.isfinite(p_value) else np.nan,
                        "q_value": np.nan,
                        "n_substrates": int(n_substrates),
                        "n_background_sites": int(n_background),
                        "evidence_threshold": np.nan,
                        "evidence_threshold_operator": "",
                        "evidence_threshold_description": (
                            "not thresholded; explicit kinase-substrate membership"
                        ),
                        "min_substrates": int(self.min_substrates),
                        "ranking_direction": str(self.ranking_direction),
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
                "condition",
                "z_score",
                "enrichment_score",
                "p_value",
                "q_value",
                "n_substrates",
                "n_background_sites",
                "evidence_threshold",
                "evidence_threshold_operator",
                "evidence_threshold_description",
                "min_substrates",
                "ranking_direction",
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
                condition_index=condition_index,
            )

        target_counts = _build_target_counts(
            aligned_membership=aligned_membership,
            kinase_index=kinase_index,
        )
        thresholded_substrate_counts = target_counts.rename("n_substrates")
        target_table = _build_target_table(aligned_membership=aligned_membership)
        summary = ActivityMethodSummary(
            kinases_evaluated=int(len(kinase_index)),
            kinase_condition_pairs_evaluated=int(
                len(kinase_index) * len(condition_index)
            ),
            kinase_condition_pairs_computed=counts[SSGSEA_STATUS_COMPUTED],
            kinase_condition_pairs_insufficient_substrates=counts[
                SSGSEA_STATUS_INSUFFICIENT_SUBSTRATES
            ],
            kinase_condition_pairs_invalid_background_variance=(
                counts[SSGSEA_STATUS_NO_FINITE_BACKGROUND_VALUES]
                + counts[SSGSEA_STATUS_INSUFFICIENT_BACKGROUND_SITES]
            ),
            kinase_condition_pairs_no_finite_background_values=counts[
                SSGSEA_STATUS_NO_FINITE_BACKGROUND_VALUES
            ],
            kinase_condition_pairs_no_finite_substrate_values=counts[
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


def _rank_sites(
    *,
    site_labels: np.ndarray,
    values: np.ndarray,
    finite_positions: np.ndarray,
    ranking_direction: str,
) -> np.ndarray:
    if finite_positions.size == 0:
        return np.asarray([], dtype=object)
    finite_values = values[finite_positions]
    if ranking_direction == SSGSEA_RANKING_DIRECTION_ASCENDING:
        order = np.argsort(finite_values, kind="mergesort")
    else:
        order = np.argsort(-finite_values, kind="mergesort")
    return site_labels[finite_positions[order]]


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


def _permutation_p_value(
    *,
    observed_score: float,
    n_background: int,
    n_substrates: int,
    permutation_count: int,
    rng: np.random.Generator,
) -> float:
    extreme_count = 0
    observed_abs = abs(float(observed_score))
    for _ in range(int(permutation_count)):
        hit_positions = rng.choice(
            int(n_background),
            size=int(n_substrates),
            replace=False,
        )
        hit_mask = np.zeros(int(n_background), dtype=bool)
        hit_mask[hit_positions] = True
        score = _score_from_hit_mask(hit_mask)
        if abs(float(score)) >= observed_abs:
            extreme_count += 1
    return float((extreme_count + 1) / (int(permutation_count) + 1))


def _apply_q_values(
    *,
    statistics_table: pd.DataFrame,
    q_value_matrix: pd.DataFrame,
    condition_index: pd.Index,
) -> None:
    for condition_name in condition_index:
        condition_mask = statistics_table.loc[:, "condition"].astype(str) == str(
            condition_name
        )
        computed_mask = (
            statistics_table.loc[:, "computability_status"] == SSGSEA_STATUS_COMPUTED
        )
        selected = condition_mask & computed_mask
        if not bool(selected.any()):
            continue
        condition_p_values = statistics_table.loc[selected, "p_value"].astype(float)
        q_values = benjamini_hochberg_q_values(condition_p_values)
        statistics_table.loc[selected, "q_value"] = q_values.to_numpy(
            dtype=float,
            copy=False,
        )
        condition_rows = statistics_table.loc[selected, "kinase"].astype(str)
        for kinase_name, q_value in zip(
            condition_rows.tolist(),
            q_values.to_numpy(dtype=float, copy=False).tolist(),
            strict=True,
        ):
            q_value_matrix.at[str(kinase_name), str(condition_name)] = float(q_value)


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
    "SSGSEA_STATUS_COMPUTED",
    "SSGSEA_STATUS_INSUFFICIENT_BACKGROUND_SITES",
    "SSGSEA_STATUS_INSUFFICIENT_SUBSTRATES",
    "SSGSEA_STATUS_NO_FINITE_BACKGROUND_VALUES",
    "SSGSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES",
    "SsgseaSubstrateEnrichmentActivityMethod",
]
