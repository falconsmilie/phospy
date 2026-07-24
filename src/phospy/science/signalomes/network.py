"""Score-derived kinase association services for signalome network outputs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from phospy.errors.workflows import WorkflowStageError
from phospy.science.configs import (
    SIGNALOME_KINASE_NETWORK_POLICIES,
    SIGNALOME_KINASE_NETWORK_POLICY_ABSOLUTE_THRESHOLD,
    SIGNALOME_KINASE_NETWORK_POLICY_POSITIVE_ONLY,
    SIGNALOME_KINASE_NETWORK_POLICY_SIGNED,
    SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_DEFAULT,
    SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_FLOOR,
    SignalomeKinaseNetworkPolicy,
)
from phospy.science.signalomes.constants import (
    CORRELATION_COLUMN,
    CORRELATION_REASON_COLUMN,
    CORRELATION_STATUS_COLUMN,
    DEGREE_COLUMN,
    KINASE_COLUMN,
    N_SUBSTRATES_COLUMN,
    SOURCE_KINASE_COLUMN,
    TARGET_KINASE_COLUMN,
    VALID_OBSERVATIONS_COLUMN,
)
from phospy.science.signalomes.models import (
    SIGNALOME_CORRELATION_STATUS_CONSTANT_PROFILE,
    SIGNALOME_CORRELATION_STATUS_FINITE,
    SIGNALOME_CORRELATION_STATUS_INSUFFICIENT_OBSERVATIONS,
    SIGNALOME_CORRELATION_STATUS_MISSING_VALUES,
    SIGNALOME_CORRELATION_STATUS_NON_FINITE_VALUES,
    SIGNALOME_CORRELATION_STATUS_UNDEFINED,
    SignalomeNetworkCorrelationDiagnostics,
)


def build_kinase_network(
    *,
    downstream_score_matrix: pd.DataFrame,
    kinase_order: Sequence[str],
    kinase_substrates: Mapping[str, Sequence[str]],
    threshold: float,
    network_policy: SignalomeKinaseNetworkPolicy = (
        SIGNALOME_KINASE_NETWORK_POLICY_SIGNED
    ),
    min_paired_observations: int = (
        SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_DEFAULT
    ),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build deterministic score-profile correlation edge and node tables."""

    edges, nodes, _, _ = build_kinase_network_with_diagnostics(
        downstream_score_matrix=downstream_score_matrix,
        kinase_order=kinase_order,
        kinase_substrates=kinase_substrates,
        threshold=threshold,
        network_policy=network_policy,
        min_paired_observations=min_paired_observations,
    )
    return edges, nodes


def build_kinase_network_with_diagnostics(
    *,
    downstream_score_matrix: pd.DataFrame,
    kinase_order: Sequence[str],
    kinase_substrates: Mapping[str, Sequence[str]],
    threshold: float,
    network_policy: SignalomeKinaseNetworkPolicy = (
        SIGNALOME_KINASE_NETWORK_POLICY_SIGNED
    ),
    min_paired_observations: int = (
        SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_DEFAULT
    ),
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    SignalomeNetworkCorrelationDiagnostics,
]:
    """Build exploratory network tables plus correlation traceability diagnostics."""

    _validate_min_paired_observations(min_paired_observations)
    kinase_index = pd.Index(
        [str(kinase) for kinase in kinase_order], name=KINASE_COLUMN
    )
    kinase_index = pd.Index(
        list(dict.fromkeys(kinase_index.tolist())), name=KINASE_COLUMN
    )
    if kinase_index.empty:
        raise WorkflowStageError("kinase network requires at least one kinase")
    available_kinases = set(downstream_score_matrix.columns.astype(str).tolist())
    missing_kinases: list[str] = [
        str(kinase) for kinase in kinase_index if str(kinase) not in available_kinases
    ]
    if missing_kinases:
        preview = ", ".join(missing_kinases[:3])
        suffix = "..." if len(missing_kinases) > 3 else ""
        raise WorkflowStageError(
            "downstream score matrix is missing kinases required for signalome network: "
            f"{preview}{suffix}"
        )

    aligned_scores = _precondition_network_scores(
        downstream_score_matrix=downstream_score_matrix,
        kinase_index=kinase_index,
    )
    candidate_correlations = _build_candidate_correlations(
        aligned_scores=aligned_scores,
        kinase_index=kinase_index,
        min_paired_observations=min_paired_observations,
    )
    finite_candidates = candidate_correlations.loc[
        candidate_correlations.loc[:, CORRELATION_STATUS_COLUMN].eq(
            SIGNALOME_CORRELATION_STATUS_FINITE
        )
        & candidate_correlations.loc[:, CORRELATION_COLUMN].notna(),
        [
            SOURCE_KINASE_COLUMN,
            TARGET_KINASE_COLUMN,
            CORRELATION_COLUMN,
            VALID_OBSERVATIONS_COLUMN,
        ],
    ].reset_index(drop=True)
    pair_correlations = finite_candidates.loc[:, CORRELATION_COLUMN].to_numpy(
        dtype=float,
        copy=False,
    )
    edge_correlations, edge_mask = _resolve_network_edges_by_policy(
        pair_correlations=pair_correlations,
        threshold=float(threshold),
        network_policy=network_policy,
    )
    selected = finite_candidates.loc[edge_mask, :].copy()
    selected.loc[:, CORRELATION_COLUMN] = edge_correlations[edge_mask]
    edges = (
        selected.loc[
            :,
            [
                SOURCE_KINASE_COLUMN,
                TARGET_KINASE_COLUMN,
                CORRELATION_COLUMN,
                VALID_OBSERVATIONS_COLUMN,
            ],
        ]
        if not selected.empty
        else _empty_edges_table()
    )
    edges = edges.astype(
        {
            SOURCE_KINASE_COLUMN: str,
            TARGET_KINASE_COLUMN: str,
            CORRELATION_COLUMN: float,
            VALID_OBSERVATIONS_COLUMN: "int64",
        }
    )
    edges = edges.sort_values(
        [SOURCE_KINASE_COLUMN, TARGET_KINASE_COLUMN],
        ascending=[True, True],
        kind="stable",
    ).reset_index(drop=True)

    degree_values = pd.Series(0, index=kinase_index.copy(), dtype="int64")
    if not edges.empty:
        edge_degrees = pd.concat(
            [
                edges.loc[:, SOURCE_KINASE_COLUMN],
                edges.loc[:, TARGET_KINASE_COLUMN],
            ],
            axis=0,
        ).value_counts()
        degree_values.loc[edge_degrees.index.astype(str)] = edge_degrees.to_numpy(
            dtype=np.int64,
            copy=False,
        )
    node_substrates = np.asarray(
        [
            len(tuple(kinase_substrates.get(str(kinase), ())))
            for kinase in kinase_index.to_numpy(dtype=object, copy=False)
        ],
        dtype=np.int64,
    )
    nodes = pd.DataFrame(
        {
            DEGREE_COLUMN: degree_values.to_numpy(dtype=np.int64, copy=False),
            N_SUBSTRATES_COLUMN: node_substrates,
        },
        index=kinase_index.copy(),
    )
    nodes.index.name = KINASE_COLUMN
    nodes = nodes.astype({DEGREE_COLUMN: "int64", N_SUBSTRATES_COLUMN: "int64"})
    diagnostics = _build_network_correlation_diagnostics(
        candidate_correlations=candidate_correlations,
        edges=edges,
    )
    return edges, nodes, candidate_correlations, diagnostics


def _build_candidate_correlations(
    *,
    aligned_scores: pd.DataFrame,
    kinase_index: pd.Index,
    min_paired_observations: int,
) -> pd.DataFrame:
    if len(kinase_index) < 2:
        return _empty_candidate_correlation_table()
    rows: list[dict[str, object]] = []
    score_values = {
        str(kinase): aligned_scores.loc[:, kinase].to_numpy(dtype=float, copy=False)
        for kinase in kinase_index
    }
    kinase_names = kinase_index.astype(str).tolist()
    for source_position, source_kinase in enumerate(kinase_names[:-1]):
        source_values = score_values[source_kinase]
        for target_kinase in kinase_names[source_position + 1 :]:
            correlation, status, valid_observations, reason = (
                _classify_pair_correlation(
                    source_values=source_values,
                    target_values=score_values[target_kinase],
                    min_paired_observations=min_paired_observations,
                )
            )
            rows.append(
                {
                    SOURCE_KINASE_COLUMN: source_kinase,
                    TARGET_KINASE_COLUMN: target_kinase,
                    CORRELATION_COLUMN: correlation,
                    CORRELATION_STATUS_COLUMN: status,
                    VALID_OBSERVATIONS_COLUMN: valid_observations,
                    CORRELATION_REASON_COLUMN: reason,
                }
            )
    candidates = pd.DataFrame.from_records(rows)
    if candidates.empty:
        return _empty_candidate_correlation_table()
    candidates = candidates.astype(
        {
            SOURCE_KINASE_COLUMN: str,
            TARGET_KINASE_COLUMN: str,
            CORRELATION_COLUMN: float,
            CORRELATION_STATUS_COLUMN: str,
            VALID_OBSERVATIONS_COLUMN: "int64",
        }
    )
    # Normalise missing reasons to NaN so CSV/parquet round-trips preserve equality.
    reasons = candidates.loc[:, CORRELATION_REASON_COLUMN]
    candidates.loc[:, CORRELATION_REASON_COLUMN] = reasons.where(
        ~reasons.isna(),
        np.nan,
    )
    return candidates.loc[
        :,
        [
            SOURCE_KINASE_COLUMN,
            TARGET_KINASE_COLUMN,
            CORRELATION_COLUMN,
            CORRELATION_STATUS_COLUMN,
            VALID_OBSERVATIONS_COLUMN,
            CORRELATION_REASON_COLUMN,
        ],
    ]


def _classify_pair_correlation(
    *,
    source_values: np.ndarray,
    target_values: np.ndarray,
    min_paired_observations: int = (
        SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_DEFAULT
    ),
) -> tuple[float, str, int, str | None]:
    source_array = np.asarray(source_values, dtype=float)
    target_array = np.asarray(target_values, dtype=float)
    finite_pair_mask = np.isfinite(source_array) & np.isfinite(target_array)
    valid_observations = int(finite_pair_mask.sum())
    non_finite_values_present = (
        ((~np.isfinite(source_array)) & (~np.isnan(source_array)))
        | ((~np.isfinite(target_array)) & (~np.isnan(target_array)))
    ).any()
    missing_values_present = (
        np.isnan(source_array).any() or np.isnan(target_array).any()
    )
    if non_finite_values_present:
        return (
            float("nan"),
            SIGNALOME_CORRELATION_STATUS_NON_FINITE_VALUES,
            valid_observations,
            "input contains non-finite values",
        )
    if valid_observations < int(min_paired_observations):
        if missing_values_present:
            return (
                float("nan"),
                SIGNALOME_CORRELATION_STATUS_MISSING_VALUES,
                valid_observations,
                "missing values reduced paired observations below minimum",
            )
        return (
            float("nan"),
            SIGNALOME_CORRELATION_STATUS_INSUFFICIENT_OBSERVATIONS,
            valid_observations,
            "fewer than the minimum paired finite observations",
        )
    valid_source = source_array[finite_pair_mask]
    valid_target = target_array[finite_pair_mask]
    source_variance = float(np.var(valid_source, ddof=0))
    target_variance = float(np.var(valid_target, ddof=0))
    if source_variance == 0.0 or target_variance == 0.0:
        return (
            float("nan"),
            SIGNALOME_CORRELATION_STATUS_CONSTANT_PROFILE,
            valid_observations,
            "zero variance among paired finite observations",
        )
    correlation = float(np.corrcoef(valid_source, valid_target)[0, 1])
    if np.isfinite(correlation):
        return (
            float(np.clip(correlation, -1.0, 1.0)),
            SIGNALOME_CORRELATION_STATUS_FINITE,
            valid_observations,
            None,
        )
    if missing_values_present:
        return (
            float("nan"),
            SIGNALOME_CORRELATION_STATUS_MISSING_VALUES,
            valid_observations,
            "missing values prevented stable correlation estimation",
        )
    return (
        float("nan"),
        SIGNALOME_CORRELATION_STATUS_UNDEFINED,
        valid_observations,
        "correlation calculation returned non-finite output",
    )


def _validate_min_paired_observations(min_paired_observations: int) -> None:
    try:
        resolved = int(min_paired_observations)
    except (TypeError, ValueError) as exc:
        raise WorkflowStageError(
            "signalome network min_paired_observations must be an integer"
        ) from exc
    if resolved < SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_FLOOR:
        raise WorkflowStageError(
            "signalome network min_paired_observations must be at least "
            f"{SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_FLOOR}; "
            f"got {resolved}. Legacy threshold 2 cannot be used for new "
            "signalome network execution; set "
            "config.output.network_min_paired_finite_observations to at least 3."
        )


def _build_network_correlation_diagnostics(
    *,
    candidate_correlations: pd.DataFrame,
    edges: pd.DataFrame,
) -> SignalomeNetworkCorrelationDiagnostics:
    if candidate_correlations.empty:
        return SignalomeNetworkCorrelationDiagnostics(
            total_candidate_correlations=0,
            finite_correlations=0,
            undefined_correlations=0,
            constant_profile_correlations=0,
            insufficient_observation_correlations=0,
            missing_value_correlations=0,
            non_finite_value_correlations=0,
            edges_created=int(edges.shape[0]),
            edges_skipped_non_finite_correlation=0,
            edges_skipped_below_threshold=0,
            edges_skipped_insufficient_paired_observations=0,
            edges_skipped_constant_profile=0,
            edges_skipped_missing_score=0,
            edges_skipped_non_finite_score=0,
            edges_skipped_undefined_correlation=0,
        )
    status_counts = (
        candidate_correlations.loc[:, CORRELATION_STATUS_COLUMN]
        .astype(str)
        .value_counts()
        .to_dict()
    )
    total_candidates = int(candidate_correlations.shape[0])
    finite_correlations = int(status_counts.get(SIGNALOME_CORRELATION_STATUS_FINITE, 0))
    undefined_correlations = int(total_candidates - finite_correlations)
    edges_created = int(edges.shape[0])
    constant_profile_correlations = int(
        status_counts.get(SIGNALOME_CORRELATION_STATUS_CONSTANT_PROFILE, 0)
    )
    insufficient_observation_correlations = int(
        status_counts.get(SIGNALOME_CORRELATION_STATUS_INSUFFICIENT_OBSERVATIONS, 0)
    )
    missing_value_correlations = int(
        status_counts.get(SIGNALOME_CORRELATION_STATUS_MISSING_VALUES, 0)
    )
    non_finite_value_correlations = int(
        status_counts.get(SIGNALOME_CORRELATION_STATUS_NON_FINITE_VALUES, 0)
    )
    undefined_status_correlations = int(
        status_counts.get(SIGNALOME_CORRELATION_STATUS_UNDEFINED, 0)
    )
    return SignalomeNetworkCorrelationDiagnostics(
        total_candidate_correlations=total_candidates,
        finite_correlations=finite_correlations,
        undefined_correlations=undefined_correlations,
        constant_profile_correlations=constant_profile_correlations,
        insufficient_observation_correlations=insufficient_observation_correlations,
        missing_value_correlations=missing_value_correlations,
        non_finite_value_correlations=non_finite_value_correlations,
        edges_created=edges_created,
        edges_skipped_non_finite_correlation=undefined_correlations,
        edges_skipped_below_threshold=max(
            0,
            finite_correlations - edges_created,
        ),
        edges_skipped_insufficient_paired_observations=(
            insufficient_observation_correlations
        ),
        edges_skipped_constant_profile=constant_profile_correlations,
        edges_skipped_missing_score=missing_value_correlations,
        edges_skipped_non_finite_score=non_finite_value_correlations,
        edges_skipped_undefined_correlation=undefined_status_correlations,
    )


def _resolve_network_edges_by_policy(
    *,
    pair_correlations: np.ndarray,
    threshold: float,
    network_policy: SignalomeKinaseNetworkPolicy,
) -> tuple[np.ndarray, np.ndarray]:
    if network_policy == SIGNALOME_KINASE_NETWORK_POLICY_POSITIVE_ONLY:
        edge_values = pair_correlations
        edge_mask = pair_correlations >= float(threshold)
        return edge_values, edge_mask
    if network_policy == SIGNALOME_KINASE_NETWORK_POLICY_ABSOLUTE_THRESHOLD:
        edge_values = np.abs(pair_correlations)
        edge_mask = edge_values >= float(threshold)
        return edge_values, edge_mask
    if network_policy == SIGNALOME_KINASE_NETWORK_POLICY_SIGNED:
        edge_values = pair_correlations
        edge_mask = np.abs(pair_correlations) >= float(threshold)
        return edge_values, edge_mask
    allowed = ", ".join(sorted(SIGNALOME_KINASE_NETWORK_POLICIES))
    raise WorkflowStageError(
        f"unsupported network_policy '{network_policy}'; expected one of: {allowed}"
    )


def _precondition_network_scores(
    *,
    downstream_score_matrix: pd.DataFrame,
    kinase_index: pd.Index,
) -> pd.DataFrame:
    aligned_scores = downstream_score_matrix.loc[:, kinase_index].astype(float)
    supported_row_mask = (
        aligned_scores.notna().any(axis=1).to_numpy(dtype=bool, copy=False)
    )
    if supported_row_mask.all():
        return aligned_scores
    return aligned_scores.iloc[supported_row_mask, :]


def _empty_edges_table() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            SOURCE_KINASE_COLUMN,
            TARGET_KINASE_COLUMN,
            CORRELATION_COLUMN,
            VALID_OBSERVATIONS_COLUMN,
        ]
    ).astype(
        {
            SOURCE_KINASE_COLUMN: str,
            TARGET_KINASE_COLUMN: str,
            CORRELATION_COLUMN: float,
            VALID_OBSERVATIONS_COLUMN: "int64",
        }
    )


def _empty_candidate_correlation_table() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            SOURCE_KINASE_COLUMN,
            TARGET_KINASE_COLUMN,
            CORRELATION_COLUMN,
            CORRELATION_STATUS_COLUMN,
            VALID_OBSERVATIONS_COLUMN,
            CORRELATION_REASON_COLUMN,
        ]
    ).astype(
        {
            SOURCE_KINASE_COLUMN: str,
            TARGET_KINASE_COLUMN: str,
            CORRELATION_COLUMN: float,
            CORRELATION_STATUS_COLUMN: str,
            VALID_OBSERVATIONS_COLUMN: "int64",
        }
    )


__all__ = ["build_kinase_network", "build_kinase_network_with_diagnostics"]
