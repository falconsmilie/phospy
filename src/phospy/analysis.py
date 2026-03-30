from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .activities import (
    build_kinase_target_table,
    compute_ksea_scores,
    compute_weighted_kinase_activity,
    count_predicted_targets,
)
from .io import load_pred_mat
from .validation.compatibility import validate_pred_mat_overlap
from .validation.tables import PredMatSchema, SiteMatrixSchema


@dataclass(slots=True)
class KinaseActivityResult:
    weighted_activity: pd.DataFrame
    ksea_scores: pd.DataFrame
    ksea_counts: pd.Series
    target_counts: pd.Series
    target_table: pd.DataFrame


def load_validated_pred_mat(pred_mat_path: str | Path) -> pd.DataFrame:
    """Load and validate a kinase prediction matrix from disk."""

    pred_mat = load_pred_mat(pred_mat_path)
    return PredMatSchema.validate(pred_mat, context="pred_mat")


def analyze_kinase_activity(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    threshold: float = 0.6,
    min_substrates: int = 3,
    top_n_substrates: int = 20,
) -> KinaseActivityResult:
    """Compute downstream kinase summaries from a validated prediction matrix."""

    validated_pred_mat = PredMatSchema.validate(pred_mat, context="pred_mat")
    validated_matrix = SiteMatrixSchema.validate(
        phospho_matrix,
        context="phospho_matrix",
    )
    validate_pred_mat_overlap(validated_pred_mat, validated_matrix)

    weighted_activity = compute_weighted_kinase_activity(
        pred_mat=validated_pred_mat,
        phospho_matrix=validated_matrix,
        top_n_substrates=top_n_substrates,
        min_substrates=min_substrates,
    )
    ksea_scores, ksea_counts = compute_ksea_scores(
        pred_mat=validated_pred_mat,
        phospho_matrix=validated_matrix,
        threshold=threshold,
        min_substrates=min_substrates,
    )
    target_counts = count_predicted_targets(validated_pred_mat, threshold=threshold)
    target_table = build_kinase_target_table(validated_pred_mat, threshold=threshold)

    return KinaseActivityResult(
        weighted_activity=weighted_activity,
        ksea_scores=ksea_scores,
        ksea_counts=ksea_counts,
        target_counts=target_counts,
        target_table=target_table,
    )


def load_and_analyze_kinase_activity(
    pred_mat_path: str | Path,
    phospho_matrix: pd.DataFrame,
    threshold: float = 0.6,
    min_substrates: int = 3,
    top_n_substrates: int = 20,
) -> KinaseActivityResult:
    """Load a prediction matrix from disk and compute downstream kinase summaries."""

    pred_mat = load_validated_pred_mat(pred_mat_path)
    return analyze_kinase_activity(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=threshold,
        min_substrates=min_substrates,
        top_n_substrates=top_n_substrates,
    )


def write_kinase_activity_outputs(
    result: KinaseActivityResult,
    outdir: str | Path,
) -> None:
    from .writers import KinaseActivityWriter

    KinaseActivityWriter.write(result, outdir)
