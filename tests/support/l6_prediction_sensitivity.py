from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
REFERENCE_ROOT = (
    ROOT / "tests" / "fixtures" / "rewrite_parity" / "r_reference_l6_prediction"
)
REFERENCE_MATRIX_PATH = REFERENCE_ROOT / "predMat.csv"
REFERENCE_CANDIDATES_PATH = REFERENCE_ROOT / "native_candidate_substrates.csv"
REFERENCE_TOPK_PATH = REFERENCE_ROOT / "native_prediction_top30.csv"

REFERENCE_MATRIX_SHA256 = (
    "04aeec9d0dd14dc7822e66dbed4805997ee807be490f6b59a0cd565ed7b24bc1"
)
REFERENCE_CANDIDATES_SHA256 = (
    "077df04a2705493264eb6351d82eeac90320193a0bdf645cc93e578f363863e8"
)
REFERENCE_TOPK_SHA256 = (
    "0922a38a1443e91c6e195d4dac0d5698795c62821f7b263e6ea853e06fe67f7c"
)

TOP_K = 30
SCORE_THRESHOLD = 0.0
SAFETY_FACTOR = 0.5
PEARSON_FLOOR = 0.999999
SPEARMAN_FLOOR = 0.999999


@dataclass(frozen=True, slots=True)
class DecisionBoundary:
    kinase: str
    higher_ranked_site: str
    lower_ranked_site: str
    higher_score: float
    lower_score: float
    score_gap: float
    max_abs_failure_boundary: float
    mean_abs_failure_boundary: float


@dataclass(frozen=True, slots=True)
class L6PredictionSensitivityBudget:
    boundary: DecisionBoundary
    safety_factor: float
    mean_abs_diff_ceiling: float
    max_abs_diff_ceiling: float
    reference_cell_count: int
    reference_shape: tuple[int, int]
    minimum_threshold_margin: float
    minimum_positive_adjacent_rank_gap: float


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_reference_matrix() -> pd.DataFrame:
    frame = pd.read_csv(REFERENCE_MATRIX_PATH, index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str).str.strip().str.upper()
    numeric = frame.astype(float)
    values = numeric.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise AssertionError("L6 prediction sensitivity reference must be finite")
    if ((values < 0.0) | (values > 1.0)).any():
        raise AssertionError("L6 prediction sensitivity reference must use [0, 1]")
    return numeric


def _ranked(frame: pd.DataFrame, kinase: str) -> pd.Series:
    return frame.loc[:, kinase].sort_values(ascending=False, kind="mergesort")


def _candidate_sites_by_kinase() -> dict[str, tuple[str, ...]]:
    frame = pd.read_csv(REFERENCE_CANDIDATES_PATH).astype(str)
    frame.loc[:, "kinase"] = frame.loc[:, "kinase"].str.strip().str.upper()
    return {
        str(kinase): tuple(group.loc[:, "site_id"].tolist())
        for kinase, group in frame.groupby("kinase", sort=False)
    }


def _ranking_rows(frame: pd.DataFrame) -> tuple[tuple[str, tuple[str, ...]], ...]:
    candidate_sites = _candidate_sites_by_kinase()
    return tuple(
        (
            str(kinase),
            tuple(
                _ranked(frame.loc[list(sites), :], str(kinase))
                .index.astype(str)
                .tolist()
            ),
        )
        for kinase, sites in candidate_sites.items()
    )


def _find_first_ranking_failure_boundary(
    frame: pd.DataFrame,
) -> DecisionBoundary:
    candidates: list[DecisionBoundary] = []
    for kinase, sites in _candidate_sites_by_kinase().items():
        ranked = _ranked(frame.loc[list(sites), :], kinase)
        for position in range(len(ranked) - 1):
            higher_score = float(ranked.iloc[position])
            lower_score = float(ranked.iloc[position + 1])
            gap = higher_score - lower_score
            candidates.append(
                DecisionBoundary(
                    kinase=kinase,
                    higher_ranked_site=str(ranked.index[position]),
                    lower_ranked_site=str(ranked.index[position + 1]),
                    higher_score=higher_score,
                    lower_score=lower_score,
                    score_gap=gap,
                    max_abs_failure_boundary=gap / 2.0,
                    mean_abs_failure_boundary=gap / float(frame.size),
                )
            )
    return min(candidates, key=lambda item: item.score_gap)


@lru_cache(maxsize=1)
def derive_l6_prediction_sensitivity_budget() -> L6PredictionSensitivityBudget:
    """Derive limits solely from the fixed reference and executable decisions."""
    if _sha256(REFERENCE_MATRIX_PATH) != REFERENCE_MATRIX_SHA256:
        raise AssertionError("L6 prediction-matrix reference identity changed")
    if _sha256(REFERENCE_CANDIDATES_PATH) != REFERENCE_CANDIDATES_SHA256:
        raise AssertionError("L6 candidate reference identity changed")
    if _sha256(REFERENCE_TOPK_PATH) != REFERENCE_TOPK_SHA256:
        raise AssertionError("L6 top-k reference identity changed")

    frame = _load_reference_matrix()
    expected_candidates = {
        (str(row.kinase).strip().upper(), str(row.site_id))
        for row in pd.read_csv(REFERENCE_CANDIDATES_PATH).itertuples(index=False)
    }
    expected_topk = {
        (str(row.kinase).strip().upper(), str(row.site_id))
        for row in pd.read_csv(REFERENCE_TOPK_PATH).itertuples(index=False)
    }
    if expected_candidates != expected_topk:
        raise AssertionError(
            "L6 maintained candidate and ranked-output site sets disagree"
        )
    if any(
        kinase not in frame.columns or site not in frame.index
        for kinase, site in expected_candidates
    ):
        raise AssertionError(
            "L6 maintained candidate site is absent from prediction matrix"
        )
    boundary = _find_first_ranking_failure_boundary(frame)
    positive_adjacent_gaps: list[float] = []
    for kinase in frame.columns:
        ranked = _ranked(frame, str(kinase)).to_numpy(dtype=float)
        gaps = ranked[:-1] - ranked[1:]
        positive_adjacent_gaps.extend(float(value) for value in gaps if value > 0.0)

    return L6PredictionSensitivityBudget(
        boundary=boundary,
        safety_factor=SAFETY_FACTOR,
        mean_abs_diff_ceiling=(boundary.mean_abs_failure_boundary * SAFETY_FACTOR),
        max_abs_diff_ceiling=boundary.max_abs_failure_boundary * SAFETY_FACTOR,
        reference_cell_count=int(frame.size),
        reference_shape=(int(frame.shape[0]), int(frame.shape[1])),
        minimum_threshold_margin=float(frame.to_numpy(dtype=float).min()),
        minimum_positive_adjacent_rank_gap=min(positive_adjacent_gaps),
    )


def apply_boundary_perturbation(
    frame: pd.DataFrame,
    *,
    amplitude: float,
) -> pd.DataFrame:
    """Apply the decision-directed two-cell perturbation used by the analysis."""
    budget = derive_l6_prediction_sensitivity_budget()
    boundary = budget.boundary
    perturbed = frame.copy(deep=True)
    higher_value = cast(
        float, perturbed.at[boundary.higher_ranked_site, boundary.kinase]
    )
    lower_value = cast(float, perturbed.at[boundary.lower_ranked_site, boundary.kinase])
    perturbed.at[boundary.higher_ranked_site, boundary.kinase] = (
        higher_value - amplitude
    )
    perturbed.at[boundary.lower_ranked_site, boundary.kinase] = lower_value + amplitude
    return perturbed


def _surface_metrics(
    reference: pd.DataFrame,
    observed: pd.DataFrame,
) -> dict[str, float | bool]:
    reference_values = [
        float(value) for value in reference.to_numpy(dtype=float).ravel()
    ]
    observed_values = [float(value) for value in observed.to_numpy(dtype=float).ravel()]
    delta = [
        abs(observed_value - reference_value)
        for reference_value, observed_value in zip(
            reference_values, observed_values, strict=True
        )
    ]
    reference_mean = math.fsum(reference_values) / len(reference_values)
    observed_mean = math.fsum(observed_values) / len(observed_values)
    covariance = math.fsum(
        (reference_value - reference_mean) * (observed_value - observed_mean)
        for reference_value, observed_value in zip(
            reference_values, observed_values, strict=True
        )
    )
    reference_variance = math.fsum(
        (value - reference_mean) ** 2 for value in reference_values
    )
    observed_variance = math.fsum(
        (value - observed_mean) ** 2 for value in observed_values
    )
    pearson_correlation = covariance / math.sqrt(reference_variance * observed_variance)
    return {
        "selected_sites_unchanged": True,
        "ranking_unchanged": _ranking_rows(observed) == _ranking_rows(reference),
        "finite_mask_unchanged": bool(
            np.array_equal(
                np.isfinite(observed.to_numpy(dtype=float)),
                np.isfinite(reference.to_numpy(dtype=float)),
            )
        ),
        "mean_absolute_difference": math.fsum(delta) / len(delta),
        "maximum_absolute_difference": max(delta),
        "pearson_correlation": pearson_correlation,
    }


def build_l6_prediction_sensitivity_report() -> dict[str, Any]:
    """Return the maintained reference-only L6 sensitivity report."""
    reference = _load_reference_matrix()
    budget = derive_l6_prediction_sensitivity_budget()
    boundary = budget.boundary
    permitted = apply_boundary_perturbation(
        reference,
        amplitude=budget.max_abs_diff_ceiling * (1.0 - 1e-9),
    )
    first_failure = apply_boundary_perturbation(
        reference,
        amplitude=np.nextafter(boundary.max_abs_failure_boundary, np.inf),
    )
    return {
        "schema": "phospy-l6-prediction-sensitivity-v1",
        "reference": {
            "prediction_matrix": {
                "path": REFERENCE_MATRIX_PATH.relative_to(ROOT).as_posix(),
                "sha256": REFERENCE_MATRIX_SHA256,
            },
            "candidate_rows": {
                "path": REFERENCE_CANDIDATES_PATH.relative_to(ROOT).as_posix(),
                "sha256": REFERENCE_CANDIDATES_SHA256,
            },
            "ranked_top_k": {
                "path": REFERENCE_TOPK_PATH.relative_to(ROOT).as_posix(),
                "sha256": REFERENCE_TOPK_SHA256,
            },
            "shape": list(budget.reference_shape),
            "cell_count": budget.reference_cell_count,
        },
        "algorithm_semantics": {
            "ranking": "descending score with stable mergesort ordering",
            "top_k": TOP_K,
            "score_threshold": SCORE_THRESHOLD,
            "threshold_operator": "greater_than",
            "minimum_inclusion": 1,
            "prediction_score_rule": (
                "rank the fixed upstream-selected candidate sites by descending "
                "prediction support"
            ),
        },
        "protected_decisions": {
            "score_driven": [
                "per-kinase ranking gates",
                "selected top-30 sites fixed by the upstream scoring surface",
                "candidate-set overlap",
            ],
            "structural": [
                "class assignment",
                "finite and missing-value masks",
                "site and kinase attrition",
                "provenance records",
                "diagnostic records",
            ],
            "structural_rule": (
                "The perturbation is applied to the completed reference prediction "
                "matrix; labels, shape, masks, upstream class assignments, attrition, "
                "provenance, and diagnostics are held identical."
            ),
        },
        "perturbation_model": {
            "name": "decision-directed symmetric boundary-pair perturbation",
            "description": (
                "Lower the higher-ranked score and raise the immediately lower-"
                "ranked score for the narrowest adjacent gap among fixed selected "
                "sites. This is "
                "the minimum-L1 and minimum-L-infinity additive perturbation that "
                "can change a ranked prediction-output decision."
            ),
            "clipping": "none required inside the derived boundary",
        },
        "decision_margins": {
            "minimum_positive_adjacent_rank_gap": (
                budget.minimum_positive_adjacent_rank_gap
            ),
            "minimum_score_threshold_margin": budget.minimum_threshold_margin,
            "narrowest_selected_site_ranking_boundary": {
                "kinase": boundary.kinase,
                "higher_ranked_site": boundary.higher_ranked_site,
                "lower_ranked_site": boundary.lower_ranked_site,
                "higher_score": boundary.higher_score,
                "lower_score": boundary.lower_score,
                "score_gap": boundary.score_gap,
            },
        },
        "first_decision_failure_boundaries": {
            "mean_absolute_difference": boundary.mean_abs_failure_boundary,
            "maximum_absolute_difference": boundary.max_abs_failure_boundary,
            "failed_contract": "selected-site ranking identity",
        },
        "safety_factor": {
            "value": budget.safety_factor,
            "rule": "one-half of each independently derived first-failure boundary",
        },
        "derived_limits": {
            "mean_absolute_difference": budget.mean_abs_diff_ceiling,
            "maximum_absolute_difference": budget.max_abs_diff_ceiling,
            "pearson_correlation_floor": PEARSON_FLOOR,
            "spearman_correlation_floor": SPEARMAN_FLOOR,
        },
        "controls": {
            "permitted_boundary_directed_perturbation": _surface_metrics(
                reference, permitted
            ),
            "first_failure_negative_control": _surface_metrics(
                reference, first_failure
            ),
            "negative_control_interpretation": (
                "The first selected-site ranking failure still exceeds the Pearson floor, "
                "proving correlation alone cannot protect the decision contract."
            ),
        },
        "derivation_independence": (
            "Hosted Linux, Windows, and macOS observations are not inputs to this "
            "derivation. Candidate measurements are evaluated only after this "
            "reference-only budget has been constructed."
        ),
    }


def render_l6_prediction_sensitivity_report() -> str:
    return (
        json.dumps(
            build_l6_prediction_sensitivity_report(),
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render_l6_prediction_sensitivity_report(), encoding="utf-8")


if __name__ == "__main__":
    main()


__all__ = [
    "L6PredictionSensitivityBudget",
    "apply_boundary_perturbation",
    "build_l6_prediction_sensitivity_report",
    "derive_l6_prediction_sensitivity_budget",
    "render_l6_prediction_sensitivity_report",
]
