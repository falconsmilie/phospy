from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..types import PredictionSvmMode
from .trace_replay import PredictionSamplingTrace


def make_prediction_random_generators(
    rng: np.random.Generator,
) -> tuple[np.random.Generator, np.random.Generator]:
    """Create independent RNG streams for prediction sampling steps."""

    return (
        np.random.default_rng(int(rng.integers(0, 2**31 - 1))),
        np.random.default_rng(int(rng.integers(0, 2**31 - 1))),
    )


def normalize_probabilities(values: np.ndarray) -> np.ndarray | None:
    total = float(np.nansum(values))
    if total <= 0.0 or not np.isfinite(total):
        return None
    return values / total


def transform_resampling_probabilities(
    values: np.ndarray,
    *,
    svm_mode: PredictionSvmMode,
) -> np.ndarray:
    """Adjust class resampling weights before normalization.

    In the native non-replayed path, scikit-learn's probability calibration can
    be slightly peakier than the e1071 path used by PhosR. Flattening the
    resampling weights a little reduces learner-side jitter in the adaptive
    sampling loop without changing candidate selection or the replay override
    seam.
    """

    weights = np.asarray(values, dtype=float)
    if svm_mode == "default":
        return np.power(weights, 0.8)
    return weights


def coerce_sampling_trace(
    sampling_trace: PredictionSamplingTrace | str | Path | None,
) -> PredictionSamplingTrace | None:
    if sampling_trace is None:
        return None
    if isinstance(sampling_trace, PredictionSamplingTrace):
        return sampling_trace
    return PredictionSamplingTrace.from_trace_directory(sampling_trace)


def resolve_sampled_site_positions(
    *,
    available_sites: pd.Index,
    sampled_sites: list[str],
    expected_size: int,
    context: str,
) -> np.ndarray:
    sampled_site_list = validate_override_sites(
        available_sites=available_sites,
        sampled_sites=sampled_sites,
        expected_size=expected_size,
        context=context,
    )
    position_lookup: dict[str, int] = {}
    for position, site in enumerate(available_sites.astype(str).tolist()):
        position_lookup.setdefault(site, position)
    return np.asarray([position_lookup[site] for site in sampled_site_list], dtype=int)


def validate_override_sites(
    *,
    available_sites: pd.Index,
    sampled_sites: list[str],
    expected_size: int,
    context: str,
) -> list[str]:
    sampled_site_list = [str(site) for site in sampled_sites]
    if len(sampled_site_list) != expected_size:
        msg = (
            f"Sampling override for {context} has {len(sampled_site_list)} rows; "
            f"expected {expected_size}"
        )
        raise ValueError(msg)
    available_site_set = set(available_sites.astype(str).tolist())
    invalid_sites = [
        site for site in sampled_site_list if site not in available_site_set
    ]
    if invalid_sites:
        unique_invalid_sites = sorted(set(invalid_sites))
        msg = (
            f"Sampling override for {context} contains sites outside the "
            f"available training rows: {', '.join(unique_invalid_sites)}"
        )
        raise ValueError(msg)
    return sampled_site_list


__all__ = [
    "coerce_sampling_trace",
    "make_prediction_random_generators",
    "normalize_probabilities",
    "resolve_sampled_site_positions",
    "transform_resampling_probabilities",
    "validate_override_sites",
]
