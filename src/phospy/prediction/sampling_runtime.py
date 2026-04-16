from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

from ..internal.types import PREDICTION_SVM_MODE_DEFAULT
from .policies import PredictionSamplingPolicy, resolve_prediction_sampling_policy
from .trace_replay import PredictionSamplingTrace

_PREDICTION_STREAM_SEED_HIGH_EXCLUSIVE = (2**31) - 1
_PREDICTION_STREAM_SEED_DIGEST_SIZE_BYTES = 16
_DEFAULT_RESAMPLING_WEIGHT_EXPONENT = 0.8


def make_prediction_random_generators(
    rng: np.random.Generator,
) -> tuple[np.random.Generator, np.random.Generator]:
    """Create independent RNG streams for prediction sampling steps."""

    return (
        np.random.default_rng(
            int(rng.integers(0, _PREDICTION_STREAM_SEED_HIGH_EXCLUSIVE))
        ),
        np.random.default_rng(
            int(rng.integers(0, _PREDICTION_STREAM_SEED_HIGH_EXCLUSIVE))
        ),
    )


def _derive_prediction_stream_seed(
    *,
    random_state: int,
    kinase: str,
    stream_name: str,
) -> int:
    seed_material = f"{random_state}:{kinase}:{stream_name}".encode()
    digest = hashlib.blake2b(
        seed_material,
        digest_size=_PREDICTION_STREAM_SEED_DIGEST_SIZE_BYTES,
    ).digest()
    return int.from_bytes(digest, "little")


def make_kinase_prediction_random_generators(
    *,
    random_state: int | None,
    kinase: str,
) -> tuple[np.random.Generator, np.random.Generator]:
    """Create deterministic per-kinase RNG streams when a base seed is provided."""

    if random_state is None:
        return make_prediction_random_generators(np.random.default_rng())

    return (
        np.random.default_rng(
            _derive_prediction_stream_seed(
                random_state=random_state,
                kinase=kinase,
                stream_name="negative_sampling",
            )
        ),
        np.random.default_rng(
            _derive_prediction_stream_seed(
                random_state=random_state,
                kinase=kinase,
                stream_name="resampling",
            )
        ),
    )


def normalize_probabilities(values: np.ndarray) -> np.ndarray | None:
    """Normalize finite probability-like weights to sum to one."""

    total = float(np.nansum(values))
    if total <= 0.0 or not np.isfinite(total):
        return None
    return values / total


def transform_resampling_probabilities(
    values: np.ndarray,
    *,
    sampling_policy: PredictionSamplingPolicy | None = None,
    svm_mode: str | None = None,
) -> np.ndarray:
    """Adjust class resampling weights before normalization.

    In the native non-replayed path, scikit-learn's probability calibration can
    be slightly peakier than the e1071 path used by PhosR. Flattening the
    resampling weights a little reduces learner-side jitter in the adaptive
    sampling loop without changing candidate selection or the replay override
    seam.
    """

    resolved_policy = sampling_policy
    if resolved_policy is None:
        if svm_mode is None:
            msg = "Either sampling_policy or svm_mode must be provided"
            raise ValueError(msg)
        resolved_policy = resolve_prediction_sampling_policy(svm_mode)

    weights = np.asarray(values, dtype=float)
    if resolved_policy.resampling_weight_mode == PREDICTION_SVM_MODE_DEFAULT:
        return np.power(weights, _DEFAULT_RESAMPLING_WEIGHT_EXPONENT)
    return weights


def coerce_sampling_trace(
    sampling_trace: PredictionSamplingTrace | str | Path | None,
) -> PredictionSamplingTrace | None:
    """Resolve a sampling-trace input into a replay object."""

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
    """Map sampled site identifiers to positional indices in available rows."""

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
    """Validate replay override site IDs against the available sampling pool."""

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
    "make_kinase_prediction_random_generators",
    "make_prediction_random_generators",
    "normalize_probabilities",
    "resolve_sampled_site_positions",
    "transform_resampling_probabilities",
    "validate_override_sites",
]
