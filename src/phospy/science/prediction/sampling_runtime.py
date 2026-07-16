"""Runtime helpers for adaptive sampling RNG and resampling weights."""

from __future__ import annotations

import hashlib

import numpy as np

from phospy.contracts.configs.prediction import KinaseAdaptivePolicy
from phospy.science.prediction.policies import (
    PredictionSamplingPolicy,
    resolve_prediction_sampling_policy,
)

_PREDICTION_STREAM_SEED_HIGH_EXCLUSIVE = (2**31) - 1
_PREDICTION_STREAM_SEED_DIGEST_SIZE_BYTES = 16
_DEFAULT_RESAMPLING_WEIGHT_EXPONENT = 0.8


class PredictionSamplingRandomSource:
    """Resolve per-kinase RNG streams from the selected sampling policy."""

    def __init__(
        self,
        *,
        policy: PredictionSamplingPolicy,
        random_state: int,
    ) -> None:
        self.policy = policy
        self.random_state = random_state
        self._global_rng: np.random.Generator | None = None
        if self.policy.seed_strategy == "global_parity":
            self._global_rng = np.random.default_rng(random_state)

    def generators_for_kinase(
        self,
        *,
        kinase: str,
    ) -> tuple[np.random.Generator, np.random.Generator]:
        if self.policy.seed_strategy == "stable_by_kinase":
            return make_kinase_prediction_random_generators(
                random_state=self.random_state,
                kinase=kinase,
            )

        if self._global_rng is None:
            self._global_rng = np.random.default_rng(self.random_state)
        return make_prediction_random_generators(self._global_rng)


def make_prediction_random_generators(
    rng: np.random.Generator,
) -> tuple[np.random.Generator, np.random.Generator]:
    """Create independent RNG streams for negatives and resampling."""

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
    random_state: int,
    kinase: str,
) -> tuple[np.random.Generator, np.random.Generator]:
    """Create deterministic per-kinase RNG streams from a base seed."""

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
    """Normalize finite, non-negative weights to a probability vector."""

    total = float(np.nansum(values))
    if total <= 0.0 or not np.isfinite(total):
        return None
    return values / total


def transform_resampling_probabilities(
    values: np.ndarray,
    *,
    sampling_policy: PredictionSamplingPolicy | None = None,
    adaptive_policy: KinaseAdaptivePolicy | None = None,
) -> np.ndarray:
    """Adjust class resampling weights before normalization."""

    resolved_policy = sampling_policy
    if resolved_policy is None:
        if adaptive_policy is None:
            msg = "Either sampling_policy or adaptive_policy must be provided"
            raise ValueError(msg)
        resolved_policy = resolve_prediction_sampling_policy(adaptive_policy)

    weights = np.asarray(values, dtype=float)
    if resolved_policy.resampling_weight_mode == "default":
        return np.power(weights, _DEFAULT_RESAMPLING_WEIGHT_EXPONENT)
    return weights


__all__ = [
    "PredictionSamplingRandomSource",
    "make_kinase_prediction_random_generators",
    "make_prediction_random_generators",
    "normalize_probabilities",
    "transform_resampling_probabilities",
]
