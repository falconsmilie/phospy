from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

from ..datasets.schema import DatasetSchema
from ..internal.constants import ComparisonSpec
from ..internal.defaults import (
    DEFAULT_KINASE_ACTIVITY_MIN_SUBSTRATES,
    DEFAULT_KINASE_ACTIVITY_THRESHOLD,
    DEFAULT_KINASE_ACTIVITY_TOP_N_SUBSTRATES,
    DEFAULT_PREDICTION_ENSEMBLE_SIZE,
    DEFAULT_PREDICTION_INCLUSION,
    DEFAULT_PREDICTION_N_ITERATIONS,
    DEFAULT_PREDICTION_SCORE_THRESHOLD,
    DEFAULT_PREDICTION_TOP,
)
from ..internal.types import (
    SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    SIGNALOME_KINASE_NETWORK_POLICY_POSITIVE_ONLY,
    PredictionSvmMode,
    SignalomeAssignmentPolicy,
    SignalomeKinaseNetworkPolicy,
)
from ..profiles import KinaseProfilePolicy
from ..signalomes import SignalomeModuleSelectionPolicy
from ..validation.values.enums import (
    validate_kinase_network_policy,
    validate_signalome_assignment_policy,
    validate_svm_mode,
)
from ..validation.values.numeric import (
    validate_fraction,
    validate_positive_int,
)

__all__ = [
    "DatasetLoadOptions",
    "KinaseActivityConfig",
    "PredictionRunConfig",
    "SignalomeRunConfig",
]


@dataclass(frozen=True, slots=True)
class DatasetLoadOptions:
    """Dataset loading options shared by high-level public workflows."""

    phospho_encoding: str | None = None
    schema: DatasetSchema = field(default_factory=DatasetSchema)
    comparisons: tuple[ComparisonSpec, ...] | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.schema, DatasetSchema):
            msg = "schema must be a DatasetSchema instance"
            raise TypeError(msg)

    @classmethod
    def from_value(cls, value: object) -> DatasetLoadOptions:
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if isinstance(value, Mapping):
            data = dict(value)
            comparisons = data.get("comparisons")
            if comparisons is not None:
                data["comparisons"] = tuple(comparisons)
            return cls(**data)
        msg = "dataset_options must be a DatasetLoadOptions or mapping"
        raise TypeError(msg)


@dataclass(frozen=True, slots=True)
class PredictionRunConfig:
    """Prediction and kinase-workflow tuning policy for public workflows."""

    min_substrates: int = 1
    min_motif_size: int = 1
    allow_profile_only_fallback: bool = False
    ensemble_size: int = DEFAULT_PREDICTION_ENSEMBLE_SIZE
    top: int = DEFAULT_PREDICTION_TOP
    score_threshold: float = DEFAULT_PREDICTION_SCORE_THRESHOLD
    inclusion: int = DEFAULT_PREDICTION_INCLUSION
    n_iterations: int = DEFAULT_PREDICTION_N_ITERATIONS
    random_state: int | None = None
    svm_mode: PredictionSvmMode | None = None
    profile_policy: KinaseProfilePolicy = field(default_factory=KinaseProfilePolicy)

    def __post_init__(self) -> None:
        validate_positive_int(self.min_substrates, name="min_substrates")
        validate_positive_int(self.min_motif_size, name="min_motif_size")
        validate_positive_int(self.ensemble_size, name="ensemble_size")
        validate_positive_int(self.top, name="top")
        validate_fraction(self.score_threshold, name="score_threshold")
        validate_positive_int(self.inclusion, name="inclusion")
        validate_positive_int(self.n_iterations, name="n_iterations")
        if self.svm_mode is not None:
            validate_svm_mode(self.svm_mode)
        KinaseProfilePolicy.from_value(self.profile_policy)

    @classmethod
    def from_value(cls, value: object) -> PredictionRunConfig:
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if isinstance(value, Mapping):
            return cls(**dict(value))
        msg = "prediction_config must be a PredictionRunConfig or mapping"
        raise TypeError(msg)


@dataclass(frozen=True, slots=True)
class KinaseActivityConfig:
    """Downstream kinase activity analysis options."""

    threshold: float = DEFAULT_KINASE_ACTIVITY_THRESHOLD
    min_substrates: int = DEFAULT_KINASE_ACTIVITY_MIN_SUBSTRATES
    top_n_substrates: int = DEFAULT_KINASE_ACTIVITY_TOP_N_SUBSTRATES

    def __post_init__(self) -> None:
        validate_fraction(self.threshold, name="threshold")
        validate_positive_int(self.min_substrates, name="min_substrates")
        validate_positive_int(self.top_n_substrates, name="top_n_substrates")

    @classmethod
    def from_value(cls, value: object) -> KinaseActivityConfig:
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if isinstance(value, Mapping):
            return cls(**dict(value))
        msg = "activity_config must be a KinaseActivityConfig or mapping"
        raise TypeError(msg)


@dataclass(frozen=True, slots=True)
class SignalomeRunConfig:
    """Signalome execution policy for public signalome workflows."""

    kinase_network_threshold: float = 0.9
    kinase_network_policy: SignalomeKinaseNetworkPolicy = (
        SIGNALOME_KINASE_NETWORK_POLICY_POSITIVE_ONLY
    )
    assignment_policy: SignalomeAssignmentPolicy = (
        SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY
    )
    signalome_cutoff: float = 0.5
    module_count: int | None = None
    min_kinase_module_share_percent: float = 1.0
    module_selection_policy: SignalomeModuleSelectionPolicy = field(
        default_factory=SignalomeModuleSelectionPolicy
    )

    def __post_init__(self) -> None:
        validate_fraction(
            self.kinase_network_threshold, name="kinase_network_threshold"
        )
        validate_kinase_network_policy(self.kinase_network_policy)
        validate_signalome_assignment_policy(self.assignment_policy)
        validate_fraction(self.signalome_cutoff, name="signalome_cutoff")
        if self.module_count is not None:
            validate_positive_int(self.module_count, name="module_count")
        if (
            not isinstance(self.min_kinase_module_share_percent, (int, float))
            or self.min_kinase_module_share_percent < 0
        ):
            msg = "min_kinase_module_share_percent must be a non-negative number"
            raise TypeError(msg)
        SignalomeModuleSelectionPolicy.from_value(self.module_selection_policy)

    @classmethod
    def from_value(cls, value: object) -> SignalomeRunConfig:
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if isinstance(value, Mapping):
            return cls(**dict(value))
        msg = "config must be a SignalomeRunConfig or mapping"
        raise TypeError(msg)
