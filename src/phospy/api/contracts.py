from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from typing import Literal, TypeAlias

from ..datasets.schema import DatasetSchema
from ..internal.constants import (
    SIMPLE_KINASE_WORKFLOW_BUNDLE_FORMAT as _SIMPLE_KINASE_WORKFLOW_BUNDLE_FORMAT,
)
from ..internal.constants import (
    SIMPLE_KINASE_WORKFLOW_RESULT_TYPE as _SIMPLE_KINASE_WORKFLOW_RESULT_TYPE,
)
from ..internal.constants import (
    ComparisonSpec,
)
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
from ..prediction.profiles import KinaseProfilePolicy
from ..preprocessing import CorePreprocessingConfig
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
    "SimpleKinaseWorkflowBundleMetadata",
    "SimpleKinaseWorkflowConfigSnapshot",
    "SignalomeRunConfig",
    "WorkflowOutputInventoryItem",
]

BundleValueType: TypeAlias = Literal["dataframe", "series"]
SIMPLE_KINASE_WORKFLOW_BUNDLE_FORMAT: str = _SIMPLE_KINASE_WORKFLOW_BUNDLE_FORMAT
SIMPLE_KINASE_WORKFLOW_RESULT_TYPE: str = _SIMPLE_KINASE_WORKFLOW_RESULT_TYPE


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


@dataclass(frozen=True, slots=True)
class SimpleKinaseWorkflowConfigSnapshot:
    """Serializable workflow-configuration snapshot for output bundles."""

    dataset_options: dict[str, object]
    preprocessing_config: dict[str, object]
    prediction_config: dict[str, object]
    activity_config: dict[str, object]

    @classmethod
    def from_workflow_inputs(
        cls,
        *,
        dataset_options: DatasetLoadOptions | Mapping[str, object] | None = None,
        preprocessing_config: CorePreprocessingConfig | None = None,
        prediction_config: PredictionRunConfig | Mapping[str, object] | None = None,
        activity_config: KinaseActivityConfig | Mapping[str, object] | None = None,
    ) -> SimpleKinaseWorkflowConfigSnapshot:
        resolved_dataset_options = DatasetLoadOptions.from_value(dataset_options)
        resolved_preprocessing = (
            CorePreprocessingConfig()
            if preprocessing_config is None
            else preprocessing_config
        )
        resolved_prediction = PredictionRunConfig.from_value(prediction_config)
        resolved_activity = KinaseActivityConfig.from_value(activity_config)
        return cls(
            dataset_options={
                "phospho_encoding": resolved_dataset_options.phospho_encoding,
                "comparisons": (
                    [
                        [left, right]
                        for left, right in resolved_dataset_options.comparisons
                    ]
                    if resolved_dataset_options.comparisons is not None
                    else None
                ),
                "schema": {
                    "total_cols": list(resolved_dataset_options.schema.total_cols),
                    "phospho_cols": list(resolved_dataset_options.schema.phospho_cols),
                    "corrected_cols": list(
                        resolved_dataset_options.schema.corrected_cols
                    ),
                },
            },
            preprocessing_config=asdict(resolved_preprocessing),
            prediction_config=asdict(resolved_prediction),
            activity_config=asdict(resolved_activity),
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "dataset_options": dict(self.dataset_options),
            "preprocessing_config": dict(self.preprocessing_config),
            "prediction_config": dict(self.prediction_config),
            "activity_config": dict(self.activity_config),
        }


@dataclass(frozen=True, slots=True)
class WorkflowOutputInventoryItem:
    """One persisted output table entry in a workflow bundle manifest."""

    table_id: str
    path: str
    value_type: BundleValueType

    def __post_init__(self) -> None:
        if not self.table_id.strip():
            msg = "table_id must not be empty"
            raise ValueError(msg)
        if not self.path.strip():
            msg = "path must not be empty"
            raise ValueError(msg)
        if self.value_type not in ("dataframe", "series"):
            msg = f"Unsupported value_type {self.value_type!r}."
            raise ValueError(msg)

    @classmethod
    def from_mapping(cls, value: Mapping[str, object]) -> WorkflowOutputInventoryItem:
        table_id = value.get("table_id")
        path = value.get("path")
        value_type = value.get("value_type")
        if not isinstance(table_id, str):
            msg = "Output inventory item 'table_id' must be a string."
            raise TypeError(msg)
        if not isinstance(path, str):
            msg = "Output inventory item 'path' must be a string."
            raise TypeError(msg)
        if not isinstance(value_type, str):
            msg = "Output inventory item 'value_type' must be a string."
            raise TypeError(msg)
        return cls(table_id=table_id, path=path, value_type=value_type)

    def to_dict(self) -> dict[str, str]:
        return {
            "table_id": self.table_id,
            "path": self.path,
            "value_type": self.value_type,
        }


@dataclass(frozen=True, slots=True)
class SimpleKinaseWorkflowBundleMetadata:
    """Minimal manifest metadata for a saved workflow output bundle."""

    workflow_type: str
    bundle_format: str
    generated_at_utc: str
    package_version: str
    config_snapshot: dict[str, object]
    reference_identity: dict[str, object]
    output_inventory: tuple[WorkflowOutputInventoryItem, ...]

    def __post_init__(self) -> None:
        if not self.workflow_type.strip():
            msg = "workflow_type must not be empty"
            raise ValueError(msg)
        if not self.bundle_format.strip():
            msg = "bundle_format must not be empty"
            raise ValueError(msg)
        if not isinstance(self.config_snapshot, Mapping):
            msg = "config_snapshot must be a mapping"
            raise TypeError(msg)
        if not isinstance(self.reference_identity, Mapping):
            msg = "reference_identity must be a mapping"
            raise TypeError(msg)
        object.__setattr__(self, "config_snapshot", dict(self.config_snapshot))
        object.__setattr__(self, "reference_identity", dict(self.reference_identity))

    @classmethod
    def from_mapping(
        cls,
        value: Mapping[str, object],
    ) -> SimpleKinaseWorkflowBundleMetadata:
        workflow_type = value.get("workflow_type")
        bundle_format = value.get("bundle_format")
        generated_at_utc = value.get("generated_at_utc")
        package_version = value.get("package_version")
        config_snapshot = value.get("config_snapshot")
        reference_identity = value.get("reference_identity")
        output_inventory_raw = value.get("output_inventory")

        if not isinstance(workflow_type, str):
            msg = "Manifest 'workflow_type' must be a string."
            raise TypeError(msg)
        if not isinstance(bundle_format, str):
            msg = "Manifest 'bundle_format' must be a string."
            raise TypeError(msg)
        if not isinstance(generated_at_utc, str):
            msg = "Manifest 'generated_at_utc' must be a string."
            raise TypeError(msg)
        if not isinstance(package_version, str):
            msg = "Manifest 'package_version' must be a string."
            raise TypeError(msg)
        if not isinstance(config_snapshot, Mapping):
            msg = "Manifest 'config_snapshot' must be a mapping."
            raise TypeError(msg)
        if not isinstance(reference_identity, Mapping):
            msg = "Manifest 'reference_identity' must be a mapping."
            raise TypeError(msg)
        if not isinstance(output_inventory_raw, list):
            msg = "Manifest 'output_inventory' must be a list."
            raise TypeError(msg)

        inventory = tuple(
            WorkflowOutputInventoryItem.from_mapping(item)
            for item in output_inventory_raw
            if isinstance(item, Mapping)
        )
        if len(inventory) != len(output_inventory_raw):
            msg = "Manifest 'output_inventory' entries must be mappings."
            raise TypeError(msg)

        return cls(
            workflow_type=workflow_type,
            bundle_format=bundle_format,
            generated_at_utc=generated_at_utc,
            package_version=package_version,
            config_snapshot=dict(config_snapshot),
            reference_identity=dict(reference_identity),
            output_inventory=inventory,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "workflow_type": self.workflow_type,
            "bundle_format": self.bundle_format,
            "generated_at_utc": self.generated_at_utc,
            "package_version": self.package_version,
            "config_snapshot": dict(self.config_snapshot),
            "reference_identity": dict(self.reference_identity),
            "output_inventory": [item.to_dict() for item in self.output_inventory],
        }
