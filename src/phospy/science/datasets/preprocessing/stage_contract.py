"""Stage-owned preprocessing contract primitives."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from phospy.errors.build import DatasetBuildError
from phospy.provenance.models import DeterminismKind
from phospy.science.datasets.preprocessing.models import (
    PREPROCESSING_STATE_TABLE_KEYS,
    PreprocessingPlan,
    PreprocessingStage,
    PreprocessingStateTableKey,
)

if TYPE_CHECKING:
    from phospy.science.datasets.preprocessing.stages.batch_correction import (
        BatchCorrectionAdequacyValidatorProtocol,
        BatchDesignMetadataValidatorProtocol,
        SpsRuvStyleBatchCorrectionRunner,
    )

_ParameterSerializer = Callable[[PreprocessingPlan], dict[str, object]]
_OperationResolver = Callable[[PreprocessingPlan], str]
_PlanValidator = Callable[[PreprocessingPlan], None]
_RandomSeedResolver = Callable[[Mapping[str, object], str], int | None]
_DeterminismDeclaration = (
    DeterminismKind | str | Callable[[PreprocessingPlan], DeterminismKind | str]
)


def _always_include(_plan: PreprocessingPlan) -> bool:
    return True


def _no_plan_validation(_plan: PreprocessingPlan) -> None:
    return None


def _default_random_seed_resolver(
    diagnostics: Mapping[str, object],
    _stage_key: str,
) -> int | None:
    value = diagnostics.get("random_seed")
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        return None
    return int(value)


@dataclass(frozen=True, slots=True)
class InterpretedPreprocessingStageContract:
    """Resolved contract payload for one stage in one preprocessing plan."""

    stage: str
    operation: str
    parameters: dict[str, object]
    consumed_input_tables: tuple[PreprocessingStateTableKey, ...]
    produced_output_tables: tuple[PreprocessingStateTableKey, ...]
    determinism_kind: DeterminismKind
    backend: str | None = None


@dataclass(frozen=True, slots=True)
class PreprocessingStageFactoryContext:
    """Composition-time collaborators for preprocessing stage factories."""

    batch_correction_runner: SpsRuvStyleBatchCorrectionRunner | None = None
    batch_correction_metadata_validator: BatchDesignMetadataValidatorProtocol | None = (
        None
    )
    batch_correction_adequacy_validator: (
        BatchCorrectionAdequacyValidatorProtocol | None
    ) = None


_StageFactory = Callable[[PreprocessingStageFactoryContext], PreprocessingStage]


@dataclass(frozen=True, slots=True)
class PreprocessingStageContract:
    """Single shared contract for a preprocessing stage."""

    stage_key: str
    display_label: str
    operation_name: _OperationResolver
    serialize_parameters: _ParameterSerializer
    consumed_input_tables: tuple[PreprocessingStateTableKey, ...]
    produced_output_tables: tuple[PreprocessingStateTableKey, ...]
    stage_factory: _StageFactory | None = None
    provenance_stage: str | None = None
    backend: str | None = None
    determinism_kind: _DeterminismDeclaration = DeterminismKind.DETERMINISTIC
    include_in_builder_provenance: bool = True
    include_when: Callable[[PreprocessingPlan], bool] = field(default=_always_include)
    diagnostics_metadata: Mapping[str, object] = field(default_factory=dict)
    validate_plan: _PlanValidator = field(default=_no_plan_validation)
    resolve_random_seed: _RandomSeedResolver = field(
        default=_default_random_seed_resolver
    )

    def __post_init__(self) -> None:
        normalized_stage_key = str(self.stage_key).strip()
        object.__setattr__(self, "stage_key", normalized_stage_key)
        object.__setattr__(
            self,
            "consumed_input_tables",
            _normalize_stage_table_keys(
                stage_key=normalized_stage_key,
                table_keys=self.consumed_input_tables,
                role="consumed_input_tables",
            ),
        )
        if not callable(self.determinism_kind):
            object.__setattr__(
                self,
                "determinism_kind",
                _normalize_determinism_kind(
                    self.determinism_kind,
                    stage_key=normalized_stage_key,
                ),
            )
        object.__setattr__(
            self,
            "produced_output_tables",
            _normalize_stage_table_keys(
                stage_key=normalized_stage_key,
                table_keys=self.produced_output_tables,
                role="produced_output_tables",
            ),
        )

    @property
    def provenance_stage_key(self) -> str:
        if self.provenance_stage is None:
            return self.stage_key
        normalized = self.provenance_stage.strip()
        if not normalized:
            return self.stage_key
        return normalized

    def interpret(
        self,
        plan: PreprocessingPlan,
    ) -> InterpretedPreprocessingStageContract:
        self.validate_plan(plan)
        return InterpretedPreprocessingStageContract(
            stage=self.provenance_stage_key,
            operation=self.operation_name(plan),
            parameters=self.serialize_parameters(plan),
            consumed_input_tables=self.consumed_input_tables,
            produced_output_tables=self.produced_output_tables,
            determinism_kind=self.resolve_determinism_kind(plan),
            backend=self.backend,
        )

    def resolve_determinism_kind(self, plan: PreprocessingPlan) -> DeterminismKind:
        declaration = self.determinism_kind
        resolved = declaration(plan) if callable(declaration) else declaration
        return _normalize_determinism_kind(
            resolved,
            stage_key=self.stage_key,
        )


def normalize_stage_table_keys(
    *,
    stage_key: str,
    table_keys: tuple[PreprocessingStateTableKey, ...],
    role: str,
) -> tuple[PreprocessingStateTableKey, ...]:
    """Public helper used by registry validation paths."""

    return _normalize_stage_table_keys(
        stage_key=stage_key,
        table_keys=table_keys,
        role=role,
    )


def _normalize_stage_table_keys(
    *,
    stage_key: str,
    table_keys: tuple[PreprocessingStateTableKey, ...],
    role: str,
) -> tuple[PreprocessingStateTableKey, ...]:
    normalized: list[PreprocessingStateTableKey] = []
    for index, table_key in enumerate(table_keys):
        if isinstance(table_key, PreprocessingStateTableKey):
            normalized.append(table_key)
            continue
        if not isinstance(table_key, str):
            raise DatasetBuildError(
                "dataset preprocessing stage metadata contains non-string table key: "
                f"stage={stage_key or '<empty>'!r}, field={role!r}[{index}], "
                f"got {table_key!r} ({type(table_key).__name__})"
            )
        try:
            normalized.append(PreprocessingStateTableKey(table_key))
        except ValueError as exc:
            supported = ", ".join(key.value for key in PREPROCESSING_STATE_TABLE_KEYS)
            raise DatasetBuildError(
                "dataset preprocessing stage metadata contains unknown table key: "
                f"stage={stage_key or '<empty>'!r}, field={role!r}[{index}], "
                f"table={table_key!r}, supported tables: {supported}"
            ) from exc
    return tuple(normalized)


def _normalize_determinism_kind(
    value: DeterminismKind | str,
    *,
    stage_key: str,
) -> DeterminismKind:
    if isinstance(value, DeterminismKind):
        return value
    if isinstance(value, str):
        normalized = value.strip()
        aliases = {
            "pure": DeterminismKind.DETERMINISTIC,
            "external_dependency": DeterminismKind.EXTERNALLY_NONDETERMINISTIC,
        }
        alias = aliases.get(normalized)
        if alias is not None:
            return alias
        try:
            return DeterminismKind(normalized)
        except ValueError as exc:
            supported = ", ".join(item.value for item in DeterminismKind)
            raise DatasetBuildError(
                "dataset preprocessing stage metadata contains unknown "
                "determinism kind: "
                f"stage={stage_key or '<empty>'!r}, determinism_kind={value!r}; "
                f"supported determinism kinds: {supported}"
            ) from exc
    raise DatasetBuildError(
        "dataset preprocessing stage metadata contains invalid determinism kind: "
        f"stage={stage_key or '<empty>'!r}, got {value!r} "
        f"({type(value).__name__})"
    )


__all__ = [
    "DeterminismKind",
    "InterpretedPreprocessingStageContract",
    "PreprocessingStageContract",
    "normalize_stage_table_keys",
]
