"""Private preprocessing stage trace assembly service."""

from __future__ import annotations

from phospy.errors.build import DatasetBuildError
from phospy.provenance.models import (
    PREPROCESSING_EXTERNAL_NONDETERMINISM_CAVEAT_CODE,
    PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3,
    DeterminismKind,
    ReproducibilityCaveat,
)
from phospy.science.datasets.preprocessing.diagnostics_normalization import (
    _coerce_int,
    _NormalizedStageDiagnostics,
)
from phospy.science.datasets.preprocessing.fingerprints import _StageFingerprintBundle
from phospy.science.datasets.preprocessing.models import (
    PreprocessingStageExecution,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.stage_contract import (
    InterpretedPreprocessingStageContract,
    PreprocessingStageContract,
)
from phospy.science.transformations.models import IntensityTransformationEvent


class _StageTraceBuilder:
    """Build immutable trace records from already executed stage outputs."""

    def run(
        self,
        *,
        stage_key: str,
        contract: PreprocessingStageContract,
        interpreted_contract: InterpretedPreprocessingStageContract,
        previous: PreprocessingState,
        current: PreprocessingState,
        stage_result: PreprocessingStageResult,
        diagnostics: _NormalizedStageDiagnostics,
        fingerprints: _StageFingerprintBundle,
        intensity_transformation_event: IntensityTransformationEvent | None,
    ) -> PreprocessingStageExecution:
        random_seed = _resolve_random_seed(
            stage_key=stage_key,
            value=contract.resolve_random_seed(diagnostics.diagnostics, stage_key),
        )
        determinism = _resolve_stage_determinism(
            stage_key=stage_key,
            determinism_kind=interpreted_contract.determinism_kind,
            random_seed=random_seed,
        )
        reproducibility_caveats = _build_stage_reproducibility_caveats(
            stage=interpreted_contract.stage,
            operation=interpreted_contract.operation,
            backend=interpreted_contract.backend,
            determinism=determinism,
        )
        return PreprocessingStageExecution(
            stage=interpreted_contract.stage,
            operation=interpreted_contract.operation,
            parameters=interpreted_contract.parameters,
            input_shape=(
                int(previous.phospho.shape[0]),
                int(previous.phospho.shape[1]),
            ),
            output_shape=(
                int(current.phospho.shape[0]),
                int(current.phospho.shape[1]),
            ),
            input_hash=fingerprints.input_hash,
            output_hash=fingerprints.output_hash,
            phospho_input_hash=fingerprints.phospho_input_hash,
            phospho_output_hash=fingerprints.phospho_output_hash,
            dropped_row_ids=tuple(diagnostics.dropped_row_ids),
            dropped_row_count=int(diagnostics.dropped_row_count),
            schema_version=PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3,
            consumed_input_tables=fingerprints.consumed_input_tables,
            produced_output_tables=fingerprints.produced_output_tables,
            backend=interpreted_contract.backend,
            random_seed=random_seed,
            determinism=determinism,
            reproducibility_caveats=reproducibility_caveats,
            is_deterministic=determinism is DeterminismKind.DETERMINISTIC,
            imputed_cell_count=int(diagnostics.imputed_cell_count),
            imputed_row_ids=tuple(diagnostics.imputed_row_ids),
            notes=diagnostics.notes,
            diagnostics=dict(diagnostics.diagnostics),
            batch_correction_provenance=stage_result.batch_correction_provenance,
            intensity_transformation_event=intensity_transformation_event,
        )


def _resolve_random_seed(
    *,
    stage_key: str,
    value: object,
) -> int | None:
    if value is None:
        return None
    return _coerce_int(
        value,
        stage_key=stage_key,
        field_name="diagnostics.random_seed",
        default=0,
    )


def _resolve_stage_determinism(
    *,
    stage_key: str,
    determinism_kind: DeterminismKind,
    random_seed: int | None,
) -> DeterminismKind:
    if determinism_kind is DeterminismKind.SEEDED_STOCHASTIC and random_seed is None:
        raise DatasetBuildError(
            "dataset preprocessing stage declared seeded stochastic execution "
            "but did not record an explicit random seed: "
            f"stage={stage_key!r}, determinism_kind={determinism_kind.value!r}"
        )
    return determinism_kind


def _build_stage_reproducibility_caveats(
    *,
    stage: str,
    operation: str,
    backend: str | None,
    determinism: DeterminismKind,
) -> tuple[ReproducibilityCaveat, ...]:
    if determinism is not DeterminismKind.EXTERNALLY_NONDETERMINISTIC:
        return ()
    return (
        ReproducibilityCaveat(
            code=PREPROCESSING_EXTERNAL_NONDETERMINISM_CAVEAT_CODE,
            severity="warning",
            message=(
                "Preprocessing stage declares externally nondeterministic "
                "execution; exact reproduction requires the external system, "
                "runtime state, and inputs used to produce the recorded output."
            ),
            details={
                "stage": stage,
                "operation": operation,
                "backend": backend,
                "determinism_kind": determinism.value,
            },
        ),
    )
