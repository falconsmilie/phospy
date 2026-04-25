"""Internal dataset preprocessing pipeline orchestration."""

from __future__ import annotations

from phospy.datasets.preprocessing.models import (
    PreprocessingStage,
    PreprocessingStageExecution,
    PreprocessingState,
)
from phospy.datasets.preprocessing.stages.comparisons import ComparisonsStage
from phospy.datasets.preprocessing.stages.intensity_transform import (
    IntensityTransformStage,
)
from phospy.datasets.preprocessing.stages.missing_data import MissingDataStage
from phospy.datasets.preprocessing.stages.normalisation import NormalisationStage
from phospy.datasets.preprocessing.stages.site_matrix import SiteMatrixStage
from phospy.datasets.preprocessing.stages.total_protein_correction import (
    TotalProteinCorrectionStage,
)
from phospy.errors.build import DatasetBuildError


class PreprocessingPipeline:
    """Apply ordered preprocessing stages for interpreted dataset input."""

    def __init__(
        self,
        *,
        stage_registry: tuple[PreprocessingStage, ...] | None = None,
    ) -> None:
        stages = stage_registry or (
            MissingDataStage(),
            TotalProteinCorrectionStage(),
            SiteMatrixStage(),
            IntensityTransformStage(),
            NormalisationStage(),
            ComparisonsStage(),
        )
        self._stages_by_key = {stage.stage_key: stage for stage in stages}
        if len(self._stages_by_key) != len(stages):
            raise DatasetBuildError(
                "dataset preprocessing stage registry contains duplicate stage keys"
            )

    def run(self, state: PreprocessingState) -> PreprocessingState:
        final_state, _ = self.run_with_trace(state)
        return final_state

    def run_with_trace(
        self,
        state: PreprocessingState,
    ) -> tuple[PreprocessingState, tuple[PreprocessingStageExecution, ...]]:
        current = state
        trace: list[PreprocessingStageExecution] = []
        for stage_key in current.plan.stage_order:
            stage = self._stages_by_key.get(stage_key)
            if stage is None:
                raise DatasetBuildError(
                    "dataset preprocessing plan references an unsupported stage: "
                    f"{stage_key}"
                )
            input_rows = int(len(current.phospho.index))
            current = stage.run(current)
            output_rows = int(len(current.phospho.index))
            trace.append(
                PreprocessingStageExecution(
                    stage=stage_key,
                    input_rows=input_rows,
                    output_rows=output_rows,
                )
            )
        return current, tuple(trace)


__all__ = ["PreprocessingPipeline"]
