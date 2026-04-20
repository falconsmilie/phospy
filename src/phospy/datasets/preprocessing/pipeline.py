"""Internal dataset preprocessing pipeline orchestration."""

from __future__ import annotations

from phospy.datasets.preprocessing.models import PreprocessingStage, PreprocessingState
from phospy.datasets.preprocessing.stages.comparisons import ComparisonsStage
from phospy.datasets.preprocessing.stages.missing_data import MissingDataStage
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
            ComparisonsStage(),
        )
        self._stages_by_key = {stage.stage_key: stage for stage in stages}
        if len(self._stages_by_key) != len(stages):
            raise DatasetBuildError(
                "dataset preprocessing stage registry contains duplicate stage keys"
            )

    def run(self, state: PreprocessingState) -> PreprocessingState:
        current = state
        for stage_key in current.plan.stage_order:
            stage = self._stages_by_key.get(stage_key)
            if stage is None:
                raise DatasetBuildError(
                    "dataset preprocessing plan references an unsupported stage: "
                    f"{stage_key}"
                )
            current = stage.run(current)
        return current


__all__ = ["PreprocessingPipeline"]
