"""Shared workflow-config validation."""

from __future__ import annotations

from phospy.contracts.configs import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    SignalomeConfig,
)
from phospy.errors.validation import WorkflowValidationError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.transformations.models import QuantitativeMeaning

_MIXED_QUANTITATIVE_MEANING = QuantitativeMeaning.MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE.value


class KinaseWorkflowConfigValidator:
    """Validate kinase workflow config objects and local invariants."""

    def run(
        self,
        *,
        scoring_config: object,
        prediction_config: object,
        activity_config: object | None,
    ) -> tuple[
        KinaseScoringConfig,
        KinasePredictionConfig,
        KinaseActivityConfig | None,
    ]:
        validated_scoring = self._validate_scoring(scoring_config)
        validated_prediction = self._validate_prediction(prediction_config)
        validated_activity = self._validate_activity(activity_config)
        return validated_scoring, validated_prediction, validated_activity

    @staticmethod
    def _validate_scoring(config: object) -> KinaseScoringConfig:
        if not isinstance(config, KinaseScoringConfig):
            raise WorkflowValidationError(
                "kinase workflow request scoring_config must be KinaseScoringConfig"
            )
        return config

    @staticmethod
    def _validate_prediction(config: object) -> KinasePredictionConfig:
        if not isinstance(config, KinasePredictionConfig):
            raise WorkflowValidationError(
                "kinase workflow request prediction_config must be KinasePredictionConfig"
            )
        return config

    @staticmethod
    def _validate_activity(config: object | None) -> KinaseActivityConfig | None:
        if config is None:
            return None
        if not isinstance(config, KinaseActivityConfig):
            raise WorkflowValidationError(
                "kinase workflow request activity_config must be KinaseActivityConfig or None"
            )
        return config


class SignalomeConfigValidator:
    """Validate signalome workflow config objects and local invariants."""

    def run(self, config: object) -> SignalomeConfig:
        if not isinstance(config, SignalomeConfig):
            raise WorkflowValidationError(
                "signalome workflow request config must be SignalomeConfig"
            )
        return config


def reject_mixed_total_protein_quantitative_meaning(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    allow_mixed: bool,
    context: str,
) -> None:
    quantity = dataset.intensity_scale_state.quantity
    if quantity is None or quantity.value != _MIXED_QUANTITATIVE_MEANING:
        return
    if allow_mixed:
        return
    diagnostics = dataset.processing_state.total_protein_correction.diagnostics
    uncorrected_row_count = None
    unmatched_policy = None
    if diagnostics is not None:
        uncorrected_row_count = diagnostics.get("uncorrected_row_count")
        unmatched_policy = diagnostics.get("unmatched_policy")
    raise WorkflowValidationError(
        f"{context} received a dataset with mixed total-protein quantitative meaning "
        f"({_MIXED_QUANTITATIVE_MEANING}). "
        f"uncorrected_rows={uncorrected_row_count!r}, "
        f"unmatched_policy={unmatched_policy!r}. "
        "This usually happens when total-protein correction is applied with "
        "unmatched_policy='allow_uncorrected'. "
        "Recommended actions: set unmatched_policy='error' or complete the "
        "phosphosite-to-total mapping. "
        "If you intentionally want this mixed dataset, set the workflow "
        "allow_mixed_total_protein_quantitative_meaning option to True."
    )


__all__ = [
    "KinaseWorkflowConfigValidator",
    "SignalomeConfigValidator",
    "reject_mixed_total_protein_quantitative_meaning",
]
