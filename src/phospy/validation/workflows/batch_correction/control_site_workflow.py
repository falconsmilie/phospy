"""Control-site validation wrapper for batch-correction workflows."""

from __future__ import annotations

from phospy.science.datasets.preprocessing.control_sites import ControlSiteMapping
from phospy.validation.workflows.batch_correction.control_sites import (
    ControlSiteEligibilityValidator,
)
from phospy.validation.workflows.batch_correction.protocols import (
    BatchCorrectionWorkflowRequestProtocol,
)

_MIN_ELIGIBLE_CONTROLS = 2


class BatchCorrectionWorkflowControlSiteValidator:
    """Validate caller control-site mappings for supported workflow execution."""

    def __init__(
        self,
        *,
        eligibility_validator: ControlSiteEligibilityValidator | None = None,
    ) -> None:
        self._eligibility_validator = (
            eligibility_validator or ControlSiteEligibilityValidator()
        )

    def run(
        self, *, request: BatchCorrectionWorkflowRequestProtocol
    ) -> ControlSiteMapping:
        config = request.config
        return self._eligibility_validator.run(
            control_set=request.control_site_set,
            method=config.method.value,
            min_eligible_controls=_MIN_ELIGIBLE_CONTROLS,
            phospho=request.phospho,
            site_metadata=request.site_metadata,
            dataset_organism=request.dataset_organism,
            n_unwanted_factors=config.n_unwanted_factors,
            control_site_source_type=config.control_site_source.value,
            supports_weights=True,
            supports_groups=False,
            supports_weighted_groups=False,
        )


__all__ = ["BatchCorrectionWorkflowControlSiteValidator"]
