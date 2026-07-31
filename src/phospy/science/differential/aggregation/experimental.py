"""Withdrawn compatibility route for peptide-to-site estimate combination.

This module is retained only so older import routes fail closed with a clear
scientific explanation. It must not execute post-hoc peptide-to-site
differential estimate combination until a future ADR defines coherent combined
effect/inference and executable mapping semantics.
"""

from __future__ import annotations

from phospy.errors.input import PhosPyInputError
from phospy.science.differential.aggregation.models import (
    PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS,
)

WITHDRAWN_POSTHOC_PEPTIDE_TO_SITE_AGGREGATION_MESSAGE = (
    "post-hoc peptide-to-site differential estimate combination has been "
    "withdrawn from public support. The lane fails closed because coherent "
    "combined effect/inference semantics and executable peptide-to-site mapping "
    "semantics are not yet implemented. Resolve peptide evidence at "
    "sample-intensity level before DifferentialAnalysisWorkflow. Future public "
    "support requires an ADR-backed model with executable mapping semantics and "
    "a coherent combined estimand and inferential result."
)

EXPERIMENTAL_INTERNAL_API = True
EXPERIMENTAL_INTERNAL_REASON = (
    "The post-hoc peptide-to-site differential estimate-combination lane was "
    "withdrawn from public support because coherent combined effect/inference "
    "and executable mapping semantics are not yet implemented."
)


class PeptideToSiteAggregator:
    """Unsupported shell retained only for fail-closed compatibility."""

    experimental_internal_api: bool = True
    scientific_support_status: str = PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS

    def __init__(
        self,
        *,
        executor: object | None = None,
    ) -> None:
        self._executor = executor

    def run(
        self,
        estimates: object,
        *,
        config: object | None = None,
        contrast_name: str = "aggregated",
    ) -> object:
        raise PhosPyInputError(WITHDRAWN_POSTHOC_PEPTIDE_TO_SITE_AGGREGATION_MESSAGE)

    def run_estimates(
        self,
        *,
        estimates: object,
        config: object | None = None,
        contrast_name: str = "aggregated",
    ) -> object:
        return self.run(
            estimates,
            config=config,
            contrast_name=contrast_name,
        )

    def run_table(
        self,
        *,
        estimate_table: object | None = None,
        estimates: object | None = None,
        peptide_differential_table: object | None = None,
        evidence: object | None = None,
        config: object | None = None,
        contrast_name: str = "aggregated",
    ) -> object:
        raise PhosPyInputError(WITHDRAWN_POSTHOC_PEPTIDE_TO_SITE_AGGREGATION_MESSAGE)


__all__ = [
    "EXPERIMENTAL_INTERNAL_API",
    "EXPERIMENTAL_INTERNAL_REASON",
    "PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS",
    "PeptideToSiteAggregator",
    "WITHDRAWN_POSTHOC_PEPTIDE_TO_SITE_AGGREGATION_MESSAGE",
]
