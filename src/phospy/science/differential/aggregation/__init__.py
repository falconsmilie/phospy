"""Withdrawn peptide-to-site differential post-hoc aggregation boundary.

The preferred supported peptide-to-site route resolves peptide evidence at
sample-intensity level during dataset building, then runs the core differential
workflow. The former post-hoc peptide differential estimate-combination lane is
not part of the supported public facade.
"""

from phospy.science.differential.aggregation.models import (
    PEPTIDE_TO_SITE_AGGREGATION_LEVEL_SAMPLE_INTENSITY,
    PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS,
)

__all__ = [
    "PEPTIDE_TO_SITE_AGGREGATION_LEVEL_SAMPLE_INTENSITY",
    "PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS",
]
