# pyright: reportUnsupportedDunderAll=false
"""Compatibility route for dataset-resolution contracts."""

from __future__ import annotations

from phospy.science.evidence.dataset_resolution.models import *  # noqa: F403
from phospy.science.evidence.dataset_resolution.models import __all__ as _MODEL_ALL
from phospy.science.evidence.dataset_resolution.policies import (
    build_multi_site_handling_config_for_dataset_policy,
    build_peptide_to_site_aggregation_policy,
    validate_dataset_multi_site_policy,
)

__all__ = [
    *_MODEL_ALL,
    "build_multi_site_handling_config_for_dataset_policy",
    "build_peptide_to_site_aggregation_policy",
    "validate_dataset_multi_site_policy",
]
