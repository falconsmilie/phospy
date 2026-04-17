"""Kinase activity analysis domain.

This package owns downstream activity scoring built on prediction outputs and
compatible phospho-derived inputs. It should remain focused on activity
analysis behaviour and result models, not public workflow orchestration.
"""

from .analysis import KinaseActivityAnalyzer
from .results import KinaseActivityResult
from .scoring import (
    build_kinase_target_table,
    compute_ksea_scores,
    compute_weighted_kinase_activity,
    count_predicted_targets,
)

__all__ = [
    "KinaseActivityAnalyzer",
    "KinaseActivityResult",
    "build_kinase_target_table",
    "compute_ksea_scores",
    "compute_weighted_kinase_activity",
    "count_predicted_targets",
]
