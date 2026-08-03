# pyright: reportUnsupportedDunderAll=false
# ruff: noqa: F401
"""Advanced supported result inspection models.

Primary workflow result containers remain stable under ``phospy.api.results``.
The names exported here are narrower result diagnostics or compatibility type
aliases that callers should import deliberately.
"""

from __future__ import annotations

from phospy._api_inventory import ADVANCED_RESULT_API
from phospy.contracts.results import DifferentialModelDiagnostics
from phospy.contracts.results.kinase import (
    KinaseEligibilityReport,
    KinaseWorkflowAttritionProvenance,
    KinaseWorkflowCaveat,
    KinaseWorkflowPreprocessingAttritionSummary,
    KinaseWorkflowScoringAttritionSummary,
    KinaseWorkflowSiteAttritionSummary,
)

__all__ = ADVANCED_RESULT_API
