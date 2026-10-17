"""Narrow inputs consumed by kinase activity-stage coordination."""

from __future__ import annotations

from typing import Protocol

import pandas as pd

from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.kinase.activity.config import ResolvedKinaseActivityExecutionConfig
from phospy.kinase.references.models import KinaseLibraryResource
from phospy.references.models import ReferenceBundle


class ActivityStageSiteUniverses(Protocol):
    """Site-universe view required by activity execution."""

    @property
    def predicted_membership_sites(self) -> pd.Index: ...


class ActivityStageRequest(Protocol):
    """Resolved workflow facts required by the activity stage."""

    @property
    def dataset(self) -> AnalysisReadyPhosphoDataset: ...

    @property
    def activity_phospho_matrix(self) -> pd.DataFrame: ...

    @property
    def ksea_background_phospho_matrix(self) -> pd.DataFrame: ...

    @property
    def ssgsea_effect_matrix(self) -> pd.DataFrame: ...

    @property
    def scoring_phospho_matrix(self) -> pd.DataFrame: ...

    @property
    def reference_membership_map(self) -> pd.DataFrame: ...

    @property
    def site_identity_map(self) -> pd.DataFrame | None: ...

    @property
    def site_universes(self) -> ActivityStageSiteUniverses | None: ...

    @property
    def references(self) -> ReferenceBundle: ...

    @property
    def kinase_library_resource(self) -> KinaseLibraryResource | None: ...


class ActivityStageExecutionConfig(Protocol):
    """Resolved workflow options required by activity coordination."""

    @property
    def activity(self) -> ResolvedKinaseActivityExecutionConfig | None: ...

    @property
    def prediction_mode(self) -> object: ...

    @property
    def prediction_top_k(self) -> int: ...

    @property
    def prediction_deterministic_max_selected_kinases(self) -> int: ...


class ActivityDownstreamScoreSource(Protocol):
    """Closed scoring-source value consumed by membership coordination."""

    @property
    def value(self) -> str: ...


__all__ = [
    "ActivityDownstreamScoreSource",
    "ActivityStageExecutionConfig",
    "ActivityStageRequest",
    "ActivityStageSiteUniverses",
]
