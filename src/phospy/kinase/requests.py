"""Public Kinase workflow request contract."""

from __future__ import annotations

from dataclasses import dataclass, field

from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.kinase.activity.config import KinaseActivityConfig
from phospy.kinase.prediction.config import KinasePredictionConfig
from phospy.kinase.references.models import KinaseLibraryResource
from phospy.kinase.scoring.config import (
    KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE,
    KinaseReferenceDisplayAmbiguityPolicy,
    KinaseScoringConfig,
    KinaseSiteSequenceConflictPolicy,
)
from phospy.references.models import ReferenceBundle, ReferencePreset


@dataclass(frozen=True, slots=True, eq=False)
class KinaseWorkflowRequest:
    """Request for the public kinase workflow.

    Construction stores the payload only. Dataset compatibility, reference
    compatibility, scoring/prediction/activity config coherence, localisation
    requirements, and reference-projection policies are validated by
    ``KinaseWorkflow.run(...)`` before interpretation or execution.

    `site_sequence_conflict_policy` controls how dataset and reference
    site-sequence disagreements are handled during interpretation.

    `reference_display_ambiguity_policy` controls whether display-level
    kinase-substrate reference rows may project to more than one dataset
    `site_key`. The default is conservative and rejects ambiguous projection.

    Python equality and hashing are identity-based. The workflow owns semantic
    request interpretation; this payload does not define content equality.
    """

    __hash__ = object.__hash__

    dataset: AnalysisReadyPhosphoDataset
    references: ReferencePreset | ReferenceBundle = ReferencePreset.AUTO
    scoring_config: KinaseScoringConfig | None = None
    prediction_config: KinasePredictionConfig = field(
        default_factory=KinasePredictionConfig
    )
    activity_config: KinaseActivityConfig | None = None
    site_sequence_conflict_policy: KinaseSiteSequenceConflictPolicy = (
        KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE
    )
    reference_display_ambiguity_policy: KinaseReferenceDisplayAmbiguityPolicy = (
        KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR
    )
    kinase_library_resource: KinaseLibraryResource | None = None


__all__ = ["KinaseWorkflowRequest"]
