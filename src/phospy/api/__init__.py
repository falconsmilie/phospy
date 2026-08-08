# pyright: reportUnsupportedDunderAll=false
# ruff: noqa: F401
"""Stable public facade for the supported PhosPy API.

Stable public API:
``phospy.api`` is intentionally smaller than the implementation package tree.
It owns the stable user-facing contract: builders, workflow requests, workflow
classes, primary result objects, reference bundle entrypoints, enums needed by
default examples, and common exception types.

Advanced supported API:
Advanced supported objects live under ``phospy.advanced``. Historical advanced
imports from this module remain temporarily available through a compatibility
adapter that emits ``PhosPyDeprecationWarning`` with the replacement import.

Internal / experimental API:
Internal and experimental objects remain unsupported and are not exported from
this facade.
"""

from __future__ import annotations

from phospy._api_inventory import (
    ADVANCED_PUBLIC_API as _ADVANCED_SUPPORTED_API,
)
from phospy._api_inventory import (
    INTERNAL_EXPERIMENTAL_API as _INTERNAL_EXPERIMENTAL_API,
)
from phospy._api_inventory import (
    STABLE_PUBLIC_API as _STABLE_PUBLIC_API,
)
from phospy.api._compat import deprecated_advanced_export
from phospy.api.builders import AnalysisReadyDatasetBuilder, ReferenceBundleBuilder
from phospy.api.configs import (
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
    EnrichmentConfig,
)
from phospy.api.datasets import AnalysisReadyPhosphoDataset
from phospy.api.enums import Organism, ReferencePreset
from phospy.api.requests import (
    BatchCovariate,
    CategoricalCovariate,
    ContinuousCovariate,
    Contrast,
    ContrastMatrix,
    DatasetBuildRequest,
    DesignMatrix,
    DifferentialAnalysisRequest,
    EnrichmentIdentifierSetProvenance,
    EnrichmentIdentifierSetSourceType,
    EnrichmentSet,
    EnrichmentSetCollection,
    EnrichmentWorkflowRequest,
    ExperimentalDesign,
    FixedEffectCovariate,
    GeneSetCollection,
    KinaseWorkflowRequest,
    PhosphositeImporter,
    PhosphositeImportRequest,
    PtmSetCollection,
    SampleDesignRecord,
    SignalomeWorkflowRequest,
    all_pairwise_contrasts,
    contrasts_vs_control,
)
from phospy.api.results import (
    DifferentialAnalysisResult,
    EnrichmentResultRecord,
    EnrichmentWorkflowResult,
    KinaseActivityResult,
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
    PhosphositeImportResult,
    ResultCaveat,
    SignalomeWorkflowResult,
)
from phospy.api.workflows import (
    DifferentialAnalysisWorkflow,
    EnrichmentWorkflow,
    KinaseWorkflow,
    SignalomeWorkflow,
)
from phospy.errors import (
    ContractValidationError,
    PhosPyError,
    PhosPyInputError,
    PhosPyReferenceError,
    PhosPyValidationError,
    PhosPyWorkflowError,
    ReferenceCompatibilityError,
    ReferenceResolutionError,
    SignalomeScaleError,
    UnsupportedInputFormatError,
    WorkflowBoundaryError,
    WorkflowValidationError,
)
from phospy.science.references.models import (
    ReferenceBundle,
    ReferenceBundleBuildRequest,
)
from phospy.science.transformations.models import IntensityScaleKind

__all__ = _STABLE_PUBLIC_API


def __getattr__(name: str) -> object:
    if name in _ADVANCED_SUPPORTED_API:
        return deprecated_advanced_export(name, old_module=__name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
