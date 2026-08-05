"""Compatibility import route for transformation-domain models.

Owned implementations are split by responsibility under this package. This
module preserves the historical ``phospy.science.transformations.models`` import
path by re-exporting the same class and constant objects.
"""

from phospy.science.transformations.policy import (
    CALLER_DECLARABLE_QUANTITATIVE_MEANINGS,
    IDENTITY_INTENSITY_SCALE_ESTABLISHER,
    QUANTITATIVE_MEANING_LEGACY_UNVERIFIED_CAVEAT_CODE,
    QUANTITATIVE_MEANING_OPERATION_CALLER_DECLARATION,
    QUANTITATIVE_MEANING_OPERATION_LEGACY_BUNDLE_MIGRATION,
    QUANTITATIVE_MEANING_OPERATION_SCALE_CONTRACT_INFERENCE,
    QUANTITATIVE_MEANING_OPERATION_TOTAL_PROTEIN_SUBTRACT_LOG_TOTAL,
    QUANTITATIVE_MEANING_PROVENANCE_SCHEMA_VERSION_V1,
    QUANTITATIVE_MEANING_USER_DECLARED_CAVEAT_CODE,
    DeclaredIntensityScaleDiagnosticPolicy,
    IntensityScaleEstablishmentMode,
    IntensityScaleEstablishmentSource,
    IntensityScaleEvidenceLevel,
    IntensityScaleKind,
    QuantitativeMeaning,
    QuantitativeMeaningEvidenceMode,
    QuantitativeMeaningScaleRule,
    caller_declarable_quantitative_meaning_values,
    default_quantitative_meaning_for_scale_kind,
    is_caller_declarable_quantitative_meaning,
)
from phospy.science.transformations.provenance import (
    IntensityScaleEstablishmentProvenance,
    QuantitativeMeaningTransitionProvenance,
)
from phospy.science.transformations.scale_state import (
    IntensityScaleState,
    establish_intensity_scale_state,
)
from phospy.science.transformations.scale_values import (
    IntensityTransformationEvent,
    MatrixIntensityScaleState,
)

__all__ = [
    "CALLER_DECLARABLE_QUANTITATIVE_MEANINGS",
    "DeclaredIntensityScaleDiagnosticPolicy",
    "IDENTITY_INTENSITY_SCALE_ESTABLISHER",
    "IntensityScaleEstablishmentMode",
    "IntensityScaleEstablishmentProvenance",
    "IntensityScaleEstablishmentSource",
    "IntensityScaleEvidenceLevel",
    "IntensityScaleKind",
    "IntensityScaleState",
    "IntensityTransformationEvent",
    "MatrixIntensityScaleState",
    "QUANTITATIVE_MEANING_LEGACY_UNVERIFIED_CAVEAT_CODE",
    "QUANTITATIVE_MEANING_OPERATION_CALLER_DECLARATION",
    "QUANTITATIVE_MEANING_OPERATION_LEGACY_BUNDLE_MIGRATION",
    "QUANTITATIVE_MEANING_OPERATION_SCALE_CONTRACT_INFERENCE",
    "QUANTITATIVE_MEANING_OPERATION_TOTAL_PROTEIN_SUBTRACT_LOG_TOTAL",
    "QUANTITATIVE_MEANING_PROVENANCE_SCHEMA_VERSION_V1",
    "QUANTITATIVE_MEANING_USER_DECLARED_CAVEAT_CODE",
    "QuantitativeMeaning",
    "QuantitativeMeaningEvidenceMode",
    "QuantitativeMeaningScaleRule",
    "QuantitativeMeaningTransitionProvenance",
    "caller_declarable_quantitative_meaning_values",
    "default_quantitative_meaning_for_scale_kind",
    "establish_intensity_scale_state",
    "is_caller_declarable_quantitative_meaning",
]
