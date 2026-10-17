"""Published composition boundary for the core preprocessing workflow."""

from phospy.preprocessing.batch_correction.contracts.stage import (
    BatchCorrectionAdequacyValidatorProtocol,
    BatchDesignMetadataValidatorProtocol,
    SpsRuvStyleBatchCorrectionRunner,
)
from phospy.preprocessing.execution.pipeline import PreprocessingPipeline
from phospy.preprocessing.provenance.composition import PreprocessingProvenanceAdapter
from phospy.preprocessing.provenance.models import (
    DeterminismKind,
    PreprocessingStageProvenance,
)
from phospy.preprocessing.provenance.policies import (
    PreprocessingStageOrderPolicy,
    build_duplicate_site_resolution_policy,
)
from phospy.preprocessing.results.report_rows import compose_stage_owned_report_tables
from phospy.preprocessing.validation.configuration import (
    DatasetPreprocessingConfigValidator,
)
from phospy.preprocessing.validation.group_coverage import (
    GroupCoverageFilterMetadataValidator,
)
from phospy.preprocessing.validation.quantitative_evidence import (
    AdditivePreprocessingScaleGuard,
)

__all__ = [
    "AdditivePreprocessingScaleGuard",
    "BatchCorrectionAdequacyValidatorProtocol",
    "BatchDesignMetadataValidatorProtocol",
    "DatasetPreprocessingConfigValidator",
    "DeterminismKind",
    "GroupCoverageFilterMetadataValidator",
    "PreprocessingPipeline",
    "PreprocessingProvenanceAdapter",
    "PreprocessingStageOrderPolicy",
    "PreprocessingStageProvenance",
    "SpsRuvStyleBatchCorrectionRunner",
    "build_duplicate_site_resolution_policy",
    "compose_stage_owned_report_tables",
]
