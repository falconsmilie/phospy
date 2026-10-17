"""Science-layer batch-correction numerical executors."""

from phospy.preprocessing.batch_correction.execution.ruv_iii import (
    RUV_III_ALGORITHM_ID,
    RUV_III_METHOD,
    RuvIIIDiagnostics,
    RuvIIIKernel,
    RuvIIIReplicateStructure,
    RuvIIIResult,
    run_ruv_iii,
)
from phospy.preprocessing.batch_correction.execution.ruv_iii_executor import (
    RUV_III_STYLE_ALGORITHM_DESCRIPTION,
    RUV_III_STYLE_EXECUTOR_ID,
    RUV_III_STYLE_METHOD,
    RuvIIIStyleExecutor,
    RuvIIIStyleExecutorDiagnostics,
    RuvIIIStyleExecutorResult,
)
from phospy.preprocessing.batch_correction.execution.sps_ruv import (
    SPS_RUV_STYLE_ALGORITHM_DESCRIPTION,
    SPS_RUV_STYLE_BATCH_TERM_ROLE,
    SPS_RUV_STYLE_EXECUTOR_ID,
    SPS_RUV_STYLE_PROTECTED_TERM_ROLE,
    SPS_RUV_STYLE_REPLICATE_METADATA_ROLE,
    BatchCorrectionExecutor,
    DeterministicSpsRuvStyleExecutor,
    SpsRuvStyleExecutor,
    SpsRuvStyleExecutorDiagnostics,
    SpsRuvStyleExecutorResult,
)

__all__ = [
    "BatchCorrectionExecutor",
    "SPS_RUV_STYLE_ALGORITHM_DESCRIPTION",
    "SPS_RUV_STYLE_BATCH_TERM_ROLE",
    "SPS_RUV_STYLE_EXECUTOR_ID",
    "SPS_RUV_STYLE_PROTECTED_TERM_ROLE",
    "SPS_RUV_STYLE_REPLICATE_METADATA_ROLE",
    "RUV_III_ALGORITHM_ID",
    "RUV_III_METHOD",
    "RUV_III_STYLE_ALGORITHM_DESCRIPTION",
    "RUV_III_STYLE_EXECUTOR_ID",
    "RUV_III_STYLE_METHOD",
    "DeterministicSpsRuvStyleExecutor",
    "SpsRuvStyleExecutor",
    "SpsRuvStyleExecutorDiagnostics",
    "SpsRuvStyleExecutorResult",
    "RuvIIIDiagnostics",
    "RuvIIIKernel",
    "RuvIIIStyleExecutor",
    "RuvIIIStyleExecutorDiagnostics",
    "RuvIIIStyleExecutorResult",
    "RuvIIIReplicateStructure",
    "RuvIIIResult",
    "run_ruv_iii",
]
