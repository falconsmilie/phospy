from .activities import (
    build_kinase_target_table,
    compute_ksea_scores,
    compute_weighted_kinase_activity,
    count_predicted_targets,
)
from .analysis import KinaseActivityAnalyzer, KinaseActivityResult
from .dataset import CoreProcessingResult, PhosphoDataset, SiteMatrixResult
from .matrices import build_site_matrix
from .pipeline import CoreOutputs, PhosRPipeline, run_core_pipeline
from .prediction import (
    KinasePredictionResult,
    KinasePredictor,
    build_candidate_substrate_list,
)
from .preprocessing import (
    add_pairwise_comparisons,
    collapse_duplicate_genes,
    correct_phospho_to_protein,
    filter_min_observed,
    replace_sentinel_with_nan,
)
from .profiles import (
    KinaseProfileBuilder,
    KinaseProfileResult,
    build_kinase_substrate_profiles,
)
from .scoring import (
    KinaseScorer,
    KinaseScoringResult,
    KinaseSubstrateScoreResult,
    combine_profile_and_motif_scores,
    kinase_substrate_score,
)

__all__ = [
    "add_pairwise_comparisons",
    "build_kinase_target_table",
    "build_site_matrix",
    "collapse_duplicate_genes",
    "compute_ksea_scores",
    "compute_weighted_kinase_activity",
    "correct_phospho_to_protein",
    "count_predicted_targets",
    "CoreOutputs",
    "CoreProcessingResult",
    "filter_min_observed",
    "KinaseActivityAnalyzer",
    "KinasePredictionResult",
    "KinaseProfileBuilder",
    "KinaseProfileResult",
    "KinasePredictor",
    "KinaseActivityResult",
    "KinaseScorer",
    "KinaseScoringResult",
    "KinaseSubstrateScoreResult",
    "PhosphoDataset",
    "build_candidate_substrate_list",
    "build_kinase_substrate_profiles",
    "combine_profile_and_motif_scores",
    "kinase_substrate_score",
    "PhosRPipeline",
    "replace_sentinel_with_nan",
    "run_core_pipeline",
    "SiteMatrixResult",
]
