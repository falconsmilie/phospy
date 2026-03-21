from .activities import (
    build_kinase_target_table,
    compute_ksea_scores,
    compute_weighted_kinase_activity,
    count_predicted_targets,
)
from .analysis import KinaseActivityAnalyzer, KinaseActivityResult
from .dataset import CoreProcessingResult, PhosphoDataset, SiteMatrixResult
from .matrices import build_site_matrix
from .motifs import (
    KinaseMotifScorer,
    MotifScoringResult,
    create_frequency_matrix,
    frequency_scoring,
    score_phosphosite_motifs,
)
from .pipeline import CoreOutputs, PhosRPipeline, run_core_pipeline
from .preprocessing import (
    add_pairwise_comparisons,
    collapse_duplicate_genes,
    correct_phospho_to_protein,
    filter_min_observed,
    replace_sentinel_with_nan,
)
from .scoring import KinaseScorer, KinaseScoringResult, combine_profile_and_motif_scores

__all__ = [
    "add_pairwise_comparisons",
    "build_kinase_target_table",
    "build_site_matrix",
    "collapse_duplicate_genes",
    "compute_ksea_scores",
    "compute_weighted_kinase_activity",
    "create_frequency_matrix",
    "correct_phospho_to_protein",
    "count_predicted_targets",
    "CoreOutputs",
    "CoreProcessingResult",
    "filter_min_observed",
    "frequency_scoring",
    "KinaseActivityAnalyzer",
    "KinaseActivityResult",
    "KinaseMotifScorer",
    "KinaseScorer",
    "KinaseScoringResult",
    "MotifScoringResult",
    "PhosphoDataset",
    "combine_profile_and_motif_scores",
    "PhosRPipeline",
    "replace_sentinel_with_nan",
    "run_core_pipeline",
    "score_phosphosite_motifs",
    "SiteMatrixResult",
]
