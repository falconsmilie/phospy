"""Phosphoproteomic preprocessing domain.

This package owns phosphoproteomic preprocessing behaviour, including stepwise
filters, protein correction, site-matrix construction, dataset-bound core
processing, and the stable preprocessing building blocks used by public
workflow entry points.
"""

from .core import (
    CorePreprocessingConfig,
    CoreProcessingResult,
    CoreProcessor,
    resolve_core_preprocessing_config,
)
from .dataset import DatasetPreprocessing
from .protein_correction import ProteinCorrectionResult, ProteinCorrectionSummary
from .services import PhosphoPreprocessor, ProteinCorrectionService, TotalPreprocessor
from .site_matrix import SiteMatrixBuilder, SiteMatrixResult
from .steps import (
    CoverageFilterResult,
    CoverageFilterSummary,
    LocalizationFilterResult,
    LocalizationFilterSummary,
    add_pairwise_comparisons,
    collapse_duplicate_genes,
    correct_phospho_to_protein,
    filter_localized_sites,
    filter_min_observed,
    filter_sites_by_coverage,
    replace_sentinel_with_nan,
)

__all__ = [
    "CorePreprocessingConfig",
    "CoreProcessingResult",
    "CoreProcessor",
    "CoverageFilterResult",
    "CoverageFilterSummary",
    "DatasetPreprocessing",
    "LocalizationFilterResult",
    "LocalizationFilterSummary",
    "PhosphoPreprocessor",
    "ProteinCorrectionResult",
    "ProteinCorrectionService",
    "ProteinCorrectionSummary",
    "SiteMatrixBuilder",
    "SiteMatrixResult",
    "TotalPreprocessor",
    "add_pairwise_comparisons",
    "collapse_duplicate_genes",
    "correct_phospho_to_protein",
    "filter_localized_sites",
    "filter_min_observed",
    "filter_sites_by_coverage",
    "replace_sentinel_with_nan",
    "resolve_core_preprocessing_config",
]
