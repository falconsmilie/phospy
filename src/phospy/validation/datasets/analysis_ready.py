"""Adapter for analysis-ready dataset model-boundary validation."""

from __future__ import annotations

import pandas as pd

from phospy.provenance.models import RunProvenance, TrustedDatasetConstructionAssertions
from phospy.science.datasets.models import (
    AnalysisReadyPhosphoDataset,
    DatasetPreprocessingReport,
)
from phospy.science.datasets.processing_state import DatasetProcessingState
from phospy.science.references.models import Organism
from phospy.science.transformations.models import IntensityScaleState


class AnalysisReadyDatasetModelBoundaryValidator:
    """Validate by delegating to `AnalysisReadyPhosphoDataset` construction.

    This adapter does not own analysis-ready invariants. The model constructor
    remains the single authoritative owner for strict dataset-boundary checks,
    and trusted callers must supply typed construction assertions.
    """

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        sample_metadata: pd.DataFrame | None,
        total: pd.DataFrame | None,
        comparisons: pd.DataFrame | None,
        organism: Organism | None,
        imputation_observation_mask: pd.DataFrame | None = None,
        preprocessing_report: DatasetPreprocessingReport | None = None,
        provenance: RunProvenance | None = None,
        trusted_construction_assertions: TrustedDatasetConstructionAssertions
        | None = None,
    ) -> AnalysisReadyPhosphoDataset:
        return AnalysisReadyPhosphoDataset.from_trusted_tables(
            phospho=phospho,
            site_metadata=site_metadata,
            intensity_scale_state=intensity_scale_state,
            processing_state=processing_state,
            sample_metadata=sample_metadata,
            total=total,
            comparisons=comparisons,
            imputation_observation_mask=imputation_observation_mask,
            organism=organism,
            preprocessing_report=preprocessing_report,
            provenance=provenance,
            trusted_construction_assertions=trusted_construction_assertions,
        )
