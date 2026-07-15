"""Derived analysis-ready quantitative dataset models."""

from __future__ import annotations

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.provenance.derived_quantitative import (
    DerivedQuantitativeDataProvenance,
)
from phospy.provenance.models import RunProvenance
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.datasets.processing_state import DatasetProcessingState
from phospy.science.references.models import Organism
from phospy.science.transformations.models import IntensityScaleState


class DerivedAnalysisReadyPhosphoDataset(AnalysisReadyPhosphoDataset):
    """Analysis-ready dataset whose quantitative matrices were derived from a parent."""

    __slots__ = ("derived_lineage",)

    def __init__(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        derived_lineage: DerivedQuantitativeDataProvenance,
        provenance: RunProvenance,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        imputation_observation_mask: pd.DataFrame | None = None,
        organism: Organism | None = None,
        allow_opaque_site_values: bool = False,
        _assume_owned: bool = False,
    ) -> None:
        if not isinstance(derived_lineage, DerivedQuantitativeDataProvenance):
            raise PhosPyInputError(
                "derived dataset requires DerivedQuantitativeDataProvenance"
            )
        if not isinstance(provenance, RunProvenance):
            raise PhosPyInputError("derived dataset requires RunProvenance")
        super().__init__(
            phospho=phospho,
            site_metadata=site_metadata,
            intensity_scale_state=intensity_scale_state,
            processing_state=processing_state,
            sample_metadata=sample_metadata,
            total=total,
            comparisons=comparisons,
            imputation_observation_mask=imputation_observation_mask,
            organism=organism,
            preprocessing_report=None,
            protein_aware_preparation=None,
            provenance=provenance,
            allow_opaque_site_values=allow_opaque_site_values,
            _assume_owned=_assume_owned,
        )
        object.__setattr__(self, "derived_lineage", derived_lineage)

    @classmethod
    def from_owned_derived_tables(
        cls,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        derived_lineage: DerivedQuantitativeDataProvenance,
        provenance: RunProvenance,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        imputation_observation_mask: pd.DataFrame | None = None,
        organism: Organism | None = None,
        allow_opaque_site_values: bool = False,
    ) -> DerivedAnalysisReadyPhosphoDataset:
        """Construct from already-owned derived tables and typed lineage."""

        return cls(
            phospho=phospho,
            site_metadata=site_metadata,
            intensity_scale_state=intensity_scale_state,
            processing_state=processing_state,
            sample_metadata=sample_metadata,
            total=total,
            comparisons=comparisons,
            imputation_observation_mask=imputation_observation_mask,
            organism=organism,
            provenance=provenance,
            derived_lineage=derived_lineage,
            allow_opaque_site_values=allow_opaque_site_values,
            _assume_owned=True,
        )


__all__ = ["DerivedAnalysisReadyPhosphoDataset"]
