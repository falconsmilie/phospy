"""Dataset-owned defensive internal view for workflow collaborators."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from phospy.science.datasets.internal_frame_store import DatasetInternalFrameStore
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.datasets.preprocessing.protein_aware_models import (
    ProteinAwarePreparationReport,
    ProteinAwarePreparationResult,
)


class ProteinAwarePreparationInternalView:
    """Narrow mutation-isolated access to a dataset-owned preparation sidecar."""

    __slots__ = ("_preparation",)

    def __init__(self, preparation: ProteinAwarePreparationResult) -> None:
        self._preparation = preparation

    @property
    def report(self) -> ProteinAwarePreparationReport:
        return self._preparation.report

    @property
    def matched_pairs(self) -> pd.DataFrame:
        return self._preparation.matched_pairs_dataframe()

    @property
    def protein_covariate_matrix(self) -> pd.DataFrame:
        return self._preparation.protein_covariate_matrix_dataframe()

    @property
    def site_eligibility(self) -> pd.DataFrame:
        return self._preparation.site_eligibility_table


class DatasetInternalView:
    """Narrow dataset-domain access to dataset-owned immutable snapshots."""

    __slots__ = ("_dataset", "_frame_store")

    def __init__(self, dataset: AnalysisReadyPhosphoDataset) -> None:
        self._dataset = dataset
        self._frame_store: DatasetInternalFrameStore = (
            dataset._internal_frame_store_for_current_frames()  # pyright: ignore[reportPrivateUsage] - DatasetInternalView is the dataset-domain internal access boundary.
        )

    @property
    def phospho(self) -> pd.DataFrame:
        return self._frame_store.phospho_frame()

    @property
    def site_metadata(self) -> pd.DataFrame:
        return self._frame_store.site_metadata_frame()

    @property
    def sample_metadata(self) -> pd.DataFrame | None:
        return self._frame_store.sample_metadata_frame()

    @property
    def total(self) -> pd.DataFrame | None:
        return self._frame_store.total_frame()

    @property
    def comparisons(self) -> pd.DataFrame | None:
        return self._frame_store.comparisons_frame()

    @property
    def protein_aware_preparation(
        self,
    ) -> ProteinAwarePreparationInternalView | None:
        preparation = self._dataset.protein_aware_preparation
        if preparation is None:
            return None
        return ProteinAwarePreparationInternalView(preparation)

    def imputation_observation_summary(
        self,
        *,
        feature_ids: Sequence[object],
        sample_ids: Sequence[object],
    ) -> pd.DataFrame | None:
        return self._dataset._imputation_observation_summary_frame(
            feature_ids=feature_ids,
            sample_ids=sample_ids,
        )

    def aggregate_imputation_observation_mask(
        self,
        *,
        sample_groups: Sequence[tuple[object, Sequence[object]]],
    ) -> pd.DataFrame | None:
        return self._dataset._aggregated_imputation_observation_mask_frame(
            sample_groups=sample_groups,
        )


__all__ = ["DatasetInternalView", "ProteinAwarePreparationInternalView"]
