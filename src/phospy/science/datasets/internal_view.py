"""Dataset-owned internal frame view for trusted workflow collaborators."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from phospy.science.datasets.models import AnalysisReadyPhosphoDataset


class DatasetInternalView:
    """Narrow borrowed-frame contract for internal workflow access."""

    __slots__ = ("_dataset",)

    def __init__(self, dataset: AnalysisReadyPhosphoDataset) -> None:
        self._dataset = dataset

    @property
    def phospho(self) -> pd.DataFrame:
        return self._dataset._borrow_phospho_frame()

    @property
    def site_metadata(self) -> pd.DataFrame:
        return self._dataset._borrow_site_metadata_frame()

    @property
    def sample_metadata(self) -> pd.DataFrame | None:
        return self._dataset._borrow_sample_metadata_frame()

    @property
    def total(self) -> pd.DataFrame | None:
        return self._dataset._borrow_total_frame()

    @property
    def comparisons(self) -> pd.DataFrame | None:
        return self._dataset._borrow_comparisons_frame()

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


__all__ = ["DatasetInternalView"]
