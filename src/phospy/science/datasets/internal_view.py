"""Dataset-owned internal frame view for trusted workflow collaborators."""

from __future__ import annotations

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

    @property
    def imputation_observed_mask(self) -> pd.DataFrame | None:
        return self._dataset._borrow_imputation_observed_mask_frame()


__all__ = ["DatasetInternalView"]
