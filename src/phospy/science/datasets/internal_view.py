"""Dataset-owned defensive internal view for workflow collaborators."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from phospy.frames.ownership import (
    ImmutableDataFrameSnapshot,
    immutable_dataframe_snapshot,
    immutable_optional_dataframe_snapshot,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset


class DatasetInternalView:
    """Narrow dataset-domain access to workflow-scoped immutable snapshots."""

    __slots__ = (
        "_comparisons_snapshot",
        "_dataset",
        "_phospho_snapshot",
        "_sample_metadata_snapshot",
        "_site_metadata_snapshot",
        "_total_snapshot",
    )

    def __init__(self, dataset: AnalysisReadyPhosphoDataset) -> None:
        self._dataset = dataset
        self._phospho_snapshot: ImmutableDataFrameSnapshot | None = None
        self._site_metadata_snapshot: ImmutableDataFrameSnapshot | None = None
        self._sample_metadata_snapshot: ImmutableDataFrameSnapshot | None = None
        self._total_snapshot: ImmutableDataFrameSnapshot | None = None
        self._comparisons_snapshot: ImmutableDataFrameSnapshot | None = None

    @property
    def phospho(self) -> pd.DataFrame:
        snapshot = self._phospho_snapshot
        if snapshot is None:
            snapshot = immutable_dataframe_snapshot(
                self._dataset._phospho,  # pyright: ignore[reportPrivateUsage] - DatasetInternalView is the dataset-domain internal access boundary.
                field_name="dataset.phospho internal snapshot",
            )
            self._phospho_snapshot = snapshot
        return snapshot.dataframe()

    @property
    def site_metadata(self) -> pd.DataFrame:
        snapshot = self._site_metadata_snapshot
        if snapshot is None:
            snapshot = immutable_dataframe_snapshot(
                self._dataset._site_metadata,  # pyright: ignore[reportPrivateUsage] - DatasetInternalView is the dataset-domain internal access boundary.
                field_name="dataset.site_metadata internal snapshot",
            )
            self._site_metadata_snapshot = snapshot
        return snapshot.dataframe()

    @property
    def sample_metadata(self) -> pd.DataFrame | None:
        if self._sample_metadata_snapshot is None:
            self._sample_metadata_snapshot = immutable_optional_dataframe_snapshot(
                self._dataset._sample_metadata,  # pyright: ignore[reportPrivateUsage] - DatasetInternalView is the dataset-domain internal access boundary.
                field_name="dataset.sample_metadata internal snapshot",
            )
        snapshot = self._sample_metadata_snapshot
        if snapshot is None:
            return None
        return snapshot.dataframe()

    @property
    def total(self) -> pd.DataFrame | None:
        if self._total_snapshot is None:
            self._total_snapshot = immutable_optional_dataframe_snapshot(
                self._dataset._total,  # pyright: ignore[reportPrivateUsage] - DatasetInternalView is the dataset-domain internal access boundary.
                field_name="dataset.total internal snapshot",
            )
        snapshot = self._total_snapshot
        if snapshot is None:
            return None
        return snapshot.dataframe()

    @property
    def comparisons(self) -> pd.DataFrame | None:
        if self._comparisons_snapshot is None:
            self._comparisons_snapshot = immutable_optional_dataframe_snapshot(
                self._dataset._comparisons,  # pyright: ignore[reportPrivateUsage] - DatasetInternalView is the dataset-domain internal access boundary.
                field_name="dataset.comparisons internal snapshot",
            )
        snapshot = self._comparisons_snapshot
        if snapshot is None:
            return None
        return snapshot.dataframe()

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
