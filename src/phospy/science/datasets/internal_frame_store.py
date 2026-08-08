"""Dataset-owned immutable frame snapshots for internal workflow reads."""

from __future__ import annotations

from threading import Lock

import pandas as pd

from phospy.frames.ownership import (
    ImmutableDataFrameSnapshot,
    immutable_dataframe_snapshot,
    immutable_optional_dataframe_snapshot,
)


class DatasetInternalFrameStore:
    """Lazy immutable snapshot store owned by one analysis-ready dataset.

    The store keeps references to the dataset-owned frames but never returns
    them. First internal access to each frame builds one owner-detached immutable
    snapshot; later workflow views receive workflow-local wrappers over that
    stored snapshot. Shareable NumPy-backed columns use immutable buffers, while
    unshareable columns are copied per wrapper.
    """

    __slots__ = (
        "_comparisons",
        "_comparisons_snapshot",
        "_lock",
        "_phospho",
        "_phospho_snapshot",
        "_sample_metadata",
        "_sample_metadata_snapshot",
        "_site_metadata",
        "_site_metadata_snapshot",
        "_total",
        "_total_snapshot",
    )

    def __init__(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        total: pd.DataFrame | None,
        comparisons: pd.DataFrame | None,
    ) -> None:
        self._phospho = phospho
        self._site_metadata = site_metadata
        self._sample_metadata = sample_metadata
        self._total = total
        self._comparisons = comparisons
        self._phospho_snapshot: ImmutableDataFrameSnapshot | None = None
        self._site_metadata_snapshot: ImmutableDataFrameSnapshot | None = None
        self._sample_metadata_snapshot: ImmutableDataFrameSnapshot | None = None
        self._total_snapshot: ImmutableDataFrameSnapshot | None = None
        self._comparisons_snapshot: ImmutableDataFrameSnapshot | None = None
        self._lock = Lock()

    @classmethod
    def from_frames(
        cls,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        total: pd.DataFrame | None,
        comparisons: pd.DataFrame | None,
    ) -> DatasetInternalFrameStore:
        return cls(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            total=total,
            comparisons=comparisons,
        )

    def phospho_frame(self) -> pd.DataFrame:
        return self._required_snapshot(
            snapshot_attr="_phospho_snapshot",
            owner=self._phospho,
            field_name="dataset.phospho internal snapshot",
        ).dataframe()

    def site_metadata_frame(self) -> pd.DataFrame:
        return self._required_snapshot(
            snapshot_attr="_site_metadata_snapshot",
            owner=self._site_metadata,
            field_name="dataset.site_metadata internal snapshot",
        ).dataframe()

    def sample_metadata_frame(self) -> pd.DataFrame | None:
        snapshot = self._optional_snapshot(
            snapshot_attr="_sample_metadata_snapshot",
            owner=self._sample_metadata,
            field_name="dataset.sample_metadata internal snapshot",
        )
        if snapshot is None:
            return None
        return snapshot.dataframe()

    def total_frame(self) -> pd.DataFrame | None:
        snapshot = self._optional_snapshot(
            snapshot_attr="_total_snapshot",
            owner=self._total,
            field_name="dataset.total internal snapshot",
        )
        if snapshot is None:
            return None
        return snapshot.dataframe()

    def comparisons_frame(self) -> pd.DataFrame | None:
        snapshot = self._optional_snapshot(
            snapshot_attr="_comparisons_snapshot",
            owner=self._comparisons,
            field_name="dataset.comparisons internal snapshot",
        )
        if snapshot is None:
            return None
        return snapshot.dataframe()

    def is_current_for(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        total: pd.DataFrame | None,
        comparisons: pd.DataFrame | None,
    ) -> bool:
        """Return whether this store still represents the dataset-owned frames."""

        return (
            self._phospho is phospho
            and self._site_metadata is site_metadata
            and self._sample_metadata is sample_metadata
            and self._total is total
            and self._comparisons is comparisons
        )

    def _required_snapshot(
        self,
        *,
        snapshot_attr: str,
        owner: pd.DataFrame,
        field_name: str,
    ) -> ImmutableDataFrameSnapshot:
        snapshot = getattr(self, snapshot_attr)
        if snapshot is not None:
            return snapshot
        with self._lock:
            snapshot = getattr(self, snapshot_attr)
            if snapshot is None:
                snapshot = immutable_dataframe_snapshot(
                    owner,
                    field_name=field_name,
                )
                setattr(self, snapshot_attr, snapshot)
            return snapshot

    def _optional_snapshot(
        self,
        *,
        snapshot_attr: str,
        owner: pd.DataFrame | None,
        field_name: str,
    ) -> ImmutableDataFrameSnapshot | None:
        if owner is None:
            return None
        snapshot = getattr(self, snapshot_attr)
        if snapshot is not None:
            return snapshot
        with self._lock:
            snapshot = getattr(self, snapshot_attr)
            if snapshot is None:
                snapshot = immutable_optional_dataframe_snapshot(
                    owner,
                    field_name=field_name,
                )
                setattr(self, snapshot_attr, snapshot)
            return snapshot


__all__: list[str] = []
