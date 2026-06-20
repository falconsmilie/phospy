"""Unsafe AnalysisReadyPhosphoDataset corruption helpers for tests only.

These helpers deliberately bypass the AnalysisReadyPhosphoDataset constructor.
They exist only for validator tests that must prove malformed post-boundary
objects are rejected. Production dataset validation remains owned by
AnalysisReadyPhosphoDataset; do not use these helpers to model ordinary dataset
construction.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

from phospy import AnalysisReadyPhosphoDataset
from phospy.science.transformations.models import IntensityScaleState


def unsafe_corrupt_dataset_to_display_index(
    dataset: AnalysisReadyPhosphoDataset,
) -> None:
    site_metadata = dataset._borrow_site_metadata_frame()
    display_index = pd.Index(
        site_metadata.loc[:, "display_id"].astype(str).tolist(),
        name="site_key",
    )
    object.__setattr__(
        dataset,
        "_phospho",
        dataset._borrow_phospho_frame().set_axis(display_index.copy(), axis="index"),
    )
    object.__setattr__(
        dataset,
        "_site_metadata",
        site_metadata.set_axis(display_index.copy(), axis="index"),
    )


def unsafe_replace_dataset_site_metadata(
    dataset: AnalysisReadyPhosphoDataset,
    site_metadata: pd.DataFrame,
) -> None:
    object.__setattr__(dataset, "_site_metadata", site_metadata.copy(deep=True))


def unsafe_drop_dataset_site_metadata_columns(
    dataset: AnalysisReadyPhosphoDataset,
    columns: str | Sequence[str],
) -> None:
    column_names = [columns] if isinstance(columns, str) else list(columns)
    unsafe_replace_dataset_site_metadata(
        dataset,
        dataset._borrow_site_metadata_frame().drop(columns=column_names),
    )


def unsafe_set_dataset_site_metadata_columns(
    dataset: AnalysisReadyPhosphoDataset,
    columns: Mapping[str, object],
) -> None:
    site_metadata = dataset._borrow_site_metadata_frame()
    for column_name, values in columns.items():
        site_metadata.loc[:, column_name] = values
    unsafe_replace_dataset_site_metadata(dataset, site_metadata)


def unsafe_set_dataset_site_metadata_index(
    dataset: AnalysisReadyPhosphoDataset,
    index: pd.Index,
    *,
    update_site_key_column: bool = False,
) -> None:
    site_metadata = dataset._borrow_site_metadata_frame().set_axis(
        index.copy(), axis="index"
    )
    if update_site_key_column:
        site_metadata.loc[:, "site_key"] = index.astype(str).tolist()
    unsafe_replace_dataset_site_metadata(dataset, site_metadata)


def unsafe_reverse_dataset_site_metadata_index(
    dataset: AnalysisReadyPhosphoDataset,
) -> None:
    site_metadata = dataset._borrow_site_metadata_frame()
    reversed_index = pd.Index(
        list(reversed(site_metadata.index.astype(str).tolist())),
        name=site_metadata.index.name,
    )
    unsafe_set_dataset_site_metadata_index(dataset, reversed_index)


def unsafe_replace_dataset_intensity_scale_state(
    dataset: AnalysisReadyPhosphoDataset,
    intensity_scale_state: IntensityScaleState,
) -> None:
    object.__setattr__(dataset, "intensity_scale_state", intensity_scale_state)


__all__ = [
    "unsafe_corrupt_dataset_to_display_index",
    "unsafe_drop_dataset_site_metadata_columns",
    "unsafe_replace_dataset_intensity_scale_state",
    "unsafe_replace_dataset_site_metadata",
    "unsafe_reverse_dataset_site_metadata_index",
    "unsafe_set_dataset_site_metadata_columns",
    "unsafe_set_dataset_site_metadata_index",
]
