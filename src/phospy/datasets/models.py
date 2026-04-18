"""Dataset domain models."""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from phospy.errors.validation import DatasetValidationError
from phospy.references.models import Organism
from phospy.site_ids import canonicalize_site_index
from phospy.transformations.models import TransformationState
from phospy.validation.datasets.analysis_ready import AnalysisReadyDatasetValidator

_DATASET_VALIDATOR = AnalysisReadyDatasetValidator()


@dataclass(frozen=True, slots=True)
class AnalysisReadyPhosphoDataset:
    """Public analysis-ready dataset contract."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame | None = None
    total: pd.DataFrame | None = None
    organism: Organism | None = None
    transformation_state: TransformationState = field(
        default_factory=TransformationState.raw
    )

    def __post_init__(self) -> None:
        phospho = _copy_frame(self.phospho)
        site_metadata = _copy_frame(self.site_metadata)
        sample_metadata = _copy_optional_frame(self.sample_metadata)
        total = _copy_optional_frame(self.total)
        if isinstance(phospho, pd.DataFrame):
            phospho.index = canonicalize_site_index(
                phospho.index,
                field_name="dataset.phospho.index",
                error_type=DatasetValidationError,
            )
        if isinstance(site_metadata, pd.DataFrame):
            site_metadata.index = canonicalize_site_index(
                site_metadata.index,
                field_name="dataset.site_metadata.index",
                error_type=DatasetValidationError,
            )
        _DATASET_VALIDATOR.run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            total=total,
            organism=self.organism,
            transformation_state=self.transformation_state,
        )
        object.__setattr__(self, "phospho", phospho)
        object.__setattr__(self, "site_metadata", site_metadata)
        object.__setattr__(self, "sample_metadata", sample_metadata)
        object.__setattr__(self, "total", total)


def _copy_frame(value: object) -> object:
    if not isinstance(value, pd.DataFrame):
        return value
    return value.copy(deep=True)


def _copy_optional_frame(value: object | None) -> object | None:
    if value is None:
        return None
    return _copy_frame(value)
