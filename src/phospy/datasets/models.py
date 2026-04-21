"""Dataset domain models."""

from __future__ import annotations

from dataclasses import InitVar, dataclass

import pandas as pd

from phospy._frame_ownership import own_dataframe, own_optional_dataframe
from phospy.errors.validation import DatasetValidationError
from phospy.references.models import Organism
from phospy.transformations.models import TransformationState
from phospy.validation.datasets.analysis_ready import AnalysisReadyDatasetValidator
from phospy.validation.transformations.state import TransformationStateValidator

_DATASET_VALIDATOR = AnalysisReadyDatasetValidator()
_TRANSFORMATION_STATE_VALIDATOR = TransformationStateValidator()


@dataclass(frozen=True, slots=True)
class AnalysisReadyPhosphoDataset:
    """Public analysis-ready dataset contract.

    `phospho` stores the quantitative matrix after builder preprocessing policy
    has been applied. When total/protein correction is enabled in the builder
    lane, corrected values are represented directly in this matrix. When
    site-matrix construction is enabled in the builder lane, this matrix already
    reflects the constructed site-matrix-ready rows. Intermediate site-matrix
    artefacts remain private to preprocessing internals.
    """

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    transformation_state: TransformationState
    sample_metadata: pd.DataFrame | None = None
    total: pd.DataFrame | None = None
    organism: Organism | None = None
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        phospho = own_dataframe(
            self.phospho,
            field_name="dataset.phospho",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        site_metadata = own_dataframe(
            self.site_metadata,
            field_name="dataset.site_metadata",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        sample_metadata = own_optional_dataframe(
            self.sample_metadata,
            field_name="dataset.sample_metadata",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        total = own_optional_dataframe(
            self.total,
            field_name="dataset.total",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        _DATASET_VALIDATOR.run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            total=total,
            organism=self.organism,
        )
        _TRANSFORMATION_STATE_VALIDATOR.run(
            transformation_state=self.transformation_state,
            has_total_matrix=total is not None,
            require_established=True,
        )
        object.__setattr__(self, "phospho", phospho)
        object.__setattr__(self, "site_metadata", site_metadata)
        object.__setattr__(self, "sample_metadata", sample_metadata)
        object.__setattr__(self, "total", total)

    @classmethod
    def _from_owned(
        cls,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        transformation_state: TransformationState,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        organism: Organism | None = None,
    ) -> AnalysisReadyPhosphoDataset:
        return cls(
            phospho=phospho,
            site_metadata=site_metadata,
            transformation_state=transformation_state,
            sample_metadata=sample_metadata,
            total=total,
            organism=organism,
            _assume_owned=True,
        )
