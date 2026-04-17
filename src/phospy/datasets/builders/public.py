"""Public dataset builder shell."""

from __future__ import annotations

import pandas as pd

from phospy.api.requests import DatasetBuildRequest
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.build import PhosPyBuildError
from phospy.errors.input import PhosPyInputError
from phospy.transformations.models import TransformationState


class AnalysisReadyDatasetBuilder:
    """Build ``AnalysisReadyPhosphoDataset`` from a structured request."""

    def run(self, request: DatasetBuildRequest) -> AnalysisReadyPhosphoDataset:
        phospho = self._resolve_frame(request.phospho, "phospho")
        site_metadata = self._resolve_optional_frame(
            request.site_metadata,
            "site_metadata",
            default=pd.DataFrame(index=phospho.index),
        )
        sample_metadata = self._resolve_optional_frame(
            request.sample_metadata,
            "sample_metadata",
            default=None,
        )
        total = self._resolve_optional_frame(
            request.total,
            "total",
            default=None,
        )
        transformation_state = request.transformation_state or TransformationState()
        if phospho.empty:
            raise PhosPyBuildError("phospho input cannot be empty for dataset build")
        return AnalysisReadyPhosphoDataset(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            total=total,
            organism=request.organism,
            transformation_state=transformation_state,
        )

    @staticmethod
    def _resolve_frame(value: object, field_name: str) -> pd.DataFrame:
        if isinstance(value, pd.DataFrame):
            return value
        raise PhosPyInputError(
            f"{field_name} must be provided as a pandas DataFrame in this rewrite phase"
        )

    @classmethod
    def _resolve_optional_frame(
        cls,
        value: object | None,
        field_name: str,
        default: pd.DataFrame | None,
    ) -> pd.DataFrame | None:
        if value is None:
            return default
        return cls._resolve_frame(value, field_name)
