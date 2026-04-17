"""Public dataset builder shell."""

from __future__ import annotations

import pandas as pd

from phospy.api.requests import DatasetBuildRequest
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.build import DatasetBuildError
from phospy.errors.input import PhosPyInputError, UnsupportedInputFormatError
from phospy.errors.validation import PhosPyValidationError
from phospy.transformations.transformers.identity import IdentityTransformer


class AnalysisReadyDatasetBuilder:
    """Build ``AnalysisReadyPhosphoDataset`` from a structured request."""

    def run(self, request: DatasetBuildRequest) -> AnalysisReadyPhosphoDataset:
        try:
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
            if request.transformation_state is None:
                transformed = IdentityTransformer().run(phospho=phospho, total=total)
                phospho = transformed.phospho
                total = transformed.total
                transformation_state = transformed.state
            else:
                transformation_state = request.transformation_state
            return AnalysisReadyPhosphoDataset(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                total=total,
                organism=request.organism,
                transformation_state=transformation_state,
            )
        except (PhosPyInputError, PhosPyValidationError):
            raise
        except Exception as exc:  # pragma: no cover - defensive boundary translation
            raise DatasetBuildError(
                "failed to construct AnalysisReadyPhosphoDataset from build request"
            ) from exc

    @staticmethod
    def _resolve_frame(value: object, field_name: str) -> pd.DataFrame:
        if isinstance(value, pd.DataFrame):
            return value
        raise UnsupportedInputFormatError(
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
