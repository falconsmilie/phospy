"""Internal executor for the dataset builder path."""

from __future__ import annotations

from phospy.datasets.builders.contracts import InterpretedDatasetBuildRequest
from phospy.datasets.builders.transformation_resolver import (
    DatasetTransformationResolver,
)
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.build import DatasetBuildError
from phospy.errors.input import PhosPyInputError
from phospy.errors.transformations import PhosPyTransformationError
from phospy.errors.validation import PhosPyValidationError
from phospy.transformations.contracts import Transformer


class DatasetBuildExecutor:
    """Construct `AnalysisReadyPhosphoDataset` from interpreted builder input."""

    def __init__(
        self,
        *,
        transformer: Transformer | None = None,
        transformation_resolver: DatasetTransformationResolver | None = None,
    ) -> None:
        self._transformation_resolver = (
            transformation_resolver
            or DatasetTransformationResolver(transformer=transformer)
        )

    def run(
        self, request: InterpretedDatasetBuildRequest
    ) -> AnalysisReadyPhosphoDataset:
        try:
            resolved = self._transformation_resolver.run(
                phospho=request.phospho,
                total=request.total,
                transformation_state=request.transformation_state,
            )
            return AnalysisReadyPhosphoDataset(
                phospho=resolved.phospho,
                site_metadata=request.site_metadata,
                sample_metadata=request.sample_metadata,
                total=resolved.total,
                organism=request.organism,
                transformation_state=resolved.transformation_state,
            )
        except (
            PhosPyInputError,
            PhosPyTransformationError,
            PhosPyValidationError,
            DatasetBuildError,
        ):
            raise
        except Exception as exc:  # pragma: no cover - defensive boundary translation
            raise DatasetBuildError(
                "failed to construct AnalysisReadyPhosphoDataset from interpreted input"
            ) from exc
