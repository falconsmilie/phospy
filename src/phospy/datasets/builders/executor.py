"""Internal executor for the dataset builder path."""

from __future__ import annotations

from phospy.datasets.builders.contracts import InterpretedDatasetBuildRequest
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.build import DatasetBuildError
from phospy.errors.input import PhosPyInputError
from phospy.errors.transformations import PhosPyTransformationError
from phospy.errors.validation import PhosPyValidationError
from phospy.transformations.contracts import Transformer
from phospy.transformations.transformers.identity import IdentityTransformer


class DatasetBuildExecutor:
    """Construct `AnalysisReadyPhosphoDataset` from interpreted builder input."""

    def __init__(self, *, transformer: Transformer | None = None) -> None:
        self._transformer = transformer or IdentityTransformer()

    def run(
        self, request: InterpretedDatasetBuildRequest
    ) -> AnalysisReadyPhosphoDataset:
        try:
            phospho = request.phospho
            total = request.total
            transformation_state = request.transformation_state
            if transformation_state is None:
                transformed = self._transformer.run(phospho=phospho, total=total)
                phospho = transformed.phospho
                total = transformed.total
                transformation_state = transformed.state
            return AnalysisReadyPhosphoDataset(
                phospho=phospho,
                site_metadata=request.site_metadata,
                sample_metadata=request.sample_metadata,
                total=total,
                organism=request.organism,
                transformation_state=transformation_state,
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
