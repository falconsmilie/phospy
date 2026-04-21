"""Internal executor for the dataset builder path.

The public builder lane stays intentionally narrow: establish supported
transformation state after applying explicit builder preprocessing policy.
"""

from __future__ import annotations

from phospy.datasets.builders.contracts import (
    DatasetPreprocessorContract,
    InterpretedDatasetBuildRequest,
)
from phospy.datasets.builders.preprocessing import DatasetPreprocessor
from phospy.datasets.builders.transformation_resolver import (
    DatasetTransformationResolver,
)
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.build import DatasetBuildError
from phospy.errors.input import PhosPyInputError
from phospy.errors.transformations import (
    PhosPyTransformationError,
    TransformationStateEstablishmentError,
)
from phospy.errors.validation import PhosPyValidationError
from phospy.transformations.contracts import Transformer
from phospy.transformations.transformers import IdentityTransformer


class DatasetBuildExecutor:
    """Construct `AnalysisReadyPhosphoDataset` from interpreted builder input.

    Default policy uses the identity transformer, which is a pass-through
    establisher for already-prepared quantitative matrices after internal
    preprocessing stages (including optional site-matrix construction).
    """

    def __init__(
        self,
        *,
        transformer: Transformer | None = None,
        transformation_resolver: DatasetTransformationResolver | None = None,
        preprocessor: DatasetPreprocessorContract | None = None,
    ) -> None:
        self._transformation_resolver = (
            transformation_resolver
            or DatasetTransformationResolver(
                transformer=transformer or IdentityTransformer()
            )
        )
        self._preprocessor = preprocessor or DatasetPreprocessor()

    def run(
        self, request: InterpretedDatasetBuildRequest
    ) -> AnalysisReadyPhosphoDataset:
        try:
            preprocessed = self._preprocessor.run(
                phospho=request.phospho,
                site_metadata=request.site_metadata,
                sample_metadata=request.sample_metadata,
                total=request.total,
                plan=request.preprocessing_plan,
            )
            resolved = self._transformation_resolver.run(
                phospho=preprocessed.phospho,
                total=preprocessed.total,
            )
            if not resolved.transformation_state.is_established:
                raise TransformationStateEstablishmentError(
                    "transformation resolver returned a non-established "
                    "transformation state; this violates the dataset boundary "
                    "contract"
                )
            return AnalysisReadyPhosphoDataset._from_owned(
                phospho=resolved.phospho,
                site_metadata=preprocessed.site_metadata,
                sample_metadata=preprocessed.sample_metadata,
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
