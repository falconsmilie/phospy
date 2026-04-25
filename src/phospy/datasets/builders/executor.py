"""Internal executor for the dataset builder path.

The public builder lane stays intentionally narrow: establish supported
transformation state after applying explicit builder preprocessing policy.
"""

from __future__ import annotations

import pandas as pd

from phospy.datasets.builders.contracts import (
    DatasetPreprocessorContract,
    InterpretedDatasetBuildRequest,
)
from phospy.datasets.builders.preprocessing import DatasetPreprocessor
from phospy.datasets.builders.transformation_resolver import (
    DatasetTransformationResolver,
)
from phospy.datasets.models import (
    AnalysisReadyPhosphoDataset,
    DatasetPreprocessingReport,
)
from phospy.errors.build import DatasetBuildError
from phospy.errors.input import PhosPyInputError
from phospy.errors.transformations import (
    PhosPyTransformationError,
    TransformationStateEstablishmentError,
)
from phospy.errors.validation import PhosPyValidationError
from phospy.transformations.contracts import Transformer
from phospy.transformations.transformers import IdentityTransformer

_ROW_COUNT_COLUMNS = ("stage", "input_rows", "output_rows", "dropped_rows")
_OPERATION_COLUMNS = (
    "step_order",
    "stage",
    "operation",
    "parameters",
    "input_rows",
    "output_rows",
    "notes",
)
_DUPLICATE_SITE_RESOLUTION_COLUMNS = (
    "site_id",
    "source_row_id",
    "retained",
    "resolution_strategy",
    "retained_reason",
    "dropped_reason",
    "observed_values",
    "mean_signal",
    "n_source_rows",
    "n_aggregated_rows",
    "source_protein_id",
    "source_gene_symbol",
    "source_site",
    "source_site_sequence",
    "metadata_conflict_detected",
)
_METADATA_CONFLICT_COLUMNS = (
    "site_id",
    "field",
    "values",
    "n_distinct_values",
    "source_row_ids",
)
_FINAL_DATASET_STAGE = "final_dataset_construction"


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
            report = _build_dataset_preprocessing_report(
                row_counts=preprocessed.preprocessing_row_counts,
                operations=preprocessed.preprocessing_operations,
                duplicate_site_resolution=preprocessed.duplicate_site_resolution,
                metadata_conflicts=preprocessed.metadata_conflicts,
                final_dataset_rows=int(len(resolved.phospho.index)),
                transformation_state_label=resolved.transformation_state.label,
            )
            return AnalysisReadyPhosphoDataset._from_owned(
                phospho=resolved.phospho,
                site_metadata=preprocessed.site_metadata,
                sample_metadata=preprocessed.sample_metadata,
                total=resolved.total,
                comparisons=preprocessed.comparisons,
                organism=request.organism,
                transformation_state=resolved.transformation_state,
                preprocessing_report=report,
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


def _build_dataset_preprocessing_report(
    *,
    row_counts: pd.DataFrame | None,
    operations: pd.DataFrame | None,
    duplicate_site_resolution: pd.DataFrame | None,
    metadata_conflicts: pd.DataFrame | None,
    final_dataset_rows: int,
    transformation_state_label: str,
) -> DatasetPreprocessingReport:
    base_row_counts = (
        pd.DataFrame.from_records([], columns=_ROW_COUNT_COLUMNS)
        if row_counts is None
        else row_counts
    )
    base_operations = (
        pd.DataFrame.from_records([], columns=_OPERATION_COLUMNS)
        if operations is None
        else operations
    )
    base_duplicate_site_resolution = (
        pd.DataFrame.from_records([], columns=_DUPLICATE_SITE_RESOLUTION_COLUMNS)
        if duplicate_site_resolution is None
        else duplicate_site_resolution
    )
    base_metadata_conflicts = (
        pd.DataFrame.from_records([], columns=_METADATA_CONFLICT_COLUMNS)
        if metadata_conflicts is None
        else metadata_conflicts
    )
    final_row_counts = pd.concat(
        [
            base_row_counts,
            pd.DataFrame.from_records(
                [
                    {
                        "stage": _FINAL_DATASET_STAGE,
                        "input_rows": final_dataset_rows,
                        "output_rows": final_dataset_rows,
                        "dropped_rows": 0,
                    }
                ],
                columns=_ROW_COUNT_COLUMNS,
            ),
        ],
        axis=0,
        ignore_index=True,
    )
    if base_operations.empty:
        final_step_order = 1
    else:
        final_step_order = int(base_operations.loc[:, "step_order"].max()) + 1
    final_operations = pd.concat(
        [
            base_operations,
            pd.DataFrame.from_records(
                [
                    {
                        "step_order": final_step_order,
                        "stage": _FINAL_DATASET_STAGE,
                        "operation": "construct_analysis_ready_dataset",
                        "parameters": {
                            "transformation_state_label": transformation_state_label
                        },
                        "input_rows": final_dataset_rows,
                        "output_rows": final_dataset_rows,
                        "notes": "analysis-ready dataset boundary construction",
                    }
                ],
                columns=_OPERATION_COLUMNS,
            ),
        ],
        axis=0,
        ignore_index=True,
    )
    return DatasetPreprocessingReport._from_owned(
        row_counts=final_row_counts,
        operations=final_operations,
        duplicate_site_resolution=base_duplicate_site_resolution,
        metadata_conflicts=base_metadata_conflicts,
    )
