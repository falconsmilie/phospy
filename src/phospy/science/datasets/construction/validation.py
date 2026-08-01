"""Analysis-ready dataset construction validation guards."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import cast

import numpy as np
import pandas as pd

from phospy.errors.validation import DatasetValidationError
from phospy.frames.ownership import own_dataframe, own_optional_dataframe
from phospy.frames.validation import (
    require_dataframe,
    require_exact_index_match,
    require_finite_numeric_dataframe,
    require_non_empty_dataframe,
    require_numeric_dataframe,
    require_unique_columns,
)
from phospy.provenance.models import (
    RunProvenance,
    TrustedDatasetConstructionAssertions,
    TrustedDatasetConstructionEvidence,
)
from phospy.science.transformations.models import IntensityScaleState
from phospy.science.transformations.state_coherence import (
    require_intensity_scale_state_coherence,
    require_quantitative_numeric_domain_coherence,
)

_NUMPY_DTYPES_WITH_NAN_SENTINELS = frozenset(("f", "c"))

_NUMPY_DTYPES_WITHOUT_MISSING_SENTINELS = frozenset(("i", "u"))


class _IntensityScaleStateValidator:
    def run(
        self,
        *,
        intensity_scale_state: IntensityScaleState,
        has_total_matrix: bool,
        require_established: bool = False,
    ) -> IntensityScaleState:
        return require_intensity_scale_state_coherence(
            intensity_scale_state=intensity_scale_state,
            has_total_matrix=has_total_matrix,
            require_established=require_established,
        )


_INTENSITY_SCALE_STATE_VALIDATOR = _IntensityScaleStateValidator()


class _QuantitativeNumericDomainValidator:
    def run(
        self,
        *,
        phospho: pd.DataFrame,
        total: pd.DataFrame | None,
        intensity_scale_state: IntensityScaleState,
        trusted_construction_assertions: TrustedDatasetConstructionAssertions
        | None = None,
    ) -> None:
        require_quantitative_numeric_domain_coherence(
            phospho=phospho,
            total=total,
            intensity_scale_state=intensity_scale_state,
            allow_numeric_semantic_domain_waiver=(
                _has_numeric_semantic_domain_waiver(trusted_construction_assertions)
            ),
            error_type=DatasetValidationError,
        )


def _has_numeric_semantic_domain_waiver(
    trusted_construction_assertions: TrustedDatasetConstructionAssertions | None,
) -> bool:
    if trusted_construction_assertions is None:
        return False
    evidence = trusted_construction_assertions.numeric_semantic_domain
    return isinstance(evidence, TrustedDatasetConstructionEvidence) and (
        evidence.is_waiver
    )


_QUANTITATIVE_NUMERIC_DOMAIN_VALIDATOR = _QuantitativeNumericDomainValidator()


@dataclass(frozen=True, slots=True)
class _OwnedDatasetFrames:
    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame | None
    total: pd.DataFrame | None
    comparisons: pd.DataFrame | None
    imputation_observation_mask: pd.DataFrame | None


def _own_dataset_frames(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    sample_metadata: pd.DataFrame | None,
    total: pd.DataFrame | None,
    comparisons: pd.DataFrame | None,
    imputation_observation_mask: pd.DataFrame | None,
    assume_owned: bool,
) -> _OwnedDatasetFrames:
    return _OwnedDatasetFrames(
        phospho=own_dataframe(
            phospho,
            field_name="dataset.phospho",
            error_type=DatasetValidationError,
            assume_owned=assume_owned,
        ),
        site_metadata=own_dataframe(
            site_metadata,
            field_name="dataset.site_metadata",
            error_type=DatasetValidationError,
            assume_owned=assume_owned,
        ),
        sample_metadata=own_optional_dataframe(
            sample_metadata,
            field_name="dataset.sample_metadata",
            error_type=DatasetValidationError,
            assume_owned=assume_owned,
        ),
        total=own_optional_dataframe(
            total,
            field_name="dataset.total",
            error_type=DatasetValidationError,
            assume_owned=assume_owned,
        ),
        comparisons=own_optional_dataframe(
            comparisons,
            field_name="dataset.comparisons",
            error_type=DatasetValidationError,
            assume_owned=assume_owned,
        ),
        imputation_observation_mask=own_optional_dataframe(
            imputation_observation_mask,
            field_name="dataset.imputation_observation_mask",
            error_type=DatasetValidationError,
            assume_owned=assume_owned,
        ),
    )


def _validate_optional_comparisons(
    *,
    comparisons: pd.DataFrame | None,
    expected_index: pd.Index,
) -> pd.DataFrame | None:
    if comparisons is None:
        return None

    comparisons_frame = require_dataframe(
        comparisons,
        field_name="dataset.comparisons",
        allow_empty=True,
        error_type=DatasetValidationError,
    )
    require_non_empty_dataframe(
        comparisons_frame,
        field_name="dataset.comparisons",
        error_type=DatasetValidationError,
    )
    require_numeric_dataframe(
        comparisons_frame,
        field_name="dataset.comparisons",
        error_type=DatasetValidationError,
    )
    require_finite_numeric_dataframe(
        comparisons_frame,
        field_name="dataset.comparisons",
        error_type=DatasetValidationError,
        allow_missing=False,
    )
    require_unique_columns(
        comparisons_frame,
        field_name="dataset.comparisons",
        error_type=DatasetValidationError,
    )
    require_exact_index_match(
        left=comparisons_frame.index,
        right=expected_index,
        left_name="dataset.comparisons.index",
        right_name="dataset.phospho.index",
        error_type=DatasetValidationError,
    )
    return comparisons_frame


_DIRECT_CONSTRUCTION_ERROR_MESSAGE = (
    "AnalysisReadyPhosphoDataset(...) direct construction is no longer supported. "
    "Use AnalysisReadyDatasetBuilder for ordinary construction, or "
    "AnalysisReadyPhosphoDataset.from_trusted_tables(...) with complete "
    "TrustedDatasetConstructionAssertions for advanced trusted reconstruction; "
    "numeric-semantic-domain conflicts additionally require the typed "
    "numeric_semantic_domain waiver."
)


def _require_builder_output_provenance(provenance: object) -> None:
    if not isinstance(provenance, RunProvenance):
        raise DatasetValidationError(
            "AnalysisReadyPhosphoDataset._from_builder_output requires "
            "builder-created RunProvenance"
        )
    construction_raw = provenance.workflow_parameters.get("construction")
    if provenance.workflow_name != "dataset_builder" or not isinstance(
        construction_raw,
        Mapping,
    ):
        raise DatasetValidationError(
            "AnalysisReadyPhosphoDataset._from_builder_output requires "
            "AnalysisReadyDatasetBuilder provenance"
        )
    construction = cast(Mapping[str, object], construction_raw)
    if construction.get("method") != "AnalysisReadyDatasetBuilder.run":
        raise DatasetValidationError(
            "AnalysisReadyPhosphoDataset._from_builder_output requires "
            "provenance.workflow_parameters['construction']['method'] to be "
            "'AnalysisReadyDatasetBuilder.run'"
        )
    if not isinstance(construction.get("processing_state_establishment"), Mapping):
        raise DatasetValidationError(
            "AnalysisReadyPhosphoDataset._from_builder_output requires "
            "builder-created processing state establishment provenance"
        )


def analysis_ready_matrix_missing_value_count(matrix: pd.DataFrame) -> int:
    if _can_use_fast_numeric_missing_value_scan(matrix):
        try:
            return _fast_numeric_missing_value_count(matrix)
        except (AttributeError, TypeError, ValueError):
            pass
    return _object_level_missing_value_count(matrix)


def _can_use_fast_numeric_missing_value_scan(matrix: pd.DataFrame) -> bool:
    for dtype in matrix.dtypes:
        if pd.api.types.is_bool_dtype(dtype):
            return False
        if not pd.api.types.is_numeric_dtype(dtype):
            return False
    return True


def _fast_numeric_missing_value_count(matrix: pd.DataFrame) -> int:
    values = matrix.to_numpy(copy=False)
    if values.dtype.kind in _NUMPY_DTYPES_WITHOUT_MISSING_SENTINELS:
        return 0
    if values.dtype.kind in _NUMPY_DTYPES_WITH_NAN_SENTINELS:
        return int(np.count_nonzero(np.isnan(values)))

    missing_mask = np.asarray(pd.isna(values), dtype=bool)
    return int(np.count_nonzero(missing_mask))


def _is_missing_value(value: object) -> bool:
    from phospy.science.datasets.processing_state import (
        is_missing_value as processing_state_is_missing_value,
    )

    return processing_state_is_missing_value(value)


def _object_level_missing_value_count(matrix: pd.DataFrame) -> int:
    return sum(
        1
        for value in matrix.to_numpy(dtype="object").ravel()
        if _is_missing_value(value)
    )
