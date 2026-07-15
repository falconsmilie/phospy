"""Dataset domain models."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from phospy.errors.validation import DatasetValidationError
from phospy.frames.ownership import (
    borrow_dataframe,
    borrow_optional_dataframe,
    export_dataframe,
    export_optional_dataframe,
    own_dataframe,
    own_optional_dataframe,
)
from phospy.provenance.models import RunProvenance, TrustedDatasetConstructionAssertions
from phospy.science.datasets.direct_construction import (
    build_direct_construction_provenance,
)
from phospy.science.datasets.imputation_metadata import (
    IMPUTATION_FEATURE_METADATA_COLUMNS as IMPUTATION_FEATURE_METADATA_COLUMNS,
)
from phospy.science.datasets.imputation_metadata import (
    IMPUTATION_OBSERVATION_SUMMARY_COLUMNS as IMPUTATION_OBSERVATION_SUMMARY_COLUMNS,
)
from phospy.science.datasets.imputation_metadata import (
    ImputationObservationMetadata,
)
from phospy.science.datasets.imputation_metadata import (
    build_imputation_observation_metadata_or_none as _build_imputation_observation_metadata_or_none,
)
from phospy.science.datasets.imputation_metadata import (
    require_boolean_observation_mask as _require_boolean_observation_mask,
)
from phospy.science.datasets.organism_coherence import (
    normalize_dataset_organism_state,
    resolve_single_dataset_organism,
)
from phospy.science.datasets.preprocessing.protein_aware_preparation import (
    ProteinAwarePreparationResult,
)
from phospy.science.datasets.processing_state import (
    DatasetPreprocessingReport,
    DatasetProcessingState,
    RuvReadinessState,
)
from phospy.science.datasets.processing_state import (
    PreprocessingSiteAttritionSummary as PreprocessingSiteAttritionSummary,
)
from phospy.science.datasets.processing_state import (
    SiteSequenceResolutionReport as SiteSequenceResolutionReport,
)
from phospy.science.datasets.processing_state import (
    is_missing_value as _is_missing_value,
)
from phospy.science.datasets.processing_state import (
    missing_data_state_claims_no_missing_values as _missing_data_state_claims_no_missing_values,
)
from phospy.science.datasets.processing_state import (
    require_instance as _require_instance,
)
from phospy.science.datasets.processing_state import (
    require_optional_instance as _require_optional_instance,
)
from phospy.science.references.models import Organism, ReferenceContext
from phospy.science.transformations.models import IntensityScaleState
from phospy.tables.datasets import (
    PhosphoIntensityMatrix,
    SampleMetadataTable,
    SiteMetadataTable,
    TotalProteinMatrix,
)
from phospy.validation.common.dataframes import (
    require_dataframe,
    require_exact_index_match,
    require_finite_numeric_dataframe,
    require_non_empty_dataframe,
    require_numeric_dataframe,
    require_unique_columns,
)
from phospy.validation.transformations.state import IntensityScaleStateValidator

_PROCESSING_STATE_COMPAT_EXPORTS = (_require_boolean_observation_mask,)
_INTENSITY_SCALE_STATE_VALIDATOR = IntensityScaleStateValidator()
_NUMPY_DTYPES_WITH_NAN_SENTINELS = frozenset(("f", "c"))
_NUMPY_DTYPES_WITHOUT_MISSING_SENTINELS = frozenset(("i", "u"))


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


def _resolve_trusted_construction_assertions(
    *,
    trusted_construction_assertions: TrustedDatasetConstructionAssertions | None,
    provenance: RunProvenance | None,
    assume_owned: bool,
) -> TrustedDatasetConstructionAssertions | None:
    if trusted_construction_assertions is not None:
        return trusted_construction_assertions
    if provenance is None or not assume_owned:
        return TrustedDatasetConstructionAssertions.missing()
    return None


def _analysis_ready_matrix_missing_value_count(matrix: pd.DataFrame) -> int:
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


def _object_level_missing_value_count(matrix: pd.DataFrame) -> int:
    return sum(
        1
        for value in matrix.to_numpy(dtype="object").ravel()
        if _is_missing_value(value)
    )


@dataclass(frozen=True, slots=True, init=False)
class AnalysisReadyPhosphoDataset:
    """Public analysis-ready dataset contract.

    Direct construction is for trusted advanced/internal use. Ordinary users
    should construct datasets through ``AnalysisReadyDatasetBuilder.run(...)``,
    which owns user-input interpretation, preprocessing, processing-state
    establishment, and construction provenance. Advanced callers who already
    own complete analysis-ready tables should prefer
    ``AnalysisReadyPhosphoDataset.from_trusted_tables(...)`` as the explicit
    trusted construction lane. The public constructor remains available for
    compatibility.

    Construction validates structural invariants, including table shape,
    alignment, analysis-ready ``site_key`` identity, processing-state
    coherence, and established transformation state. It cannot prove the
    biological correctness of user-asserted provenance or scientific claims
    supplied by a direct caller.
    Direct construction without supplied assertion metadata records explicit
    missing trusted-construction assertions; direct construction without
    supplied provenance receives a minimal direct-construction provenance marker
    that records this audit limitation.

    `phospho` stores the quantitative matrix after builder preprocessing policy
    has been applied. When total/protein correction is enabled in the builder
    lane, corrected values are represented directly in this matrix. When
    site-matrix construction is enabled in the builder lane, this matrix already
    reflects the constructed site-matrix-ready rows. Intermediate site-matrix
    artefacts remain private to preprocessing internals.
    Optional `comparisons` can carry builder-constructed dataset-level pairwise
    columns aligned to `phospho.index`.
    Site identity is strict at this boundary: rows are indexed by encoded
    `site_key` values, while `display_id` preserves the human-readable
    `GENE;SITE;` label and must be coherent with the protein-context metadata.

    Provenance in this object describes owned internal state at creation time.
    Public export helpers return defensive snapshots; mutating exports does not
    mutate this owning dataset.
    Internal `_borrow_*` accessors are reserved for dataset-domain internal
    view construction and return mutation-isolated internal frame snapshots.
    """

    intensity_scale_state: IntensityScaleState
    processing_state: DatasetProcessingState
    organism: Organism | None = None
    preprocessing_report: DatasetPreprocessingReport | None = None
    protein_aware_preparation: ProteinAwarePreparationResult | None = None
    provenance: RunProvenance | None = None
    trusted_construction_assertions: TrustedDatasetConstructionAssertions | None = None
    allow_opaque_site_values: bool = False
    _phospho: pd.DataFrame = field(init=False, repr=False)
    _site_metadata: pd.DataFrame = field(init=False, repr=False)
    _sample_metadata: pd.DataFrame | None = field(init=False, repr=False)
    _total: pd.DataFrame | None = field(init=False, repr=False)
    _comparisons: pd.DataFrame | None = field(init=False, repr=False)
    _imputation_observation_metadata: ImputationObservationMetadata | None = field(
        init=False,
        repr=False,
    )
    _allow_opaque_site_values: bool = field(init=False, repr=False, default=False)

    def __init__(
        self,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        imputation_observation_mask: pd.DataFrame | None = None,
        organism: Organism | None = None,
        preprocessing_report: DatasetPreprocessingReport | None = None,
        protein_aware_preparation: ProteinAwarePreparationResult | None = None,
        provenance: RunProvenance | None = None,
        allow_opaque_site_values: bool = False,
        trusted_construction_assertions: TrustedDatasetConstructionAssertions
        | None = None,
        _assume_owned: bool = False,
    ) -> None:
        _require_instance(
            allow_opaque_site_values,
            expected_type=bool,
            error_message="dataset.allow_opaque_site_values must be a bool",
        )
        frames = _own_dataset_frames(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            total=total,
            comparisons=comparisons,
            imputation_observation_mask=imputation_observation_mask,
            assume_owned=_assume_owned,
        )
        phospho = frames.phospho
        site_metadata = frames.site_metadata
        sample_metadata = frames.sample_metadata
        total = frames.total
        comparisons = frames.comparisons
        imputation_observation_mask = frames.imputation_observation_mask
        normalized_organism_state = normalize_dataset_organism_state(
            phospho=phospho,
            site_metadata=site_metadata,
            error_type=DatasetValidationError,
        )
        phospho = normalized_organism_state.phospho
        site_metadata = normalized_organism_state.site_metadata
        _require_instance(
            processing_state,
            expected_type=DatasetProcessingState,
            error_message=(
                "dataset.processing_state must be a DatasetProcessingState instance"
            ),
        )
        raw_missing_value_count = _analysis_ready_matrix_missing_value_count(phospho)
        if raw_missing_value_count > 0 and _missing_data_state_claims_no_missing_values(
            processing_state
        ):
            raise DatasetValidationError(
                "dataset.phospho must not contain missing values; "
                "dataset.processing_state.missing_data claims no missing values "
                "but dataset.phospho contains missing values"
            )
        phospho_table = PhosphoIntensityMatrix(
            frame=phospho,
            _assume_owned=True,
        )
        site_metadata_table = SiteMetadataTable(
            frame=site_metadata,
            expected_index=phospho_table.frame.index,
            allow_opaque_site_values=allow_opaque_site_values,
            _assume_owned=True,
        )
        sample_metadata_table = (
            None
            if sample_metadata is None
            else SampleMetadataTable(
                frame=sample_metadata,
                expected_index=phospho_table.frame.columns,
                _assume_owned=True,
            )
        )
        total_table = (
            None
            if total is None
            else TotalProteinMatrix(
                frame=total,
                expected_sample_index=phospho_table.frame.columns,
                _assume_owned=True,
            )
        )
        if processing_state.total_protein_correction.applied and total_table is None:
            raise DatasetValidationError(
                "dataset.processing_state.total_protein_correction.applied "
                "requires dataset.total"
            )
        validated_intensity_scale_state = _INTENSITY_SCALE_STATE_VALIDATOR.run(
            intensity_scale_state=intensity_scale_state,
            has_total_matrix=total_table is not None,
            require_established=True,
        )
        _require_optional_instance(
            organism,
            expected_type=Organism,
            error_message="dataset.organism must be an Organism enum value or None",
        )
        resolved_dataset_organism = resolve_single_dataset_organism(
            site_metadata=site_metadata_table.frame,
            organism=organism,
            error_type=DatasetValidationError,
        )
        _require_instance(
            processing_state.ruv_readiness,
            expected_type=RuvReadinessState,
            error_message=(
                "dataset.processing_state.ruv_readiness must be a "
                "RuvReadinessState instance"
            ),
        )
        if (
            not processing_state.ruv_readiness.enabled
            and processing_state.ruv_readiness.ready
        ):
            raise DatasetValidationError(
                "dataset.processing_state.ruv_readiness.ready must be False when "
                "ruv_readiness.enabled is False"
            )
        if processing_state.intensity_scale != validated_intensity_scale_state:
            raise DatasetValidationError(
                "dataset.processing_state.intensity_scale must match "
                "dataset.intensity_scale_state"
            )
        if not bool(processing_state.missing_data.complete_matrix):
            raise DatasetValidationError(
                "dataset.processing_state.missing_data.complete_matrix must be True "
                "at AnalysisReadyPhosphoDataset boundary"
            )
        _require_optional_instance(
            preprocessing_report,
            expected_type=DatasetPreprocessingReport,
            error_message=(
                "dataset.preprocessing_report must be DatasetPreprocessingReport "
                "or None"
            ),
        )
        _require_optional_instance(
            protein_aware_preparation,
            expected_type=ProteinAwarePreparationResult,
            error_message=(
                "dataset.protein_aware_preparation must be "
                "ProteinAwarePreparationResult or None"
            ),
        )
        _require_optional_instance(
            provenance,
            expected_type=RunProvenance,
            error_message="dataset.provenance must be RunProvenance or None",
        )
        _require_optional_instance(
            trusted_construction_assertions,
            expected_type=TrustedDatasetConstructionAssertions,
            error_message=(
                "dataset.trusted_construction_assertions must be "
                "TrustedDatasetConstructionAssertions or None"
            ),
        )
        resolved_trusted_construction_assertions = (
            _resolve_trusted_construction_assertions(
                trusted_construction_assertions=trusted_construction_assertions,
                provenance=provenance,
                assume_owned=_assume_owned,
            )
        )
        comparisons = _validate_optional_comparisons(
            comparisons=comparisons,
            expected_index=phospho_table.frame.index,
        )
        imputation_observation_metadata = (
            _build_imputation_observation_metadata_or_none(
                imputation_observation_mask=imputation_observation_mask,
                phospho_index=phospho_table.frame.index,
                sample_index=phospho_table.frame.columns,
            )
        )
        if provenance is None:
            provenance = build_direct_construction_provenance(
                phospho=phospho_table.frame,
                site_metadata=site_metadata_table.frame,
                sample_metadata=(
                    None
                    if sample_metadata_table is None
                    else sample_metadata_table.frame
                ),
                total=None if total_table is None else total_table.frame,
                comparisons=comparisons,
                imputation_observation_mask=imputation_observation_mask,
                trusted_construction_assertions=(
                    resolved_trusted_construction_assertions
                ),
            )
        object.__setattr__(
            self, "intensity_scale_state", validated_intensity_scale_state
        )
        object.__setattr__(self, "processing_state", processing_state)
        object.__setattr__(self, "organism", resolved_dataset_organism)
        object.__setattr__(self, "preprocessing_report", preprocessing_report)
        object.__setattr__(
            self,
            "trusted_construction_assertions",
            resolved_trusted_construction_assertions,
        )
        object.__setattr__(
            self,
            "protein_aware_preparation",
            protein_aware_preparation,
        )
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(self, "allow_opaque_site_values", allow_opaque_site_values)
        object.__setattr__(self, "_phospho", phospho_table.frame)
        object.__setattr__(self, "_site_metadata", site_metadata_table.frame)
        object.__setattr__(
            self,
            "_sample_metadata",
            None if sample_metadata_table is None else sample_metadata_table.frame,
        )
        object.__setattr__(
            self, "_total", None if total_table is None else total_table.frame
        )
        object.__setattr__(self, "_comparisons", comparisons)
        object.__setattr__(
            self,
            "_imputation_observation_metadata",
            imputation_observation_metadata,
        )
        object.__setattr__(self, "_allow_opaque_site_values", allow_opaque_site_values)

    @property
    def phospho(self) -> pd.DataFrame:
        return export_dataframe(self._phospho)

    @property
    def site_metadata(self) -> pd.DataFrame:
        return export_dataframe(self._site_metadata)

    @property
    def sample_metadata(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._sample_metadata)

    @property
    def total(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._total)

    @property
    def comparisons(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._comparisons)

    @property
    def imputation_observation_metadata(
        self,
    ) -> ImputationObservationMetadata | None:
        return self._imputation_observation_metadata

    @property
    def imputation_feature_metadata(self) -> pd.DataFrame | None:
        metadata = self._imputation_observation_metadata
        if metadata is None:
            return None
        return metadata.feature_summary

    @property
    def opaque_site_values_allowed(self) -> bool:
        return bool(self._allow_opaque_site_values)

    @property
    def reference_context(self) -> ReferenceContext | None:
        if self.provenance is None:
            return None
        return self.provenance.reference_context

    def _borrow_phospho_frame(self) -> pd.DataFrame:
        """Package-private phospho snapshot for DatasetInternalView."""

        return borrow_dataframe(self._phospho)

    def _borrow_site_metadata_frame(self) -> pd.DataFrame:
        """Package-private site-metadata snapshot for DatasetInternalView."""

        return borrow_dataframe(self._site_metadata)

    def _borrow_sample_metadata_frame(self) -> pd.DataFrame | None:
        """Package-private sample-metadata snapshot for DatasetInternalView."""

        return borrow_optional_dataframe(self._sample_metadata)

    def _borrow_total_frame(self) -> pd.DataFrame | None:
        """Package-private total-protein snapshot for DatasetInternalView."""

        return borrow_optional_dataframe(self._total)

    def _borrow_comparisons_frame(self) -> pd.DataFrame | None:
        """Package-private comparisons snapshot for DatasetInternalView."""

        return borrow_optional_dataframe(self._comparisons)

    def _borrow_imputation_observed_mask_frame(self) -> pd.DataFrame | None:
        """Package-private observation-mask snapshot for internal read paths."""

        metadata = self._imputation_observation_metadata
        if metadata is None:
            return None
        return metadata.observed_mask_dataframe()

    def _imputation_observation_summary_frame(
        self,
        *,
        feature_ids: Sequence[object],
        sample_ids: Sequence[object],
    ) -> pd.DataFrame | None:
        """Package-private summary of imputation observations for internals."""

        return self.imputation_observation_summary_dataframe(
            feature_ids=feature_ids,
            sample_ids=sample_ids,
        )

    @classmethod
    def from_trusted_tables(
        cls,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        imputation_observation_mask: pd.DataFrame | None = None,
        organism: Organism | None = None,
        preprocessing_report: DatasetPreprocessingReport | None = None,
        protein_aware_preparation: ProteinAwarePreparationResult | None = None,
        provenance: RunProvenance | None = None,
        trusted_construction_assertions: TrustedDatasetConstructionAssertions
        | None = None,
        allow_opaque_site_values: bool = False,
    ) -> AnalysisReadyPhosphoDataset:
        """Construct from caller-owned analysis-ready tables.

        This explicit trusted factory is for advanced/internal callers that
        already own fully prepared ``site_key``-indexed tables, complete
        ``site_sequence`` evidence, established intensity-scale state, and
        coherent processing state. It delegates to the existing constructor and
        enforces the same structural invariants as direct construction,
        including table shape, alignment, site identity, ``site_sequence``
        validation, and state coherence.

        Validation can confirm structural consistency, but it cannot prove the
        biological correctness of caller-asserted analysis-ready state,
        provenance, or scientific claims. Pass
        ``trusted_construction_assertions`` when sequence, identity,
        quantitative meaning, and reference context are asserted by the caller.
        Without supplied assertion metadata, the dataset records explicit
        missing assertion provenance and downstream workflows can surface a
        stronger trusted-construction caveat. Without supplied run provenance,
        the dataset receives the same direct trusted-construction marker used
        by the compatibility constructor.
        """

        return cls(
            phospho=phospho,
            site_metadata=site_metadata,
            intensity_scale_state=intensity_scale_state,
            processing_state=processing_state,
            sample_metadata=sample_metadata,
            total=total,
            comparisons=comparisons,
            imputation_observation_mask=imputation_observation_mask,
            organism=organism,
            preprocessing_report=preprocessing_report,
            protein_aware_preparation=protein_aware_preparation,
            provenance=provenance,
            trusted_construction_assertions=trusted_construction_assertions,
            allow_opaque_site_values=allow_opaque_site_values,
        )

    @classmethod
    def _from_owned(
        cls,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        imputation_observation_mask: pd.DataFrame | None = None,
        organism: Organism | None = None,
        preprocessing_report: DatasetPreprocessingReport | None = None,
        protein_aware_preparation: ProteinAwarePreparationResult | None = None,
        provenance: RunProvenance | None = None,
        allow_opaque_site_values: bool = False,
    ) -> AnalysisReadyPhosphoDataset:
        return cls(
            phospho=phospho,
            site_metadata=site_metadata,
            intensity_scale_state=intensity_scale_state,
            processing_state=processing_state,
            sample_metadata=sample_metadata,
            total=total,
            comparisons=comparisons,
            imputation_observation_mask=imputation_observation_mask,
            organism=organism,
            preprocessing_report=preprocessing_report,
            protein_aware_preparation=protein_aware_preparation,
            provenance=provenance,
            allow_opaque_site_values=allow_opaque_site_values,
            _assume_owned=True,
        )

    def to_dataframe(self) -> pd.DataFrame:
        """Return a phospho snapshot; mutating it does not mutate this dataset."""

        return export_dataframe(self._phospho)

    def site_metadata_dataframe(self) -> pd.DataFrame:
        """Return a site-metadata snapshot isolated from this dataset."""

        return export_dataframe(self._site_metadata)

    def sample_metadata_dataframe(self) -> pd.DataFrame | None:
        """Return an optional sample-metadata snapshot isolated from this dataset."""

        return export_optional_dataframe(self._sample_metadata)

    def total_dataframe(self) -> pd.DataFrame | None:
        """Return an optional total-protein snapshot isolated from this dataset."""

        return export_optional_dataframe(self._total)

    def comparisons_dataframe(self) -> pd.DataFrame | None:
        """Return an optional comparisons snapshot isolated from this dataset."""

        return export_optional_dataframe(self._comparisons)

    def imputation_feature_metadata_dataframe(self) -> pd.DataFrame | None:
        """Return per-feature imputation metadata isolated from this dataset."""

        metadata = self._imputation_observation_metadata
        if metadata is None:
            return None
        return metadata.feature_summary_dataframe()

    def imputation_observation_summary_dataframe(
        self,
        *,
        feature_ids: Sequence[object],
        sample_ids: Sequence[object],
    ) -> pd.DataFrame | None:
        """Return imputation observation counts for requested features/samples."""

        metadata = self._imputation_observation_metadata
        if metadata is None:
            if bool(self.processing_state.missing_data.imputed):
                raise DatasetValidationError(
                    "dataset.imputation_observation_summary requires "
                    "dataset.imputation_observation_mask because "
                    "dataset.processing_state.missing_data.imputed is True"
                )
            return None
        return metadata.feature_observation_summary_dataframe(
            feature_ids=feature_ids,
            sample_ids=sample_ids,
        )

    def _aggregated_imputation_observation_mask_frame(
        self,
        *,
        sample_groups: Sequence[tuple[object, Sequence[object]]],
    ) -> pd.DataFrame | None:
        """Package-private aggregated observation mask for dataset rebuilding."""

        metadata = self._imputation_observation_metadata
        if metadata is None:
            if bool(self.processing_state.missing_data.imputed):
                raise DatasetValidationError(
                    "dataset.imputation_observation_mask aggregation requires "
                    "dataset.imputation_observation_mask because "
                    "dataset.processing_state.missing_data.imputed is True"
                )
            return None
        return metadata.aggregated_observed_mask_dataframe(sample_groups=sample_groups)

    def imputation_observed_mask_dataframe(self) -> pd.DataFrame | None:
        """Return an optional observed-cell mask snapshot."""

        metadata = self._imputation_observation_metadata
        if metadata is None:
            return None
        return metadata.observed_mask_dataframe()
