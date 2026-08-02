"""Private service for validating analysis-ready dataset construction tables."""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from phospy.errors.validation import DatasetValidationError
from phospy.provenance.models import (
    RunProvenance,
    TrustedDatasetConstructionAssertions,
)
from phospy.science.datasets.construction.fingerprints import (
    _fingerprints_for_analysis_ready_tables,
    _require_trusted_provenance_table_fingerprints,
)
from phospy.science.datasets.construction.trusted_assertions import (
    _require_complete_trusted_assertion_metadata,
    _require_complete_trusted_assertions,
    _resolve_trusted_construction_assertions,
)
from phospy.science.datasets.construction.validation import (
    _INTENSITY_SCALE_STATE_VALIDATOR,
    _QUANTITATIVE_NUMERIC_DOMAIN_VALIDATOR,
    _own_dataset_frames,
    _require_builder_output_provenance,
    _validate_optional_comparisons,
    analysis_ready_matrix_missing_value_count,
)
from phospy.science.datasets.direct_construction import (
    build_direct_construction_provenance,
)
from phospy.science.datasets.imputation_metadata import ImputationObservationMetadata
from phospy.science.datasets.imputation_metadata import (
    build_imputation_observation_metadata_or_none as _build_imputation_observation_metadata_or_none,
)
from phospy.science.datasets.organism_coherence import (
    normalize_dataset_organism_state,
    require_dataset_provenance_organism_coherence,
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
    missing_data_state_claims_no_missing_values as _missing_data_state_claims_no_missing_values,
)
from phospy.science.datasets.processing_state import (
    require_instance as _require_instance,
)
from phospy.science.datasets.processing_state import (
    require_optional_instance as _require_optional_instance,
)
from phospy.science.references.models import Organism
from phospy.science.tables.datasets import (
    PhosphoIntensityMatrix,
    SampleMetadataTable,
    SiteMetadataTable,
    TotalProteinMatrix,
)
from phospy.science.transformations.models import IntensityScaleState

_VALIDATED_TABLES_AUTHORITY = object()
_VALIDATED_TABLES_DIRECT_CONSTRUCTION_ERROR_MESSAGE = (
    "_ValidatedAnalysisReadyTables must be produced by the private analysis-ready "
    "dataset construction service"
)


@dataclass(frozen=True, slots=True, init=False)
class _ValidatedAnalysisReadyTables:
    """Owned, validated tables ready to be installed into the dataset object."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    intensity_scale_state: IntensityScaleState
    processing_state: DatasetProcessingState
    sample_metadata: pd.DataFrame | None = None
    total: pd.DataFrame | None = None
    comparisons: pd.DataFrame | None = None
    imputation_observation_metadata: ImputationObservationMetadata | None = field(
        default=None,
        repr=False,
    )
    organism: Organism | None = None
    preprocessing_report: DatasetPreprocessingReport | None = None
    protein_aware_preparation: ProteinAwarePreparationResult | None = None
    provenance: RunProvenance | None = None
    trusted_construction_assertions: TrustedDatasetConstructionAssertions | None = None
    allow_opaque_site_values: bool = False

    def __init__(self, *args: object, **kwargs: object) -> None:
        del args, kwargs
        raise TypeError(_VALIDATED_TABLES_DIRECT_CONSTRUCTION_ERROR_MESSAGE)

    @classmethod
    def _from_validated_parts(
        cls,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        sample_metadata: pd.DataFrame | None,
        total: pd.DataFrame | None,
        comparisons: pd.DataFrame | None,
        imputation_observation_metadata: ImputationObservationMetadata | None,
        organism: Organism | None,
        preprocessing_report: DatasetPreprocessingReport | None,
        protein_aware_preparation: ProteinAwarePreparationResult | None,
        provenance: RunProvenance,
        trusted_construction_assertions: TrustedDatasetConstructionAssertions | None,
        allow_opaque_site_values: bool,
        _authority: object,
    ) -> _ValidatedAnalysisReadyTables:
        if _authority is not _VALIDATED_TABLES_AUTHORITY:
            raise TypeError(_VALIDATED_TABLES_DIRECT_CONSTRUCTION_ERROR_MESSAGE)
        tables = object.__new__(cls)
        object.__setattr__(tables, "phospho", phospho)
        object.__setattr__(tables, "site_metadata", site_metadata)
        object.__setattr__(tables, "intensity_scale_state", intensity_scale_state)
        object.__setattr__(tables, "processing_state", processing_state)
        object.__setattr__(tables, "sample_metadata", sample_metadata)
        object.__setattr__(tables, "total", total)
        object.__setattr__(tables, "comparisons", comparisons)
        object.__setattr__(
            tables,
            "imputation_observation_metadata",
            imputation_observation_metadata,
        )
        object.__setattr__(tables, "organism", organism)
        object.__setattr__(tables, "preprocessing_report", preprocessing_report)
        object.__setattr__(
            tables,
            "protein_aware_preparation",
            protein_aware_preparation,
        )
        object.__setattr__(tables, "provenance", provenance)
        object.__setattr__(
            tables,
            "trusted_construction_assertions",
            trusted_construction_assertions,
        )
        object.__setattr__(
            tables,
            "allow_opaque_site_values",
            allow_opaque_site_values,
        )
        return tables


class _AnalysisReadyDatasetConstructionService:
    """Resolve construction evidence, validate tables, and return an aggregate."""

    def from_builder_output(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        provenance: RunProvenance,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        imputation_observation_mask: pd.DataFrame | None = None,
        organism: Organism | None = None,
        preprocessing_report: DatasetPreprocessingReport | None = None,
        protein_aware_preparation: ProteinAwarePreparationResult | None = None,
        allow_opaque_site_values: bool = False,
        assume_owned: bool = True,
    ) -> _ValidatedAnalysisReadyTables:
        _require_builder_output_provenance(provenance)
        return self._validate_analysis_ready_tables(
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
            trusted_construction_assertions=None,
            allow_opaque_site_values=allow_opaque_site_values,
            assume_owned=assume_owned,
        )

    def from_trusted_tables(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        trusted_construction_assertions: TrustedDatasetConstructionAssertions,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        imputation_observation_mask: pd.DataFrame | None = None,
        organism: Organism | None = None,
        preprocessing_report: DatasetPreprocessingReport | None = None,
        protein_aware_preparation: ProteinAwarePreparationResult | None = None,
        provenance: RunProvenance | None = None,
        allow_opaque_site_values: bool = False,
    ) -> _ValidatedAnalysisReadyTables:
        _require_instance(
            trusted_construction_assertions,
            expected_type=TrustedDatasetConstructionAssertions,
            error_message=(
                "dataset.trusted_construction_assertions must be "
                "TrustedDatasetConstructionAssertions"
            ),
        )
        resolved_trusted_construction_assertions = (
            _resolve_trusted_construction_assertions(
                trusted_construction_assertions=trusted_construction_assertions,
                provenance=provenance,
                assume_owned=False,
            )
        )
        _require_complete_trusted_assertion_metadata(
            assertions=resolved_trusted_construction_assertions,
        )
        tables = self._validate_analysis_ready_tables(
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
            trusted_construction_assertions=(resolved_trusted_construction_assertions),
            allow_opaque_site_values=allow_opaque_site_values,
            assume_owned=False,
        )
        _require_complete_trusted_assertions(
            assertions=tables.trusted_construction_assertions,
            provenance=tables.provenance,
        )
        return tables

    def from_provenanced_tables(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        provenance: RunProvenance,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        imputation_observation_mask: pd.DataFrame | None = None,
        organism: Organism | None = None,
        preprocessing_report: DatasetPreprocessingReport | None = None,
        protein_aware_preparation: ProteinAwarePreparationResult | None = None,
        allow_opaque_site_values: bool = False,
        assume_owned: bool = False,
    ) -> _ValidatedAnalysisReadyTables:
        _require_instance(
            provenance,
            expected_type=RunProvenance,
            error_message="analysis-ready dataset construction requires RunProvenance",
        )
        return self._validate_analysis_ready_tables(
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
            trusted_construction_assertions=None,
            allow_opaque_site_values=allow_opaque_site_values,
            assume_owned=assume_owned,
        )

    def _validate_analysis_ready_tables(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        sample_metadata: pd.DataFrame | None,
        total: pd.DataFrame | None,
        comparisons: pd.DataFrame | None,
        imputation_observation_mask: pd.DataFrame | None,
        organism: Organism | None,
        preprocessing_report: DatasetPreprocessingReport | None,
        protein_aware_preparation: ProteinAwarePreparationResult | None,
        provenance: RunProvenance | None,
        trusted_construction_assertions: TrustedDatasetConstructionAssertions | None,
        allow_opaque_site_values: bool,
        assume_owned: bool,
    ) -> _ValidatedAnalysisReadyTables:
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
            assume_owned=assume_owned,
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
        raw_missing_value_count = analysis_ready_matrix_missing_value_count(phospho)
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
        _QUANTITATIVE_NUMERIC_DOMAIN_VALIDATOR.run(
            phospho=phospho_table.frame,
            total=None if total_table is None else total_table.frame,
            intensity_scale_state=validated_intensity_scale_state,
            trusted_construction_assertions=trusted_construction_assertions,
        )
        if provenance is None and trusted_construction_assertions is None:
            raise DatasetValidationError(
                "AnalysisReadyPhosphoDataset private construction requires "
                "builder-created provenance or complete "
                "TrustedDatasetConstructionAssertions supplied through "
                "AnalysisReadyPhosphoDataset.from_trusted_tables(...)"
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
            if trusted_construction_assertions is None:
                raise DatasetValidationError(
                    "AnalysisReadyPhosphoDataset.from_trusted_tables requires "
                    "complete TrustedDatasetConstructionAssertions before "
                    "trusted-table reconstruction provenance can be created"
                )
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
                trusted_construction_assertions=trusted_construction_assertions,
            )
        require_dataset_provenance_organism_coherence(
            organism=resolved_dataset_organism,
            provenance=provenance,
            error_type=DatasetValidationError,
        )
        actual_fingerprints = _fingerprints_for_analysis_ready_tables(
            phospho=phospho_table.frame,
            site_metadata=site_metadata_table.frame,
            sample_metadata=(
                None if sample_metadata_table is None else sample_metadata_table.frame
            ),
            total=None if total_table is None else total_table.frame,
            comparisons=comparisons,
            imputation_observation_mask=imputation_observation_mask,
        )
        _require_trusted_provenance_table_fingerprints(
            provenance=provenance,
            actual_fingerprints=actual_fingerprints,
        )
        return _ValidatedAnalysisReadyTables._from_validated_parts(
            phospho=phospho_table.frame,
            site_metadata=site_metadata_table.frame,
            intensity_scale_state=validated_intensity_scale_state,
            processing_state=processing_state,
            sample_metadata=(
                None if sample_metadata_table is None else sample_metadata_table.frame
            ),
            total=None if total_table is None else total_table.frame,
            comparisons=comparisons,
            imputation_observation_metadata=imputation_observation_metadata,
            organism=resolved_dataset_organism,
            preprocessing_report=preprocessing_report,
            protein_aware_preparation=protein_aware_preparation,
            provenance=provenance,
            trusted_construction_assertions=trusted_construction_assertions,
            allow_opaque_site_values=allow_opaque_site_values,
            _authority=_VALIDATED_TABLES_AUTHORITY,
        )


__all__: list[str] = []
