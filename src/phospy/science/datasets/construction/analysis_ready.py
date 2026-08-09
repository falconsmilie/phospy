"""Analysis-ready phosphoproteomics dataset model implementation."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import pandas as pd

from phospy.errors.validation import DatasetValidationError
from phospy.frames.comparison import dataframe_equals, optional_dataframe_equals
from phospy.frames.ownership import (
    borrow_dataframe,
    borrow_optional_dataframe,
    export_dataframe,
    export_optional_dataframe,
)
from phospy.provenance.models import (
    ReferenceContextProtocol,
    RunProvenance,
    TrustedDatasetConstructionAssertions,
)
from phospy.science.datasets.construction.service import (
    _AnalysisReadyDatasetConstructionService,
    _ValidatedAnalysisReadyTables,
)
from phospy.science.datasets.construction.validation import (
    _DIRECT_CONSTRUCTION_ERROR_MESSAGE,
)
from phospy.science.datasets.imputation_metadata import ImputationObservationMetadata
from phospy.science.datasets.internal_frame_store import DatasetInternalFrameStore
from phospy.science.datasets.preprocessing.protein_aware_models import (
    ProteinAwarePreparationResult,
)
from phospy.science.datasets.processing_state import (
    DatasetPreprocessingReport,
    DatasetProcessingState,
)
from phospy.science.references.models import Organism
from phospy.science.transformations.models import IntensityScaleState


@dataclass(frozen=True, slots=True, init=False, eq=False)
class AnalysisReadyPhosphoDataset:
    """Public analysis-ready dataset contract.

    ``AnalysisReadyPhosphoDataset`` remains a stable public result/domain type,
    but ordinary direct construction is not a supported creation path. Use
    ``AnalysisReadyDatasetBuilder.run(...)`` for ordinary dataset construction;
    it owns user-input interpretation, preprocessing, processing-state
    establishment, private validation, and construction provenance. Advanced
    callers who already own complete analysis-ready tables must use
    ``AnalysisReadyPhosphoDataset.from_trusted_tables(...)`` with complete
    ``TrustedDatasetConstructionAssertions`` and provenance-linked table
    fingerprints.

    Construction validates structural invariants, including table shape,
    alignment, analysis-ready ``site_key`` identity, processing-state
    coherence, established transformation state, and numeric-semantic
    coherence between quantitative meaning and observed numeric sign domain. It
    cannot prove the biological correctness of caller-asserted provenance or
    scientific claims.
    The ``from_trusted_tables(...)`` lane requires typed evidence or an explicit
    waiver for identity, intensity scale, quantitative meaning, aligned table
    structure, localisation, sequence, and reference context. A
    numeric-semantic-domain conflict can only be bypassed by the optional typed
    ``numeric_semantic_domain`` waiver, which remains visible in provenance.
    Any supplied provenance must fingerprint the actual represented tables.

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

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for an explicit comparison of owned
    analysis-ready table content and declared scientific state.
    """

    __hash__ = object.__hash__

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
    _internal_frame_store: DatasetInternalFrameStore = field(init=False, repr=False)
    _allow_opaque_site_values: bool = field(init=False, repr=False, default=False)

    def __init__(self, *args: object, **kwargs: object) -> None:
        del args, kwargs
        raise TypeError(_DIRECT_CONSTRUCTION_ERROR_MESSAGE)

    def _init_from_validated_tables(
        self,
        tables: _ValidatedAnalysisReadyTables,
    ) -> None:
        if not isinstance(tables, _ValidatedAnalysisReadyTables):
            raise TypeError(
                "AnalysisReadyPhosphoDataset private initialization requires "
                "_ValidatedAnalysisReadyTables"
            )
        object.__setattr__(self, "intensity_scale_state", tables.intensity_scale_state)
        object.__setattr__(self, "processing_state", tables.processing_state)
        object.__setattr__(self, "organism", tables.organism)
        object.__setattr__(self, "preprocessing_report", tables.preprocessing_report)
        object.__setattr__(
            self,
            "trusted_construction_assertions",
            tables.trusted_construction_assertions,
        )
        object.__setattr__(
            self,
            "protein_aware_preparation",
            tables.protein_aware_preparation,
        )
        object.__setattr__(self, "provenance", tables.provenance)
        object.__setattr__(
            self,
            "allow_opaque_site_values",
            tables.allow_opaque_site_values,
        )
        object.__setattr__(self, "_phospho", tables.phospho)
        object.__setattr__(self, "_site_metadata", tables.site_metadata)
        object.__setattr__(self, "_sample_metadata", tables.sample_metadata)
        object.__setattr__(self, "_total", tables.total)
        object.__setattr__(self, "_comparisons", tables.comparisons)
        object.__setattr__(
            self,
            "_imputation_observation_metadata",
            tables.imputation_observation_metadata,
        )
        object.__setattr__(
            self,
            "_internal_frame_store",
            tables.internal_frame_store,
        )
        object.__setattr__(
            self,
            "_allow_opaque_site_values",
            tables.allow_opaque_site_values,
        )

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
    def reference_context(self) -> ReferenceContextProtocol | None:
        if self.provenance is None:
            return None
        return self.provenance.reference_context

    def _internal_frame_store_for_current_frames(self) -> DatasetInternalFrameStore:
        """Return the dataset-owned internal store for current private frames.

        Normal construction installs the store once. This freshness check keeps
        package-private validator tests honest when they deliberately replace
        private frame attributes to model malformed post-boundary objects.
        """

        store = self._internal_frame_store
        if store.is_current_for(
            phospho=self._phospho,
            site_metadata=self._site_metadata,
            sample_metadata=self._sample_metadata,
            total=self._total,
            comparisons=self._comparisons,
        ):
            return store
        store = DatasetInternalFrameStore.from_frames(
            phospho=self._phospho,
            site_metadata=self._site_metadata,
            sample_metadata=self._sample_metadata,
            total=self._total,
            comparisons=self._comparisons,
        )
        object.__setattr__(self, "_internal_frame_store", store)
        return store

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
    def _construct_from_validated_tables(
        cls,
        tables: _ValidatedAnalysisReadyTables,
    ) -> AnalysisReadyPhosphoDataset:
        dataset = object.__new__(cls)
        AnalysisReadyPhosphoDataset._init_from_validated_tables(dataset, tables)
        return dataset

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
        coherent processing state. It uses the private construction service and
        enforces the same structural invariants as the builder-owned path,
        including table shape, alignment, site identity, ``site_sequence``
        validation, state coherence, and numeric-semantic coherence between
        established quantitative meaning and observed numeric sign domain.

        Validation can confirm structural consistency, but it cannot prove the
        biological correctness of caller-asserted analysis-ready state,
        provenance, or scientific claims. The primary advanced lane requires
        ``trusted_construction_assertions`` with typed evidence or an explicit
        waiver for identity, intensity scale, quantitative meaning, aligned
        structure, localisation, sequence, and reference context. Localisation
        evidence must record source, policy, and threshold; otherwise callers
        must record an explicit waiver. A numeric-semantic-domain conflict
        requires an explicit typed ``numeric_semantic_domain`` waiver; generic
        trusted construction or quantitative-meaning waivers do not bypass it.
        Supplied run provenance must fingerprint the exact analysis-ready
        tables; false table fingerprints are rejected.
        Without supplied run provenance, the dataset receives a trusted-table
        reconstruction provenance marker with the assertion fingerprint linked
        into provenance.
        """

        if trusted_construction_assertions is None:
            raise DatasetValidationError(
                "AnalysisReadyPhosphoDataset.from_trusted_tables requires "
                "trusted_construction_assertions with typed evidence or an "
                "explicit waiver for identity, intensity scale, quantitative "
                "meaning, aligned structure, localisation, sequence, and "
                "reference context"
            )
        tables = _AnalysisReadyDatasetConstructionService().from_trusted_tables(
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
        return cls._construct_from_validated_tables(tables)

    @classmethod
    def _from_builder_output(
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
        provenance: RunProvenance,
        allow_opaque_site_values: bool = False,
    ) -> AnalysisReadyPhosphoDataset:
        tables = _AnalysisReadyDatasetConstructionService().from_builder_output(
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
        )
        return cls._construct_from_validated_tables(tables)

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

    def scientifically_equals(
        self,
        other: object,
        *,
        include_provenance: bool = True,
    ) -> bool:
        """Return ``True`` when another dataset has the same scientific content.

        The comparison covers the owned analysis-ready tables, imputation
        observation metadata, declared analysis state, organism, opaque-site
        policy, and trusted construction assertions. Run provenance is included
        by default because it defines the audit identity of the represented
        tables; pass ``include_provenance=False`` to compare only table/state
        content.

        Preprocessing and protein-aware report sidecars are not part of this
        dataset-level content contract; compare those sidecars with their own
        named methods when that report content matters.
        """

        if not isinstance(other, AnalysisReadyPhosphoDataset):
            return False
        same_content = (
            dataframe_equals(self._phospho, other._phospho)
            and dataframe_equals(self._site_metadata, other._site_metadata)
            and optional_dataframe_equals(self._sample_metadata, other._sample_metadata)
            and optional_dataframe_equals(self._total, other._total)
            and optional_dataframe_equals(self._comparisons, other._comparisons)
            and _optional_imputation_observation_metadata_equals(
                self._imputation_observation_metadata,
                other._imputation_observation_metadata,
            )
            and self.intensity_scale_state == other.intensity_scale_state
            and self.processing_state == other.processing_state
            and self.organism == other.organism
            and self.trusted_construction_assertions
            == other.trusted_construction_assertions
            and self.allow_opaque_site_values == other.allow_opaque_site_values
            and self._allow_opaque_site_values == other._allow_opaque_site_values
        )
        if not same_content:
            return False
        if include_provenance and self.provenance != other.provenance:
            return False
        return True


def _optional_imputation_observation_metadata_equals(
    left: ImputationObservationMetadata | None,
    right: ImputationObservationMetadata | None,
) -> bool:
    if left is None or right is None:
        return left is right
    return left.scientifically_equals(right)
