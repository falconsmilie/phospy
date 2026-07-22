"""Explicit test-only AnalysisReadyPhosphoDataset construction lanes."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset
from phospy.api import DatasetBuildRequest
from phospy.contracts.configs import DatasetPreprocessingConfig
from phospy.provenance.models import (
    JsonValue,
    RunProvenance,
    TrustedDatasetConstructionAssertions,
    TrustedDatasetConstructionEvidence,
)
from phospy.science.datasets.models import DatasetPreprocessingReport
from phospy.science.datasets.preprocessing.protein_aware_preparation import (
    ProteinAwarePreparationResult,
)
from phospy.science.datasets.processing_state import DatasetProcessingState
from phospy.science.references.models import Organism
from phospy.science.transformations.models import (
    IntensityScaleKind,
    IntensityScaleState,
    QuantitativeMeaning,
)


def complete_trusted_dataset_construction_assertions_for_tests(
    *,
    phospho: pd.DataFrame | None = None,
    site_metadata: pd.DataFrame | None = None,
    sample_metadata: pd.DataFrame | None = None,
    total: pd.DataFrame | None = None,
    intensity_scale_state: IntensityScaleState | None = None,
    processing_state: DatasetProcessingState | None = None,
    asserted_by: str = "tests.support.analysis_ready_dataset_factories",
) -> TrustedDatasetConstructionAssertions:
    """Return deterministic complete trusted assertions for table-focused tests."""

    return TrustedDatasetConstructionAssertions(
        identity=TrustedDatasetConstructionEvidence.evidence(
            source="test fixture protein-scoped site_key tables",
            policy="require_site_key_identity_columns",
            details={
                "phospho_index_name": _index_name(phospho),
                "site_metadata_index_name": _index_name(site_metadata),
                "site_key_column_present": _has_column(site_metadata, "site_key"),
                "display_id_column_present": _has_column(site_metadata, "display_id"),
            },
        ),
        intensity_scale=TrustedDatasetConstructionEvidence.evidence(
            source="test fixture established IntensityScaleState",
            policy="require_established_intensity_scale_state",
            details={
                "label": None
                if intensity_scale_state is None
                else str(intensity_scale_state.label),
                "processing_state_matches": (
                    processing_state is not None
                    and intensity_scale_state is not None
                    and processing_state.intensity_scale == intensity_scale_state
                ),
            },
        ),
        quantitative_meaning=TrustedDatasetConstructionEvidence.evidence(
            source="test fixture analysis-ready quantitative matrix",
            policy="trusted_test_fixture_quantitative_meaning",
            details={
                "quantity": None
                if intensity_scale_state is None
                else _enum_value(intensity_scale_state.quantity)
            },
        ),
        aligned_structure=TrustedDatasetConstructionEvidence.evidence(
            source="test fixture aligned DataFrames",
            policy="analysis_ready_private_initializer_alignment_checks",
            details=_alignment_details(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                total=total,
            ),
        ),
        localisation=_localisation_assertion(site_metadata),
        sequence=TrustedDatasetConstructionEvidence.evidence(
            source="test fixture site_metadata.site_sequence",
            policy="require_site_sequence_column",
            details={
                "site_sequence_column_present": _has_column(
                    site_metadata,
                    "site_sequence",
                )
            },
        ),
        reference_context=TrustedDatasetConstructionEvidence.waiver(
            reason="test fixture does not model external reference-context metadata",
            policy="unit_test_reference_context_not_under_test",
        ),
        asserted_by=asserted_by,
        assertion_source="test-only trusted table fixture",
    )


def trusted_analysis_ready_dataset_from_tables(
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
    trusted_construction_assertions: TrustedDatasetConstructionAssertions | None = None,
) -> AnalysisReadyPhosphoDataset:
    """Construct through the public trusted-table lane for downstream tests."""

    assertions = (
        trusted_construction_assertions
        if trusted_construction_assertions is not None
        else complete_trusted_dataset_construction_assertions_for_tests(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            total=total,
            intensity_scale_state=intensity_scale_state,
            processing_state=processing_state,
        )
    )
    return AnalysisReadyPhosphoDataset.from_trusted_tables(
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
        trusted_construction_assertions=assertions,
    )


def builder_backed_analysis_ready_dataset_from_tables(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    sample_metadata: pd.DataFrame | None = None,
    total: pd.DataFrame | None = None,
    organism: Organism | None = None,
    preprocessing_config: DatasetPreprocessingConfig | None = None,
    input_intensity_scale: IntensityScaleKind | str | None = "linear",
    quantitative_meaning: QuantitativeMeaning | str | None = None,
    allow_opaque_site_values: bool = False,
) -> AnalysisReadyPhosphoDataset:
    """Construct through AnalysisReadyDatasetBuilder for boundary-focused tests."""

    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            total=total,
            organism=organism,
            preprocessing_config=(
                DatasetPreprocessingConfig()
                if preprocessing_config is None
                else preprocessing_config
            ),
            input_intensity_scale=input_intensity_scale,
            quantitative_meaning=quantitative_meaning,
            allow_opaque_site_values=allow_opaque_site_values,
        )
    )


def _alignment_details(
    *,
    phospho: pd.DataFrame | None,
    site_metadata: pd.DataFrame | None,
    sample_metadata: pd.DataFrame | None,
    total: pd.DataFrame | None,
) -> Mapping[str, JsonValue]:
    details: dict[str, JsonValue] = {
        "phospho_present": phospho is not None,
        "site_metadata_present": site_metadata is not None,
        "sample_metadata_present": sample_metadata is not None,
        "total_present": total is not None,
    }
    if isinstance(phospho, pd.DataFrame):
        details["phospho_rows"] = int(phospho.shape[0])
        details["phospho_columns"] = int(phospho.shape[1])
    if isinstance(phospho, pd.DataFrame) and isinstance(site_metadata, pd.DataFrame):
        details["site_metadata_index_matches_phospho"] = bool(
            site_metadata.index.equals(phospho.index)
        )
    if isinstance(phospho, pd.DataFrame) and isinstance(sample_metadata, pd.DataFrame):
        details["sample_metadata_index_matches_phospho_columns"] = bool(
            sample_metadata.index.equals(phospho.columns)
        )
    if isinstance(phospho, pd.DataFrame) and isinstance(total, pd.DataFrame):
        details["total_columns_match_phospho_columns"] = bool(
            total.columns.equals(phospho.columns)
        )
    return details


def _localisation_assertion(
    site_metadata: pd.DataFrame | None,
) -> TrustedDatasetConstructionEvidence:
    if isinstance(site_metadata, pd.DataFrame):
        for column_name in ("localisation_confidence", "localisation_probability"):
            if column_name in site_metadata.columns:
                return TrustedDatasetConstructionEvidence.evidence(
                    source=f"test fixture site_metadata.{column_name}",
                    policy="recorded_localisation_column",
                    threshold=0.0,
                    details={
                        "column": column_name,
                        "non_missing_count": int(
                            site_metadata.loc[:, column_name].count()
                        ),
                    },
                )
    return TrustedDatasetConstructionEvidence.waiver(
        reason="test fixture does not model localisation-confidence evidence",
        policy="unit_test_localisation_not_under_test",
    )


def _has_column(frame: pd.DataFrame | None, column: str) -> bool:
    return isinstance(frame, pd.DataFrame) and column in frame.columns


def _index_name(frame: pd.DataFrame | None) -> str | None:
    if not isinstance(frame, pd.DataFrame):
        return None
    if frame.index.name is None:
        return None
    return str(frame.index.name)


def _enum_value(value: object) -> str | None:
    if value is None:
        return None
    return str(getattr(value, "value", value))


__all__ = [
    "builder_backed_analysis_ready_dataset_from_tables",
    "complete_trusted_dataset_construction_assertions_for_tests",
    "trusted_analysis_ready_dataset_from_tables",
]
