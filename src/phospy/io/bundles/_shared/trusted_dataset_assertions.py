"""Trusted construction assertions for bundle dataset reconstruction."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from phospy.provenance.models import (
    JsonValue,
    RunProvenance,
    TrustedDatasetConstructionAssertions,
    TrustedDatasetConstructionEvidence,
)
from phospy.science.datasets.processing_state import DatasetProcessingState
from phospy.science.transformations.models import IntensityScaleState
from phospy.science.transformations.state_coherence import observe_numeric_domain


def build_bundle_reconstruction_assertions(
    *,
    bundle_kind: str,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    intensity_scale_state: IntensityScaleState,
    processing_state: DatasetProcessingState,
    sample_metadata: pd.DataFrame | None = None,
    total: pd.DataFrame | None = None,
    provenance: RunProvenance | None = None,
) -> TrustedDatasetConstructionAssertions:
    """Return complete assertions for reconstructing a saved bundle dataset."""

    return TrustedDatasetConstructionAssertions(
        identity=TrustedDatasetConstructionEvidence.evidence(
            source=f"{bundle_kind} bundle dataset tables",
            policy="require_site_key_identity_columns",
            details={
                "phospho_index_name": _optional_text(phospho.index.name),
                "site_metadata_index_name": _optional_text(site_metadata.index.name),
                "site_key_column_present": "site_key" in site_metadata.columns,
                "display_id_column_present": "display_id" in site_metadata.columns,
            },
        ),
        intensity_scale=TrustedDatasetConstructionEvidence.evidence(
            source=f"{bundle_kind} bundle manifest intensity_scale_state",
            policy="require_established_intensity_scale_state",
            details={
                "label": str(intensity_scale_state.label),
                "quantity": _enum_value(intensity_scale_state.quantity),
                "processing_state_matches": (
                    processing_state.intensity_scale == intensity_scale_state
                ),
            },
        ),
        quantitative_meaning=TrustedDatasetConstructionEvidence.evidence(
            source=f"{bundle_kind} bundle manifest intensity_scale_state",
            policy="require_serialized_quantitative_meaning",
            details={
                "quantity": _enum_value(intensity_scale_state.quantity),
                "has_quantitative_meaning_provenance": (
                    intensity_scale_state.quantitative_meaning_provenance is not None
                ),
            },
        ),
        aligned_structure=TrustedDatasetConstructionEvidence.evidence(
            source=f"{bundle_kind} bundle dataset tables",
            policy="analysis_ready_private_initializer_alignment_checks",
            details=_alignment_details(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                total=total,
            ),
        ),
        localisation=_localisation_assertion(
            bundle_kind=bundle_kind,
            site_metadata=site_metadata,
        ),
        sequence=TrustedDatasetConstructionEvidence.evidence(
            source=f"{bundle_kind} bundle dataset site_metadata.site_sequence",
            policy="require_site_sequence_column",
            details={
                "site_sequence_column_present": "site_sequence"
                in site_metadata.columns,
                "row_count": int(site_metadata.shape[0]),
            },
        ),
        reference_context=_reference_context_assertion(
            bundle_kind=bundle_kind,
            provenance=provenance,
        ),
        numeric_semantic_domain=TrustedDatasetConstructionEvidence.evidence(
            source=f"{bundle_kind} bundle dataset phospho matrix",
            policy="analysis_ready_numeric_semantic_domain_preserved",
            details=_numeric_domain_details(phospho),
        ),
        asserted_by="phospy.io.bundles",
        assertion_source=f"{bundle_kind} bundle reconstruction",
    )


def _alignment_details(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    sample_metadata: pd.DataFrame | None,
    total: pd.DataFrame | None,
) -> Mapping[str, JsonValue]:
    details: dict[str, JsonValue] = {
        "phospho_rows": int(phospho.shape[0]),
        "phospho_columns": int(phospho.shape[1]),
        "site_metadata_rows": int(site_metadata.shape[0]),
        "site_metadata_index_matches_phospho": bool(
            site_metadata.index.equals(phospho.index)
        ),
        "sample_metadata_present": sample_metadata is not None,
        "sample_metadata_index_matches_phospho_columns": (
            sample_metadata is not None
            and sample_metadata.index.equals(phospho.columns)
        ),
        "total_present": total is not None,
    }
    if total is not None:
        details["total_columns_match_phospho_columns"] = bool(
            total.columns.equals(phospho.columns)
        )
    return details


def _localisation_assertion(
    *,
    bundle_kind: str,
    site_metadata: pd.DataFrame,
) -> TrustedDatasetConstructionEvidence:
    for column_name in ("localisation_confidence", "localisation_probability"):
        if column_name in site_metadata.columns:
            return TrustedDatasetConstructionEvidence.evidence(
                source=f"{bundle_kind} bundle dataset site_metadata.{column_name}",
                policy="recorded_localisation_column",
                threshold=0.0,
                details={
                    "column": column_name,
                    "non_missing_count": int(site_metadata.loc[:, column_name].count()),
                },
            )
    return TrustedDatasetConstructionEvidence.waiver(
        reason=(
            "bundle dataset site_metadata does not include localisation confidence "
            "metadata; reconstruction preserves the saved analysis-ready state"
        ),
        policy="bundle_reconstruction_preserves_saved_analysis_ready_state",
    )


def _reference_context_assertion(
    *,
    bundle_kind: str,
    provenance: RunProvenance | None,
) -> TrustedDatasetConstructionEvidence:
    if provenance is not None and provenance.reference_context is not None:
        return TrustedDatasetConstructionEvidence.evidence(
            source=f"{bundle_kind} bundle provenance.reference_context",
            policy="restore_serialized_reference_context",
        )
    return TrustedDatasetConstructionEvidence.waiver(
        reason=(
            "bundle provenance does not record dataset reference_context; "
            "reconstruction preserves the saved analysis-ready tables"
        ),
        policy="bundle_reconstruction_reference_context_unavailable",
    )


def _numeric_domain_details(phospho: pd.DataFrame) -> Mapping[str, JsonValue]:
    observation = observe_numeric_domain(phospho, table_name="dataset.phospho")
    return {
        "table": observation.table_name,
        "observed_numeric_domain": observation.observed_domain.value,
        "value_count": observation.value_count,
        "negative_count": observation.negative_count,
        "zero_count": observation.zero_count,
        "positive_count": observation.positive_count,
        "min": observation.minimum,
        "max": observation.maximum,
    }


def _enum_value(value: object) -> str | None:
    if value is None:
        return None
    return str(getattr(value, "value", value))


def _optional_text(value: object) -> str | None:
    if value is None:
        return None
    return str(value)


__all__ = ["build_bundle_reconstruction_assertions"]
