"""Reconstruct typed kinase workflow models from decoded bundle sections."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import replace
from pathlib import Path
from typing import NoReturn, TypeGuard

import pandas as pd

from phospy.contracts.kinase_reference_projection import (
    KinaseReferenceProjectionSummary,
)
from phospy.contracts.result_caveats import result_caveats_from_payloads
from phospy.contracts.results import (
    KinaseWorkflowAttritionProvenance,
    KinaseWorkflowCaveat,
    KinaseWorkflowResult,
)
from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import PhosPyValidationError
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.io.bundles._kinase.manifest import KinaseManifestSections
from phospy.io.bundles._shared.intensity_scale_state import (
    intensity_scale_state_from_payload,
)
from phospy.io.bundles._shared.organisms import (
    parse_optional_organism,
    parse_required_organism,
)
from phospy.io.bundles._shared.primitives import (
    require_mapping,
    validate_json_safe_mapping,
)
from phospy.io.bundles._shared.processing_state import (
    processing_state_from_payload,
)
from phospy.io.bundles._shared.tables import (
    read_optional_series,
    read_optional_table,
    read_required_table,
)
from phospy.io.bundles._shared.trusted_dataset_assertions import (
    build_bundle_reconstruction_assertions,
)
from phospy.provenance.models import RunProvenance
from phospy.provenance.serialization import from_payload as provenance_from_payload
from phospy.science.activities.membership import ActivityMembershipSelection
from phospy.science.activities.models import (
    KSEA_ZSCORE_ACTIVITY_METHOD,
    ActivityMethodMetadata,
    ActivityMethodSummary,
    KinaseActivityResult,
)
from phospy.science.activities.semantics import (
    ActivityInputSemantics,
    ActivityProfileAxis,
    ActivityProfileMetadata,
    ActivityQuantitativeSemantics,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.datasets.processing_state import DatasetProcessingState
from phospy.science.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.science.references.models import ReferenceBundle
from phospy.science.transformations.models import IntensityScaleState

_LEGACY_KINASE_BUNDLE_SCHEMA_ERROR = (
    "Legacy kinase bundle schemas are no longer supported. Regenerate the bundle "
    "with the current PhosPy version."
)


def reconstruct_kinase_result(
    *,
    bundle_root: Path,
    sections: KinaseManifestSections,
) -> KinaseWorkflowResult:
    """Rebuild a KinaseWorkflowResult from already-validated manifest sections."""

    provenance = _parse_bundle_provenance(sections.provenance_payload)
    provenance = _validate_kinase_reference_projection_provenance(provenance)
    processing_state_payload = require_mapping(
        sections.dataset_metadata.get("processing_state"),
        field_name="bundle manifest.dataset.metadata.processing_state",
    )
    processing_state = _parse_bundle_processing_state(processing_state_payload)
    intensity_scale_payload = require_mapping(
        sections.dataset_metadata.get("intensity_scale_state"),
        field_name="bundle manifest.dataset.metadata.intensity_scale_state",
    )
    site_metadata = read_required_table(
        bundle_root=bundle_root,
        tables=sections.dataset_tables,
        table_key="site_metadata",
        field_name="bundle manifest.dataset.tables.site_metadata",
    )
    site_metadata = _normalise_site_metadata_bundle_table(site_metadata)
    phospho = read_required_table(
        bundle_root=bundle_root,
        tables=sections.dataset_tables,
        table_key="phospho",
        field_name="bundle manifest.dataset.tables.phospho",
    )
    sample_metadata = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.dataset_tables,
        table_key="sample_metadata",
        field_name="bundle manifest.dataset.tables.sample_metadata",
    )
    total = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.dataset_tables,
        table_key="total",
        field_name="bundle manifest.dataset.tables.total",
    )
    intensity_scale_state = _parse_bundle_intensity_scale_state(intensity_scale_payload)
    dataset = AnalysisReadyPhosphoDataset.from_trusted_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        total=total,
        organism=parse_optional_organism(
            sections.dataset_metadata.get("organism"),
            field_name="bundle manifest.dataset.metadata.organism",
        ),
        intensity_scale_state=intensity_scale_state,
        processing_state=processing_state,
        trusted_construction_assertions=build_bundle_reconstruction_assertions(
            bundle_kind="kinase_workflow_result",
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            total=total,
            intensity_scale_state=intensity_scale_state,
            processing_state=processing_state,
            provenance=provenance,
        ),
    )

    references = ReferenceBundle.from_trusted_owned(
        organism=parse_required_organism(
            sections.references_metadata.get("organism"),
            field_name="bundle manifest.resolved_references.metadata.organism",
        ),
        kinase_substrate_map=read_required_table(
            bundle_root=bundle_root,
            tables=sections.reference_tables,
            table_key="kinase_substrate_map",
            field_name="bundle manifest.resolved_references.tables.kinase_substrate_map",
        ),
        site_sequences=read_required_table(
            bundle_root=bundle_root,
            tables=sections.reference_tables,
            table_key="site_sequences",
            field_name="bundle manifest.resolved_references.tables.site_sequences",
        ),
    )

    scoring_result = KinaseScoringResult.from_trusted_owned(
        profile_scores=read_required_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="profile_scores",
            field_name="bundle manifest.outputs.scoring.tables.profile_scores",
        ),
        motif_scores=read_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="motif_scores",
            field_name="bundle manifest.outputs.scoring.tables.motif_scores",
        ),
        rank_weighted_fusion_scores=read_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="rank_weighted_fusion_scores",
            field_name=(
                "bundle manifest.outputs.scoring.tables.rank_weighted_fusion_scores"
            ),
        ),
        kinase_library_motif_scores=_read_absent_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="kinase_library_motif_scores",
            field_name=(
                "bundle manifest.outputs.scoring.tables.kinase_library_motif_scores"
            ),
        ),
        combined_profile_motif_scores=_read_absent_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="combined_profile_motif_scores",
            field_name=(
                "bundle manifest.outputs.scoring.tables.combined_profile_motif_scores"
            ),
        ),
        score_fusion_weights=read_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="score_fusion_weights",
            field_name="bundle manifest.outputs.scoring.tables.score_fusion_weights",
        ),
        kinase_library_site_diagnostics=_read_absent_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="kinase_library_site_diagnostics",
            field_name=(
                "bundle manifest.outputs.scoring.tables.kinase_library_site_diagnostics"
            ),
        ),
        kinase_library_kinase_diagnostics=_read_absent_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="kinase_library_kinase_diagnostics",
            field_name=(
                "bundle manifest.outputs.scoring.tables."
                "kinase_library_kinase_diagnostics"
            ),
        ),
        profile_self_inclusion_policy=_profile_self_inclusion_policy_from_provenance(
            provenance
        ),
    )

    prediction_result = KinasePredictionResult.from_trusted_owned(
        pred_mat=read_required_table(
            bundle_root=bundle_root,
            tables=sections.prediction_tables,
            table_key="pred_mat",
            field_name="bundle manifest.outputs.prediction.tables.pred_mat",
        ),
        substrate_list=read_optional_table(
            bundle_root=bundle_root,
            tables=sections.prediction_tables,
            table_key="substrate_list",
            field_name="bundle manifest.outputs.prediction.tables.substrate_list",
        ),
    )
    substrate_contributions = _read_absent_optional_table(
        bundle_root=bundle_root,
        tables=sections.scoring_tables,
        table_key="substrate_contributions",
        field_name="bundle manifest.outputs.scoring.tables.substrate_contributions",
    )

    weighted_activity = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.activity_tables,
        table_key="weighted_activity",
        field_name="bundle manifest.outputs.activity.tables.weighted_activity",
    )
    thresholded_substrate_mean_activity = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.activity_tables,
        table_key="thresholded_substrate_mean_activity",
        field_name=(
            "bundle manifest.outputs.activity.tables."
            "thresholded_substrate_mean_activity"
        ),
    )
    thresholded_substrate_counts = read_optional_series(
        bundle_root=bundle_root,
        tables=sections.activity_tables,
        table_key="thresholded_substrate_counts",
        field_name="bundle manifest.outputs.activity.tables.thresholded_substrate_counts",
        series_name="n_substrates",
    )
    activity_substrate_counts = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.activity_tables,
        table_key="activity_substrate_counts",
        field_name="bundle manifest.outputs.activity.tables.activity_substrate_counts",
    )
    target_counts = read_optional_series(
        bundle_root=bundle_root,
        tables=sections.activity_tables,
        table_key="target_counts",
        field_name="bundle manifest.outputs.activity.tables.target_counts",
        series_name="n_targets",
    )
    target_table = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.activity_tables,
        table_key="target_table",
        field_name="bundle manifest.outputs.activity.tables.target_table",
    )
    statistics_table = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.activity_tables,
        table_key="statistics_table",
        field_name="bundle manifest.outputs.activity.tables.statistics_table",
    )
    statistics_table = _normalise_activity_statistics_bundle_table(statistics_table)

    if sections.activity_enabled:
        if (
            weighted_activity is None
            or thresholded_substrate_mean_activity is None
            or thresholded_substrate_counts is None
            or target_counts is None
            or target_table is None
        ):
            raise PhosPyInputError(
                "bundle manifest outputs.activity.tables are incomplete for enabled activity outputs"
            )
        if sections.activity_method_metadata is None:
            raise PhosPyInputError(
                "bundle manifest.outputs.activity.method is required when activity is enabled"
            )
        try:
            activity_method = ActivityMethodMetadata.from_payload(
                sections.activity_method_metadata
            )
        except ValueError as exc:
            raise PhosPyInputError(
                f"bundle manifest.outputs.activity.method is invalid: {exc}"
            ) from exc
        activity_method_summary = None
        if sections.activity_method_summary is not None:
            try:
                activity_method_summary = ActivityMethodSummary.from_payload(
                    sections.activity_method_summary
                )
            except (TypeError, ValueError, PhosPyInputError) as exc:
                raise PhosPyInputError(
                    "bundle manifest.outputs.activity.summary is invalid: "
                    f"{exc}; regenerate the bundle from the original "
                    "KinaseActivityResult"
                ) from exc
        input_semantics = _parse_activity_input_semantics(
            sections.activity_input_semantics
        )
        profile_metadata = _parse_activity_profile_metadata(
            sections.activity_profile_metadata
        )
        membership_selection = _parse_activity_membership_selection(
            sections.activity_membership_selection,
            activity_method=activity_method,
            weighted_activity=weighted_activity,
        )
        _validate_activity_semantic_metadata(
            input_semantics=input_semantics,
            profile_metadata=profile_metadata,
            activity_matrix=weighted_activity,
        )
        _validate_activity_provenance_agreement(
            provenance=provenance,
            input_semantics=input_semantics,
        )
        try:
            activity_result = KinaseActivityResult.from_trusted_owned(
                weighted_activity=weighted_activity,
                thresholded_substrate_mean_activity=(
                    thresholded_substrate_mean_activity
                ),
                thresholded_substrate_counts=thresholded_substrate_counts,
                activity_substrate_counts=activity_substrate_counts,
                target_counts=target_counts,
                target_table=target_table,
                statistics_table=statistics_table,
                method_summary=activity_method_summary,
                activity_method=activity_method,
                input_semantics=input_semantics,
                profile_metadata=profile_metadata,
                membership_selection=membership_selection,
            )
        except (WorkflowBoundaryError, PhosPyValidationError, ValueError) as exc:
            raise PhosPyInputError(
                "bundle manifest.outputs.activity semantic metadata is "
                "inconsistent with activity tables: "
                f"{exc}; correct the manifest or regenerate the bundle from the "
                "original KinaseActivityResult"
            ) from exc
    else:
        if (
            weighted_activity is not None
            or thresholded_substrate_mean_activity is not None
            or thresholded_substrate_counts is not None
            or activity_substrate_counts is not None
            or target_counts is not None
            or target_table is not None
            or statistics_table is not None
        ):
            raise PhosPyInputError(
                "bundle manifest outputs.activity.enabled=false must not declare populated activity tables"
            )
        if sections.activity_method_metadata is not None:
            raise PhosPyInputError(
                "bundle manifest outputs.activity.enabled=false must not declare activity method metadata"
            )
        if sections.activity_method_summary is not None:
            raise PhosPyInputError(
                "bundle manifest outputs.activity.enabled=false must not declare activity method summary metadata"
            )
        if sections.activity_input_semantics is not None:
            raise PhosPyInputError(
                "bundle manifest outputs.activity.enabled=false must not declare "
                "activity input_semantics; remove the semantic payload or "
                "regenerate the bundle"
            )
        if sections.activity_profile_metadata is not None:
            raise PhosPyInputError(
                "bundle manifest outputs.activity.enabled=false must not declare "
                "activity profile_metadata; remove the semantic payload or "
                "regenerate the bundle"
            )
        activity_result = None

    return KinaseWorkflowResult.from_trusted_owned(
        dataset=dataset,
        references=references,
        scoring_result=scoring_result,
        prediction_result=prediction_result,
        activity_result=activity_result,
        provenance=provenance,
        substrate_contributions=substrate_contributions,
        attrition_provenance=_kinase_attrition_provenance_from_provenance(provenance),
        caveats=(
            result_caveats_from_payloads(sections.caveats_payload)
            or _kinase_caveats_from_provenance(provenance)
        ),
    )


def _parse_activity_input_semantics(
    payload: Mapping[str, object] | None,
) -> ActivityInputSemantics:
    field_name = "bundle manifest.outputs.activity.input_semantics"
    if payload is None:
        raise PhosPyInputError(
            f"{field_name} is required when activity is enabled; regenerate the "
            "bundle from the original KinaseActivityResult"
        )
    try:
        return ActivityInputSemantics.from_payload(payload)
    except (TypeError, ValueError, PhosPyValidationError, WorkflowBoundaryError) as exc:
        raise PhosPyInputError(
            f"{field_name} is invalid: {exc}; correct the manifest or regenerate "
            "the bundle from the original KinaseActivityResult"
        ) from exc


def _parse_activity_profile_metadata(
    payload: Mapping[str, object] | None,
) -> ActivityProfileMetadata:
    field_name = "bundle manifest.outputs.activity.profile_metadata"
    if payload is None:
        raise PhosPyInputError(
            f"{field_name} is required when activity is enabled; regenerate the "
            "bundle from the original KinaseActivityResult"
        )
    try:
        return ActivityProfileMetadata.from_payload(payload)
    except (TypeError, ValueError, PhosPyValidationError, WorkflowBoundaryError) as exc:
        raise PhosPyInputError(
            f"{field_name} is invalid: {exc}; correct the manifest or regenerate "
            "the bundle from the original KinaseActivityResult"
        ) from exc


def _parse_activity_membership_selection(
    payload: Mapping[str, object] | None,
    *,
    activity_method: ActivityMethodMetadata,
    weighted_activity: pd.DataFrame,
) -> ActivityMembershipSelection | None:
    if payload is None:
        if (
            str(activity_method.activity_method_id)
            != KSEA_ZSCORE_ACTIVITY_METHOD.activity_method_id
        ):
            return None
        return ActivityMembershipSelection.missing(
            selected_kinase_universe=weighted_activity.index.astype(str).tolist(),
            selected_substrate_universe=(),
        )
    try:
        return ActivityMembershipSelection.from_payload(payload)
    except (TypeError, ValueError, WorkflowBoundaryError) as exc:
        raise PhosPyInputError(
            "bundle manifest.outputs.activity.membership_selection is invalid: "
            f"{exc}; correct the manifest or regenerate the bundle from the "
            "original KinaseActivityResult"
        ) from exc


def _validate_activity_semantic_metadata(
    *,
    input_semantics: ActivityInputSemantics,
    profile_metadata: ActivityProfileMetadata,
    activity_matrix: pd.DataFrame,
) -> None:
    profile_ids = tuple(str(column) for column in activity_matrix.columns)
    if profile_metadata.axis is not input_semantics.profile_axis:
        _raise_activity_semantic_manifest_error(
            "bundle manifest.outputs.activity.profile_metadata.axis",
            "must match bundle manifest.outputs.activity.input_semantics.profile_axis",
        )
    _require_exact_manifest_labels(
        observed=profile_metadata.profile_ids,
        expected=profile_ids,
        field_name="bundle manifest.outputs.activity.profile_metadata.profile_ids",
        expected_label="activity/weighted_activity table columns",
    )

    axis = _activity_profile_axis(input_semantics)
    if axis is ActivityProfileAxis.SAMPLE:
        _require_exact_manifest_labels(
            observed=profile_metadata.sample_ids,
            expected=profile_ids,
            field_name="bundle manifest.outputs.activity.profile_metadata.sample_ids",
            expected_label="activity/weighted_activity table columns",
        )
        _require_empty_manifest_labels(
            profile_metadata.condition_ids,
            field_name=(
                "bundle manifest.outputs.activity.profile_metadata.condition_ids"
            ),
            axis=axis,
        )
        _require_empty_manifest_labels(
            profile_metadata.contrast_ids,
            field_name=(
                "bundle manifest.outputs.activity.profile_metadata.contrast_ids"
            ),
            axis=axis,
        )
        _require_no_aggregation_metadata(profile_metadata, axis=axis)
        return

    if axis is ActivityProfileAxis.CONDITION_SUMMARY:
        _require_empty_manifest_labels(
            profile_metadata.sample_ids,
            field_name="bundle manifest.outputs.activity.profile_metadata.sample_ids",
            axis=axis,
        )
        _require_exact_manifest_labels(
            observed=profile_metadata.condition_ids,
            expected=profile_ids,
            field_name=(
                "bundle manifest.outputs.activity.profile_metadata.condition_ids"
            ),
            expected_label="activity/weighted_activity table columns",
        )
        _require_empty_manifest_labels(
            profile_metadata.contrast_ids,
            field_name=(
                "bundle manifest.outputs.activity.profile_metadata.contrast_ids"
            ),
            axis=axis,
        )
        aggregation = profile_metadata.aggregation_metadata
        if aggregation is None:
            _raise_activity_semantic_manifest_error(
                "bundle manifest.outputs.activity.profile_metadata."
                "aggregation_metadata",
                "must be a valid ActivityAggregationMetadata object for "
                "condition-summary activity semantics",
            )
        aggregation_profile_ids = tuple(
            record.profile_id for record in aggregation.records
        )
        if len(aggregation_profile_ids) != len(set(aggregation_profile_ids)):
            _raise_activity_semantic_manifest_error(
                "bundle manifest.outputs.activity.profile_metadata."
                "aggregation_metadata.records",
                "must not contain duplicate profile_id values",
            )
        _require_exact_manifest_labels(
            observed=aggregation_profile_ids,
            expected=profile_ids,
            field_name=(
                "bundle manifest.outputs.activity.profile_metadata."
                "aggregation_metadata.records[].profile_id"
            ),
            expected_label="activity/weighted_activity table columns",
        )
        return

    if axis is ActivityProfileAxis.CONTRAST:
        _require_empty_manifest_labels(
            profile_metadata.sample_ids,
            field_name="bundle manifest.outputs.activity.profile_metadata.sample_ids",
            axis=axis,
        )
        _require_empty_manifest_labels(
            profile_metadata.condition_ids,
            field_name=(
                "bundle manifest.outputs.activity.profile_metadata.condition_ids"
            ),
            axis=axis,
        )
        _require_exact_manifest_labels(
            observed=profile_metadata.contrast_ids,
            expected=profile_ids,
            field_name="bundle manifest.outputs.activity.profile_metadata.contrast_ids",
            expected_label="activity/weighted_activity table columns",
        )
        _require_no_aggregation_metadata(profile_metadata, axis=axis)
        return

    if axis is ActivityProfileAxis.EFFECT:
        _require_empty_manifest_labels(
            profile_metadata.sample_ids,
            field_name="bundle manifest.outputs.activity.profile_metadata.sample_ids",
            axis=axis,
        )
        _require_empty_manifest_labels(
            profile_metadata.condition_ids,
            field_name=(
                "bundle manifest.outputs.activity.profile_metadata.condition_ids"
            ),
            axis=axis,
        )
        _require_empty_manifest_labels(
            profile_metadata.contrast_ids,
            field_name=(
                "bundle manifest.outputs.activity.profile_metadata.contrast_ids"
            ),
            axis=axis,
        )
        _require_no_aggregation_metadata(profile_metadata, axis=axis)


def _validate_activity_provenance_agreement(
    *,
    provenance: RunProvenance,
    input_semantics: ActivityInputSemantics,
) -> None:
    workflow_parameters = provenance.workflow_parameters
    activity_config = _optional_json_mapping(
        workflow_parameters.get("activity_config"),
        field_name="bundle manifest.provenance.workflow_parameters.activity_config",
    )
    if activity_config is None:
        return
    method_input_contract = _optional_json_mapping(
        activity_config.get("method_input_contract"),
        field_name=(
            "bundle manifest.provenance.workflow_parameters.activity_config."
            "method_input_contract"
        ),
    )
    if method_input_contract is None:
        return
    expected_axis = _activity_profile_axis(input_semantics).value
    expected_quantity = _activity_quantitative_semantics(input_semantics).value
    _require_optional_provenance_semantic_agreement(
        method_input_contract,
        key="resolved_activity_profile_axis",
        expected=expected_axis,
        manifest_field=(
            "bundle manifest.outputs.activity.input_semantics.profile_axis"
        ),
    )
    _require_optional_provenance_semantic_agreement(
        method_input_contract,
        key="resolved_activity_quantitative_semantics",
        expected=expected_quantity,
        manifest_field=(
            "bundle manifest.outputs.activity.input_semantics.quantitative_semantics"
        ),
    )


def _require_optional_provenance_semantic_agreement(
    payload: Mapping[str, object],
    *,
    key: str,
    expected: str,
    manifest_field: str,
) -> None:
    if key not in payload or payload.get(key) is None:
        return
    provenance_field = (
        "bundle manifest.provenance.workflow_parameters.activity_config."
        f"method_input_contract.{key}"
    )
    observed = payload.get(key)
    if not isinstance(observed, str):
        _raise_activity_semantic_manifest_error(
            provenance_field,
            f"must be a string matching {manifest_field}; regenerate the bundle",
        )
    if observed != expected:
        _raise_activity_semantic_manifest_error(
            provenance_field,
            f"must agree with {manifest_field}; expected {expected!r}, "
            f"got {observed!r}",
        )


def _require_exact_manifest_labels(
    *,
    observed: tuple[str, ...],
    expected: tuple[str, ...],
    field_name: str,
    expected_label: str,
) -> None:
    observed_values = tuple(str(value) for value in observed)
    expected_values = tuple(str(value) for value in expected)
    if observed_values == expected_values:
        return
    _raise_activity_semantic_manifest_error(
        field_name,
        f"must exactly match {expected_label} in order; "
        f"expected={expected_values!r}, got={observed_values!r}",
    )


def _require_empty_manifest_labels(
    observed: tuple[str, ...],
    *,
    field_name: str,
    axis: ActivityProfileAxis,
) -> None:
    if not observed:
        return
    _raise_activity_semantic_manifest_error(
        field_name,
        f"must be empty when profile_metadata.axis is {axis.value!r}; "
        f"got={tuple(str(value) for value in observed)!r}",
    )


def _require_no_aggregation_metadata(
    profile_metadata: ActivityProfileMetadata,
    *,
    axis: ActivityProfileAxis,
) -> None:
    if profile_metadata.aggregation_metadata is None:
        return
    _raise_activity_semantic_manifest_error(
        "bundle manifest.outputs.activity.profile_metadata.aggregation_metadata",
        f"must be null when profile_metadata.axis is {axis.value!r}",
    )


def _raise_activity_semantic_manifest_error(
    field_name: str,
    message: str,
) -> NoReturn:
    raise PhosPyInputError(
        f"{field_name} {message}; correct the manifest or regenerate the bundle "
        "from the original KinaseActivityResult"
    )


def _normalise_site_metadata_bundle_table(table: pd.DataFrame) -> pd.DataFrame:
    if "site_key" in table.columns:
        return table
    if "site_key.1" not in table.columns:
        return table
    normalised = table.copy(deep=True)
    normalised = normalised.rename(columns={"site_key.1": "site_key"})
    return normalised


def _normalise_activity_statistics_bundle_table(
    table: pd.DataFrame | None,
) -> pd.DataFrame | None:
    if table is None:
        return None
    string_columns = (
        "kinase",
        "condition",
        "profile_id",
        "evidence_threshold_operator",
        "evidence_threshold_description",
        "computability_status",
        "reason",
        "significance_status",
        "inferential_status",
        "inferential_reason",
        "membership_source_category",
        "membership_selection_method",
    )
    normalized = table.copy(deep=True)
    for column_name in string_columns:
        if column_name in normalized.columns:
            normalized[column_name] = normalized[column_name].fillna("").astype(str)
    return normalized


def _parse_bundle_provenance(payload: Mapping[str, object]) -> RunProvenance:
    try:
        return provenance_from_payload(payload)
    except PhosPyInputError as exc:
        _raise_legacy_bundle_schema(exc)


def _validate_kinase_reference_projection_provenance(
    provenance: RunProvenance,
) -> RunProvenance:
    workflow_parameters = provenance.workflow_parameters
    if "reference_projection_summary" not in workflow_parameters:
        _raise_projection_manifest_error(
            "bundle manifest.provenance.workflow_parameters.reference_projection_summary",
            "is required for current kinase bundle reconstruction; regenerate the "
            "bundle with the current PhosPy version",
        )
    raw_summary = workflow_parameters.get("reference_projection_summary")
    if raw_summary is None:
        _raise_projection_manifest_error(
            "bundle manifest.provenance.workflow_parameters.reference_projection_summary",
            "must be a current schema projection-summary object; null does not "
            "prove zero reference attrition",
        )
    summary_field_name = (
        "bundle manifest.provenance.workflow_parameters.reference_projection_summary"
    )
    raw_summary = _projection_mapping_payload(
        raw_summary,
        field_name=summary_field_name,
    )
    try:
        projection_summary = KinaseReferenceProjectionSummary.from_payload(raw_summary)
    except WorkflowBoundaryError as exc:
        _raise_projection_manifest_error(
            "bundle manifest.provenance.workflow_parameters.reference_projection_summary",
            f"is invalid: {exc}",
        )
    if "universe_attrition" not in workflow_parameters:
        _raise_projection_manifest_error(
            "bundle manifest.provenance.workflow_parameters.universe_attrition",
            "is required for current kinase bundle reconstruction",
        )
    raw_universe_attrition = workflow_parameters.get("universe_attrition")
    raw_universe_attrition = _projection_mapping_payload(
        raw_universe_attrition,
        field_name="bundle manifest.provenance.workflow_parameters.universe_attrition",
    )
    canonical_reference_attrition = _validate_reference_attrition_agreement(
        projection_summary=projection_summary,
        universe_attrition=raw_universe_attrition,
    )
    canonical_universe_attrition = dict(raw_universe_attrition)
    canonical_universe_attrition["reference_attrition"] = canonical_reference_attrition
    canonical_workflow_parameters = dict(workflow_parameters)
    canonical_workflow_parameters["reference_projection_summary"] = (
        projection_summary.to_payload()
    )
    canonical_workflow_parameters["universe_attrition"] = canonical_universe_attrition
    return replace(provenance, workflow_parameters=canonical_workflow_parameters)


def _validate_reference_attrition_agreement(
    *,
    projection_summary: KinaseReferenceProjectionSummary,
    universe_attrition: Mapping[str, object],
) -> list[dict[str, object]]:
    raw_reference_attrition = universe_attrition.get("reference_attrition")
    field_name = (
        "bundle manifest.provenance.workflow_parameters.universe_attrition."
        "reference_attrition"
    )
    if not _is_sequence_payload(raw_reference_attrition):
        _raise_projection_manifest_error(field_name, "must be an array")
    reference_attrition = list(raw_reference_attrition)
    expected_records = _canonical_reference_attrition_records(projection_summary)
    if len(reference_attrition) != len(expected_records):
        _raise_projection_manifest_error(
            field_name,
            "must contain exactly the reference-attrition record implied by "
            "reference_projection_summary; expected "
            f"{len(expected_records)} record(s), got {len(reference_attrition)}",
        )
    for index, (observed, expected) in enumerate(
        zip(reference_attrition, expected_records, strict=True)
    ):
        record_field = f"{field_name}[{index}]"
        observed = _projection_mapping_payload(observed, field_name=record_field)
        _require_reference_attrition_record_matches(
            observed=observed,
            expected=expected,
            field_name=record_field,
        )
    return expected_records


def _canonical_reference_attrition_records(
    summary: KinaseReferenceProjectionSummary,
) -> list[dict[str, object]]:
    if summary.unmatched_source_substrate_identifier_count == 0:
        return []
    return [
        {
            "attrition_type": "reference_attrition",
            "stage": "reference_projection_to_dataset_site_key",
            "reason": (
                "source_reference_substrate_identifier_has_no_dataset_site_key_"
                "or_display_id_match"
            ),
            "input_universe": "source_reference_substrate_identifiers",
            "output_universe": (
                "source_reference_substrate_identifiers_with_dataset_projection"
            ),
            "input_identifier_namespace": summary.source_identifier_namespace,
            "output_identifier_namespace": summary.source_identifier_namespace,
            "projected_output_identifier_namespace": (
                summary.output_identifier_namespace
            ),
            "input_identity_semantics": summary.source_identity_semantics,
            "output_identity_semantics": (
                "source reference substrate identifiers that have at least one "
                "dataset projection; these are not dataset site_key rows"
            ),
            "projected_output_identity_semantics": summary.output_identity_semantics,
            "input_sites": int(summary.unique_source_substrate_identifier_count),
            "output_sites": int(summary.matched_source_substrate_identifier_count),
            "removed_sites": int(summary.unmatched_source_substrate_identifier_count),
            "input_identifier_count": int(
                summary.unique_source_substrate_identifier_count
            ),
            "output_identifier_count": int(
                summary.matched_source_substrate_identifier_count
            ),
            "removed_identifier_count": int(
                summary.unmatched_source_substrate_identifier_count
            ),
            "examples": list(summary.unmatched_source_substrate_identifier_examples),
            "removed_identifier_examples": list(
                summary.unmatched_source_substrate_identifier_examples
            ),
            "projected_dataset_site_key_count": int(
                summary.projected_dataset_site_key_count
            ),
            "one_to_many_display_reference_match_count": int(
                summary.one_to_many_display_reference_match_count
            ),
            "one_to_many_display_reference_site_key_rows": int(
                summary.one_to_many_display_reference_site_key_rows
            ),
            "one_to_many_projection_diagnostics": (
                "display_reference_matching.one_to_many_display_reference_matches"
            ),
            "interpreter_version": summary.interpreter_version,
        }
    ]


def _require_reference_attrition_record_matches(
    *,
    observed: Mapping[str, object],
    expected: Mapping[str, object],
    field_name: str,
) -> None:
    unsupported = sorted(str(key) for key in observed if str(key) not in expected)
    if unsupported:
        _raise_projection_manifest_error(
            field_name,
            "contains unsupported field(s): " + ", ".join(unsupported),
        )
    missing = sorted(key for key in expected if key not in observed)
    if missing:
        _raise_projection_manifest_error(
            field_name,
            "is missing required field(s): " + ", ".join(missing),
        )
    for key, expected_value in expected.items():
        observed_value = observed.get(key)
        if observed_value != expected_value:
            _raise_projection_manifest_error(
                f"{field_name}.{key}",
                "must agree with reference_projection_summary; "
                f"expected={expected_value!r}, got={observed_value!r}",
            )


def _is_sequence_payload(value: object) -> TypeGuard[Sequence[object]]:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


def _raise_projection_manifest_error(field_name: str, message: str) -> NoReturn:
    raise PhosPyInputError(
        f"{field_name} {message}; correct the manifest or regenerate the bundle "
        "from the original KinaseWorkflowResult"
    )


def _profile_self_inclusion_policy_from_provenance(
    provenance: RunProvenance,
) -> str:
    workflow_parameters = provenance.workflow_parameters
    scoring_config = _optional_json_mapping(
        workflow_parameters.get("scoring_config"),
        field_name="bundle manifest.provenance.workflow_parameters.scoring_config",
    )
    if scoring_config is None:
        return "allow"
    policy = scoring_config.get("profile_self_inclusion_policy")
    return policy if isinstance(policy, str) else "allow"


def _kinase_caveats_from_provenance(
    provenance: RunProvenance,
) -> tuple[KinaseWorkflowCaveat, ...]:
    workflow_parameters = provenance.workflow_parameters
    scoring_diagnostics = _optional_json_mapping(
        workflow_parameters.get("scoring_diagnostics"),
        field_name=(
            "bundle manifest.provenance.workflow_parameters.scoring_diagnostics"
        ),
    )
    if scoring_diagnostics is None:
        return ()
    violation_field_name = (
        "bundle manifest.provenance.workflow_parameters.scoring_diagnostics."
        "attrition_policy_violations"
    )
    raw_violations = scoring_diagnostics.get("attrition_policy_violations")
    if _is_object_list(raw_violations):
        violations = _mapping_list_payload(
            raw_violations,
            field_name=violation_field_name,
        )
    else:
        attrition_provenance = _optional_json_mapping(
            workflow_parameters.get("attrition_provenance"),
            field_name=(
                "bundle manifest.provenance.workflow_parameters.attrition_provenance"
            ),
        )
        violations = (
            []
            if attrition_provenance is None
            else _mapping_list_payload(
                attrition_provenance.get("policy_violations"),
                field_name=(
                    "bundle manifest.provenance.workflow_parameters."
                    "attrition_provenance.policy_violations"
                ),
            )
        )
    caveats: list[KinaseWorkflowCaveat] = []
    for violation in violations:
        raw_message = violation.get("message")
        if not isinstance(raw_message, str) or raw_message.strip() == "":
            continue
        raw_code = violation.get("code")
        code = (
            raw_code
            if isinstance(raw_code, str) and raw_code.strip() != ""
            else "kinase_attrition_policy_violation"
        )
        caveats.append(
            KinaseWorkflowCaveat(
                code=code,
                severity="warning",
                message=raw_message,
                details=dict(violation),
            )
        )
    return tuple(caveats)


def _kinase_attrition_provenance_from_provenance(
    provenance: RunProvenance,
) -> KinaseWorkflowAttritionProvenance | None:
    workflow_parameters = provenance.workflow_parameters
    raw_payload = workflow_parameters.get("attrition_provenance")
    payload = _optional_json_mapping(
        raw_payload,
        field_name="bundle manifest.provenance.workflow_parameters.attrition_provenance",
    )
    if payload is not None:
        return _kinase_attrition_provenance_from_payload(payload)
    scoring_diagnostics = _optional_json_mapping(
        workflow_parameters.get("scoring_diagnostics"),
        field_name=(
            "bundle manifest.provenance.workflow_parameters.scoring_diagnostics"
        ),
    )
    scoring_config = _optional_json_mapping(
        workflow_parameters.get("scoring_config"),
        field_name="bundle manifest.provenance.workflow_parameters.scoring_config",
    )
    if scoring_diagnostics is None or scoring_config is None:
        return None
    metrics = _optional_json_mapping(
        scoring_diagnostics.get("attrition_metrics"),
        field_name=(
            "bundle manifest.provenance.workflow_parameters.scoring_diagnostics."
            "attrition_metrics"
        ),
    )
    policy = _optional_json_mapping(
        scoring_config.get("attrition_policy"),
        field_name=(
            "bundle manifest.provenance.workflow_parameters.scoring_config."
            "attrition_policy"
        ),
    )
    if metrics is None or policy is None:
        return None
    raw_violations = scoring_diagnostics.get("attrition_policy_violations", [])
    violations = _mapping_list_payload(
        raw_violations,
        field_name=(
            "bundle manifest.provenance.workflow_parameters.scoring_diagnostics."
            "attrition_policy_violations"
        ),
    )
    outcome = "passed"
    if violations:
        outcome = "failed" if policy.get("on_violation") == "error" else "warned"
    return KinaseWorkflowAttritionProvenance(
        metrics=metrics,
        policy=policy,
        policy_outcome=outcome,
        policy_violations=tuple(violations),
        warning_messages=tuple(
            message
            for item in violations
            if isinstance((message := item.get("message")), str)
            and message.strip() != ""
        ),
    )


def _kinase_attrition_provenance_from_payload(
    payload: Mapping[str, object],
) -> KinaseWorkflowAttritionProvenance | None:
    metrics = _optional_json_mapping(
        payload.get("metrics"),
        field_name=(
            "bundle manifest.provenance.workflow_parameters.attrition_provenance."
            "metrics"
        ),
    )
    policy = _optional_json_mapping(
        payload.get("policy"),
        field_name=(
            "bundle manifest.provenance.workflow_parameters.attrition_provenance.policy"
        ),
    )
    policy_outcome = payload.get("policy_outcome")
    if metrics is None or policy is None:
        return None
    if not isinstance(policy_outcome, str):
        return None
    raw_violations = payload.get("policy_violations", [])
    violations = _mapping_list_payload(
        raw_violations,
        field_name=(
            "bundle manifest.provenance.workflow_parameters.attrition_provenance."
            "policy_violations"
        ),
    )
    raw_warnings = payload.get("warning_messages", [])
    warnings = _text_list_payload(raw_warnings)
    return KinaseWorkflowAttritionProvenance(
        metrics=metrics,
        policy=policy,
        policy_outcome=policy_outcome,
        policy_violations=tuple(violations),
        warning_messages=tuple(item for item in warnings if item.strip() != ""),
    )


def _parse_bundle_processing_state(
    payload: Mapping[str, object],
) -> DatasetProcessingState:
    try:
        return processing_state_from_payload(payload)
    except PhosPyInputError as exc:
        _raise_legacy_bundle_schema(exc)


def _parse_bundle_intensity_scale_state(
    payload: Mapping[str, object],
) -> IntensityScaleState:
    try:
        return intensity_scale_state_from_payload(
            payload,
            legacy_quantitative_meaning_policy="migrate_unverified",
        )
    except PhosPyInputError as exc:
        _raise_legacy_bundle_schema(exc)


def _raise_legacy_bundle_schema(exc: PhosPyInputError) -> NoReturn:
    raise PhosPyInputError(f"{_LEGACY_KINASE_BUNDLE_SCHEMA_ERROR} {exc}") from exc


def _read_absent_optional_table(
    *,
    bundle_root: Path,
    tables: Mapping[str, object],
    table_key: str,
    field_name: str,
) -> pd.DataFrame | None:
    if table_key not in tables:
        return None
    return read_optional_table(
        bundle_root=bundle_root,
        tables=tables,
        table_key=table_key,
        field_name=field_name,
    )


def _activity_profile_axis(
    input_semantics: ActivityInputSemantics,
) -> ActivityProfileAxis:
    axis = input_semantics.profile_axis
    if isinstance(axis, ActivityProfileAxis):
        return axis
    _raise_activity_semantic_manifest_error(
        "bundle manifest.outputs.activity.input_semantics.profile_axis",
        "must be a normalized ActivityProfileAxis",
    )


def _activity_quantitative_semantics(
    input_semantics: ActivityInputSemantics,
) -> ActivityQuantitativeSemantics:
    quantity = input_semantics.quantitative_semantics
    if isinstance(quantity, ActivityQuantitativeSemantics):
        return quantity
    _raise_activity_semantic_manifest_error(
        "bundle manifest.outputs.activity.input_semantics.quantitative_semantics",
        "must be a normalized ActivityQuantitativeSemantics",
    )


def _projection_mapping_payload(
    value: object,
    *,
    field_name: str,
) -> Mapping[str, object]:
    mapping = _optional_json_mapping(value, field_name=field_name)
    if mapping is not None:
        return mapping
    _raise_projection_manifest_error(field_name, "must be an object")


def _optional_json_mapping(
    value: object,
    *,
    field_name: str,
) -> Mapping[str, object] | None:
    if not _is_object_mapping(value):
        return None
    return validate_json_safe_mapping(value, field_name=field_name)


def _mapping_list_payload(
    value: object,
    *,
    field_name: str,
) -> list[Mapping[str, object]]:
    if not _is_object_list(value):
        return []
    mappings: list[Mapping[str, object]] = []
    for index, item in enumerate(value):
        mapping = _optional_json_mapping(item, field_name=f"{field_name}[{index}]")
        if mapping is not None:
            mappings.append(mapping)
    return mappings


def _text_list_payload(value: object) -> list[str]:
    if not _is_object_list(value):
        return []
    return [item for item in value if isinstance(item, str)]


def _is_object_mapping(value: object) -> TypeGuard[Mapping[object, object]]:
    return isinstance(value, Mapping)


def _is_object_list(value: object) -> TypeGuard[list[object]]:
    return isinstance(value, list)
