"""Reconstruct typed kinase workflow models from decoded bundle sections."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import NoReturn

from phospy.api.results import (
    KinaseWorkflowAttritionProvenance,
    KinaseWorkflowCaveat,
    KinaseWorkflowResult,
)
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._kinase.manifest import KinaseManifestSections
from phospy.io.bundles._shared.intensity_scale_state import (
    intensity_scale_state_from_payload,
)
from phospy.io.bundles._shared.organisms import (
    parse_optional_organism,
    parse_required_organism,
)
from phospy.io.bundles._shared.primitives import require_mapping
from phospy.io.bundles._shared.processing_state import (
    processing_state_from_payload,
)
from phospy.io.bundles._shared.tables import (
    read_optional_series,
    read_optional_table,
    read_required_table,
)
from phospy.provenance.models import RunProvenance
from phospy.provenance.serialization import from_payload as provenance_from_payload
from phospy.science.activities.models import (
    ActivityMethodMetadata,
    ActivityMethodSummary,
    KinaseActivityResult,
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
    dataset = AnalysisReadyPhosphoDataset._from_owned(
        phospho=read_required_table(
            bundle_root=bundle_root,
            tables=sections.dataset_tables,
            table_key="phospho",
            field_name="bundle manifest.dataset.tables.phospho",
        ),
        site_metadata=site_metadata,
        sample_metadata=read_optional_table(
            bundle_root=bundle_root,
            tables=sections.dataset_tables,
            table_key="sample_metadata",
            field_name="bundle manifest.dataset.tables.sample_metadata",
        ),
        total=read_optional_table(
            bundle_root=bundle_root,
            tables=sections.dataset_tables,
            table_key="total",
            field_name="bundle manifest.dataset.tables.total",
        ),
        organism=parse_optional_organism(
            sections.dataset_metadata.get("organism"),
            field_name="bundle manifest.dataset.metadata.organism",
        ),
        intensity_scale_state=_parse_bundle_intensity_scale_state(
            intensity_scale_payload
        ),
        processing_state=processing_state,
    )

    references = ReferenceBundle(
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

    scoring_result = KinaseScoringResult(
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
    )

    prediction_result = KinasePredictionResult(
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
            activity_method_summary = ActivityMethodSummary.from_payload(
                sections.activity_method_summary
            )
        activity_result = KinaseActivityResult(
            weighted_activity=weighted_activity,
            thresholded_substrate_mean_activity=thresholded_substrate_mean_activity,
            thresholded_substrate_counts=thresholded_substrate_counts,
            activity_substrate_counts=activity_substrate_counts,
            target_counts=target_counts,
            target_table=target_table,
            statistics_table=statistics_table,
            method_summary=activity_method_summary,
            activity_method=activity_method,
        )
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
        activity_result = None

    return KinaseWorkflowResult(
        dataset=dataset,
        references=references,
        scoring_result=scoring_result,
        prediction_result=prediction_result,
        activity_result=activity_result,
        provenance=provenance,
        substrate_contributions=substrate_contributions,
        attrition_provenance=_kinase_attrition_provenance_from_provenance(provenance),
        caveats=_kinase_caveats_from_provenance(provenance),
        _assume_owned=True,
    )


def _normalise_site_metadata_bundle_table(table):
    if "site_key" in table.columns:
        return table
    if "site_key.1" not in table.columns:
        return table
    normalised = table.copy(deep=True)
    normalised = normalised.rename(columns={"site_key.1": "site_key"})
    return normalised


def _parse_bundle_provenance(payload: Mapping[str, object]) -> RunProvenance:
    try:
        return provenance_from_payload(payload)
    except PhosPyInputError as exc:
        _raise_legacy_bundle_schema(exc)


def _kinase_caveats_from_provenance(
    provenance: RunProvenance,
) -> tuple[KinaseWorkflowCaveat, ...]:
    workflow_parameters = provenance.workflow_parameters
    if not isinstance(workflow_parameters, Mapping):
        return ()
    scoring_diagnostics = workflow_parameters.get("scoring_diagnostics")
    if not isinstance(scoring_diagnostics, Mapping):
        return ()
    raw_violations = scoring_diagnostics.get("attrition_policy_violations")
    if not isinstance(raw_violations, list):
        attrition_provenance = workflow_parameters.get("attrition_provenance")
        if isinstance(attrition_provenance, Mapping):
            raw_violations = attrition_provenance.get("policy_violations")
    if not isinstance(raw_violations, list):
        return ()
    caveats: list[KinaseWorkflowCaveat] = []
    for raw_violation in raw_violations:
        if not isinstance(raw_violation, Mapping):
            continue
        raw_message = raw_violation.get("message")
        if not isinstance(raw_message, str) or raw_message.strip() == "":
            continue
        raw_code = raw_violation.get("code")
        code = (
            raw_code
            if isinstance(raw_code, str) and raw_code.strip() != ""
            else "kinase_attrition_policy_violation"
        )
        caveats.append(
            KinaseWorkflowCaveat(
                code=code,
                message=raw_message,
                details=dict(raw_violation),
            )
        )
    return tuple(caveats)


def _kinase_attrition_provenance_from_provenance(
    provenance: RunProvenance,
) -> KinaseWorkflowAttritionProvenance | None:
    workflow_parameters = provenance.workflow_parameters
    if not isinstance(workflow_parameters, Mapping):
        return None
    raw_payload = workflow_parameters.get("attrition_provenance")
    if isinstance(raw_payload, Mapping):
        return _kinase_attrition_provenance_from_payload(raw_payload)
    scoring_diagnostics = workflow_parameters.get("scoring_diagnostics")
    scoring_config = workflow_parameters.get("scoring_config")
    if not isinstance(scoring_diagnostics, Mapping) or not isinstance(
        scoring_config, Mapping
    ):
        return None
    metrics = scoring_diagnostics.get("attrition_metrics")
    policy = scoring_config.get("attrition_policy")
    if not isinstance(metrics, Mapping) or not isinstance(policy, Mapping):
        return None
    raw_violations = scoring_diagnostics.get("attrition_policy_violations", [])
    violations = raw_violations if isinstance(raw_violations, list) else []
    outcome = "passed"
    if violations:
        outcome = "failed" if policy.get("on_violation") == "error" else "warned"
    return KinaseWorkflowAttritionProvenance(
        metrics=metrics,
        policy=policy,
        policy_outcome=outcome,
        policy_violations=tuple(
            item for item in violations if isinstance(item, Mapping)
        ),
        warning_messages=tuple(
            str(item.get("message"))
            for item in violations
            if isinstance(item, Mapping)
            and isinstance(item.get("message"), str)
            and str(item.get("message")).strip() != ""
        ),
    )


def _kinase_attrition_provenance_from_payload(
    payload: Mapping[str, object],
) -> KinaseWorkflowAttritionProvenance | None:
    metrics = payload.get("metrics")
    policy = payload.get("policy")
    policy_outcome = payload.get("policy_outcome")
    if not isinstance(metrics, Mapping) or not isinstance(policy, Mapping):
        return None
    if not isinstance(policy_outcome, str):
        return None
    raw_violations = payload.get("policy_violations", [])
    violations = raw_violations if isinstance(raw_violations, list) else []
    raw_warnings = payload.get("warning_messages", [])
    warnings = raw_warnings if isinstance(raw_warnings, list) else []
    return KinaseWorkflowAttritionProvenance(
        metrics=metrics,
        policy=policy,
        policy_outcome=policy_outcome,
        policy_violations=tuple(
            item for item in violations if isinstance(item, Mapping)
        ),
        warning_messages=tuple(
            str(item)
            for item in warnings
            if isinstance(item, str) and item.strip() != ""
        ),
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
        return intensity_scale_state_from_payload(payload)
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
):
    if table_key not in tables:
        return None
    return read_optional_table(
        bundle_root=bundle_root,
        tables=tables,
        table_key=table_key,
        field_name=field_name,
    )
