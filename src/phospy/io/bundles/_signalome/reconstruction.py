"""Reconstruct typed signalome workflow models from decoded bundle sections."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import NoReturn, TypeGuard

import pandas as pd

from phospy.contracts.result_caveats import result_caveats_from_payloads
from phospy.contracts.results import (
    KinaseWorkflowAttritionProvenance,
    KinaseWorkflowCaveat,
    KinaseWorkflowResult,
    SignalomeWorkflowResult,
)
from phospy.errors.input import PhosPyInputError
from phospy.errors.workflows import WorkflowStageError
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
from phospy.io.bundles._signalome.diagnostics import (
    signalome_alignment_diagnostics_from_payload,
    signalome_clustering_preparation_diagnostics_from_payload,
    signalome_module_selection_diagnostics_from_payload,
    signalome_network_correlation_diagnostics_from_payload,
    signalome_score_preconditioning_diagnostics_from_payload,
)
from phospy.io.bundles._signalome.manifest import SignalomeManifestSections
from phospy.io.bundles._signalome.tables import (
    migrate_signalome_protein_group_id_column,
    normalize_module_assignments_table,
)
from phospy.provenance.models import RunProvenance
from phospy.provenance.serialization import from_payload as provenance_from_payload
from phospy.science.activities.models import KinaseActivityResult
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.datasets.processing_state import DatasetProcessingState
from phospy.science.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.science.references.models import ReferenceBundle
from phospy.science.signalomes._result_validation import (
    validate_signalome_result_site_level_identity,
)
from phospy.science.signalomes.constants import (
    LEGACY_PROTEIN_GROUP_ID_COLUMN,
    PROTEIN_GROUP_ID_COLUMN,
)
from phospy.science.signalomes.context import SITE_MEMBERSHIP_EXCLUDED_REASON_COLUMN
from phospy.science.signalomes.models import (
    KinaseNetwork,
    SignalomeAlignmentDiagnostics,
    SignalomeAssignments,
    SignalomeClusteringPreparationDiagnostics,
    SignalomeModules,
    SignalomeModuleSelectionDiagnostics,
    SignalomeNetworkCorrelationDiagnostics,
    SignalomeScorePreconditioningDiagnostics,
    default_signalome_alignment_diagnostics,
)
from phospy.science.transformations.models import IntensityScaleState

_LEGACY_SIGNALOME_BUNDLE_SCHEMA_ERROR = (
    "Legacy signalome bundle schemas are no longer supported. Regenerate the bundle "
    "with the current PhosPy version."
)
_LEGACY_SIGNALOME_DIAGNOSTIC_FIELDS = frozenset({"tree_engine", "tree_engine_version"})


@dataclass(frozen=True, slots=True)
class _BundleProvenances:
    signalome: RunProvenance
    upstream_kinase: RunProvenance


@dataclass(frozen=True, slots=True)
class _SignalomeDiagnostics:
    module_selection: SignalomeModuleSelectionDiagnostics
    clustering_preparation: SignalomeClusteringPreparationDiagnostics
    score_preconditioning: SignalomeScorePreconditioningDiagnostics
    alignment: SignalomeAlignmentDiagnostics
    network_correlation: SignalomeNetworkCorrelationDiagnostics


@dataclass(frozen=True, slots=True)
class _SignalomeOptionalTables:
    candidate_correlations: pd.DataFrame | None
    site_membership: pd.DataFrame | None
    protein_site_context: pd.DataFrame | None
    expanded_signalome: pd.DataFrame | None


def _normalize_site_metadata_for_dataset_contract(
    site_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Repair CSV drift for valid site_key indexes."""

    normalized = site_metadata.copy(deep=True)
    if "site_key" not in normalized.columns and "site_key.1" in normalized.columns:
        normalized = normalized.rename(columns={"site_key.1": "site_key"})
    if "site_key" not in normalized.columns and normalized.index.name == "site_key":
        normalized.loc[:, "site_key"] = normalized.index.astype(str).tolist()
    if (
        PROTEIN_GROUP_ID_COLUMN in normalized.columns
        and LEGACY_PROTEIN_GROUP_ID_COLUMN in normalized.columns
    ):
        current = (
            normalized.loc[:, PROTEIN_GROUP_ID_COLUMN]
            .fillna("")
            .astype(str)
            .str.strip()
        )
        legacy = (
            normalized.loc[:, LEGACY_PROTEIN_GROUP_ID_COLUMN]
            .fillna("")
            .astype(str)
            .str.strip()
        )
        if bool(current.ne(legacy).any()):
            raise PhosPyInputError(
                "bundle dataset site_metadata has conflicting Signalome grouping "
                "columns protein_group_id and legacy protein_id"
            )
    return normalized


def _normalize_optional_string_columns(
    table: pd.DataFrame,
    *,
    columns: tuple[str, ...],
) -> pd.DataFrame:
    normalized = table.copy(deep=True)
    for column_name in columns:
        if column_name in normalized.columns:
            column_index = _unique_column_position(
                normalized,
                column_name=column_name,
                field_name="bundle signalome table",
            )
            series = (
                normalized.loc[:, column_name].astype(object).fillna("").astype(str)
            )
            normalized = normalized.drop(columns=[column_name])
            normalized.insert(column_index, column_name, series)
    return normalized


def _unique_column_position(
    table: pd.DataFrame,
    *,
    column_name: str,
    field_name: str,
) -> int:
    position = table.columns.get_loc(column_name)
    if isinstance(position, int):
        return position
    raise PhosPyInputError(f"{field_name}.{column_name} must be a unique column")


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


def reconstruct_signalome_result(
    *,
    bundle_root: Path,
    sections: SignalomeManifestSections,
) -> SignalomeWorkflowResult:
    """Rebuild a SignalomeWorkflowResult from validated manifest sections."""

    provenances = _parse_bundle_provenances(sections)
    dataset = _reconstruct_dataset(
        bundle_root=bundle_root,
        sections=sections,
        provenance=provenances.signalome,
    )
    references = _reconstruct_references(bundle_root=bundle_root, sections=sections)
    kinase_result = _reconstruct_kinase_result(
        bundle_root=bundle_root,
        sections=sections,
        dataset=dataset,
        references=references,
        provenance=provenances.upstream_kinase,
    )
    diagnostics = _reconstruct_signalome_diagnostics(sections)
    optional_tables = _read_signalome_optional_tables(
        bundle_root=bundle_root,
        sections=sections,
    )
    module_assignments = normalize_module_assignments_table(
        read_required_table(
            bundle_root=bundle_root,
            tables=sections.signalome_tables,
            table_key="module_assignments",
            field_name="bundle manifest.signalome_outputs.tables.module_assignments",
        )
    )
    signalome_modules = read_required_table(
        bundle_root=bundle_root,
        tables=sections.signalome_tables,
        table_key="signalome_modules",
        field_name="bundle manifest.signalome_outputs.tables.signalome_modules",
    )
    network_edges = read_required_table(
        bundle_root=bundle_root,
        tables=sections.signalome_tables,
        table_key="kinase_network_edges",
        field_name="bundle manifest.signalome_outputs.tables.kinase_network_edges",
    )
    network_nodes = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.signalome_tables,
        table_key="kinase_network_nodes",
        field_name="bundle manifest.signalome_outputs.tables.kinase_network_nodes",
    )
    result = SignalomeWorkflowResult.from_trusted_owned(
        dataset=dataset,
        kinase_result=kinase_result,
        module_assignments=SignalomeAssignments.from_trusted_owned(
            table=module_assignments
        ),
        signalome_modules=SignalomeModules.from_trusted_owned(table=signalome_modules),
        kinase_network=KinaseNetwork.from_trusted_owned(
            edges=network_edges,
            nodes=network_nodes,
            candidate_correlations=optional_tables.candidate_correlations,
            correlation_diagnostics=diagnostics.network_correlation,
        ),
        module_selection_diagnostics=diagnostics.module_selection,
        clustering_preparation_diagnostics=diagnostics.clustering_preparation,
        score_preconditioning_diagnostics=diagnostics.score_preconditioning,
        alignment_diagnostics=diagnostics.alignment,
        expanded_signalome=optional_tables.expanded_signalome,
        site_membership=optional_tables.site_membership,
        protein_site_context=optional_tables.protein_site_context,
        provenance=provenances.signalome,
        caveats=result_caveats_from_payloads(sections.signalome_caveats_payload),
    )
    try:
        validate_signalome_result_site_level_identity(
            module_assignments=result.module_assignments.table,
            expanded_signalome=result.expanded_signalome,
            site_membership=result.site_membership,
            site_metadata=DatasetInternalView(dataset).site_metadata,
        )
    except WorkflowStageError as exc:
        raise PhosPyInputError(
            f"bundle signalome result identity validation failed: {exc}"
        ) from exc
    return result


def _parse_bundle_provenances(
    sections: SignalomeManifestSections,
) -> _BundleProvenances:
    signalome_provenance = _parse_bundle_provenance(sections.provenance_payload)
    _reject_legacy_signalome_diagnostic_fields(
        signalome_provenance.workflow_parameters,
        field_path="bundle manifest.provenance.workflow_parameters",
    )
    upstream_raw = signalome_provenance.workflow_parameters.get(
        "upstream_kinase_provenance"
    )
    upstream_payload = _optional_json_mapping(
        upstream_raw,
        field_name=(
            "bundle manifest provenance.workflow_parameters.upstream_kinase_provenance"
        ),
    )
    if upstream_payload is None:
        raise PhosPyInputError(
            "bundle manifest provenance.workflow_parameters.upstream_kinase_provenance "
            "is required for signalome bundles; regenerate this bundle with the "
            "current PhosPy version"
        )
    return _BundleProvenances(
        signalome=signalome_provenance,
        upstream_kinase=_parse_bundle_provenance(upstream_payload),
    )


def _reconstruct_dataset(
    *,
    bundle_root: Path,
    sections: SignalomeManifestSections,
    provenance: RunProvenance,
) -> AnalysisReadyPhosphoDataset:
    processing_state_payload = require_mapping(
        sections.dataset_metadata.get("processing_state"),
        field_name="bundle manifest.dataset.metadata.processing_state",
    )
    processing_state = _parse_bundle_processing_state(processing_state_payload)
    intensity_scale_payload = require_mapping(
        sections.dataset_metadata.get("intensity_scale_state"),
        field_name="bundle manifest.dataset.metadata.intensity_scale_state",
    )

    dataset_site_metadata = read_required_table(
        bundle_root=bundle_root,
        tables=sections.dataset_tables,
        table_key="site_metadata",
        field_name="bundle manifest.dataset.tables.site_metadata",
    )
    dataset_site_metadata = _normalize_site_metadata_for_dataset_contract(
        dataset_site_metadata
    )
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
    return AnalysisReadyPhosphoDataset.from_trusted_tables(
        phospho=phospho,
        site_metadata=dataset_site_metadata,
        sample_metadata=sample_metadata,
        total=total,
        organism=parse_optional_organism(
            sections.dataset_metadata.get("organism"),
            field_name="bundle manifest.dataset.metadata.organism",
        ),
        intensity_scale_state=intensity_scale_state,
        processing_state=processing_state,
        trusted_construction_assertions=build_bundle_reconstruction_assertions(
            bundle_kind="signalome_workflow_result",
            phospho=phospho,
            site_metadata=dataset_site_metadata,
            sample_metadata=sample_metadata,
            total=total,
            intensity_scale_state=intensity_scale_state,
            processing_state=processing_state,
            provenance=provenance,
        ),
    )


def _reconstruct_references(
    *,
    bundle_root: Path,
    sections: SignalomeManifestSections,
) -> ReferenceBundle:
    return ReferenceBundle.from_trusted_owned(
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


def _reconstruct_kinase_result(
    *,
    bundle_root: Path,
    sections: SignalomeManifestSections,
    dataset: AnalysisReadyPhosphoDataset,
    references: ReferenceBundle,
    provenance: RunProvenance,
) -> KinaseWorkflowResult:
    return KinaseWorkflowResult.from_trusted_owned(
        dataset=dataset,
        references=references,
        scoring_result=_reconstruct_scoring_result(
            bundle_root=bundle_root,
            sections=sections,
            provenance=provenance,
        ),
        prediction_result=_reconstruct_prediction_result(
            bundle_root=bundle_root,
            sections=sections,
        ),
        activity_result=_reconstruct_activity_result(
            bundle_root=bundle_root,
            sections=sections,
        ),
        provenance=provenance,
        attrition_provenance=_kinase_attrition_provenance_from_provenance(provenance),
        caveats=(
            result_caveats_from_payloads(sections.upstream_kinase_caveats_payload)
            or _kinase_caveats_from_provenance(provenance)
        ),
    )


def _kinase_caveats_from_provenance(
    provenance: RunProvenance,
) -> tuple[KinaseWorkflowCaveat, ...]:
    workflow_parameters = provenance.workflow_parameters
    scoring_diagnostics = _optional_json_mapping(
        workflow_parameters.get("scoring_diagnostics"),
        field_name=(
            "bundle manifest.upstream_kinase_provenance.workflow_parameters."
            "scoring_diagnostics"
        ),
    )
    if scoring_diagnostics is None:
        return ()
    violation_field_name = (
        "bundle manifest.upstream_kinase_provenance.workflow_parameters."
        "scoring_diagnostics.attrition_policy_violations"
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
                "bundle manifest.upstream_kinase_provenance.workflow_parameters."
                "attrition_provenance"
            ),
        )
        violations = (
            []
            if attrition_provenance is None
            else _mapping_list_payload(
                attrition_provenance.get("policy_violations"),
                field_name=(
                    "bundle manifest.upstream_kinase_provenance.workflow_parameters."
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
        field_name=(
            "bundle manifest.upstream_kinase_provenance.workflow_parameters."
            "attrition_provenance"
        ),
    )
    if payload is not None:
        return _kinase_attrition_provenance_from_payload(payload)
    scoring_diagnostics = _optional_json_mapping(
        workflow_parameters.get("scoring_diagnostics"),
        field_name=(
            "bundle manifest.upstream_kinase_provenance.workflow_parameters."
            "scoring_diagnostics"
        ),
    )
    scoring_config = _optional_json_mapping(
        workflow_parameters.get("scoring_config"),
        field_name=(
            "bundle manifest.upstream_kinase_provenance.workflow_parameters."
            "scoring_config"
        ),
    )
    if scoring_diagnostics is None or scoring_config is None:
        return None
    metrics = _optional_json_mapping(
        scoring_diagnostics.get("attrition_metrics"),
        field_name=(
            "bundle manifest.upstream_kinase_provenance.workflow_parameters."
            "scoring_diagnostics.attrition_metrics"
        ),
    )
    policy = _optional_json_mapping(
        scoring_config.get("attrition_policy"),
        field_name=(
            "bundle manifest.upstream_kinase_provenance.workflow_parameters."
            "scoring_config.attrition_policy"
        ),
    )
    if metrics is None or policy is None:
        return None
    raw_violations = scoring_diagnostics.get("attrition_policy_violations", [])
    violations = _mapping_list_payload(
        raw_violations,
        field_name=(
            "bundle manifest.upstream_kinase_provenance.workflow_parameters."
            "scoring_diagnostics.attrition_policy_violations"
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
            "bundle manifest.upstream_kinase_provenance.workflow_parameters."
            "attrition_provenance.metrics"
        ),
    )
    policy = _optional_json_mapping(
        payload.get("policy"),
        field_name=(
            "bundle manifest.upstream_kinase_provenance.workflow_parameters."
            "attrition_provenance.policy"
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
            "bundle manifest.upstream_kinase_provenance.workflow_parameters."
            "attrition_provenance.policy_violations"
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


def _reconstruct_scoring_result(
    *,
    bundle_root: Path,
    sections: SignalomeManifestSections,
    provenance: RunProvenance,
) -> KinaseScoringResult:
    return KinaseScoringResult.from_trusted_owned(
        profile_scores=read_required_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="profile_scores",
            field_name=(
                "bundle manifest.upstream_kinase_outputs.scoring.tables.profile_scores"
            ),
        ),
        motif_scores=read_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="motif_scores",
            field_name=(
                "bundle manifest.upstream_kinase_outputs.scoring.tables.motif_scores"
            ),
        ),
        rank_weighted_fusion_scores=read_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="rank_weighted_fusion_scores",
            field_name=(
                "bundle manifest.upstream_kinase_outputs.scoring.tables."
                "rank_weighted_fusion_scores"
            ),
        ),
        kinase_library_motif_scores=_read_absent_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="kinase_library_motif_scores",
            field_name=(
                "bundle manifest.upstream_kinase_outputs.scoring.tables."
                "kinase_library_motif_scores"
            ),
        ),
        combined_profile_motif_scores=_read_absent_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="combined_profile_motif_scores",
            field_name=(
                "bundle manifest.upstream_kinase_outputs.scoring.tables."
                "combined_profile_motif_scores"
            ),
        ),
        score_fusion_weights=read_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="score_fusion_weights",
            field_name=(
                "bundle manifest.upstream_kinase_outputs.scoring.tables."
                "score_fusion_weights"
            ),
        ),
        kinase_library_site_diagnostics=_read_absent_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="kinase_library_site_diagnostics",
            field_name=(
                "bundle manifest.upstream_kinase_outputs.scoring.tables."
                "kinase_library_site_diagnostics"
            ),
        ),
        kinase_library_kinase_diagnostics=_read_absent_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="kinase_library_kinase_diagnostics",
            field_name=(
                "bundle manifest.upstream_kinase_outputs.scoring.tables."
                "kinase_library_kinase_diagnostics"
            ),
        ),
        profile_self_inclusion_policy=_profile_self_inclusion_policy_from_provenance(
            provenance
        ),
    )


def _reconstruct_prediction_result(
    *,
    bundle_root: Path,
    sections: SignalomeManifestSections,
) -> KinasePredictionResult:
    return KinasePredictionResult.from_trusted_owned(
        pred_mat=read_required_table(
            bundle_root=bundle_root,
            tables=sections.prediction_tables,
            table_key="pred_mat",
            field_name=(
                "bundle manifest.upstream_kinase_outputs.prediction.tables.pred_mat"
            ),
        ),
        substrate_list=read_optional_table(
            bundle_root=bundle_root,
            tables=sections.prediction_tables,
            table_key="substrate_list",
            field_name=(
                "bundle manifest.upstream_kinase_outputs.prediction.tables."
                "substrate_list"
            ),
        ),
    )


def _reconstruct_activity_result(
    *,
    bundle_root: Path,
    sections: SignalomeManifestSections,
) -> KinaseActivityResult | None:
    weighted_activity = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.activity_tables,
        table_key="weighted_activity",
        field_name="bundle manifest.upstream_kinase_outputs.activity.tables.weighted_activity",
    )
    thresholded_substrate_mean_activity = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.activity_tables,
        table_key="thresholded_substrate_mean_activity",
        field_name=(
            "bundle manifest.upstream_kinase_outputs.activity.tables."
            "thresholded_substrate_mean_activity"
        ),
    )
    thresholded_substrate_counts = read_optional_series(
        bundle_root=bundle_root,
        tables=sections.activity_tables,
        table_key="thresholded_substrate_counts",
        field_name=(
            "bundle manifest.upstream_kinase_outputs.activity.tables."
            "thresholded_substrate_counts"
        ),
        series_name="n_substrates",
    )
    activity_substrate_counts = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.activity_tables,
        table_key="activity_substrate_counts",
        field_name=(
            "bundle manifest.upstream_kinase_outputs.activity.tables."
            "activity_substrate_counts"
        ),
    )
    target_counts = read_optional_series(
        bundle_root=bundle_root,
        tables=sections.activity_tables,
        table_key="target_counts",
        field_name="bundle manifest.upstream_kinase_outputs.activity.tables.target_counts",
        series_name="n_targets",
    )
    target_table = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.activity_tables,
        table_key="target_table",
        field_name="bundle manifest.upstream_kinase_outputs.activity.tables.target_table",
    )
    if sections.upstream_activity_enabled:
        if (
            weighted_activity is None
            or thresholded_substrate_mean_activity is None
            or thresholded_substrate_counts is None
            or target_counts is None
            or target_table is None
        ):
            raise PhosPyInputError(
                "bundle manifest upstream_kinase_outputs.activity.tables are incomplete for enabled activity outputs"
            )
        return KinaseActivityResult.from_trusted_owned(
            weighted_activity=weighted_activity,
            thresholded_substrate_mean_activity=thresholded_substrate_mean_activity,
            thresholded_substrate_counts=thresholded_substrate_counts,
            activity_substrate_counts=activity_substrate_counts,
            target_counts=target_counts,
            target_table=target_table,
        )

    if (
        weighted_activity is not None
        or thresholded_substrate_mean_activity is not None
        or thresholded_substrate_counts is not None
        or activity_substrate_counts is not None
        or target_counts is not None
        or target_table is not None
    ):
        raise PhosPyInputError(
            "bundle manifest upstream_kinase_outputs.activity.enabled=false must "
            "not declare populated activity tables"
        )
    return None


def _reconstruct_signalome_diagnostics(
    sections: SignalomeManifestSections,
) -> _SignalomeDiagnostics:
    module_selection_diagnostics = signalome_module_selection_diagnostics_from_payload(
        sections.signalome_metadata.get("module_selection_diagnostics"),
        scope="bundle manifest.signalome_outputs.metadata",
    )
    clustering_preparation_diagnostics = (
        signalome_clustering_preparation_diagnostics_from_payload(
            sections.signalome_metadata.get("clustering_preparation_diagnostics"),
            scope="bundle manifest.signalome_outputs.metadata",
        )
    )
    score_preconditioning_diagnostics = (
        signalome_score_preconditioning_diagnostics_from_payload(
            sections.signalome_metadata.get("score_preconditioning_diagnostics"),
            scope="bundle manifest.signalome_outputs.metadata",
        )
    )
    alignment_payload = sections.signalome_metadata.get("alignment_diagnostics")
    alignment_diagnostics: SignalomeAlignmentDiagnostics
    if alignment_payload is None:
        alignment_diagnostics = default_signalome_alignment_diagnostics()
    else:
        alignment_diagnostics = signalome_alignment_diagnostics_from_payload(
            alignment_payload,
            scope="bundle manifest.signalome_outputs.metadata",
        )
    network_correlation_diagnostics = (
        signalome_network_correlation_diagnostics_from_payload(
            sections.signalome_metadata.get("network_correlation_diagnostics"),
            scope="bundle manifest.signalome_outputs.metadata",
        )
    )
    return _SignalomeDiagnostics(
        module_selection=module_selection_diagnostics,
        clustering_preparation=clustering_preparation_diagnostics,
        score_preconditioning=score_preconditioning_diagnostics,
        alignment=alignment_diagnostics,
        network_correlation=network_correlation_diagnostics,
    )


def _read_signalome_optional_tables(
    *,
    bundle_root: Path,
    sections: SignalomeManifestSections,
) -> _SignalomeOptionalTables:
    candidate_correlations = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.signalome_tables,
        table_key="kinase_network_candidate_correlations",
        field_name="bundle manifest.signalome_outputs.tables.kinase_network_candidate_correlations",
    )
    return _SignalomeOptionalTables(
        candidate_correlations=candidate_correlations,
        site_membership=_read_site_membership_table(
            bundle_root=bundle_root,
            sections=sections,
        ),
        protein_site_context=_read_protein_site_context_table(
            bundle_root=bundle_root,
            sections=sections,
        ),
        expanded_signalome=_read_expanded_signalome_table(
            bundle_root=bundle_root,
            sections=sections,
        ),
    )


def _read_site_membership_table(
    *,
    bundle_root: Path,
    sections: SignalomeManifestSections,
) -> pd.DataFrame | None:
    site_membership = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.signalome_tables,
        table_key="site_membership",
        field_name="bundle manifest.signalome_outputs.tables.site_membership",
    )
    if (
        site_membership is not None
        and SITE_MEMBERSHIP_EXCLUDED_REASON_COLUMN in site_membership.columns
    ):
        site_membership = site_membership.copy(deep=True)
        site_membership.loc[:, SITE_MEMBERSHIP_EXCLUDED_REASON_COLUMN] = (
            site_membership.loc[:, SITE_MEMBERSHIP_EXCLUDED_REASON_COLUMN]
            .fillna("")
            .astype(str)
        )
    if site_membership is not None:
        site_membership = migrate_signalome_protein_group_id_column(
            site_membership,
            field_name="bundle manifest.signalome_outputs.tables.site_membership",
        )
        site_membership = _normalize_optional_string_columns(
            site_membership,
            columns=(
                "site_key",
                "display_id",
                "site_id",
                "gene_symbol",
                "site",
                "protein_group_id",
                "protein_accession",
                "isoform_id",
                "top_kinase",
                "excluded_reason",
            ),
        )
    return site_membership


def _read_protein_site_context_table(
    *,
    bundle_root: Path,
    sections: SignalomeManifestSections,
) -> pd.DataFrame | None:
    protein_site_context = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.signalome_tables,
        table_key="protein_site_context",
        field_name="bundle manifest.signalome_outputs.tables.protein_site_context",
    )
    if protein_site_context is not None:
        protein_site_context = migrate_signalome_protein_group_id_column(
            protein_site_context,
            field_name=(
                "bundle manifest.signalome_outputs.tables.protein_site_context"
            ),
        )
        protein_site_context = _normalize_optional_string_columns(
            protein_site_context,
            columns=(
                "protein_group_id",
                "gene_symbol",
                "site",
                "protein_accession",
                "isoform_id",
                "site_ids",
                "site_keys",
                "display_ids",
                "site_clusters",
                "top_kinases_by_site",
                "module_ids_by_site",
                "site_key_to_display_id",
            ),
        )
    return protein_site_context


def _read_expanded_signalome_table(
    *,
    bundle_root: Path,
    sections: SignalomeManifestSections,
) -> pd.DataFrame | None:
    expanded_signalome = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.signalome_tables,
        table_key="expanded_signalome",
        field_name="bundle manifest.signalome_outputs.tables.expanded_signalome",
    )
    if expanded_signalome is not None:
        expanded_signalome = migrate_signalome_protein_group_id_column(
            expanded_signalome,
            field_name="bundle manifest.signalome_outputs.tables.expanded_signalome",
        )
        expanded_signalome = _normalize_optional_string_columns(
            expanded_signalome,
            columns=(
                "site_key",
                "display_id",
                "site_id",
                "gene_symbol",
                "site",
                "protein_group_id",
                "protein_accession",
                "isoform_id",
                "top_kinase",
                "expanded_signalome_kinase",
                "expanded_signalome_row_kind",
                "expanded_signalome_assignment_policy",
                "expanded_signalome_linked_kinases",
                "expanded_signalome_regulated_module_ids",
                "expanded_signalome_support_kinases",
            ),
        )
    return expanded_signalome


def _parse_bundle_provenance(payload: Mapping[str, object]) -> RunProvenance:
    try:
        return provenance_from_payload(payload)
    except PhosPyInputError as exc:
        _raise_legacy_bundle_schema(exc)


def _profile_self_inclusion_policy_from_provenance(
    provenance: RunProvenance,
) -> str:
    workflow_parameters = provenance.workflow_parameters
    scoring_config = _optional_json_mapping(
        workflow_parameters.get("scoring_config"),
        field_name=(
            "bundle manifest.upstream_kinase_provenance.workflow_parameters."
            "scoring_config"
        ),
    )
    if scoring_config is None:
        return "allow"
    policy = scoring_config.get("profile_self_inclusion_policy")
    return policy if isinstance(policy, str) else "allow"


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


def _reject_legacy_signalome_diagnostic_fields(
    value: object,
    *,
    field_path: str,
) -> None:
    if _is_object_mapping(value):
        present = sorted(
            key for key in _LEGACY_SIGNALOME_DIAGNOSTIC_FIELDS if key in value
        )
        if present:
            fields = ", ".join(present)
            _raise_legacy_bundle_schema(
                PhosPyInputError(
                    f"{field_path} contains legacy signalome diagnostic field(s): "
                    f"{fields}"
                )
            )
        for key, item in value.items():
            _reject_legacy_signalome_diagnostic_fields(
                item,
                field_path=f"{field_path}.{str(key)}",
            )
        return
    if _is_object_sequence(value):
        for position, item in enumerate(value):
            _reject_legacy_signalome_diagnostic_fields(
                item,
                field_path=f"{field_path}[{position}]",
            )


def _raise_legacy_bundle_schema(exc: PhosPyInputError) -> NoReturn:
    raise PhosPyInputError(f"{_LEGACY_SIGNALOME_BUNDLE_SCHEMA_ERROR} {exc}") from exc


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


def _is_object_sequence(value: object) -> TypeGuard[Sequence[object]]:
    return isinstance(value, Sequence) and not isinstance(
        value,
        (str, bytes, bytearray),
    )
