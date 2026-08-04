"""Reconstruct typed signalome workflow models from decoded bundle sections."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import NoReturn

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
from phospy.io.bundles._shared.primitives import require_mapping
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
            column_index = normalized.columns.get_loc(column_name)
            series = (
                normalized.loc[:, column_name].astype(object).fillna("").astype(str)
            )
            normalized = normalized.drop(columns=[column_name])
            normalized.insert(column_index, column_name, series)  # pyright: ignore[reportArgumentType] - pandas-stubs cannot narrow get_loc to int for unique columns.
    return normalized


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
    result = SignalomeWorkflowResult._from_owned(
        dataset=dataset,
        kinase_result=kinase_result,
        module_assignments=SignalomeAssignments._from_owned(table=module_assignments),
        signalome_modules=SignalomeModules._from_owned(table=signalome_modules),
        kinase_network=KinaseNetwork._from_owned(
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
    if not isinstance(upstream_raw, Mapping):
        raise PhosPyInputError(
            "bundle manifest provenance.workflow_parameters.upstream_kinase_provenance "
            "is required for signalome bundles; regenerate this bundle with the "
            "current PhosPy version"
        )
    return _BundleProvenances(
        signalome=signalome_provenance,
        upstream_kinase=_parse_bundle_provenance(upstream_raw),
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
    return ReferenceBundle._from_owned(
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
    return KinaseWorkflowResult._from_owned(
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
                severity="warning",
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


def _reconstruct_scoring_result(
    *,
    bundle_root: Path,
    sections: SignalomeManifestSections,
    provenance: RunProvenance,
) -> KinaseScoringResult:
    return KinaseScoringResult._from_owned(
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
    return KinasePredictionResult._from_owned(
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
        return KinaseActivityResult._from_owned(
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
    if not isinstance(workflow_parameters, Mapping):
        return "allow"
    scoring_config = workflow_parameters.get("scoring_config")
    if not isinstance(scoring_config, Mapping):
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
    if isinstance(value, Mapping):
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
    if isinstance(value, (list, tuple)):
        for position, item in enumerate(value):
            _reject_legacy_signalome_diagnostic_fields(
                item,
                field_path=f"{field_path}[{position}]",
            )


def _raise_legacy_bundle_schema(exc: PhosPyInputError) -> NoReturn:
    raise PhosPyInputError(f"{_LEGACY_SIGNALOME_BUNDLE_SCHEMA_ERROR} {exc}") from exc
