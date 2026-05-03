"""Reconstruct typed signalome workflow models from decoded bundle sections."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from phospy.activities.models import KinaseActivityResult
from phospy.api.results import KinaseWorkflowResult, SignalomeWorkflowResult
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.input import PhosPyInputError
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
from phospy.io.bundles._signalome.compatibility import (
    normalize_module_assignments_table,
    signalome_alignment_diagnostics_from_payload,
    signalome_module_selection_diagnostics_from_payload,
    signalome_network_correlation_diagnostics_from_payload,
    signalome_score_preconditioning_diagnostics_from_payload,
)
from phospy.io.bundles._signalome.manifest import SignalomeManifestSections
from phospy.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.provenance.serialization import from_payload as provenance_from_payload
from phospy.references.models import ReferenceBundle
from phospy.signalomes.context import SITE_MEMBERSHIP_EXCLUDED_REASON_COLUMN
from phospy.signalomes.models import (
    KinaseNetwork,
    SignalomeAlignmentDiagnostics,
    SignalomeAssignments,
    SignalomeModules,
    default_signalome_alignment_diagnostics,
)


def reconstruct_signalome_result(
    *,
    bundle_root: Path,
    sections: SignalomeManifestSections,
) -> SignalomeWorkflowResult:
    """Rebuild a SignalomeWorkflowResult from validated manifest sections."""

    signalome_provenance = provenance_from_payload(sections.provenance_payload)
    upstream_raw = signalome_provenance.workflow_parameters.get(
        "upstream_kinase_provenance"
    )
    if not isinstance(upstream_raw, Mapping):
        raise PhosPyInputError(
            "bundle manifest provenance.workflow_parameters.upstream_kinase_provenance "
            "is required for signalome bundles; regenerate this bundle with the "
            "current PhosPy version"
        )
    upstream_kinase_provenance = provenance_from_payload(upstream_raw)
    processing_state_payload = require_mapping(
        sections.dataset_metadata.get("processing_state"),
        field_name="bundle manifest.dataset.metadata.processing_state",
    )
    processing_state = processing_state_from_payload(processing_state_payload)
    intensity_scale_payload = require_mapping(
        sections.dataset_metadata.get("intensity_scale_state"),
        field_name="bundle manifest.dataset.metadata.intensity_scale_state",
    )

    dataset = AnalysisReadyPhosphoDataset(
        phospho=read_required_table(
            bundle_root=bundle_root,
            tables=sections.dataset_tables,
            table_key="phospho",
            field_name="bundle manifest.dataset.tables.phospho",
        ),
        site_metadata=read_required_table(
            bundle_root=bundle_root,
            tables=sections.dataset_tables,
            table_key="site_metadata",
            field_name="bundle manifest.dataset.tables.site_metadata",
        ),
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
        intensity_scale_state=intensity_scale_state_from_payload(
            intensity_scale_payload,
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
        score_fusion_weights=read_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="score_fusion_weights",
            field_name=(
                "bundle manifest.upstream_kinase_outputs.scoring.tables."
                "score_fusion_weights"
            ),
        ),
    )

    prediction_result = KinasePredictionResult(
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
        activity_result = KinaseActivityResult(
            weighted_activity=weighted_activity,
            thresholded_substrate_mean_activity=thresholded_substrate_mean_activity,
            thresholded_substrate_counts=thresholded_substrate_counts,
            activity_substrate_counts=activity_substrate_counts,
            target_counts=target_counts,
            target_table=target_table,
        )
    else:
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
        activity_result = None

    kinase_result = KinaseWorkflowResult(
        dataset=dataset,
        references=references,
        scoring_result=scoring_result,
        prediction_result=prediction_result,
        activity_result=activity_result,
        provenance=upstream_kinase_provenance,
    )
    module_selection_diagnostics = signalome_module_selection_diagnostics_from_payload(
        sections.signalome_metadata.get("module_selection_diagnostics"),
        scope="bundle manifest.signalome_outputs.metadata",
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
    candidate_correlations = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.signalome_tables,
        table_key="kinase_network_candidate_correlations",
        field_name="bundle manifest.signalome_outputs.tables.kinase_network_candidate_correlations",
    )
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
    protein_site_context = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.signalome_tables,
        table_key="protein_site_context",
        field_name="bundle manifest.signalome_outputs.tables.protein_site_context",
    )

    return SignalomeWorkflowResult(
        dataset=dataset,
        kinase_result=kinase_result,
        module_assignments=SignalomeAssignments(
            table=normalize_module_assignments_table(
                read_required_table(
                    bundle_root=bundle_root,
                    tables=sections.signalome_tables,
                    table_key="module_assignments",
                    field_name="bundle manifest.signalome_outputs.tables.module_assignments",
                )
            )
        ),
        signalome_modules=SignalomeModules(
            table=read_required_table(
                bundle_root=bundle_root,
                tables=sections.signalome_tables,
                table_key="signalome_modules",
                field_name="bundle manifest.signalome_outputs.tables.signalome_modules",
            )
        ),
        kinase_network=KinaseNetwork(
            edges=read_required_table(
                bundle_root=bundle_root,
                tables=sections.signalome_tables,
                table_key="kinase_network_edges",
                field_name="bundle manifest.signalome_outputs.tables.kinase_network_edges",
            ),
            nodes=read_optional_table(
                bundle_root=bundle_root,
                tables=sections.signalome_tables,
                table_key="kinase_network_nodes",
                field_name="bundle manifest.signalome_outputs.tables.kinase_network_nodes",
            ),
            candidate_correlations=candidate_correlations,
            correlation_diagnostics=network_correlation_diagnostics,
        ),
        module_selection_diagnostics=module_selection_diagnostics,
        score_preconditioning_diagnostics=score_preconditioning_diagnostics,
        alignment_diagnostics=alignment_diagnostics,
        expanded_signalome=read_optional_table(
            bundle_root=bundle_root,
            tables=sections.signalome_tables,
            table_key="expanded_signalome",
            field_name="bundle manifest.signalome_outputs.tables.expanded_signalome",
        ),
        site_membership=site_membership,
        protein_site_context=protein_site_context,
        provenance=signalome_provenance,
    )
