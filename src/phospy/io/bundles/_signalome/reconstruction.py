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
    signalome_module_selection_diagnostics_from_payload_with_compatibility_support,
    signalome_score_preconditioning_diagnostics_from_payload_with_compatibility_support,
)
from phospy.io.bundles._signalome.manifest import SignalomeManifestSections
from phospy.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.provenance.serialization import from_payload as provenance_from_payload
from phospy.references.models import ReferenceBundle
from phospy.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeModules,
)


def reconstruct_signalome_result(
    *,
    bundle_root: Path,
    sections: SignalomeManifestSections,
) -> SignalomeWorkflowResult:
    """Rebuild a SignalomeWorkflowResult from validated manifest sections."""

    signalome_provenance = (
        None
        if sections.provenance_payload is None
        else provenance_from_payload(sections.provenance_payload)
    )
    upstream_kinase_provenance = None
    if signalome_provenance is not None:
        upstream_raw = signalome_provenance.workflow_parameters.get(
            "upstream_kinase_provenance"
        )
        if isinstance(upstream_raw, Mapping):
            upstream_kinase_provenance = provenance_from_payload(upstream_raw)

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
            require_mapping(
                sections.dataset_metadata.get("intensity_scale_state"),
                field_name="bundle manifest.dataset.metadata.intensity_scale_state",
            )
        ),
        processing_state=processing_state_from_payload(
            require_mapping(
                sections.dataset_metadata.get("processing_state"),
                field_name="bundle manifest.dataset.metadata.processing_state",
            )
        ),
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
        combined_scores=read_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="combined_scores",
            field_name=(
                "bundle manifest.upstream_kinase_outputs.scoring.tables.combined_scores"
            ),
        ),
        weights=read_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="weights",
            field_name="bundle manifest.upstream_kinase_outputs.scoring.tables.weights",
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
    ksea_scores = read_optional_table(
        bundle_root=bundle_root,
        tables=sections.activity_tables,
        table_key="ksea_scores",
        field_name="bundle manifest.upstream_kinase_outputs.activity.tables.ksea_scores",
    )
    ksea_counts = read_optional_series(
        bundle_root=bundle_root,
        tables=sections.activity_tables,
        table_key="ksea_counts",
        field_name="bundle manifest.upstream_kinase_outputs.activity.tables.ksea_counts",
        series_name="n_substrates",
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
    if weighted_activity is None:
        activity_result = None
    else:
        if (
            ksea_scores is None
            or ksea_counts is None
            or target_counts is None
            or target_table is None
        ):
            raise PhosPyInputError(
                "bundle manifest upstream_kinase_outputs.activity.tables are incomplete for enabled activity outputs"
            )
        activity_result = KinaseActivityResult(
            weighted_activity=weighted_activity,
            ksea_scores=ksea_scores,
            ksea_counts=ksea_counts,
            target_counts=target_counts,
            target_table=target_table,
        )

    kinase_result = KinaseWorkflowResult(
        dataset=dataset,
        references=references,
        scoring_result=scoring_result,
        prediction_result=prediction_result,
        activity_result=activity_result,
        provenance=upstream_kinase_provenance,
    )
    module_selection_diagnostics = (
        signalome_module_selection_diagnostics_from_payload_with_compatibility_support(
            sections.signalome_metadata.get("module_selection_diagnostics"),
            scope="bundle manifest.signalome_outputs.metadata",
        )
    )
    score_preconditioning_diagnostics = signalome_score_preconditioning_diagnostics_from_payload_with_compatibility_support(
        sections.signalome_metadata.get("score_preconditioning_diagnostics"),
        scope="bundle manifest.signalome_outputs.metadata",
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
        ),
        module_selection_diagnostics=module_selection_diagnostics,
        score_preconditioning_diagnostics=score_preconditioning_diagnostics,
        expanded_signalome=read_optional_table(
            bundle_root=bundle_root,
            tables=sections.signalome_tables,
            table_key="expanded_signalome",
            field_name="bundle manifest.signalome_outputs.tables.expanded_signalome",
        ),
        provenance=signalome_provenance,
    )
