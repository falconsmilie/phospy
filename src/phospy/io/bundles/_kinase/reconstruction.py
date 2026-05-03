"""Reconstruct typed kinase workflow models from decoded bundle sections."""

from __future__ import annotations

from pathlib import Path

from phospy.activities.models import (
    ActivityMethodMetadata,
    ActivityMethodSummary,
    KinaseActivityResult,
)
from phospy.api.results import KinaseWorkflowResult
from phospy.datasets.models import AnalysisReadyPhosphoDataset
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
from phospy.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.provenance.serialization import from_payload as provenance_from_payload
from phospy.references.models import ReferenceBundle


def reconstruct_kinase_result(
    *,
    bundle_root: Path,
    sections: KinaseManifestSections,
) -> KinaseWorkflowResult:
    """Rebuild a KinaseWorkflowResult from already-validated manifest sections."""

    provenance = provenance_from_payload(sections.provenance_payload)
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
        score_fusion_weights=read_optional_table(
            bundle_root=bundle_root,
            tables=sections.scoring_tables,
            table_key="score_fusion_weights",
            field_name="bundle manifest.outputs.scoring.tables.score_fusion_weights",
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
    )
