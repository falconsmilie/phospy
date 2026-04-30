"""Filesystem publishing helpers for dataset and workflow outputs."""

from __future__ import annotations

import json
from pathlib import Path

from phospy.api.results import KinaseWorkflowResult, SignalomeWorkflowResult
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.processing_state import processing_state_to_payload
from phospy.io.readers.tables import table_suffix_for_format, write_table
from phospy.provenance.serialization import to_payload as provenance_to_payload


def publish_dataset(
    dataset: AnalysisReadyPhosphoDataset,
    output_root: Path,
    *,
    output_format: str = "csv",
) -> dict[str, Path]:
    """Publish an analysis-ready dataset into an output directory."""

    suffix = table_suffix_for_format(output_format)
    dataset_dir = Path(output_root) / "dataset"
    written: dict[str, Path] = {}

    phospho_path = dataset_dir / f"phospho{suffix}"
    write_table(dataset.phospho, phospho_path)
    written["dataset.phospho"] = phospho_path

    site_metadata_path = dataset_dir / f"site_metadata{suffix}"
    write_table(dataset.site_metadata, site_metadata_path)
    written["dataset.site_metadata"] = site_metadata_path

    if dataset.sample_metadata is not None:
        sample_metadata_path = dataset_dir / f"sample_metadata{suffix}"
        write_table(dataset.sample_metadata, sample_metadata_path)
        written["dataset.sample_metadata"] = sample_metadata_path

    if dataset.total is not None:
        total_path = dataset_dir / f"total{suffix}"
        write_table(dataset.total, total_path)
        written["dataset.total"] = total_path

    manifest_path = dataset_dir / "manifest.json"
    _write_manifest(
        manifest_path,
        {
            "organism": None if dataset.organism is None else dataset.organism.value,
            "intensity_scale": dataset.intensity_scale_state.label,
            "quantitative_meaning": dataset.intensity_scale_state.quantity.value,
            "processing_state": processing_state_to_payload(dataset.processing_state),
            "output_format": output_format,
            "provenance": (
                None
                if dataset.provenance is None
                else provenance_to_payload(dataset.provenance)
            ),
        },
    )
    written["dataset.manifest"] = manifest_path
    return written


def publish_kinase_workflow(
    result: KinaseWorkflowResult,
    output_root: Path,
    *,
    output_format: str = "csv",
) -> dict[str, Path]:
    """Publish kinase workflow outputs into an output directory."""

    suffix = table_suffix_for_format(output_format)
    written = publish_dataset(result.dataset, output_root, output_format=output_format)
    workflow_dir = Path(output_root) / "kinase"

    scoring_dir = workflow_dir / "scoring"
    profile_scores_path = scoring_dir / f"profile_scores{suffix}"
    write_table(result.scoring_result.profile_scores, profile_scores_path)
    written["kinase.scoring.profile_scores"] = profile_scores_path

    if result.scoring_result.motif_scores is not None:
        motif_scores_path = scoring_dir / f"motif_scores{suffix}"
        write_table(result.scoring_result.motif_scores, motif_scores_path)
        written["kinase.scoring.motif_scores"] = motif_scores_path

    if result.scoring_result.rank_weighted_fusion_scores is not None:
        rank_weighted_fusion_scores_path = (
            scoring_dir / f"rank_weighted_fusion_scores{suffix}"
        )
        write_table(
            result.scoring_result.rank_weighted_fusion_scores,
            rank_weighted_fusion_scores_path,
        )
        written["kinase.scoring.rank_weighted_fusion_scores"] = (
            rank_weighted_fusion_scores_path
        )

    if result.scoring_result.score_fusion_weights is not None:
        score_fusion_weights_path = scoring_dir / f"score_fusion_weights{suffix}"
        write_table(
            result.scoring_result.score_fusion_weights, score_fusion_weights_path
        )
        written["kinase.scoring.score_fusion_weights"] = score_fusion_weights_path

    prediction_dir = workflow_dir / "prediction"
    pred_mat_path = prediction_dir / f"pred_mat{suffix}"
    write_table(result.prediction_result.pred_mat, pred_mat_path)
    written["kinase.prediction.pred_mat"] = pred_mat_path

    if result.prediction_result.substrate_list is not None:
        substrate_list_path = prediction_dir / f"substrate_list{suffix}"
        write_table(result.prediction_result.substrate_list, substrate_list_path)
        written["kinase.prediction.substrate_list"] = substrate_list_path

    if result.activity_result is not None:
        activity_dir = workflow_dir / "activity"
        weighted_activity_path = activity_dir / f"weighted_activity{suffix}"
        write_table(result.activity_result.weighted_activity, weighted_activity_path)
        written["kinase.activity.weighted_activity"] = weighted_activity_path

        thresholded_substrate_mean_activity_path = (
            activity_dir / f"thresholded_substrate_mean_activity{suffix}"
        )
        write_table(
            result.activity_result.thresholded_substrate_mean_activity,
            thresholded_substrate_mean_activity_path,
        )
        written["kinase.activity.thresholded_substrate_mean_activity"] = (
            thresholded_substrate_mean_activity_path
        )

        thresholded_substrate_counts_path = (
            activity_dir / f"thresholded_substrate_counts{suffix}"
        )
        write_table(
            result.activity_result.thresholded_substrate_counts.to_frame(
                name="n_substrates"
            ),
            thresholded_substrate_counts_path,
        )
        written["kinase.activity.thresholded_substrate_counts"] = (
            thresholded_substrate_counts_path
        )

        target_counts_path = activity_dir / f"target_counts{suffix}"
        write_table(
            result.activity_result.target_counts.to_frame(name="n_targets"),
            target_counts_path,
        )
        written["kinase.activity.target_counts"] = target_counts_path

        target_table_path = activity_dir / f"target_table{suffix}"
        write_table(result.activity_result.target_table, target_table_path)
        written["kinase.activity.target_table"] = target_table_path

    references_dir = workflow_dir / "references"
    kinase_substrate_map_path = references_dir / f"kinase_substrate_map{suffix}"
    write_table(result.references.kinase_substrate_map, kinase_substrate_map_path)
    written["kinase.references.kinase_substrate_map"] = kinase_substrate_map_path

    site_sequences_path = references_dir / f"site_sequences{suffix}"
    write_table(result.references.site_sequences, site_sequences_path)
    written["kinase.references.site_sequences"] = site_sequences_path

    manifest_path = workflow_dir / "manifest.json"
    _write_manifest(
        manifest_path,
        {
            "reference_organism": result.references.organism.value,
            "activity_enabled": result.activity_result is not None,
            "output_format": output_format,
            "provenance": (
                None
                if result.provenance is None
                else provenance_to_payload(result.provenance)
            ),
        },
    )
    written["kinase.manifest"] = manifest_path
    return written


def publish_signalome_workflow(
    result: SignalomeWorkflowResult,
    output_root: Path,
    *,
    output_format: str = "csv",
) -> dict[str, Path]:
    """Publish signalome workflow outputs into the supported output layout."""

    suffix = table_suffix_for_format(output_format)
    written = publish_kinase_workflow(
        result.kinase_result,
        output_root,
        output_format=output_format,
    )
    workflow_dir = Path(output_root) / "signalome"

    module_assignments_path = workflow_dir / f"module_assignments{suffix}"
    write_table(result.module_assignments.table, module_assignments_path)
    written["signalome.module_assignments"] = module_assignments_path

    signalome_modules_path = workflow_dir / f"signalome_modules{suffix}"
    write_table(result.signalome_modules.table, signalome_modules_path)
    written["signalome.signalome_modules"] = signalome_modules_path

    kinase_network_edges_path = workflow_dir / f"kinase_network_edges{suffix}"
    write_table(result.kinase_network.edges, kinase_network_edges_path)
    written["signalome.kinase_network.edges"] = kinase_network_edges_path

    if result.kinase_network.nodes is not None:
        kinase_network_nodes_path = workflow_dir / f"kinase_network_nodes{suffix}"
        write_table(result.kinase_network.nodes, kinase_network_nodes_path)
        written["signalome.kinase_network.nodes"] = kinase_network_nodes_path

    if result.kinase_network.candidate_correlations is not None:
        candidate_correlations_path = (
            workflow_dir / f"kinase_network_candidate_correlations{suffix}"
        )
        write_table(
            result.kinase_network.candidate_correlations,
            candidate_correlations_path,
        )
        written["signalome.kinase_network.candidate_correlations"] = (
            candidate_correlations_path
        )

    if result.expanded_signalome is not None:
        expanded_signalome_path = workflow_dir / f"expanded_signalome{suffix}"
        write_table(result.expanded_signalome, expanded_signalome_path)
        written["signalome.expanded_signalome"] = expanded_signalome_path

    manifest_path = workflow_dir / "manifest.json"
    _write_manifest(
        manifest_path,
        {
            "reference_organism": result.kinase_result.references.organism.value,
            "expanded_signalome_present": result.expanded_signalome is not None,
            "kinase_network_nodes_present": result.kinase_network.nodes is not None,
            "module_selection_strategy": result.module_selection_diagnostics.strategy,
            "selected_module_count": int(
                result.module_selection_diagnostics.selected_module_count
            ),
            "used_automatic_module_selection": bool(
                result.module_selection_diagnostics.used_automatic_selection
            ),
            "score_preconditioning_diagnostics": {
                "policy": result.score_preconditioning_diagnostics.policy,
                "input_row_count": int(
                    result.score_preconditioning_diagnostics.input_row_count
                ),
                "dropped_all_missing_row_count": int(
                    result.score_preconditioning_diagnostics.dropped_all_missing_row_count
                ),
                "retained_row_count": int(
                    result.score_preconditioning_diagnostics.retained_row_count
                ),
            },
            "alignment_diagnostics": {
                "dataset_sites": {
                    "provided_count": int(
                        result.alignment_diagnostics.dataset_sites.provided_count
                    ),
                    "retained_count": int(
                        result.alignment_diagnostics.dataset_sites.retained_count
                    ),
                    "dropped_count": int(
                        result.alignment_diagnostics.dataset_sites.dropped_count
                    ),
                    "dropped_reasons": dict(
                        result.alignment_diagnostics.dataset_sites.dropped_reasons
                    ),
                },
                "prediction_score_sites": {
                    "provided_count": int(
                        result.alignment_diagnostics.prediction_score_sites.provided_count
                    ),
                    "retained_count": int(
                        result.alignment_diagnostics.prediction_score_sites.retained_count
                    ),
                    "dropped_count": int(
                        result.alignment_diagnostics.prediction_score_sites.dropped_count
                    ),
                    "dropped_reasons": dict(
                        result.alignment_diagnostics.prediction_score_sites.dropped_reasons
                    ),
                },
                "downstream_score_sites": {
                    "provided_count": int(
                        result.alignment_diagnostics.downstream_score_sites.provided_count
                    ),
                    "retained_count": int(
                        result.alignment_diagnostics.downstream_score_sites.retained_count
                    ),
                    "dropped_count": int(
                        result.alignment_diagnostics.downstream_score_sites.dropped_count
                    ),
                    "dropped_reasons": dict(
                        result.alignment_diagnostics.downstream_score_sites.dropped_reasons
                    ),
                },
                "kinases": {
                    "provided_count": int(
                        result.alignment_diagnostics.kinases.provided_count
                    ),
                    "retained_count": int(
                        result.alignment_diagnostics.kinases.retained_count
                    ),
                    "dropped_count": int(
                        result.alignment_diagnostics.kinases.dropped_count
                    ),
                    "dropped_reasons": dict(
                        result.alignment_diagnostics.kinases.dropped_reasons
                    ),
                },
                "protein_identifiers": {
                    "provided_count": int(
                        result.alignment_diagnostics.protein_identifiers.provided_count
                    ),
                    "retained_count": int(
                        result.alignment_diagnostics.protein_identifiers.retained_count
                    ),
                    "dropped_count": int(
                        result.alignment_diagnostics.protein_identifiers.dropped_count
                    ),
                    "dropped_reasons": dict(
                        result.alignment_diagnostics.protein_identifiers.dropped_reasons
                    ),
                },
            },
            "network_correlation_diagnostics": {
                "total_candidate_correlations": int(
                    result.kinase_network.correlation_diagnostics.total_candidate_correlations
                ),
                "finite_correlations": int(
                    result.kinase_network.correlation_diagnostics.finite_correlations
                ),
                "undefined_correlations": int(
                    result.kinase_network.correlation_diagnostics.undefined_correlations
                ),
                "constant_profile_correlations": int(
                    result.kinase_network.correlation_diagnostics.constant_profile_correlations
                ),
                "insufficient_observation_correlations": int(
                    result.kinase_network.correlation_diagnostics.insufficient_observation_correlations
                ),
                "missing_value_correlations": int(
                    result.kinase_network.correlation_diagnostics.missing_value_correlations
                ),
                "non_finite_value_correlations": int(
                    result.kinase_network.correlation_diagnostics.non_finite_value_correlations
                ),
                "edges_created": int(
                    result.kinase_network.correlation_diagnostics.edges_created
                ),
                "edges_skipped_non_finite_correlation": int(
                    result.kinase_network.correlation_diagnostics.edges_skipped_non_finite_correlation
                ),
            },
            "output_format": output_format,
            "provenance": (
                None
                if result.provenance is None
                else provenance_to_payload(result.provenance)
            ),
        },
    )
    written["signalome.manifest"] = manifest_path
    return written


def _write_manifest(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True),
            encoding="utf-8",
        )
    except OSError as exc:
        raise PhosPyInputError(f"failed to write manifest '{path}': {exc}") from exc
