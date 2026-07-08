"""Filesystem publishing helpers for dataset and workflow outputs."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from phospy.api.results import KinaseWorkflowResult, SignalomeWorkflowResult
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.processing_state import processing_state_to_payload
from phospy.io.readers.tables import table_suffix_for_format, write_table
from phospy.provenance.serialization import to_payload as provenance_to_payload
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.signalomes.models import SignalomeAlignmentInputDiagnostics


def publish_dataset(
    dataset: AnalysisReadyPhosphoDataset,
    output_root: Path,
    *,
    output_format: str = "csv",
) -> dict[str, Path]:
    """Publish an analysis-ready dataset into an output directory."""

    quantitative_meaning = dataset.intensity_scale_state.quantity
    if quantitative_meaning is None:
        raise PhosPyInputError(
            "dataset.intensity_scale_state.quantity must be established before publishing"
        )
    suffix = table_suffix_for_format(output_format)
    dataset_dir = Path(output_root) / "dataset"
    written: dict[str, Path] = {}

    _write_required_output_table(
        dataset.phospho,
        dataset_dir / f"phospho{suffix}",
        written=written,
        written_key="dataset.phospho",
    )
    _write_required_output_table(
        dataset.site_metadata,
        dataset_dir / f"site_metadata{suffix}",
        written=written,
        written_key="dataset.site_metadata",
    )
    _write_optional_output_table(
        dataset.sample_metadata,
        dataset_dir / f"sample_metadata{suffix}",
        written=written,
        written_key="dataset.sample_metadata",
    )
    _write_optional_output_table(
        dataset.total,
        dataset_dir / f"total{suffix}",
        written=written,
        written_key="dataset.total",
    )

    manifest_path = dataset_dir / "manifest.json"
    _write_manifest(
        manifest_path,
        _dataset_manifest_payload(
            dataset,
            output_format=output_format,
            quantitative_meaning=quantitative_meaning.value,
        ),
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

    _write_kinase_scoring_outputs(result, workflow_dir, suffix, written)
    _write_kinase_prediction_outputs(result, workflow_dir, suffix, written)
    _write_kinase_activity_outputs(result, workflow_dir, suffix, written)
    _write_kinase_reference_outputs(result, workflow_dir, suffix, written)

    manifest_path = workflow_dir / "manifest.json"
    _write_manifest(
        manifest_path,
        _kinase_manifest_payload(result, output_format=output_format),
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

    _write_signalome_outputs(result, workflow_dir, suffix, written)

    manifest_path = workflow_dir / "manifest.json"
    _write_manifest(
        manifest_path,
        _signalome_manifest_payload(result, output_format=output_format),
    )
    written["signalome.manifest"] = manifest_path
    return written


def _write_required_output_table(
    table: pd.DataFrame,
    path: Path,
    *,
    written: dict[str, Path],
    written_key: str,
) -> None:
    write_table(table, path)
    written[written_key] = path


def _write_optional_output_table(
    table: pd.DataFrame | None,
    path: Path,
    *,
    written: dict[str, Path],
    written_key: str,
) -> None:
    if table is None:
        return
    _write_required_output_table(
        table,
        path,
        written=written,
        written_key=written_key,
    )


def _dataset_manifest_payload(
    dataset: AnalysisReadyPhosphoDataset,
    *,
    output_format: str,
    quantitative_meaning: str,
) -> dict[str, object]:
    return {
        "organism": None if dataset.organism is None else dataset.organism.value,
        "intensity_scale": dataset.intensity_scale_state.label,
        "quantitative_meaning": quantitative_meaning,
        "processing_state": processing_state_to_payload(dataset.processing_state),
        "output_format": output_format,
        "provenance": (
            None
            if dataset.provenance is None
            else provenance_to_payload(dataset.provenance)
        ),
    }


def _write_kinase_scoring_outputs(
    result: KinaseWorkflowResult,
    workflow_dir: Path,
    suffix: str,
    written: dict[str, Path],
) -> None:
    scoring_dir = workflow_dir / "scoring"
    _write_required_output_table(
        result.scoring_result.profile_scores,
        scoring_dir / f"profile_scores{suffix}",
        written=written,
        written_key="kinase.scoring.profile_scores",
    )
    _write_optional_output_table(
        result.scoring_result.profile_score_diagnostics,
        scoring_dir / f"profile_score_diagnostics{suffix}",
        written=written,
        written_key="kinase.scoring.profile_score_diagnostics",
    )
    _write_optional_output_table(
        result.scoring_result.motif_scores,
        scoring_dir / f"motif_scores{suffix}",
        written=written,
        written_key="kinase.scoring.motif_scores",
    )
    _write_optional_output_table(
        result.scoring_result.rank_weighted_fusion_scores,
        scoring_dir / f"rank_weighted_fusion_scores{suffix}",
        written=written,
        written_key="kinase.scoring.rank_weighted_fusion_scores",
    )
    _write_optional_output_table(
        result.scoring_result.kinase_library_motif_scores,
        scoring_dir / f"kinase_library_motif_scores{suffix}",
        written=written,
        written_key="kinase.scoring.kinase_library_motif_scores",
    )
    _write_optional_output_table(
        result.scoring_result.combined_profile_motif_scores,
        scoring_dir / f"combined_profile_motif_scores{suffix}",
        written=written,
        written_key="kinase.scoring.combined_profile_motif_scores",
    )
    _write_optional_output_table(
        result.scoring_result.kinase_library_site_diagnostics,
        scoring_dir / f"kinase_library_site_diagnostics{suffix}",
        written=written,
        written_key="kinase.scoring.kinase_library_site_diagnostics",
    )
    _write_optional_output_table(
        result.scoring_result.kinase_library_kinase_diagnostics,
        scoring_dir / f"kinase_library_kinase_diagnostics{suffix}",
        written=written,
        written_key="kinase.scoring.kinase_library_kinase_diagnostics",
    )
    _write_optional_output_table(
        result.scoring_result.score_fusion_weights,
        scoring_dir / f"score_fusion_weights{suffix}",
        written=written,
        written_key="kinase.scoring.score_fusion_weights",
    )
    _write_optional_output_table(
        result.substrate_contributions,
        scoring_dir / f"substrate_contributions{suffix}",
        written=written,
        written_key="kinase.scoring.substrate_contributions",
    )


def _write_kinase_prediction_outputs(
    result: KinaseWorkflowResult,
    workflow_dir: Path,
    suffix: str,
    written: dict[str, Path],
) -> None:
    prediction_dir = workflow_dir / "prediction"
    _write_required_output_table(
        result.prediction_result.pred_mat,
        prediction_dir / f"pred_mat{suffix}",
        written=written,
        written_key="kinase.prediction.pred_mat",
    )
    _write_optional_output_table(
        result.prediction_result.substrate_list,
        prediction_dir / f"substrate_list{suffix}",
        written=written,
        written_key="kinase.prediction.substrate_list",
    )


def _write_kinase_activity_outputs(
    result: KinaseWorkflowResult,
    workflow_dir: Path,
    suffix: str,
    written: dict[str, Path],
) -> None:
    activity_result = result.activity_result
    if activity_result is None:
        return

    activity_dir = workflow_dir / "activity"
    _write_required_output_table(
        activity_result.activity_matrix,
        activity_dir / f"weighted_activity{suffix}",
        written=written,
        written_key="kinase.activity.weighted_activity",
    )
    _write_required_output_table(
        activity_result.thresholded_substrate_mean_activity,
        activity_dir / f"thresholded_substrate_mean_activity{suffix}",
        written=written,
        written_key="kinase.activity.thresholded_substrate_mean_activity",
    )
    _write_required_output_table(
        activity_result.thresholded_substrate_counts.to_frame(name="n_substrates"),
        activity_dir / f"thresholded_substrate_counts{suffix}",
        written=written,
        written_key="kinase.activity.thresholded_substrate_counts",
    )
    _write_optional_output_table(
        activity_result.activity_substrate_counts,
        activity_dir / f"activity_substrate_counts{suffix}",
        written=written,
        written_key="kinase.activity.activity_substrate_counts",
    )
    _write_required_output_table(
        activity_result.target_counts.to_frame(name="n_targets"),
        activity_dir / f"target_counts{suffix}",
        written=written,
        written_key="kinase.activity.target_counts",
    )
    _write_required_output_table(
        activity_result.target_table,
        activity_dir / f"target_table{suffix}",
        written=written,
        written_key="kinase.activity.target_table",
    )
    _write_optional_output_table(
        activity_result.statistics_table,
        activity_dir / f"statistics_table{suffix}",
        written=written,
        written_key="kinase.activity.statistics_table",
    )


def _write_kinase_reference_outputs(
    result: KinaseWorkflowResult,
    workflow_dir: Path,
    suffix: str,
    written: dict[str, Path],
) -> None:
    references_dir = workflow_dir / "references"
    _write_required_output_table(
        result.references.kinase_substrate_map,
        references_dir / f"kinase_substrate_map{suffix}",
        written=written,
        written_key="kinase.references.kinase_substrate_map",
    )
    _write_required_output_table(
        result.references.site_sequences,
        references_dir / f"site_sequences{suffix}",
        written=written,
        written_key="kinase.references.site_sequences",
    )


def _kinase_manifest_payload(
    result: KinaseWorkflowResult,
    *,
    output_format: str,
) -> dict[str, object]:
    activity_result = result.activity_result
    return {
        "reference_organism": result.references.organism.value,
        "scoring_mode": result.scoring_result.scoring_mode,
        "score_source": result.scoring_result.score_source,
        "score_scale": result.scoring_result.score_scale,
        "profile_self_inclusion_policy": str(
            result.scoring_result.profile_self_inclusion_policy
        ),
        "activity_enabled": activity_result is not None,
        "activity_method": (
            None
            if activity_result is None
            else activity_result.activity_method.to_payload()
        ),
        "activity_method_summary": (
            None
            if activity_result is None or activity_result.method_summary is None
            else activity_result.method_summary.to_payload()
        ),
        "activity_threshold_membership_diagnostics": (
            None
            if activity_result is None
            or activity_result.threshold_membership_diagnostics is None
            else activity_result.threshold_membership_diagnostics.to_payload()
        ),
        "activity_count_field_semantics": (
            None if activity_result is None else activity_result.count_field_semantics
        ),
        "site_attrition_summary": _site_attrition_summary_payload(result),
        "attrition_provenance": (
            None
            if result.attrition_provenance is None
            else result.attrition_provenance.to_payload()
        ),
        "eligibility_report": _eligibility_report_payload(result),
        "caveats": [caveat.to_payload() for caveat in result.caveats],
        "output_format": output_format,
        "provenance": (
            None
            if result.provenance is None
            else provenance_to_payload(result.provenance)
        ),
    }


def _site_attrition_summary_payload(
    result: KinaseWorkflowResult,
) -> dict[str, object] | None:
    summary = result.site_attrition_summary
    if summary is None:
        return None
    return {
        "preprocessing": {
            "input_rows": int(summary.preprocessing.input_rows),
            "rows_removed_during_preprocessing": int(
                summary.preprocessing.rows_removed_during_preprocessing
            ),
            "rows_removed_invalid_or_missing_site_identifiers": int(
                summary.preprocessing.rows_removed_invalid_or_missing_site_identifiers
            ),
            "duplicate_sites_merged_or_resolved": int(
                summary.preprocessing.duplicate_sites_merged_or_resolved
            ),
            "output_rows": int(summary.preprocessing.output_rows),
            "sequence_complete_sites": (
                None
                if summary.preprocessing.sequence_complete_sites is None
                else int(summary.preprocessing.sequence_complete_sites)
            ),
        },
        "scoring": {
            "rows_removed_invalid_or_missing_site_identifiers": int(
                summary.scoring.rows_removed_invalid_or_missing_site_identifiers
            ),
            "final_quantitative_sites_entering_scoring": int(
                summary.scoring.final_quantitative_sites_entering_scoring
            ),
            "sites_with_valid_site_sequence": int(
                summary.scoring.sites_with_valid_site_sequence
            ),
            "sites_without_usable_site_sequence": int(
                summary.scoring.sites_without_usable_site_sequence
            ),
            "sites_eligible_for_motif_scoring": int(
                summary.scoring.sites_eligible_for_motif_scoring
            ),
            "sites_with_kinase_substrate_reference_profile_evidence": int(
                summary.scoring.sites_with_kinase_substrate_reference_profile_evidence
            ),
            "sites_contributing_to_final_fused_prediction_scoring_output": int(
                summary.scoring.sites_contributing_to_final_fused_prediction_scoring_output
            ),
            "sites_contributing_to_activity_scoring": (
                None
                if summary.scoring.sites_contributing_to_activity_scoring is None
                else int(summary.scoring.sites_contributing_to_activity_scoring)
            ),
        },
    }


def _eligibility_report_payload(
    result: KinaseWorkflowResult,
) -> dict[str, object] | None:
    report = result.eligibility_report
    if report is None:
        return None
    return {
        "total_dataset_sites": int(report.total_dataset_sites),
        "sequence_complete_sites": int(report.sequence_complete_sites),
        "localisation_eligible_sites": (
            None
            if report.localisation_eligible_sites is None
            else int(report.localisation_eligible_sites)
        ),
        "reference_overlap_sites": int(report.reference_overlap_sites),
        "excluded_no_reference_match": int(report.excluded_no_reference_match),
        "excluded_low_localisation": (
            None
            if report.excluded_low_localisation is None
            else int(report.excluded_low_localisation)
        ),
        "eligible_kinases": int(report.eligible_kinases),
        "excluded_kinases_below_min_substrates": int(
            report.excluded_kinases_below_min_substrates
        ),
    }


def _write_signalome_outputs(
    result: SignalomeWorkflowResult,
    workflow_dir: Path,
    suffix: str,
    written: dict[str, Path],
) -> None:
    _write_required_output_table(
        result.module_assignments.table,
        workflow_dir / f"module_assignments{suffix}",
        written=written,
        written_key="signalome.module_assignments",
    )
    _write_required_output_table(
        result.signalome_modules.table,
        workflow_dir / f"signalome_modules{suffix}",
        written=written,
        written_key="signalome.signalome_modules",
    )
    _write_required_output_table(
        result.kinase_network.edges,
        workflow_dir / f"kinase_network_edges{suffix}",
        written=written,
        written_key="signalome.kinase_network.edges",
    )
    _write_optional_output_table(
        result.kinase_network.nodes,
        workflow_dir / f"kinase_network_nodes{suffix}",
        written=written,
        written_key="signalome.kinase_network.nodes",
    )
    _write_optional_output_table(
        result.kinase_network.candidate_correlations,
        workflow_dir / f"kinase_network_candidate_correlations{suffix}",
        written=written,
        written_key="signalome.kinase_network.candidate_correlations",
    )
    _write_optional_output_table(
        result.expanded_signalome,
        workflow_dir / f"expanded_signalome{suffix}",
        written=written,
        written_key="signalome.expanded_signalome",
    )


def _signalome_manifest_payload(
    result: SignalomeWorkflowResult,
    *,
    output_format: str,
) -> dict[str, object]:
    return {
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
        "score_preconditioning_diagnostics": _score_preconditioning_payload(result),
        "alignment_diagnostics": _alignment_diagnostics_payload(result),
        "network_correlation_diagnostics": _network_correlation_payload(result),
        "caveats": [caveat.to_payload() for caveat in result.caveats],
        "output_format": output_format,
        "provenance": (
            None
            if result.provenance is None
            else provenance_to_payload(result.provenance)
        ),
    }


def _score_preconditioning_payload(
    result: SignalomeWorkflowResult,
) -> dict[str, object]:
    diagnostics = result.score_preconditioning_diagnostics
    return {
        "policy": diagnostics.policy,
        "input_row_count": int(diagnostics.input_row_count),
        "dropped_all_missing_row_count": int(diagnostics.dropped_all_missing_row_count),
        "retained_row_count": int(diagnostics.retained_row_count),
    }


def _alignment_diagnostics_payload(
    result: SignalomeWorkflowResult,
) -> dict[str, object]:
    diagnostics = result.alignment_diagnostics
    return {
        "dataset_sites": _alignment_input_diagnostics_payload(
            diagnostics.dataset_sites
        ),
        "prediction_score_sites": _alignment_input_diagnostics_payload(
            diagnostics.prediction_score_sites
        ),
        "downstream_score_sites": _alignment_input_diagnostics_payload(
            diagnostics.downstream_score_sites
        ),
        "kinases": _alignment_input_diagnostics_payload(diagnostics.kinases),
        "protein_identifiers": _alignment_input_diagnostics_payload(
            diagnostics.protein_identifiers
        ),
    }


def _alignment_input_diagnostics_payload(
    diagnostics: SignalomeAlignmentInputDiagnostics,
) -> dict[str, object]:
    return {
        "provided_count": int(diagnostics.provided_count),
        "retained_count": int(diagnostics.retained_count),
        "dropped_count": int(diagnostics.dropped_count),
        "dropped_reasons": dict(diagnostics.dropped_reasons),
    }


def _network_correlation_payload(
    result: SignalomeWorkflowResult,
) -> dict[str, object]:
    diagnostics = result.kinase_network.correlation_diagnostics
    return {
        "total_candidate_correlations": int(diagnostics.total_candidate_correlations),
        "finite_correlations": int(diagnostics.finite_correlations),
        "undefined_correlations": int(diagnostics.undefined_correlations),
        "constant_profile_correlations": int(diagnostics.constant_profile_correlations),
        "insufficient_observation_correlations": int(
            diagnostics.insufficient_observation_correlations
        ),
        "missing_value_correlations": int(diagnostics.missing_value_correlations),
        "non_finite_value_correlations": int(diagnostics.non_finite_value_correlations),
        "edges_created": int(diagnostics.edges_created),
        "edges_skipped_non_finite_correlation": int(
            diagnostics.edges_skipped_non_finite_correlation
        ),
    }


def _write_manifest(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True),
            encoding="utf-8",
        )
    except OSError as exc:
        raise PhosPyInputError(f"failed to write manifest '{path}': {exc}") from exc
