"""Signalome bundle writing orchestration."""

from __future__ import annotations

from pathlib import Path

from phospy.contracts.results import SignalomeWorkflowResult
from phospy.io.bundles._shared.json_files import write_json
from phospy.io.bundles._shared.tables import (
    write_bundle_table,
    write_optional_bundle_table,
)
from phospy.io.bundles._signalome.constants import (
    CONFIG_SNAPSHOT_RELATIVE_PATH,
    MANIFEST_FILENAME,
)
from phospy.io.bundles._signalome.manifest import build_manifest
from phospy.io.bundles._signalome.snapshots import SignalomeWorkflowConfigSnapshot
from phospy.io.readers.tables import table_suffix_for_format


def save_signalome_workflow_bundle(
    result: SignalomeWorkflowResult,
    output_root: Path,
    *,
    config_snapshot: SignalomeWorkflowConfigSnapshot,
    output_format: str = "csv",
) -> dict[str, Path]:
    """Write a reproducible signalome output bundle and return written paths."""

    bundle_root = Path(output_root)
    suffix = table_suffix_for_format(output_format)
    normalized_format = output_format.strip().lower()
    written: dict[str, Path] = {}

    dataset_tables = _write_dataset_tables(result, bundle_root, suffix, written)
    reference_tables = _write_reference_tables(result, bundle_root, suffix, written)
    scoring_tables = _write_scoring_tables(result, bundle_root, suffix, written)
    prediction_tables = _write_prediction_tables(result, bundle_root, suffix, written)
    activity_tables = _write_activity_tables(result, bundle_root, suffix, written)
    signalome_tables = _write_signalome_tables(result, bundle_root, suffix, written)

    config_path = bundle_root / Path(CONFIG_SNAPSHOT_RELATIVE_PATH)
    write_json(config_path, config_snapshot.to_payload(), label="config snapshot")
    written["config_snapshot"] = config_path

    manifest = build_manifest(
        result=result,
        table_format=normalized_format,
        dataset_tables=dataset_tables,
        reference_tables=reference_tables,
        scoring_tables=scoring_tables,
        prediction_tables=prediction_tables,
        activity_tables=activity_tables,
        signalome_tables=signalome_tables,
    )
    manifest_path = bundle_root / MANIFEST_FILENAME
    write_json(manifest_path, manifest, label="bundle manifest")
    written["manifest"] = manifest_path
    return written


def _write_dataset_tables(
    result: SignalomeWorkflowResult,
    bundle_root: Path,
    suffix: str,
    written: dict[str, Path],
) -> dict[str, object]:
    return {
        "phospho": write_bundle_table(
            table=result.dataset.phospho,
            bundle_root=bundle_root,
            relative_path=Path("dataset") / f"phospho{suffix}",
            written=written,
            written_key="dataset.phospho",
        ),
        "site_metadata": write_bundle_table(
            table=result.dataset.site_metadata,
            bundle_root=bundle_root,
            relative_path=Path("dataset") / f"site_metadata{suffix}",
            written=written,
            written_key="dataset.site_metadata",
        ),
        "sample_metadata": write_optional_bundle_table(
            table=result.dataset.sample_metadata,
            bundle_root=bundle_root,
            relative_path=Path("dataset") / f"sample_metadata{suffix}",
            written=written,
            written_key="dataset.sample_metadata",
        ),
        "total": write_optional_bundle_table(
            table=result.dataset.total,
            bundle_root=bundle_root,
            relative_path=Path("dataset") / f"total{suffix}",
            written=written,
            written_key="dataset.total",
        ),
    }


def _write_reference_tables(
    result: SignalomeWorkflowResult,
    bundle_root: Path,
    suffix: str,
    written: dict[str, Path],
) -> dict[str, object]:
    return {
        "kinase_substrate_map": write_bundle_table(
            table=result.kinase_result.references.kinase_substrate_map,
            bundle_root=bundle_root,
            relative_path=Path("references") / f"kinase_substrate_map{suffix}",
            written=written,
            written_key="references.kinase_substrate_map",
        ),
        "site_sequences": write_bundle_table(
            table=result.kinase_result.references.site_sequences,
            bundle_root=bundle_root,
            relative_path=Path("references") / f"site_sequences{suffix}",
            written=written,
            written_key="references.site_sequences",
        ),
    }


def _write_scoring_tables(
    result: SignalomeWorkflowResult,
    bundle_root: Path,
    suffix: str,
    written: dict[str, Path],
) -> dict[str, object]:
    scoring_tables = {
        "profile_scores": write_bundle_table(
            table=result.kinase_result.scoring_result.profile_scores,
            bundle_root=bundle_root,
            relative_path=Path("scoring") / f"profile_scores{suffix}",
            written=written,
            written_key="scoring.profile_scores",
        ),
        "motif_scores": write_optional_bundle_table(
            table=result.kinase_result.scoring_result.motif_scores,
            bundle_root=bundle_root,
            relative_path=Path("scoring") / f"motif_scores{suffix}",
            written=written,
            written_key="scoring.motif_scores",
        ),
        "rank_weighted_fusion_scores": write_optional_bundle_table(
            table=result.kinase_result.scoring_result.rank_weighted_fusion_scores,
            bundle_root=bundle_root,
            relative_path=Path("scoring") / f"rank_weighted_fusion_scores{suffix}",
            written=written,
            written_key="scoring.rank_weighted_fusion_scores",
        ),
        "kinase_library_motif_scores": write_optional_bundle_table(
            table=result.kinase_result.scoring_result.kinase_library_motif_scores,
            bundle_root=bundle_root,
            relative_path=Path("scoring") / f"kinase_library_motif_scores{suffix}",
            written=written,
            written_key="scoring.kinase_library_motif_scores",
        ),
        "combined_profile_motif_scores": write_optional_bundle_table(
            table=result.kinase_result.scoring_result.combined_profile_motif_scores,
            bundle_root=bundle_root,
            relative_path=Path("scoring") / f"combined_profile_motif_scores{suffix}",
            written=written,
            written_key="scoring.combined_profile_motif_scores",
        ),
        "score_fusion_weights": write_optional_bundle_table(
            table=result.kinase_result.scoring_result.score_fusion_weights,
            bundle_root=bundle_root,
            relative_path=Path("scoring") / f"score_fusion_weights{suffix}",
            written=written,
            written_key="scoring.score_fusion_weights",
        ),
        "kinase_library_site_diagnostics": write_optional_bundle_table(
            table=result.kinase_result.scoring_result.kinase_library_site_diagnostics,
            bundle_root=bundle_root,
            relative_path=Path("scoring") / f"kinase_library_site_diagnostics{suffix}",
            written=written,
            written_key="scoring.kinase_library_site_diagnostics",
        ),
        "kinase_library_kinase_diagnostics": write_optional_bundle_table(
            table=result.kinase_result.scoring_result.kinase_library_kinase_diagnostics,
            bundle_root=bundle_root,
            relative_path=Path("scoring")
            / f"kinase_library_kinase_diagnostics{suffix}",
            written=written,
            written_key="scoring.kinase_library_kinase_diagnostics",
        ),
    }
    return _drop_absent_kinase_library_scoring_tables(scoring_tables)


def _write_prediction_tables(
    result: SignalomeWorkflowResult,
    bundle_root: Path,
    suffix: str,
    written: dict[str, Path],
) -> dict[str, object]:
    return {
        "pred_mat": write_bundle_table(
            table=result.kinase_result.prediction_result.pred_mat,
            bundle_root=bundle_root,
            relative_path=Path("prediction") / f"pred_mat{suffix}",
            written=written,
            written_key="prediction.pred_mat",
        ),
        "substrate_list": write_optional_bundle_table(
            table=result.kinase_result.prediction_result.substrate_list,
            bundle_root=bundle_root,
            relative_path=Path("prediction") / f"substrate_list{suffix}",
            written=written,
            written_key="prediction.substrate_list",
        ),
    }


def _write_activity_tables(
    result: SignalomeWorkflowResult,
    bundle_root: Path,
    suffix: str,
    written: dict[str, Path],
) -> dict[str, object]:
    activity_result = result.kinase_result.activity_result
    return {
        "weighted_activity": write_optional_bundle_table(
            table=None if activity_result is None else activity_result.activity_matrix,
            bundle_root=bundle_root,
            relative_path=Path("activity") / f"weighted_activity{suffix}",
            written=written,
            written_key="activity.weighted_activity",
        ),
        "thresholded_substrate_mean_activity": write_optional_bundle_table(
            table=(
                None
                if activity_result is None
                else activity_result.thresholded_substrate_mean_activity
            ),
            bundle_root=bundle_root,
            relative_path=Path("activity")
            / f"thresholded_substrate_mean_activity{suffix}",
            written=written,
            written_key="activity.thresholded_substrate_mean_activity",
        ),
        "thresholded_substrate_counts": write_optional_bundle_table(
            table=(
                None
                if activity_result is None
                else activity_result.thresholded_substrate_counts.to_frame(
                    name="n_substrates"
                )
            ),
            bundle_root=bundle_root,
            relative_path=Path("activity") / f"thresholded_substrate_counts{suffix}",
            written=written,
            written_key="activity.thresholded_substrate_counts",
        ),
        "activity_substrate_counts": write_optional_bundle_table(
            table=(
                None
                if activity_result is None
                else activity_result.activity_substrate_counts
            ),
            bundle_root=bundle_root,
            relative_path=Path("activity") / f"activity_substrate_counts{suffix}",
            written=written,
            written_key="activity.activity_substrate_counts",
        ),
        "target_counts": write_optional_bundle_table(
            table=(
                None
                if activity_result is None
                else activity_result.target_counts.to_frame(name="n_targets")
            ),
            bundle_root=bundle_root,
            relative_path=Path("activity") / f"target_counts{suffix}",
            written=written,
            written_key="activity.target_counts",
        ),
        "target_table": write_optional_bundle_table(
            table=None if activity_result is None else activity_result.target_table,
            bundle_root=bundle_root,
            relative_path=Path("activity") / f"target_table{suffix}",
            written=written,
            written_key="activity.target_table",
        ),
    }


def _write_signalome_tables(
    result: SignalomeWorkflowResult,
    bundle_root: Path,
    suffix: str,
    written: dict[str, Path],
) -> dict[str, object]:
    return {
        "module_assignments": write_bundle_table(
            table=result.module_assignments.table,
            bundle_root=bundle_root,
            relative_path=Path("signalome") / f"module_assignments{suffix}",
            written=written,
            written_key="signalome.module_assignments",
        ),
        "signalome_modules": write_bundle_table(
            table=result.signalome_modules.table,
            bundle_root=bundle_root,
            relative_path=Path("signalome") / f"signalome_modules{suffix}",
            written=written,
            written_key="signalome.signalome_modules",
        ),
        "kinase_network_edges": write_bundle_table(
            table=result.kinase_network.edges,
            bundle_root=bundle_root,
            relative_path=Path("signalome") / f"kinase_network_edges{suffix}",
            written=written,
            written_key="signalome.kinase_network.edges",
        ),
        "kinase_network_nodes": write_optional_bundle_table(
            table=result.kinase_network.nodes,
            bundle_root=bundle_root,
            relative_path=Path("signalome") / f"kinase_network_nodes{suffix}",
            written=written,
            written_key="signalome.kinase_network.nodes",
        ),
        "kinase_network_candidate_correlations": write_optional_bundle_table(
            table=result.kinase_network.candidate_correlations,
            bundle_root=bundle_root,
            relative_path=Path("signalome")
            / f"kinase_network_candidate_correlations{suffix}",
            written=written,
            written_key="signalome.kinase_network.candidate_correlations",
        ),
        "expanded_signalome": write_optional_bundle_table(
            table=result.expanded_signalome,
            bundle_root=bundle_root,
            relative_path=Path("signalome") / f"expanded_signalome{suffix}",
            written=written,
            written_key="signalome.expanded_signalome",
        ),
        "site_membership": write_optional_bundle_table(
            table=result.site_membership,
            bundle_root=bundle_root,
            relative_path=Path("signalome") / f"site_membership{suffix}",
            written=written,
            written_key="signalome.site_membership",
        ),
        "protein_site_context": write_optional_bundle_table(
            table=result.protein_site_context,
            bundle_root=bundle_root,
            relative_path=Path("signalome") / f"protein_site_context{suffix}",
            written=written,
            written_key="signalome.protein_site_context",
        ),
    }


def _drop_absent_kinase_library_scoring_tables(
    scoring_tables: dict[str, object],
) -> dict[str, object]:
    """Omit absent optional Kinase Library extension tables from manifests."""

    optional_keys = (
        "kinase_library_motif_scores",
        "combined_profile_motif_scores",
        "kinase_library_site_diagnostics",
        "kinase_library_kinase_diagnostics",
    )
    return {
        key: value
        for key, value in scoring_tables.items()
        if key not in optional_keys or value is not None
    }
