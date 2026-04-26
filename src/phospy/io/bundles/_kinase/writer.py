"""Kinase bundle writing orchestration."""

from __future__ import annotations

from pathlib import Path

from phospy.api.results import KinaseWorkflowResult
from phospy.io.bundles._kinase.constants import (
    CONFIG_SNAPSHOT_RELATIVE_PATH,
    MANIFEST_FILENAME,
)
from phospy.io.bundles._kinase.manifest import build_manifest
from phospy.io.bundles._kinase.snapshots import KinaseWorkflowConfigSnapshot
from phospy.io.bundles._shared.json_files import write_json
from phospy.io.bundles._shared.tables import (
    write_bundle_table,
    write_optional_bundle_table,
)
from phospy.io.readers.tables import table_suffix_for_format


def save_kinase_workflow_bundle(
    result: KinaseWorkflowResult,
    output_root: Path,
    *,
    config_snapshot: KinaseWorkflowConfigSnapshot,
    output_format: str = "csv",
) -> dict[str, Path]:
    """Write a reproducible kinase output bundle and return written paths."""

    bundle_root = Path(output_root)
    suffix = table_suffix_for_format(output_format)
    normalized_format = output_format.strip().lower()
    written: dict[str, Path] = {}

    dataset_tables = {
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

    reference_tables = {
        "kinase_substrate_map": write_bundle_table(
            table=result.references.kinase_substrate_map,
            bundle_root=bundle_root,
            relative_path=Path("references") / f"kinase_substrate_map{suffix}",
            written=written,
            written_key="references.kinase_substrate_map",
        ),
        "site_sequences": write_bundle_table(
            table=result.references.site_sequences,
            bundle_root=bundle_root,
            relative_path=Path("references") / f"site_sequences{suffix}",
            written=written,
            written_key="references.site_sequences",
        ),
    }

    scoring_tables = {
        "profile_scores": write_bundle_table(
            table=result.scoring_result.profile_scores,
            bundle_root=bundle_root,
            relative_path=Path("scoring") / f"profile_scores{suffix}",
            written=written,
            written_key="scoring.profile_scores",
        ),
        "motif_scores": write_optional_bundle_table(
            table=result.scoring_result.motif_scores,
            bundle_root=bundle_root,
            relative_path=Path("scoring") / f"motif_scores{suffix}",
            written=written,
            written_key="scoring.motif_scores",
        ),
        "rank_weighted_fusion_scores": write_optional_bundle_table(
            table=result.scoring_result.rank_weighted_fusion_scores,
            bundle_root=bundle_root,
            relative_path=Path("scoring") / f"rank_weighted_fusion_scores{suffix}",
            written=written,
            written_key="scoring.rank_weighted_fusion_scores",
        ),
        "score_fusion_weights": write_optional_bundle_table(
            table=result.scoring_result.score_fusion_weights,
            bundle_root=bundle_root,
            relative_path=Path("scoring") / f"score_fusion_weights{suffix}",
            written=written,
            written_key="scoring.score_fusion_weights",
        ),
    }

    prediction_tables = {
        "pred_mat": write_bundle_table(
            table=result.prediction_result.pred_mat,
            bundle_root=bundle_root,
            relative_path=Path("prediction") / f"pred_mat{suffix}",
            written=written,
            written_key="prediction.pred_mat",
        ),
        "substrate_list": write_optional_bundle_table(
            table=result.prediction_result.substrate_list,
            bundle_root=bundle_root,
            relative_path=Path("prediction") / f"substrate_list{suffix}",
            written=written,
            written_key="prediction.substrate_list",
        ),
    }

    activity_tables = {
        "weighted_activity": write_optional_bundle_table(
            table=(
                None
                if result.activity_result is None
                else result.activity_result.weighted_activity
            ),
            bundle_root=bundle_root,
            relative_path=Path("activity") / f"weighted_activity{suffix}",
            written=written,
            written_key="activity.weighted_activity",
        ),
        "thresholded_substrate_mean_activity": write_optional_bundle_table(
            table=(
                None
                if result.activity_result is None
                else result.activity_result.thresholded_substrate_mean_activity
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
                if result.activity_result is None
                else result.activity_result.thresholded_substrate_counts.to_frame(
                    name="n_substrates"
                )
            ),
            bundle_root=bundle_root,
            relative_path=Path("activity") / f"thresholded_substrate_counts{suffix}",
            written=written,
            written_key="activity.thresholded_substrate_counts",
        ),
        "target_counts": write_optional_bundle_table(
            table=(
                None
                if result.activity_result is None
                else result.activity_result.target_counts.to_frame(name="n_targets")
            ),
            bundle_root=bundle_root,
            relative_path=Path("activity") / f"target_counts{suffix}",
            written=written,
            written_key="activity.target_counts",
        ),
        "target_table": write_optional_bundle_table(
            table=(
                None
                if result.activity_result is None
                else result.activity_result.target_table
            ),
            bundle_root=bundle_root,
            relative_path=Path("activity") / f"target_table{suffix}",
            written=written,
            written_key="activity.target_table",
        ),
    }

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
    )
    manifest_path = bundle_root / MANIFEST_FILENAME
    write_json(manifest_path, manifest, label="bundle manifest")
    written["manifest"] = manifest_path
    return written
