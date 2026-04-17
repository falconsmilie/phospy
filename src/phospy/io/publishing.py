"""Filesystem publishing helpers for dataset and workflow outputs."""

from __future__ import annotations

import json
from pathlib import Path

from phospy.api.results import SimpleKinaseWorkflowResult
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.input import PhosPyInputError
from phospy.io.tables import table_suffix_for_format, write_table


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
            "transformation_state": dataset.transformation_state.label,
            "output_format": output_format,
        },
    )
    written["dataset.manifest"] = manifest_path
    return written


def publish_simple_kinase_workflow(
    result: SimpleKinaseWorkflowResult,
    output_root: Path,
    *,
    output_format: str = "csv",
) -> dict[str, Path]:
    """Publish simple kinase workflow outputs into an output directory."""

    suffix = table_suffix_for_format(output_format)
    written = publish_dataset(result.dataset, output_root, output_format=output_format)
    workflow_dir = Path(output_root) / "simple_kinase"

    scoring_dir = workflow_dir / "scoring"
    profile_scores_path = scoring_dir / f"profile_scores{suffix}"
    write_table(result.scoring_result.profile_scores, profile_scores_path)
    written["simple_kinase.scoring.profile_scores"] = profile_scores_path

    if result.scoring_result.motif_scores is not None:
        motif_scores_path = scoring_dir / f"motif_scores{suffix}"
        write_table(result.scoring_result.motif_scores, motif_scores_path)
        written["simple_kinase.scoring.motif_scores"] = motif_scores_path

    if result.scoring_result.combined_scores is not None:
        combined_scores_path = scoring_dir / f"combined_scores{suffix}"
        write_table(result.scoring_result.combined_scores, combined_scores_path)
        written["simple_kinase.scoring.combined_scores"] = combined_scores_path

    if result.scoring_result.weights is not None:
        weights_path = scoring_dir / f"weights{suffix}"
        write_table(result.scoring_result.weights, weights_path)
        written["simple_kinase.scoring.weights"] = weights_path

    prediction_dir = workflow_dir / "prediction"
    pred_mat_path = prediction_dir / f"pred_mat{suffix}"
    write_table(result.prediction_result.pred_mat, pred_mat_path)
    written["simple_kinase.prediction.pred_mat"] = pred_mat_path

    if result.prediction_result.substrate_list is not None:
        substrate_list_path = prediction_dir / f"substrate_list{suffix}"
        write_table(result.prediction_result.substrate_list, substrate_list_path)
        written["simple_kinase.prediction.substrate_list"] = substrate_list_path

    if result.activity_result is not None:
        activity_dir = workflow_dir / "activity"
        activity_scores_path = activity_dir / f"activity_scores{suffix}"
        write_table(result.activity_result.activity_scores, activity_scores_path)
        written["simple_kinase.activity.activity_scores"] = activity_scores_path

    references_dir = workflow_dir / "references"
    kinase_substrate_map_path = references_dir / f"kinase_substrate_map{suffix}"
    write_table(result.references.kinase_substrate_map, kinase_substrate_map_path)
    written["simple_kinase.references.kinase_substrate_map"] = kinase_substrate_map_path

    site_sequences_path = references_dir / f"site_sequences{suffix}"
    write_table(result.references.site_sequences, site_sequences_path)
    written["simple_kinase.references.site_sequences"] = site_sequences_path

    manifest_path = workflow_dir / "manifest.json"
    _write_manifest(
        manifest_path,
        {
            "reference_organism": result.references.organism.value,
            "activity_enabled": result.activity_result is not None,
            "output_format": output_format,
        },
    )
    written["simple_kinase.manifest"] = manifest_path
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
