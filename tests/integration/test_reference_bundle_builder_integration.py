from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
from phospy.api import (
    DatasetBuildRequest,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundleBuilder,
    ReferenceBundleBuildRequest,
)

pytestmark = pytest.mark.integration


def test_local_mouse_reference_bundle_builder_output_runs_in_kinase_workflow(
    tmp_path: Path,
) -> None:
    kinase_path = tmp_path / "mouse_kinase.csv"
    sequence_path = tmp_path / "mouse_sequences.csv"
    pd.DataFrame(
        {
            "kinase": ["Map2k1", "Map2k1"],
            "site_id": ["Mapk1;S123;", "Mapk1;Y185;"],
            "organism": ["mouse", "mouse"],
        }
    ).to_csv(kinase_path, index=False)
    pd.DataFrame(
        {
            "site_id": ["Mapk1;S123;", "Mapk1;Y185;"],
            "site_sequence": [_window("S"), _window("Y")],
            "gene_symbol": ["Mapk1", "Mapk1"],
            "protein_id": ["MAPK1_MOUSE", "MAPK1_MOUSE"],
            "organism": ["mouse", "mouse"],
        }
    ).to_csv(sequence_path, index=False)

    references = ReferenceBundleBuilder().run(
        ReferenceBundleBuildRequest(
            organism=Organism.MOUSE,
            kinase_substrate_path=kinase_path,
            site_sequence_path=sequence_path,
            source_name="synthetic mouse local source",
            source_version="integration-v1",
            retrieved_at="2026-06-11",
            license="synthetic test license",
            redistribution_status="redistributable synthetic fixture",
            identifier_namespace="display_id (GENE_SYMBOL;RESIDUE;)",
        )
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=pd.DataFrame(
                {
                    "sample_a": [1.0, 0.7],
                    "sample_b": [1.2, 0.9],
                },
                index=["MAPK1;S123;", "MAPK1;Y185;"],
            ),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK1", "MAPK1"],
                    "site": ["S123", "Y185"],
                    "site_sequence": [_window("S"), _window("Y")],
                    "display_id": ["MAPK1;S123;", "MAPK1;Y185;"],
                    "organism": ["mouse", "mouse"],
                    "protein_namespace": ["protein_id", "protein_id"],
                    "protein_identifier": ["MAPK1_MOUSE", "MAPK1_MOUSE"],
                    "protein_id": ["MAPK1_MOUSE", "MAPK1_MOUSE"],
                    "localisation_confidence": [0.95, 0.95],
                },
                index=["MAPK1;S123;", "MAPK1;Y185;"],
            ),
            organism=Organism.MOUSE,
            input_intensity_scale="linear",
        )
    )

    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=1,
                deterministic_max_selected_kinases=1,
                adaptive_ensemble_runs=1,
            ),
            activity_config=None,
        )
    )

    assert result.references is references
    assert result.references.provenance is not None
    assert result.references.provenance.source_type == "local"
    assert not result.scoring_result.profile_scores.empty
    assert "MAP2K1" in result.prediction_result.pred_mat.columns


def test_kinase_workflow_code_does_not_know_local_reference_file_loading() -> None:
    workflow_root = (
        Path(__file__).resolve().parents[2] / "src" / "phospy" / "workflows" / "kinase"
    )
    searched = "\n".join(
        path.read_text(encoding="utf-8")
        for path in workflow_root.glob("*.py")
        if path.name != "__init__.py"
    )

    assert "ReferenceBundleBuilder" not in searched
    assert "ReferenceBundleBuildRequest" not in searched
    assert "kinase_substrate_path" not in searched
    assert "site_sequence_path" not in searched


def _window(center: str) -> str:
    return f"{'A' * 15}{center}{'A' * 15}"
