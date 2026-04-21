from __future__ import annotations

from pathlib import Path

import phospy
import phospy.api as public_api

ROOT = Path(__file__).resolve().parents[2]

USER_FACING_IMPORT_FILES = (
    ROOT / "README.md",
    ROOT / "docs" / "api.md",
    ROOT / "examples" / "dataset_builder_demo.py",
    ROOT / "examples" / "kinase_workflow_demo.py",
    ROOT / "examples" / "signalome_workflow_demo.py",
)

TOP_LEVEL_API_FACADE = {
    "AnalysisReadyDatasetBuilder",
    "AnalysisReadyPhosphoDataset",
    "DatasetBuildRequest",
    "DatasetComparisonBuildingConfig",
    "DatasetMissingDataConfig",
    "DatasetPreprocessingConfig",
    "DatasetSiteMatrixConfig",
    "DatasetTotalProteinCorrectionConfig",
    "KinaseActivityConfig",
    "KinasePredictionConfig",
    "KinaseScoringConfig",
    "KinaseWorkflow",
    "KinaseWorkflowRequest",
    "KinaseWorkflowResult",
    "Organism",
    "ReferenceBundle",
    "ReferencePreset",
    "SignalomeConfig",
    "SignalomeWorkflow",
    "SignalomeWorkflowRequest",
    "SignalomeWorkflowResult",
}

NON_FACADE_API_TYPES = {
    "KinaseActivityResult",
    "KinasePredictionResult",
    "KinaseScoringResult",
}


def test_top_level_facade_re_exports_curated_canonical_api_types() -> None:
    assert TOP_LEVEL_API_FACADE.issubset(set(public_api.__all__))
    assert TOP_LEVEL_API_FACADE.issubset(set(phospy.__all__))
    for exported in TOP_LEVEL_API_FACADE:
        assert getattr(phospy, exported) is getattr(public_api, exported)
    for exported in NON_FACADE_API_TYPES:
        assert exported in public_api.__all__
        assert exported not in phospy.__all__


def test_readme_and_api_guide_document_import_contract() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    api_guide = (ROOT / "docs" / "api.md").read_text(encoding="utf-8")

    canonical_phrase = (
        "`phospy.api` is the canonical namespace where public API types are defined"
    )
    primary_route_phrase = "Top-level `phospy` is the primary supported import route"

    assert canonical_phrase in readme
    assert canonical_phrase in api_guide
    assert primary_route_phrase in readme
    assert primary_route_phrase in api_guide


def test_user_facing_guides_and_examples_use_top_level_import_route() -> None:
    for file_path in USER_FACING_IMPORT_FILES:
        source = file_path.read_text(encoding="utf-8")
        assert "from phospy.api import" not in source
        assert "import phospy.api" not in source
