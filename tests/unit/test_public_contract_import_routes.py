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

TOP_LEVEL_CONVENIENCE_SURFACE = {
    "AnalysisReadyDatasetBuilder",
    "AnalysisReadyPhosphoDataset",
    "KinaseWorkflow",
    "SignalomeWorkflow",
}

API_ONLY_CONTRACT_TYPES = {
    "DatasetBuildRequest",
    "DatasetPreprocessingReport",
    "DatasetPreprocessingConfig",
    "DatasetSiteMatrixConfig",
    "KinaseWorkflowRequest",
    "KinaseWorkflowResult",
    "Organism",
    "ReferenceBundle",
    "ReferencePreset",
    "SignalomeConfig",
    "SignalomeWorkflowRequest",
    "SignalomeWorkflowResult",
    "PhosPyValidationError",
    "UnsupportedInputFormatError",
}


def test_top_level_package_exports_only_curated_convenience_surface() -> None:
    assert set(phospy.__all__) == TOP_LEVEL_CONVENIENCE_SURFACE
    assert TOP_LEVEL_CONVENIENCE_SURFACE.issubset(set(public_api.__all__))
    for exported in TOP_LEVEL_CONVENIENCE_SURFACE:
        assert getattr(phospy, exported) is getattr(public_api, exported)
    for exported in API_ONLY_CONTRACT_TYPES:
        assert exported in public_api.__all__
        assert exported not in phospy.__all__


def test_readme_and_api_guide_document_import_contract() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    api_guide = (ROOT / "docs" / "api.md").read_text(encoding="utf-8")

    canonical_phrase = (
        "`phospy.api` is the canonical namespace where public API types are defined"
    )
    convenience_phrase = "top-level `phospy` is a curated convenience surface"

    assert canonical_phrase in readme
    assert canonical_phrase in api_guide
    assert convenience_phrase in readme
    assert convenience_phrase in api_guide


def test_user_facing_guides_and_examples_use_phospy_api_for_contract_types() -> None:
    for file_path in USER_FACING_IMPORT_FILES:
        source = file_path.read_text(encoding="utf-8")
        assert "from phospy.api import" in source
        assert "from phospy import DatasetBuildRequest" not in source
        assert "from phospy import KinaseWorkflowRequest" not in source
        assert "from phospy import SignalomeWorkflowRequest" not in source


def test_public_docs_present_analysis_ready_builder_lane_as_missing_value_free() -> (
    None
):
    api_guide = (ROOT / "docs" / "api.md").read_text(encoding="utf-8")
    validation_guide = (ROOT / "docs" / "validation.md").read_text(encoding="utf-8")
    quickstart = (
        ROOT / "docs" / "getting-started" / "quickstart-first-workflow.md"
    ).read_text(encoding="utf-8")

    assert "missing-value-free `AnalysisReadyPhosphoDataset`" in api_guide
    assert (
        "`AnalysisReadyPhosphoDataset` itself is strict, missing-value-free"
        in validation_guide
    )
    assert "build an analysis-ready, missing-value-free dataset" in quickstart
