from __future__ import annotations

from pathlib import Path

import phospy
import phospy.api as public_api
import phospy.signalomes as public_signalomes
import phospy.signalomes.clustering as signalome_clustering

ROOT = Path(__file__).resolve().parents[2]

USER_FACING_IMPORT_FILES = (
    ROOT / "README.md",
    ROOT / "docs" / "guide.md",
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

QUICKSTART_CANDIDATES = (
    ROOT / "docs" / "getting-started" / "quickstart-first-workflow.md",
    ROOT / "docs" / "quickstart.md",
)


def _read_first_existing(*paths: Path) -> str:
    for path in paths:
        if path.exists():
            return path.read_text(encoding="utf-8")
    searched = ", ".join(str(path) for path in paths)
    raise FileNotFoundError(f"None of the expected files exist: {searched}")


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
    api_guide = (ROOT / "docs" / "guide.md").read_text(encoding="utf-8")

    assert "Use `phospy.api`" in readme
    assert "Use `phospy.api`" in api_guide
    assert "Use top-level `phospy`" in readme
    assert "Use top-level `phospy`" in api_guide


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
    api_guide = (ROOT / "docs" / "guide.md").read_text(encoding="utf-8")
    validation_guide = (ROOT / "docs" / "validation.md").read_text(encoding="utf-8")
    quickstart = _read_first_existing(*QUICKSTART_CANDIDATES)

    assert "must be missing-value-free" in api_guide
    assert "Analysis-Ready Dataset Boundary" in validation_guide
    assert "missing-value-free" in validation_guide
    assert "workflow" in quickstart.lower()
    assert "analysis" in quickstart.lower()


def test_build_cluster_tree_not_exported_from_signalome_public_api() -> None:
    assert "build_cluster_tree" not in public_signalomes.__all__
    assert not hasattr(public_signalomes, "build_cluster_tree")
    assert "build_cluster_tree" not in signalome_clustering.__all__
    assert not hasattr(signalome_clustering, "build_cluster_tree")


def test_public_signalome_clustering_exports_do_not_expose_exact_tree_builders() -> (
    None
):
    public_tree_builders = {
        name
        for name in signalome_clustering.__all__
        if "tree" in name and name.startswith("build")
    }
    assert public_tree_builders == {"build_cluster_labels_from_tree"}


def test_public_signalome_clustering_exports_do_not_expose_internal_helpers() -> None:
    internal_exports = {
        name for name in signalome_clustering.__all__ if name.startswith("_")
    }
    assert internal_exports == set()


def test_public_signalome_docs_do_not_reference_raw_build_cluster_tree() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    docs_markdown_files = sorted((ROOT / "docs").rglob("*.md"))
    docs_text = "\n".join(
        file_path.read_text(encoding="utf-8") for file_path in docs_markdown_files
    )
    assert "build_cluster_tree" not in readme
    assert "build_cluster_tree" not in docs_text


def test_split_module_import_routes_remain_backward_compatible() -> None:
    from phospy.io.bundles._signalome.compatibility import (
        signalome_config_from_payload,
    )
    from phospy.prediction.motif_scoring import score_phosphosite_motifs
    from phospy.tables.signalome import SignalomeSiteContext

    assert SignalomeSiteContext.__name__ == "SignalomeSiteContext"
    assert callable(signalome_config_from_payload)
    assert callable(score_phosphosite_motifs)
