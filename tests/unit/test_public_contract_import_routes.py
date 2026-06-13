from __future__ import annotations

import pytest

import phospy
import phospy.api as public_api
import phospy.science.signalomes as public_signalomes
import phospy.science.signalomes.clustering as signalome_clustering

TOP_LEVEL_CONVENIENCE_SURFACE = {
    "AnalysisReadyDatasetBuilder",
    "AnalysisReadyPhosphoDataset",
    "DifferentialAnalysisWorkflow",
    "KinaseWorkflow",
    "SignalomeWorkflow",
}

API_ONLY_CONTRACT_TYPES = {
    "DatasetBuildRequest",
    "DatasetBatchCorrectionConfig",
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


def test_differential_analysis_is_not_supported_from_phospy_api_namespace() -> None:
    with pytest.raises(ImportError):
        exec("from phospy.api import DifferentialAnalysis")


def test_differential_workflow_supported_import_routes() -> None:
    namespace: dict[str, object] = {}

    exec("from phospy import DifferentialAnalysisWorkflow", namespace)
    exec("from phospy.api import DifferentialAnalysisRequest", namespace)

    assert "DifferentialAnalysisWorkflow" in namespace
    assert "DifferentialAnalysisRequest" in namespace


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


def test_split_module_import_routes_remain_backward_compatible() -> None:
    from phospy.io.bundles._signalome.compatibility import (
        signalome_config_from_payload,
    )
    from phospy.science.prediction.motif_scoring import score_phosphosite_motifs
    from phospy.tables.signalome import SignalomeSiteContext

    assert SignalomeSiteContext.__name__ == "SignalomeSiteContext"
    assert callable(signalome_config_from_payload)
    assert callable(score_phosphosite_motifs)
