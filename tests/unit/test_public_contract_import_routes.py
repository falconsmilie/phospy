from __future__ import annotations

import pytest

import phospy
import phospy.api as public_api
import phospy.science.datasets.preprocessing.stage_registry as preprocessing_stage_registry
import phospy.science.prediction.scoring as prediction_scoring
import phospy.science.signalomes as public_signalomes
import phospy.science.signalomes.clustering as signalome_clustering
import phospy.science.signalomes.clustering.exact_python as exact_clustering
import phospy.science.signalomes.clustering.protocol as clustering_protocol
from phospy.science.datasets.preprocessing.stage_contract import (
    PreprocessingStageContract,
)

TOP_LEVEL_CONVENIENCE_SURFACE = {
    "AnalysisReadyDatasetBuilder",
    "AnalysisReadyPhosphoDataset",
    "DifferentialAnalysisWorkflow",
    "KinaseWorkflow",
    "SignalomeWorkflow",
}

API_ONLY_CONTRACT_TYPES = {
    "BatchCorrectionDiagnostics",
    "BatchCorrectionPolicy",
    "BatchCorrectionReport",
    "ControlSiteAnnotation",
    "ControlSiteSet",
    "ControlSiteSourceMetadata",
    "ControlSiteStatus",
    "CorrectionMissingnessPolicy",
    "DatasetBuildRequest",
    "DatasetBatchCorrectionConfig",
    "DatasetGroupCoverageFilterConfig",
    "DatasetPreprocessingReport",
    "DatasetPreprocessingConfig",
    "DatasetSiteMatrixConfig",
    "EnrichmentConfig",
    "EnrichmentIdentifierKind",
    "EnrichmentWorkflowRequest",
    "EnrichmentWorkflowResult",
    "GeneSetCollection",
    "IMPORTER_QUALITY_STATUS_NOT_APPLICABLE",
    "IMPORTER_QUALITY_STATUS_NOT_REPORTED",
    "IMPORTER_QUALITY_STATUS_REPORTED",
    "ImporterDetectedIntensityColumn",
    "ImporterDuplicateKeySummary",
    "ImporterFlaggedRowSummary",
    "ImporterLocalisationConfidenceSummary",
    "ImporterMissingIntensitySummary",
    "ImporterQualityCount",
    "ImporterQualityReport",
    "ImporterQualityStatus",
    "KinaseWorkflowRequest",
    "KinaseWorkflowResult",
    "Organism",
    "ProteinAwarePreparationResult",
    "ProteinAwarePreparationReport",
    "ReferenceBundle",
    "ReferencePreset",
    "PtmSetCollection",
    "SignalomeConfig",
    "SignalomeWorkflowRequest",
    "SignalomeWorkflowResult",
    "SpsRuvBatchCorrectionConfig",
    "TemporaryImputationPolicy",
    "PhosPyValidationError",
    "UnsupportedInputFormatError",
}

SEMI_PUBLIC_SCIENCE_IMPORTS = {
    (
        "phospy.science.datasets.preprocessing.stage_registry",
        "PreprocessingStageMetadata",
    ),
    ("phospy.science.signalomes.clustering.protocol", "ClusterTreeEngine"),
    ("phospy.science.signalomes.clustering.protocol", "SignalomeClusteringEngine"),
    (
        "phospy.science.signalomes.clustering.exact_python",
        "ExactPythonClusteringBackend",
    ),
    (
        "phospy.science.prediction.scoring",
        "fuse_profile_and_motif_scores_by_rank_weight",
    ),
}

PRIVATE_HELPERS_NOT_PUBLIC_EXPORTS = {
    "phospy.science.signalomes.clustering.exact_python": {
        "_WardClusterTree",
        "_build_cluster_tree",
    },
    "phospy.science.prediction.scoring": {
        "_fuse_profile_and_motif_scores_by_rank_weight",
        "_profile_only_weights",
    },
}
PRIVATE_HELPER_NAMES = {
    symbol_name
    for symbol_names in PRIVATE_HELPERS_NOT_PUBLIC_EXPORTS.values()
    for symbol_name in symbol_names
}


def _assert_from_import_fails(module_name: str, symbol_name: str) -> None:
    with pytest.raises(ImportError):
        exec(f"from {module_name} import {symbol_name}", {})


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


@pytest.mark.parametrize(
    "symbol_name",
    (
        "DATASET_BATCH_CORRECTION_METHOD_CONTROL_SITE_RUV_STYLE",
        "DATASET_BATCH_CORRECTION_METHOD_RUV_III_STYLE",
    ),
)
def test_unsupported_sps_ruv_style_method_constants_are_not_public_api(
    symbol_name: str,
) -> None:
    for module_name in (
        "phospy.api.configs",
        "phospy.api.configs.preprocessing",
        "phospy.api.configs.preprocessing.batch_correction",
    ):
        module_namespace: dict[str, object] = {}
        wildcard_namespace: dict[str, object] = {}

        exec(f"import {module_name} as module", module_namespace)
        module = module_namespace["module"]
        exec(f"from {module_name} import *", wildcard_namespace)

        assert symbol_name not in module.__all__
        assert symbol_name not in wildcard_namespace
        assert not hasattr(module, symbol_name)
        _assert_from_import_fails(module_name, symbol_name)


def test_deprecated_differential_analysis_shell_warns_and_delegates() -> None:
    from phospy.science.differential.public import DifferentialAnalysis

    class _Workflow:
        def run(self, request: object) -> str:
            return f"delegated:{request!r}"

    with pytest.warns(
        DeprecationWarning,
        match="DifferentialAnalysis is deprecated",
    ):
        shell = DifferentialAnalysis(workflow=_Workflow())  # type: ignore[arg-type]

    assert shell.run("request") == "delegated:'request'"  # type: ignore[arg-type]


def test_differential_workflow_supported_import_routes() -> None:
    namespace: dict[str, object] = {}

    exec("from phospy import DifferentialAnalysisWorkflow", namespace)
    exec("from phospy.api import DifferentialAnalysisRequest", namespace)

    assert "DifferentialAnalysisWorkflow" in namespace
    assert "DifferentialAnalysisRequest" in namespace


@pytest.mark.parametrize(("module_name", "symbol_name"), SEMI_PUBLIC_SCIENCE_IMPORTS)
def test_supported_semi_public_science_import_routes(
    module_name: str,
    symbol_name: str,
) -> None:
    namespace: dict[str, object] = {}

    exec(f"from {module_name} import {symbol_name}", namespace)

    assert symbol_name in namespace


def test_preprocessing_stage_metadata_route_remains_compatibility_alias() -> None:
    assert (
        preprocessing_stage_registry.PreprocessingStageMetadata
        is PreprocessingStageContract
    )
    assert "PreprocessingStageMetadata" in preprocessing_stage_registry.__all__


def test_clustering_protocol_route_exports_only_protocol_contracts() -> None:
    assert set(clustering_protocol.__all__) == {
        "ClusterTreeEngine",
        "SignalomeClusteringEngine",
    }


def test_rank_fusion_helper_is_semi_public_but_not_api_public() -> None:
    assert "fuse_profile_and_motif_scores_by_rank_weight" in prediction_scoring.__all__
    assert callable(prediction_scoring.fuse_profile_and_motif_scores_by_rank_weight)
    assert "fuse_profile_and_motif_scores_by_rank_weight" not in public_api.__all__
    assert not hasattr(public_api, "fuse_profile_and_motif_scores_by_rank_weight")
    assert not hasattr(phospy, "fuse_profile_and_motif_scores_by_rank_weight")


def test_semi_public_science_helpers_are_not_promoted_to_root_or_api() -> None:
    unsupported_public_promotions = {
        "PreprocessingStageMetadata",
        "ClusterTreeEngine",
        "SignalomeClusteringEngine",
        "ExactPythonClusteringBackend",
        "fuse_profile_and_motif_scores_by_rank_weight",
    }
    for symbol_name in unsupported_public_promotions:
        _assert_from_import_fails("phospy", symbol_name)
        _assert_from_import_fails("phospy.api", symbol_name)


@pytest.mark.parametrize("symbol_name", sorted(PRIVATE_HELPER_NAMES))
def test_private_helpers_are_not_supported_from_root_or_api(symbol_name: str) -> None:
    _assert_from_import_fails("phospy", symbol_name)
    _assert_from_import_fails("phospy.api", symbol_name)


@pytest.mark.parametrize(
    ("module_name", "symbol_names"),
    sorted(PRIVATE_HELPERS_NOT_PUBLIC_EXPORTS.items()),
)
def test_private_helpers_are_rejected_by_public_export_routes(
    module_name: str,
    symbol_names: set[str],
) -> None:
    module_namespace: dict[str, object] = {}
    wildcard_namespace: dict[str, object] = {}

    exec(f"import {module_name} as module", module_namespace)
    module = module_namespace["module"]
    exec(f"from {module_name} import *", wildcard_namespace)

    for symbol_name in symbol_names:
        assert symbol_name not in module.__all__
        assert symbol_name not in wildcard_namespace


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


def test_exact_python_compatibility_facade_does_not_export_internal_aliases() -> None:
    assert "ExactPythonClusteringBackend" in exact_clustering.__all__
    assert all(not symbol.startswith("_") for symbol in exact_clustering.__all__)


def test_removed_signalome_bundle_compatibility_import_fails() -> None:
    with pytest.raises(ModuleNotFoundError):
        exec("import phospy.io.bundles._signalome.compatibility", {})


def test_split_module_import_routes_remain_supported() -> None:
    from phospy.science.prediction.motif_scoring import score_phosphosite_motifs
    from phospy.tables.signalome import SignalomeSiteContext

    assert SignalomeSiteContext.__name__ == "SignalomeSiteContext"
    assert callable(score_phosphosite_motifs)
