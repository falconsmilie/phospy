from __future__ import annotations


def test_stable_config_imports_remain_from_api_configs() -> None:
    from phospy.api.configs import (
        DatasetLocalisationConfig,
        DatasetPreprocessingConfig,
        EnrichmentConfig,
    )

    assert DatasetLocalisationConfig is not None
    assert DatasetPreprocessingConfig is not None
    assert EnrichmentConfig is not None


def test_advanced_config_imports_have_advanced_owner() -> None:
    from phospy.advanced.configs import (
        DatasetBatchCorrectionConfig,
        DatasetGroupCoverageFilterConfig,
        DatasetIntensityTransformConfig,
        DatasetTotalProteinCorrectionConfig,
        KinaseActivityConfig,
        KinaseAttritionPolicy,
        KinasePredictionConfig,
        ProfileSelfInclusionPolicy,
        ReferenceContextCompatibilityPolicy,
        SignalomeConfig,
    )

    assert DatasetBatchCorrectionConfig is not None
    assert DatasetGroupCoverageFilterConfig is not None
    assert DatasetIntensityTransformConfig is not None
    assert DatasetTotalProteinCorrectionConfig is not None
    assert KinaseActivityConfig is not None
    assert KinaseAttritionPolicy is not None
    assert KinasePredictionConfig is not None
    assert ProfileSelfInclusionPolicy is not None
    assert ReferenceContextCompatibilityPolicy is not None
    assert SignalomeConfig is not None


def test_config_all_exports_public_api() -> None:
    import phospy.api.configs as configs

    assert set(configs.__all__) == {
        "DatasetLocalisationConfig",
        "DatasetPreprocessingConfig",
        "EnrichmentConfig",
    }
    assert "DatasetPreprocessingConfig" in configs.__all__
    assert "EnrichmentConfig" in configs.__all__


def test_config_refactor_does_not_create_import_cycles() -> None:
    import phospy
    import phospy.api.configs
    from phospy import AnalysisReadyPhosphoDataset
    from phospy.api.configs import DatasetPreprocessingConfig

    assert phospy is not None
    assert phospy.api.configs is not None
    assert AnalysisReadyPhosphoDataset is not None
    assert DatasetPreprocessingConfig is not None
