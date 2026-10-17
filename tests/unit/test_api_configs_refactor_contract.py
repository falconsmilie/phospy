from __future__ import annotations


def test_stable_config_imports_remain_from_api_configs() -> None:
    from phospy.api import (
        DatasetLocalisationConfig,
        DatasetPreprocessingConfig,
        EnrichmentConfig,
    )

    assert DatasetLocalisationConfig is not None
    assert DatasetPreprocessingConfig is not None
    assert EnrichmentConfig is not None


def test_advanced_config_imports_have_advanced_owner() -> None:
    from phospy.advanced import (
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
    import phospy.api as public_api

    assert {
        "DatasetLocalisationConfig",
        "DatasetPreprocessingConfig",
        "EnrichmentConfig",
    } <= set(public_api.__all__)
    assert not hasattr(public_api, "configs")


def test_config_refactor_does_not_create_import_cycles() -> None:
    import phospy
    import phospy.api
    from phospy.api import (
        AnalysisReadyPhosphoDataset,
        DatasetPreprocessingConfig,
    )

    assert phospy is not None
    assert "__getattr__" not in vars(phospy.api)
    assert AnalysisReadyPhosphoDataset is not None
    assert DatasetPreprocessingConfig is not None
