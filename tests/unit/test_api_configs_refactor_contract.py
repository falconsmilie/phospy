from __future__ import annotations


def test_public_config_imports_remain_stable() -> None:
    from phospy.api.configs import (
        DatasetIntensityTransformConfig,
        DatasetPreprocessingConfig,
        DatasetTotalProteinCorrectionConfig,
        KinaseActivityConfig,
        KinasePredictionConfig,
        SignalomeConfig,
    )

    assert DatasetPreprocessingConfig is not None
    assert DatasetIntensityTransformConfig is not None
    assert DatasetTotalProteinCorrectionConfig is not None
    assert KinaseActivityConfig is not None
    assert KinasePredictionConfig is not None
    assert SignalomeConfig is not None


def test_config_all_exports_public_api() -> None:
    import phospy.api.configs as configs

    assert "DatasetPreprocessingConfig" in configs.__all__
    assert "KinaseActivityConfig" in configs.__all__
    assert "KinasePredictionConfig" in configs.__all__
    assert "SignalomeConfig" in configs.__all__
    assert "DatasetTotalProteinCorrectionConfig" in configs.__all__


def test_config_refactor_does_not_create_import_cycles() -> None:
    import phospy
    import phospy.api.configs
    from phospy import AnalysisReadyPhosphoDataset
    from phospy.api.configs import DatasetPreprocessingConfig

    assert phospy is not None
    assert phospy.api.configs is not None
    assert AnalysisReadyPhosphoDataset is not None
    assert DatasetPreprocessingConfig is not None
