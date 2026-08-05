from __future__ import annotations

import importlib

import pytest


def test_science_package_exists() -> None:
    import phospy.science  # noqa: F401


@pytest.mark.parametrize(
    "module_name",
    [
        "phospy.science.activities",
        "phospy.science.datasets",
        "phospy.science.design",
        "phospy.science.differential",
        "phospy.science.evidence",
        "phospy.science.prediction",
        "phospy.science.references",
        "phospy.science.scoring",
        "phospy.science.sequences",
        "phospy.science.signalomes",
        "phospy.science.sites",
        "phospy.science.transformations",
    ],
)
def test_moved_science_packages_import_from_new_paths(module_name: str) -> None:
    importlib.import_module(module_name)


@pytest.mark.parametrize(
    "module_name",
    [
        "phospy.activities",
        "phospy.datasets",
        "phospy.design",
        "phospy.differential",
        "phospy.evidence",
        "phospy.prediction",
        "phospy.references",
        "phospy.scoring",
        "phospy.sequences",
        "phospy.signalomes",
        "phospy.sites",
        "phospy.transformations",
    ],
)
def test_old_root_science_packages_are_removed(module_name: str) -> None:
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module(module_name)


def test_public_api_contract_remains_stable() -> None:
    from phospy.api import (
        AnalysisReadyDatasetBuilder,
        AnalysisReadyPhosphoDataset,
        DifferentialAnalysisWorkflow,
        KinaseWorkflow,
        Organism,
        ReferenceBundle,
        ReferencePreset,
        SignalomeWorkflow,
        TechnicalReplicatePolicy,
    )

    assert AnalysisReadyDatasetBuilder is not None
    assert AnalysisReadyPhosphoDataset is not None
    assert DifferentialAnalysisWorkflow is not None
    assert KinaseWorkflow is not None
    assert Organism is not None
    assert ReferenceBundle is not None
    assert ReferencePreset is not None
    assert SignalomeWorkflow is not None
    assert TechnicalReplicatePolicy is not None


def test_root_convenience_contract_remains_stable() -> None:
    from phospy import (
        AnalysisReadyDatasetBuilder,
        AnalysisReadyPhosphoDataset,
        DifferentialAnalysisWorkflow,
        KinaseWorkflow,
        SignalomeWorkflow,
    )

    assert AnalysisReadyDatasetBuilder is not None
    assert AnalysisReadyPhosphoDataset is not None
    assert DifferentialAnalysisWorkflow is not None
    assert KinaseWorkflow is not None
    assert SignalomeWorkflow is not None


def test_representative_moved_object_ownership() -> None:
    from phospy.api import (
        AnalysisReadyPhosphoDataset,
        ReferenceBundle,
        TechnicalReplicatePolicy,
    )
    from phospy.science.scoring.policy_models import ThresholdMode
    from phospy.science.sites.identifiers import canonicalize_site_identifier
    from phospy.science.transformations.models import IntensityScaleState
    from phospy.science.transformations.scale_state import (
        IntensityScaleState as OwnedIntensityScaleState,
    )

    assert (
        AnalysisReadyPhosphoDataset.__module__
        == "phospy.science.datasets.construction.analysis_ready"
    )
    assert ReferenceBundle.__module__ == "phospy.science.references.models"
    assert (
        TechnicalReplicatePolicy.__module__
        == "phospy.science.differential.policy_models"
    )
    assert ThresholdMode.__module__ == "phospy.science.scoring.policy_models"
    assert canonicalize_site_identifier.__module__ == "phospy.science.sites.identifiers"
    assert (
        IntensityScaleState.__module__ == "phospy.science.transformations.scale_state"
    )
    assert IntensityScaleState is OwnedIntensityScaleState
