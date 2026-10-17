from __future__ import annotations

import importlib

import pytest


def test_science_package_is_retired() -> None:
    with pytest.raises(ModuleNotFoundError) as captured:
        importlib.import_module("phospy.science")

    assert captured.value.name == "phospy.science"


@pytest.mark.parametrize(
    "module_name",
    [
        "phospy.activities",
        "phospy.kinase.activity",
        "phospy.differential.design",
        "phospy.differential",
        "phospy.evidence",
        "phospy.kinase.prediction",
        "phospy.references",
        "phospy.kinase.scoring",
        "phospy.signalome",
    ],
)
def test_moved_science_packages_import_from_new_paths(module_name: str) -> None:
    importlib.import_module(module_name)


@pytest.mark.parametrize(
    "module_name",
    [
        "phospy.design",
        "phospy.prediction",
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


def test_phosphosites_owner_package_is_importable() -> None:
    import phospy.phosphosites  # noqa: F401


def test_datasets_owner_package_is_importable() -> None:
    import phospy.datasets  # noqa: F401


def test_quantitative_and_statistics_owner_packages_are_importable() -> None:
    import phospy.quantitative  # noqa: F401
    import phospy.statistics  # noqa: F401


def test_public_api_contract_remains_stable() -> None:
    from phospy.advanced import TechnicalReplicatePolicy
    from phospy.api import (
        AnalysisReadyDatasetBuilder,
        AnalysisReadyPhosphoDataset,
        DifferentialAnalysisWorkflow,
        KinaseWorkflow,
        Organism,
        ReferenceBundle,
        ReferencePreset,
        SignalomeWorkflow,
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
        KinaseWorkflow,
        SignalomeWorkflow,
    )
    from phospy.api import (
        AnalysisReadyPhosphoDataset,
        DifferentialAnalysisWorkflow,
    )

    assert AnalysisReadyDatasetBuilder is not None
    assert AnalysisReadyPhosphoDataset is not None
    assert DifferentialAnalysisWorkflow is not None
    assert KinaseWorkflow is not None
    assert SignalomeWorkflow is not None


def test_representative_moved_object_ownership() -> None:
    from phospy.advanced import TechnicalReplicatePolicy
    from phospy.api import (
        AnalysisReadyPhosphoDataset,
        ReferenceBundle,
    )
    from phospy.kinase.policies import ThresholdMode
    from phospy.phosphosites.identifiers import canonicalize_site_identifier
    from phospy.quantitative.models import IntensityScaleState
    from phospy.quantitative.scale_state import (
        IntensityScaleState as OwnedIntensityScaleState,
    )

    assert AnalysisReadyPhosphoDataset.__module__ == "phospy.datasets.models"
    assert ReferenceBundle.__module__ == "phospy.references.models"
    assert TechnicalReplicatePolicy.__module__ == "phospy.differential.policies"
    assert ThresholdMode.__module__ == "phospy.kinase.policies"
    assert canonicalize_site_identifier.__module__ == "phospy.phosphosites.identifiers"
    assert IntensityScaleState.__module__ == "phospy.quantitative.scale_state"
    assert IntensityScaleState is OwnedIntensityScaleState
