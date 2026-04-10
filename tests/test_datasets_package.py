from __future__ import annotations

from pathlib import Path

from phospy import AnalysisReadyPhosphoDataset, PhosphoDataset
from phospy.datasets import DatasetLoader, DatasetSchema, DatasetSiteMatrix


def test_public_dataset_types_are_defined_under_datasets_package() -> None:
    assert PhosphoDataset.__module__ == "phospy.datasets.models"
    assert AnalysisReadyPhosphoDataset.__module__ == "phospy.datasets.models"
    assert DatasetLoader.__module__ == "phospy.datasets.loaders"
    assert DatasetSchema.__module__ == "phospy.datasets.schema"
    assert DatasetSiteMatrix.__module__ == "phospy.datasets.builders"


def test_legacy_dataset_flat_modules_have_been_removed() -> None:
    package_root = Path(__file__).resolve().parents[1] / "src" / "phospy"
    assert not (package_root / "dataset.py").exists()
    assert not (package_root / "dataset_loader.py").exists()
    assert not (package_root / "dataset_schema.py").exists()
    assert not (package_root / "dataset_site_matrix.py").exists()
