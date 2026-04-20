from __future__ import annotations

import importlib.util
import inspect
from dataclasses import fields
from pathlib import Path
from typing import get_type_hints

import pytest
import tomllib

import phospy
import phospy.io as phospy_io
from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
)
from phospy.datasets.models import AnalysisReadyPhosphoDataset

ROOT = Path(__file__).resolve().parents[2]


def _public_methods(cls: type[object]) -> set[str]:
    return {
        name
        for name, value in inspect.getmembers(cls)
        if callable(value) and not name.startswith("_")
    }


def test_packaging_contract_excludes_legacy_package_namespace() -> None:
    pyproject = ROOT / "pyproject.toml"
    config = tomllib.loads(pyproject.read_text(encoding="utf-8"))
    package_find = config["tool"]["setuptools"]["packages"]["find"]
    assert set(package_find["include"]) == {"phospy", "phospy.*"}
    assert {"phospy_legacy", "phospy_legacy.*"}.issubset(set(package_find["exclude"]))


def test_source_tree_does_not_expose_legacy_namespace() -> None:
    assert importlib.util.find_spec("phospy_legacy") is None


def test_public_dataset_ingestion_story_is_builder_only() -> None:
    dataset_exports = {
        name for name in phospy.__all__ if "Dataset" in name or name.endswith("Builder")
    }
    dataset_exports = {name for name in dataset_exports if not name.endswith("Error")}
    assert dataset_exports == {
        "AnalysisReadyDatasetBuilder",
        "AnalysisReadyPhosphoDataset",
        "DatasetComparisonBuildingConfig",
        "DatasetBuildRequest",
        "DatasetMissingDataConfig",
        "DatasetPreprocessingConfig",
        "DatasetSiteMatrixConfig",
        "DatasetTotalProteinCorrectionConfig",
    }
    assert not hasattr(phospy, "build_dataset_from_files")
    assert not hasattr(phospy_io, "build_dataset_from_files")


def test_builder_exposes_only_run_request_contract() -> None:
    assert _public_methods(AnalysisReadyDatasetBuilder) == {"run"}
    assert not hasattr(AnalysisReadyDatasetBuilder, "execute")
    hints = get_type_hints(AnalysisReadyDatasetBuilder.run)
    assert hints["request"] is DatasetBuildRequest
    assert hints["return"] is AnalysisReadyPhosphoDataset


def test_dataset_build_request_excludes_user_declared_transformation_state() -> None:
    request_fields = {field.name for field in fields(DatasetBuildRequest)}
    assert "transformation_state" not in request_fields


def test_dataset_build_request_rejects_user_declared_transformation_state() -> None:
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        DatasetBuildRequest(
            phospho=object(),
            site_metadata=object(),
            transformation_state=object(),  # type: ignore[call-arg]
        )


def test_dataset_preprocessing_config_is_top_level_public_type() -> None:
    assert "DatasetPreprocessingConfig" in phospy.__all__
    assert DatasetPreprocessingConfig().missing_data.policy == "forbid"
