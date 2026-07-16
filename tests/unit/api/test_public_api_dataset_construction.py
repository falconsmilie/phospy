from __future__ import annotations

from pathlib import Path
from typing import get_type_hints

import phospy.api as public_api
from phospy.api import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DatasetBuildRequest,
)

ROOT = Path(__file__).resolve().parents[3]


def test_public_api_documents_builder_as_supported_construction_path() -> None:
    assert "AnalysisReadyDatasetBuilder" in public_api.__all__
    assert "AnalysisReadyPhosphoDataset" in public_api.__all__
    assert get_type_hints(AnalysisReadyDatasetBuilder.run)["request"] is (
        DatasetBuildRequest
    )

    builder_doc = AnalysisReadyDatasetBuilder.__doc__
    assert builder_doc is not None
    assert "Supported public path" in builder_doc
    assert "construction provenance" in builder_doc


def test_public_api_marks_direct_dataset_construction_advanced_trusted() -> None:
    model_doc = AnalysisReadyPhosphoDataset.__doc__
    factory_doc = AnalysisReadyPhosphoDataset.from_trusted_tables.__doc__

    assert model_doc is not None
    assert factory_doc is not None
    assert "trusted advanced/internal use" in model_doc
    assert "Ordinary users" in model_doc
    assert "AnalysisReadyDatasetBuilder.run" in model_doc
    assert "from_trusted_tables" in model_doc
    assert "not the primary advanced construction API" in model_doc
    assert "typed evidence or an explicit" in model_doc
    assert "localisation" in model_doc
    assert "cannot prove" in model_doc
    assert "biological correctness" in model_doc
    assert "minimal" in model_doc
    assert "direct-construction provenance marker" in model_doc
    assert "same structural invariants as direct construction" in factory_doc
    assert "site_sequence" in factory_doc
    assert "source, policy" in factory_doc
    assert "threshold" in factory_doc
    assert "cannot prove" in factory_doc
    assert "biological correctness" in factory_doc


def test_public_docs_examples_use_builder_path() -> None:
    documentation_paths = (
        ROOT / "README.md",
        ROOT / "docs" / "quickstart.md",
        ROOT / "docs" / "api" / "guide.md",
        ROOT / "docs" / "api" / "dataset-build-workflow.md",
        ROOT / "docs" / "api" / "dataset-builders.md",
    )
    documentation = "\n".join(
        path.read_text(encoding="utf-8") for path in documentation_paths
    )

    assert "AnalysisReadyDatasetBuilder().run(" in documentation
    assert "AnalysisReadyPhosphoDataset(" not in documentation
    assert (
        "from phospy import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset"
    ) not in documentation
    assert "advanced/trusted" in documentation
    assert "AnalysisReadyPhosphoDataset.from_trusted_tables" in documentation
    assert "minimal direct-construction provenance marker" in documentation
