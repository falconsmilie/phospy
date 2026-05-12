from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
README = ROOT / "README.md"
API_GUIDE = ROOT / "docs" / "api" / "guide.md"

README_IMPORT_SNIPPET = """from phospy import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset
from phospy import DifferentialAnalysisWorkflow, KinaseWorkflow, SignalomeWorkflow
"""

API_GUIDE_IMPORT_SNIPPET = """from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DifferentialAnalysisWorkflow,
    KinaseWorkflow,
    SignalomeWorkflow,
)
"""

API_GUIDE_API_IMPORT_SNIPPET = """from phospy.api import (
    DatasetBuildRequest,
    ExperimentalDesign,
    Contrast,
    SampleDesignRecord,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferencePreset,
    SignalomeConfig,
    SignalomeWorkflowRequest,
    UnsupportedInputFormatError,
    WorkflowValidationError,
)
"""


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_documented_public_imports_are_importable() -> None:
    namespace: dict[str, object] = {}

    exec(README_IMPORT_SNIPPET, namespace)
    exec(API_GUIDE_IMPORT_SNIPPET, namespace)
    exec(API_GUIDE_API_IMPORT_SNIPPET, namespace)

    assert "AnalysisReadyDatasetBuilder" in namespace
    assert "AnalysisReadyPhosphoDataset" in namespace
    assert "DifferentialAnalysisWorkflow" in namespace
    assert "KinaseWorkflow" in namespace
    assert "SignalomeWorkflow" in namespace
    assert "DifferentialAnalysisRequest" in namespace


def test_readme_differential_import_example_matches_supported_route() -> None:
    source = _read(README)

    assert README_IMPORT_SNIPPET in source
    assert "DifferentialAnalysis().run(" not in source


def test_api_guide_differential_import_examples_match_supported_route() -> None:
    source = _read(API_GUIDE)

    assert API_GUIDE_IMPORT_SNIPPET in source
    assert API_GUIDE_API_IMPORT_SNIPPET in source
    assert "DifferentialAnalysisWorkflow().run(differential_request)" in source
    assert "from phospy import DifferentialAnalysis," not in source
