from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
README = ROOT / "README.md"
API_GUIDE = ROOT / "docs" / "api" / "guide.md"
DATASET_WORKFLOW_DOC = ROOT / "docs" / "api" / "dataset-build-workflow.md"
WORKFLOW_DOCS_DIR = ROOT / "docs" / "api"
DIFFERENTIAL_WORKFLOW_DOC = WORKFLOW_DOCS_DIR / "differential-analysis.md"
ENRICHMENT_WORKFLOW_DOC = WORKFLOW_DOCS_DIR / "enrichment.md"
KINASE_WORKFLOW_DOC = WORKFLOW_DOCS_DIR / "kinase.md"
SIGNALOME_WORKFLOW_DOC = WORKFLOW_DOCS_DIR / "signalome.md"

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


def _normalise_whitespace(text: str) -> str:
    return " ".join(text.split())


def _iter_documentation_statements(source: str) -> tuple[str, ...]:
    statements: list[str] = []
    for block in re.split(r"\n\s*\n", source):
        stripped_block = block.strip()
        if not stripped_block:
            continue
        chunks = (
            stripped_block.splitlines()
            if stripped_block.startswith("|")
            else re.split(r"\n(?=\s*[-*]\s+)", stripped_block)
        )
        for chunk in chunks:
            normalised = _normalise_whitespace(chunk)
            if not normalised or normalised.startswith("```"):
                continue
            statements.extend(
                sentence
                for sentence in re.split(r"(?<=[.!?])\s+", normalised)
                if sentence
            )
    return tuple(statements)


def _assert_statement_contains_all(
    source: str,
    required_terms: tuple[str, ...],
    *,
    context: str,
) -> None:
    assert any(
        all(
            re.search(
                rf"(?<![a-z0-9_]){re.escape(term.lower())}(?![a-z0-9_])",
                statement.lower(),
            )
            for term in required_terms
        )
        for statement in _iter_documentation_statements(source)
    ), f"{context} must include a statement containing {required_terms}"


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


def test_readme_primary_workflow_example_is_kinase() -> None:
    source = _read(README)

    assert "## Kinase Workflow Example" in source
    assert "KinaseWorkflow().run(" in source
    assert "KinaseWorkflowRequest(" in source
    assert "DifferentialAnalysisWorkflow().run(" not in source
    assert "site_sequence" in source
    assert "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW" in source
    assert "PSGGGPGGSGRARTSSFAEPGGGGGGGGGGP" in source
    assert "FDDTPEKDSFRARSTSLNERPKSLRIARAPK" in source
    assert "DatasetLocalisationConfig(" in source
    assert 'confidence_column="localisation_confidence"' in source
    assert "min_confidence=0.75" in source


def test_readme_links_to_existing_api_workflow_docs() -> None:
    source = _read(README)

    assert "[Dataset building](docs/api/dataset-build-workflow.md)" in source
    assert "[Differential workflow](docs/api/differential-analysis.md)" in source
    assert "[Enrichment workflow](docs/api/enrichment.md)" in source
    assert "[Kinase workflow](docs/api/kinase.md)" in source
    assert "[Signalome workflow](docs/api/signalome.md)" in source

    assert DATASET_WORKFLOW_DOC.exists()
    assert DIFFERENTIAL_WORKFLOW_DOC.exists()
    assert ENRICHMENT_WORKFLOW_DOC.exists()
    assert KINASE_WORKFLOW_DOC.exists()
    assert SIGNALOME_WORKFLOW_DOC.exists()


def test_api_guide_differential_import_examples_match_supported_route() -> None:
    source = _read(API_GUIDE)

    assert API_GUIDE_IMPORT_SNIPPET in source
    assert API_GUIDE_API_IMPORT_SNIPPET in source
    assert "DifferentialAnalysisWorkflow().run(differential_request)" in source
    assert "from phospy import DifferentialAnalysis," not in source


def test_api_guide_small_working_example_includes_localisation_policy() -> None:
    source = _read(API_GUIDE)

    assert "## Small Working Example" in source
    assert "DatasetLocalisationConfig(" in source
    assert 'confidence_column="localisation_confidence"' in source
    assert "min_confidence=0.75" in source


def test_public_kinase_docs_prefer_activity_matrix() -> None:
    guide_source = _read(API_GUIDE)
    kinase_source = _read(KINASE_WORKFLOW_DOC)
    old_primary_phrase = (
        "`activity_result.activity_" + "scores` is the method-neutral primary"
    )

    _assert_statement_contains_all(
        guide_source,
        ("result models", "typed containers"),
        context="API guide result model contract",
    )
    _assert_statement_contains_all(
        kinase_source,
        ("activity_result.activity_matrix", "primary", "activity matrix"),
        context="kinase activity primary result contract",
    )
    _assert_statement_contains_all(
        kinase_source,
        ("activity_scores", "weighted_activity", "not preferred"),
        context="kinase activity compatibility alias guidance",
    )
    assert old_primary_phrase not in kinase_source
