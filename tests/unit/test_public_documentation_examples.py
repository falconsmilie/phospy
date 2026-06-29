from __future__ import annotations

import ast
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

README_IMPORT_SNIPPET = """from phospy import AnalysisReadyDatasetBuilder
from phospy import DifferentialAnalysisWorkflow, KinaseWorkflow, SignalomeWorkflow
"""

API_GUIDE_IMPORT_SNIPPET = """from phospy import (
    AnalysisReadyDatasetBuilder,
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


def _iter_python_code_blocks(source: str) -> tuple[str, ...]:
    return tuple(
        match.group("code").strip()
        for match in re.finditer(
            r"```python\s*\n(?P<code>.*?)\n```",
            source,
            flags=re.DOTALL,
        )
    )


def _parse_python_code_blocks(source: str) -> tuple[ast.Module, ...]:
    parsed_blocks: list[ast.Module] = []
    for block in _iter_python_code_blocks(source):
        try:
            parsed_blocks.append(ast.parse(block))
        except SyntaxError as exc:  # pragma: no cover - assertion context only
            raise AssertionError(f"invalid documented Python example: {exc}") from exc
    return tuple(parsed_blocks)


def _imported_names(source: str, module_name: str) -> set[str]:
    names: set[str] = set()
    for tree in _parse_python_code_blocks(source):
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module == module_name:
                names.update(alias.name for alias in node.names)
    return names


def _assert_python_imports(
    source: str,
    module_name: str,
    names: tuple[str, ...],
    *,
    context: str,
) -> None:
    imported = _imported_names(source, module_name)
    missing = sorted(set(names) - imported)
    assert not missing, f"{context} missing imports from {module_name}: {missing}"


def _assert_python_imports_absent(
    source: str,
    module_name: str,
    names: tuple[str, ...],
    *,
    context: str,
) -> None:
    imported = _imported_names(source, module_name)
    present = sorted(set(names) & imported)
    assert not present, f"{context} documents unsupported imports: {present}"


def _call_name(call: ast.Call) -> str | None:
    if isinstance(call.func, ast.Name):
        return call.func.id
    return None


def _is_workflow_run_call(call: ast.Call, workflow_name: str) -> bool:
    return (
        isinstance(call.func, ast.Attribute)
        and call.func.attr == "run"
        and isinstance(call.func.value, ast.Call)
        and isinstance(call.func.value.func, ast.Name)
        and call.func.value.func.id == workflow_name
    )


def _assert_python_call(source: str, call_name: str, *, context: str) -> None:
    assert any(
        _call_name(call) == call_name
        for tree in _parse_python_code_blocks(source)
        for call in ast.walk(tree)
        if isinstance(call, ast.Call)
    ), f"{context} must include a documented {call_name}(...) example"


def _assert_python_run_call(
    source: str,
    workflow_name: str,
    *,
    context: str,
) -> None:
    assert any(
        _is_workflow_run_call(call, workflow_name)
        for tree in _parse_python_code_blocks(source)
        for call in ast.walk(tree)
        if isinstance(call, ast.Call)
    ), f"{context} must include {workflow_name}().run(...)"


def _assert_no_python_run_call(
    source: str,
    workflow_name: str,
    *,
    context: str,
) -> None:
    assert not any(
        _is_workflow_run_call(call, workflow_name)
        for tree in _parse_python_code_blocks(source)
        for call in ast.walk(tree)
        if isinstance(call, ast.Call)
    ), f"{context} must not document {workflow_name}().run(...)"


_MISSING = object()


def _literal_value(node: ast.AST) -> object:
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
        return f"{node.value.id}.{node.attr}"
    return _MISSING


def _assert_python_call_keyword(
    source: str,
    call_name: str,
    keyword_name: str,
    expected_value: object = _MISSING,
    *,
    context: str,
) -> None:
    for tree in _parse_python_code_blocks(source):
        for call in ast.walk(tree):
            if not isinstance(call, ast.Call) or _call_name(call) != call_name:
                continue
            for keyword in call.keywords:
                if keyword.arg != keyword_name:
                    continue
                if expected_value is _MISSING:
                    return
                if _literal_value(keyword.value) == expected_value:
                    return

    expected = (
        keyword_name
        if expected_value is _MISSING
        else f"{keyword_name}={expected_value!r}"
    )
    raise AssertionError(f"{context} must document {call_name}({expected}, ...)")


def _assert_python_constant(source: str, expected: object, *, context: str) -> None:
    assert any(
        isinstance(node, ast.Constant) and node.value == expected
        for tree in _parse_python_code_blocks(source)
        for node in ast.walk(tree)
    ), f"{context} must include literal {expected!r} in a Python example"


def _assert_markdown_links_to_targets(
    source: str,
    targets: tuple[str, ...],
    *,
    context: str,
) -> None:
    linked_targets = {
        match.group("target").split("#", 1)[0]
        for match in re.finditer(r"\[[^\]]+\]\((?P<target>[^)]+)\)", source)
    }
    missing = sorted(set(targets) - linked_targets)
    assert not missing, f"{context} missing Markdown links to: {missing}"


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
    assert "DifferentialAnalysisWorkflow" in namespace
    assert "KinaseWorkflow" in namespace
    assert "SignalomeWorkflow" in namespace
    assert "DifferentialAnalysisRequest" in namespace


def test_readme_differential_import_example_matches_supported_route() -> None:
    source = _read(README)

    _assert_python_imports(
        source,
        "phospy",
        (
            "AnalysisReadyDatasetBuilder",
            "DifferentialAnalysisWorkflow",
            "KinaseWorkflow",
            "SignalomeWorkflow",
        ),
        context="README public import contract",
    )
    _assert_python_imports_absent(
        source,
        "phospy",
        ("DifferentialAnalysis",),
        context="README differential import route",
    )
    _assert_no_python_run_call(
        source,
        "DifferentialAnalysis",
        context="README differential import route",
    )


def test_readme_primary_workflow_example_is_kinase() -> None:
    source = _read(README)

    _assert_python_imports(
        source,
        "phospy",
        ("AnalysisReadyDatasetBuilder", "KinaseWorkflow"),
        context="README primary workflow example",
    )
    _assert_python_imports(
        source,
        "phospy.api",
        (
            "DatasetBuildRequest",
            "DatasetLocalisationConfig",
            "DatasetPreprocessingConfig",
            "KinaseWorkflowRequest",
            "Organism",
            "ReferencePreset",
        ),
        context="README primary workflow example",
    )
    _assert_python_run_call(
        source,
        "KinaseWorkflow",
        context="README primary workflow example",
    )
    _assert_python_call(
        source,
        "KinaseWorkflowRequest",
        context="README primary workflow example",
    )
    _assert_no_python_run_call(
        source,
        "DifferentialAnalysisWorkflow",
        context="README primary workflow example",
    )
    _assert_python_constant(
        source,
        "site_sequence",
        context="README primary workflow example",
    )
    _assert_python_call_keyword(
        source,
        "DatasetLocalisationConfig",
        "confidence_column",
        "localisation_confidence",
        context="README primary workflow example",
    )
    _assert_python_call_keyword(
        source,
        "DatasetLocalisationConfig",
        "min_confidence",
        0.75,
        context="README primary workflow example",
    )
    _assert_python_call_keyword(
        source,
        "KinaseWorkflowRequest",
        "references",
        "ReferencePreset.AUTO",
        context="README primary workflow example",
    )


def test_readme_links_to_existing_api_workflow_docs() -> None:
    source = _read(README)

    _assert_markdown_links_to_targets(
        source,
        (
            "docs/api/dataset-build-workflow.md",
            "docs/api/differential-analysis.md",
            "docs/api/enrichment.md",
            "docs/api/kinase.md",
            "docs/api/signalome.md",
        ),
        context="README workflow documentation index",
    )

    assert DATASET_WORKFLOW_DOC.exists()
    assert DIFFERENTIAL_WORKFLOW_DOC.exists()
    assert ENRICHMENT_WORKFLOW_DOC.exists()
    assert KINASE_WORKFLOW_DOC.exists()
    assert SIGNALOME_WORKFLOW_DOC.exists()


def test_api_guide_differential_import_examples_match_supported_route() -> None:
    source = _read(API_GUIDE)

    _assert_python_imports(
        source,
        "phospy",
        (
            "AnalysisReadyDatasetBuilder",
            "DifferentialAnalysisWorkflow",
            "KinaseWorkflow",
            "SignalomeWorkflow",
        ),
        context="API guide top-level import contract",
    )
    _assert_python_imports(
        source,
        "phospy.api",
        (
            "DatasetBuildRequest",
            "ExperimentalDesign",
            "Contrast",
            "SampleDesignRecord",
            "DatasetPreprocessingConfig",
            "DifferentialAnalysisRequest",
            "KinaseWorkflowRequest",
            "Organism",
            "ReferenceBundle",
            "ReferencePreset",
            "SignalomeConfig",
            "SignalomeWorkflowRequest",
            "UnsupportedInputFormatError",
            "WorkflowValidationError",
        ),
        context="API guide phospy.api import contract",
    )
    _assert_python_run_call(
        source,
        "DifferentialAnalysisWorkflow",
        context="API guide differential workflow route",
    )
    _assert_python_imports_absent(
        source,
        "phospy",
        ("DifferentialAnalysis",),
        context="API guide differential workflow route",
    )


def test_api_guide_small_working_example_includes_localisation_policy() -> None:
    source = _read(API_GUIDE)

    _assert_python_imports(
        source,
        "phospy.api",
        ("DatasetLocalisationConfig", "DatasetPreprocessingConfig"),
        context="API guide working example localisation policy",
    )
    _assert_python_call_keyword(
        source,
        "DatasetLocalisationConfig",
        "confidence_column",
        "localisation_confidence",
        context="API guide working example localisation policy",
    )
    _assert_python_call_keyword(
        source,
        "DatasetLocalisationConfig",
        "min_confidence",
        0.75,
        context="API guide working example localisation policy",
    )


def test_public_kinase_docs_prefer_activity_matrix() -> None:
    guide_source = _read(API_GUIDE)
    kinase_source = _read(KINASE_WORKFLOW_DOC)

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
