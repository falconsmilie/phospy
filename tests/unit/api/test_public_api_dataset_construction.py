from __future__ import annotations

import inspect
import re
from pathlib import Path
from typing import get_type_hints

import phospy
import phospy.api as public_api
from phospy.api import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DatasetBuildRequest,
)

ROOT = Path(__file__).resolve().parents[3]
PYTHON_CODE_BLOCK_PATTERN = re.compile(
    r"```python\s*\n(?P<code>.*?)\n```",
    flags=re.DOTALL,
)
DIRECT_DATASET_CONSTRUCTOR_CALL_PATTERN = re.compile(
    r"(?<![.\w])AnalysisReadyPhosphoDataset\s*\("
)
STALE_DATASET_CONSTRUCTION_WORDING_PATTERNS = {
    "direct construction availability": re.compile(
        r"\bdirect(?:\s+dataset)?\s+construction\b"
        r"[^.!?]{0,160}"
        r"\b(?:remains|still|continues|may\s+still|is)\b"
        r"[^.!?]{0,80}"
        r"\b(?:available|valid|supported|succeeds|works)\b",
        flags=re.IGNORECASE,
    ),
    "direct constructor warning compatibility": re.compile(
        r"(?:\b(?:direct(?:\s+public)?\s+construction|direct\s+constructor)\b|"
        r"AnalysisReadyPhosphoDataset\s*\([^)]*\))"
        r"[^.!?]{0,160}"
        r"\b(?:always\s+emits|emits|warns?|warning-based)\b"
        r"[^.!?]{0,80}"
        r"\b(?:DeprecationWarning|compatibility)\b",
        flags=re.IGNORECASE,
    ),
    "minimal direct-construction provenance": re.compile(
        r"\b(?:direct[-\s]construction\b[^.!?]{0,120}"
        r"\b(?:receives|creates|gets|manufactures|mints|generates|builds)\b|"
        r"(?:receives|creates|gets|manufactures|mints|generates|builds)\b"
        r"[^.!?]{0,120}\bdirect[-\s]construction\b)"
        r"[^.!?]{0,80}\bprovenance\b",
        flags=re.IGNORECASE,
    ),
    "advanced direct-construction lane": re.compile(
        r"\badvanced/trusted\s+direct[-\s]"
        r"(?:construction|analysis-ready|dataset)",
        flags=re.IGNORECASE,
    ),
}
NEGATED_OR_SEALED_BOUNDARY_PATTERN = re.compile(
    r"\b(?:no|not|does\s+not|do\s+not|must\s+not|cannot|can't|"
    r"unsupported|sealed|raises\s+immediately|fail(?:s|ed)?\s+immediately|"
    r"reject(?:s|ed)?)\b",
    flags=re.IGNORECASE,
)
NON_DATASET_CONSTRUCTION_CONTEXT_PATTERN = re.compile(
    r"\b(?:public\s+result|result\s+model|result\s+tables?)\b",
    flags=re.IGNORECASE,
)


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


def test_exported_dataset_signature_has_no_private_validation_controls() -> None:
    root_parameters = inspect.signature(phospy.AnalysisReadyPhosphoDataset).parameters
    api_parameters = inspect.signature(
        public_api.AnalysisReadyPhosphoDataset
    ).parameters

    assert root_parameters == api_parameters
    assert tuple(root_parameters) == ("args", "kwargs")
    assert "_emit_direct_constructor_deprecation" not in root_parameters
    assert "_enforce_trusted_table_fingerprints" not in root_parameters


def test_public_api_marks_direct_dataset_construction_sealed() -> None:
    model_doc = AnalysisReadyPhosphoDataset.__doc__
    factory_doc = AnalysisReadyPhosphoDataset.from_trusted_tables.__doc__

    assert model_doc is not None
    assert factory_doc is not None
    normalized_factory_doc = " ".join(factory_doc.split())
    assert "stable public result/domain type" in model_doc
    assert "ordinary direct construction is not a supported creation path" in model_doc
    assert "AnalysisReadyDatasetBuilder.run" in model_doc
    assert "from_trusted_tables" in model_doc
    assert "typed evidence or an explicit" in model_doc
    assert "localisation" in model_doc
    assert "cannot prove" in model_doc
    assert "biological correctness" in model_doc
    assert "same structural invariants as the builder-owned path" in (
        normalized_factory_doc
    )
    assert "site_sequence" in factory_doc
    assert "source, policy" in factory_doc
    assert "threshold" in factory_doc
    assert "cannot prove" in factory_doc
    assert "biological correctness" in factory_doc


def _maintained_public_documentation_paths() -> tuple[Path, ...]:
    return (
        ROOT / "README.md",
        ROOT / "CHANGELOG.md",
        *sorted((ROOT / "docs").rglob("*.md")),
        *sorted((ROOT / "examples").rglob("*.md")),
        *sorted((ROOT / "examples").rglob("*.py")),
    )


def _python_code_blocks(source: str) -> tuple[str, ...]:
    return tuple(
        match.group("code") for match in PYTHON_CODE_BLOCK_PATTERN.finditer(source)
    )


def _documentation_statements_without_code_blocks(source: str) -> tuple[str, ...]:
    source_without_code = re.sub(
        r"```.*?```",
        "",
        source,
        flags=re.DOTALL,
    )
    statements: list[str] = []
    for block in re.split(r"\n\s*\n", source_without_code):
        stripped = block.strip()
        if not stripped:
            continue
        chunks = stripped.splitlines() if stripped.startswith("|") else (stripped,)
        for chunk in chunks:
            normalised = " ".join(chunk.split())
            statements.append(normalised)
    return tuple(statements)


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
    direct_constructor_examples = [
        f"{path.relative_to(ROOT).as_posix()}: {block.strip()}"
        for path in documentation_paths
        for block in _python_code_blocks(path.read_text(encoding="utf-8"))
        if DIRECT_DATASET_CONSTRUCTOR_CALL_PATTERN.search(block)
    ]
    assert direct_constructor_examples == []
    assert (
        "from phospy import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset"
    ) not in documentation
    assert "advanced/trusted" in documentation
    assert "AnalysisReadyPhosphoDataset.from_trusted_tables" in documentation
    assert "raises immediately" in documentation


def test_public_docs_reject_stale_direct_dataset_construction_wording() -> None:
    stale_matches: list[str] = []
    for path in _maintained_public_documentation_paths():
        source = path.read_text(encoding="utf-8")
        for statement in _documentation_statements_without_code_blocks(source):
            if (
                NON_DATASET_CONSTRUCTION_CONTEXT_PATTERN.search(statement)
                and "AnalysisReadyPhosphoDataset" not in statement
            ):
                continue
            for label, pattern in STALE_DATASET_CONSTRUCTION_WORDING_PATTERNS.items():
                match = pattern.search(statement)
                if match is None:
                    continue
                match_context = statement[
                    max(0, match.start() - 120) : min(len(statement), match.end() + 120)
                ]
                if NEGATED_OR_SEALED_BOUNDARY_PATTERN.search(match_context):
                    continue
                stale_matches.append(
                    f"{path.relative_to(ROOT).as_posix()}: {label}: {statement}"
                )

    assert stale_matches == []
