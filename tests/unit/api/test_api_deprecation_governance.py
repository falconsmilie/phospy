from __future__ import annotations

import ast
import io
import os
import re
import subprocess
import sys
import textwrap
import tokenize
import tomllib
from dataclasses import dataclass
from pathlib import Path

import phospy.advanced as advanced_api
import phospy.api as public_api
from phospy._api_inventory import (
    ADVANCED_PUBLIC_API_BASELINE_COUNT,
    STABLE_PUBLIC_API_BASELINE_COUNT,
)
from phospy._deprecations import (
    DeprecationSourceUse,
    RetainedDeprecation,
    retained_deprecations,
)

ROOT = Path(__file__).resolve().parents[3]
EXEMPTION_MARKER = "phospy-deprecation-compat:"


@dataclass(frozen=True, slots=True)
class _RegisteredSourceUse:
    record: RetainedDeprecation
    source_use: DeprecationSourceUse


@dataclass(frozen=True, slots=True)
class _DeprecationFinding:
    identifier: str
    line_number: int
    message: str
    node: ast.AST | None


@dataclass(frozen=True, slots=True)
class _Exemption:
    identifier: str
    line_number: int


def test_python_consumer_sources_use_canonical_registered_deprecations() -> None:
    violations: list[str] = []
    for root_name in ("tests", "benchmarks", "examples"):
        root = ROOT / root_name
        if not root.exists():
            continue
        for path in sorted(root.rglob("*.py")):
            violations.extend(
                _registered_deprecation_violations(
                    path,
                    path.read_text(encoding="utf-8"),
                    python_exemptions_require_warning_capture=True,
                )
            )

    assert violations == []


def test_markdown_consumer_examples_use_canonical_registered_deprecations() -> None:
    violations: list[str] = []
    paths = (
        ROOT / "README.md",
        ROOT / "CHANGELOG.md",
        *sorted((ROOT / "docs").rglob("*.md")),
        *sorted((ROOT / "examples").rglob("*.md")),
    )

    for path in paths:
        violations.extend(
            _markdown_deprecation_violations(
                path,
                path.read_text(encoding="utf-8"),
            )
        )

    assert violations == []


def test_pytest_errors_on_unexpected_phospy_deprecation_warnings(
    tmp_path: Path,
) -> None:
    config = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    filterwarnings = config["tool"]["pytest"]["ini_options"]["filterwarnings"]

    assert "error::phospy._deprecations.PhosPyDeprecationWarning" in filterwarnings

    test_file = tmp_path / "test_uncaptured_phospy_deprecation.py"
    test_file.write_text(
        "\n".join(
            [
                "from phospy._deprecations import warn_deprecated",
                "",
                "def test_uncaptured_phospy_deprecation_fails():",
                "    warn_deprecated(",
                "        'science.differential.DifferentialAnalysis',",
                "        stacklevel=1,",
                "    )",
                "",
            ]
        ),
        encoding="utf-8",
    )
    completed = _run_pytest_file(test_file)

    output = f"{completed.stdout}\n{completed.stderr}"
    assert completed.returncode != 0, output
    assert "PhosPyDeprecationWarning" in output
    assert "DifferentialAnalysis" in output


def test_third_party_deprecation_warning_is_not_globally_escalated(
    tmp_path: Path,
) -> None:
    test_file = tmp_path / "test_third_party_deprecation.py"
    test_file.write_text(
        "\n".join(
            [
                "import warnings",
                "",
                "def test_third_party_deprecation_warning_is_not_an_error():",
                "    warnings.warn('third-party compatibility notice', DeprecationWarning)",
                "",
            ]
        ),
        encoding="utf-8",
    )
    completed = _run_pytest_file(test_file)

    output = f"{completed.stdout}\n{completed.stderr}"
    assert completed.returncode == 0, output


def test_api_surface_counts_remain_unchanged() -> None:
    assert len(public_api.__all__) == STABLE_PUBLIC_API_BASELINE_COUNT
    assert len(advanced_api.__all__) == ADVANCED_PUBLIC_API_BASELINE_COUNT


def test_deprecated_exact_value_in_ordinary_python_is_reported() -> None:
    source = """
    def test_ordinary_consumer():
        scoring_mode = "kinase_library_motif"
        assert scoring_mode
    """

    violations = _registered_deprecation_violations(
        Path("tests/unit/test_bad.py"),
        textwrap.dedent(source),
        python_exemptions_require_warning_capture=True,
    )

    assert len(violations) == 1
    assert "science.kinase.scoring_mode.kinase_library_motif" in violations[0]
    assert "exact-string-value" in violations[0]


def test_canonical_replacement_value_passes_python_source_policy() -> None:
    source = """
    def test_ordinary_consumer():
        scoring_mode = "kinase_library_contextual_motif"
        assert scoring_mode
    """

    violations = _registered_deprecation_violations(
        Path("tests/unit/test_good.py"),
        textwrap.dedent(source),
        python_exemptions_require_warning_capture=True,
    )

    assert violations == []


def test_longer_legitimate_value_does_not_match_deprecated_value_alias() -> None:
    source = """
    def test_ordinary_consumer():
        valid_mode = "kinase_library_motif_only"
        valid_policy = "kinase_library_motif_scoring_v1"
        assert valid_mode and valid_policy
    """

    violations = _registered_deprecation_violations(
        Path("tests/unit/test_good.py"),
        textwrap.dedent(source),
        python_exemptions_require_warning_capture=True,
    )

    assert violations == []


def test_deprecated_keyword_and_method_aliases_are_reported() -> None:
    source = """
    from phospy.science.datasets.preprocessing.pipeline import PreprocessingPipeline

    def test_ordinary_consumer(result):
        PreprocessingPipeline(stage_metadata_registry=())
        result.legacy_condition_statistics_table_dataframe()
    """

    violations = _registered_deprecation_violations(
        Path("tests/unit/test_bad.py"),
        textwrap.dedent(source),
        python_exemptions_require_warning_capture=True,
    )

    assert len(violations) == 2
    assert any(
        "preprocessing.pipeline.stage_metadata_registry" in item for item in violations
    )
    assert any(
        "activities.result.legacy_condition_statistics_table" in item
        for item in violations
    )


def test_deprecated_documentation_inline_and_code_block_examples_are_reported() -> None:
    source = """
    Inline bad example: `scoring_mode="kinase_library_motif"`.

    ```python
    scoring_mode = "kinase_library_motif"
    ```
    """

    violations = _markdown_deprecation_violations(
        Path("docs/api/bad.md"),
        textwrap.dedent(source),
    )

    assert len(violations) == 2
    assert all(
        "science.kinase.scoring_mode.kinase_library_motif" in item
        for item in violations
    )


def test_explicit_compatibility_use_under_warning_capture_is_accepted() -> None:
    source = """
    import pytest

    from phospy._deprecations import PhosPyDeprecationWarning
    from phospy.advanced import KinaseScoringConfig

    def test_compatibility_case():
        with pytest.warns(PhosPyDeprecationWarning, match="exploratory"):
            # phospy-deprecation-compat: contracts.kinase.KinaseScoringConfig.default
            config = KinaseScoringConfig.default()
        assert config
    """

    violations = _registered_deprecation_violations(
        Path("tests/unit/test_compat.py"),
        textwrap.dedent(source),
        python_exemptions_require_warning_capture=True,
    )

    assert violations == []


def test_uncaptured_deprecated_use_in_compatibility_test_still_fails() -> None:
    source = """
    from phospy.advanced import KinaseScoringConfig

    def test_compatibility_case():
        # phospy-deprecation-compat: contracts.kinase.KinaseScoringConfig.default
        config = KinaseScoringConfig.default()
        assert config
    """

    violations = _registered_deprecation_violations(
        Path("tests/unit/test_compat.py"),
        textwrap.dedent(source),
        python_exemptions_require_warning_capture=True,
    )

    assert len(violations) == 1
    assert "must capture and inspect PhosPyDeprecationWarning" in violations[0]


def test_unknown_and_unused_deprecation_exemptions_fail() -> None:
    unknown_source = """
    # phospy-deprecation-compat: missing.registry.identifier
    def test_noop():
        assert True
    """
    unused_source = """
    # phospy-deprecation-compat: contracts.kinase.KinaseScoringConfig.default
    def test_noop():
        assert True
    """

    unknown_violations = _registered_deprecation_violations(
        Path("tests/unit/test_unknown.py"),
        textwrap.dedent(unknown_source),
        python_exemptions_require_warning_capture=True,
    )
    unused_violations = _registered_deprecation_violations(
        Path("tests/unit/test_unused.py"),
        textwrap.dedent(unused_source),
        python_exemptions_require_warning_capture=True,
    )

    assert len(unknown_violations) == 1
    assert "unknown deprecation exemption" in unknown_violations[0]
    assert len(unused_violations) == 1
    assert "unused deprecation exemption" in unused_violations[0]


def test_existing_deprecated_api_import_route_check_uses_registry_metadata() -> None:
    source = """
    from phospy.api import KinaseScoringConfig
    """

    violations = _registered_deprecation_violations(
        Path("tests/unit/test_bad.py"),
        textwrap.dedent(source),
        python_exemptions_require_warning_capture=True,
    )

    assert len(violations) == 1
    assert "api-import:phospy.api.KinaseScoringConfig" in violations[0]
    assert "phospy.advanced" in violations[0]


def _run_pytest_file(test_file: Path) -> subprocess.CompletedProcess[str]:
    environment = dict(os.environ)
    environment["PYTHONPATH"] = os.pathsep.join(
        [str(ROOT / "src"), environment.get("PYTHONPATH", "")]
    )
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "-c",
            str(ROOT / "pyproject.toml"),
            str(test_file),
            "-q",
        ],
        cwd=ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )


def _markdown_deprecation_violations(path: Path, source: str) -> list[str]:
    findings: list[_DeprecationFinding] = []
    for line_number, block in _python_code_blocks(source):
        findings.extend(
            _python_deprecation_findings(
                path=path,
                source=block,
                line_offset=line_number - 1,
            )
        )
    if _scan_inline_code(path):
        findings.extend(_markdown_inline_deprecation_findings(path, source))
    return _violations_from_findings(
        path=path,
        findings=findings,
        exemptions=_text_exemptions(source),
        python_exemptions_require_warning_capture=False,
        parent_map={},
    )


def _registered_deprecation_violations(
    path: Path,
    source: str,
    *,
    python_exemptions_require_warning_capture: bool,
) -> list[str]:
    findings, parent_map = _python_deprecation_findings_with_parents(path, source)
    return _violations_from_findings(
        path=path,
        findings=findings,
        exemptions=_python_exemptions(source),
        python_exemptions_require_warning_capture=(
            python_exemptions_require_warning_capture
        ),
        parent_map=parent_map,
    )


def _violations_from_findings(
    *,
    path: Path,
    findings: list[_DeprecationFinding],
    exemptions: tuple[_Exemption, ...],
    python_exemptions_require_warning_capture: bool,
    parent_map: dict[ast.AST, ast.AST],
) -> list[str]:
    registry = {record.identifier: record for record in retained_deprecations()}
    used_exemptions: set[_Exemption] = set()
    violations: list[str] = []

    for exemption in exemptions:
        if exemption.identifier not in registry:
            violations.append(
                f"{path}:{exemption.line_number}: unknown deprecation exemption "
                f"{exemption.identifier!r}"
            )

    for finding in findings:
        exemption = _matching_exemption(finding, exemptions)
        if exemption is None:
            violations.append(finding.message)
            continue
        used_exemptions.add(exemption)
        if (
            python_exemptions_require_warning_capture
            and finding.node is not None
            and not _inside_inspected_phospy_deprecation_warning_capture(
                finding.node,
                parent_map,
            )
        ):
            violations.append(
                f"{path}:{finding.line_number}: {finding.identifier} exemption "
                "must capture and inspect PhosPyDeprecationWarning"
            )

    for exemption in exemptions:
        if exemption.identifier in registry and exemption not in used_exemptions:
            violations.append(
                f"{path}:{exemption.line_number}: unused deprecation exemption "
                f"{exemption.identifier!r}"
            )

    return violations


def _python_deprecation_findings_with_parents(
    path: Path,
    source: str,
) -> tuple[list[_DeprecationFinding], dict[ast.AST, ast.AST]]:
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return [], {}
    parent_map = _parent_map(tree)
    return (
        _python_deprecation_findings_from_tree(
            path=path,
            tree=tree,
            line_offset=0,
        ),
        parent_map,
    )


def _python_deprecation_findings(
    *,
    path: Path,
    source: str,
    line_offset: int,
) -> list[_DeprecationFinding]:
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return []
    return _python_deprecation_findings_from_tree(
        path=path,
        tree=tree,
        line_offset=line_offset,
    )


def _python_deprecation_findings_from_tree(
    *,
    path: Path,
    tree: ast.AST,
    line_offset: int,
) -> list[_DeprecationFinding]:
    findings: list[_DeprecationFinding] = []
    import_routes = _source_uses_by_import_route()
    exact_strings = _source_uses_by_token("exact-string-value")
    keyword_arguments = _source_uses_by_token("keyword-argument")
    class_aliases = _source_uses_by_token("class-alias")
    function_aliases = _source_uses_by_token("function-alias")
    classmethod_aliases = _source_uses_by_token("classmethod-alias")
    method_aliases = _source_uses_by_token("method-alias")
    property_aliases = _source_uses_by_token("property-alias")
    symbol_aliases = _source_uses_by_token("symbol")
    module_aliases = _module_aliases(tree)

    for node in ast.walk(tree):
        line_number = getattr(node, "lineno", 0) + line_offset
        if isinstance(node, ast.ImportFrom):
            if node.module is None:
                continue
            for alias in node.names:
                registered = import_routes.get((node.module, alias.name))
                if registered is not None:
                    findings.append(
                        _finding(
                            path=path,
                            line_number=line_number,
                            node=node,
                            registered=registered,
                        )
                    )
            continue

        if isinstance(node, ast.Attribute):
            if isinstance(node.value, ast.Name):
                registered = import_routes.get(
                    (module_aliases.get(node.value.id, ""), node.attr)
                )
                if registered is not None:
                    findings.append(
                        _finding(
                            path=path,
                            line_number=line_number,
                            node=node,
                            registered=registered,
                        )
                    )
            for registered in property_aliases.get(node.attr, ()):
                findings.append(
                    _finding(
                        path=path,
                        line_number=line_number,
                        node=node,
                        registered=registered,
                    )
                )
            continue

        if isinstance(node, ast.keyword) and node.arg is not None:
            for registered in keyword_arguments.get(node.arg, ()):
                findings.append(
                    _finding(
                        path=path,
                        line_number=line_number,
                        node=node,
                        registered=registered,
                    )
                )
            continue

        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            for registered in exact_strings.get(node.value, ()):
                findings.append(
                    _finding(
                        path=path,
                        line_number=line_number,
                        node=node,
                        registered=registered,
                    )
                )
            continue

        if isinstance(node, ast.Call):
            call_name = _call_name(node.func)
            if call_name is None:
                continue
            for registered in (
                *class_aliases.get(call_name, ()),
                *function_aliases.get(call_name, ()),
                *symbol_aliases.get(call_name, ()),
            ):
                findings.append(
                    _finding(
                        path=path,
                        line_number=line_number,
                        node=node,
                        registered=registered,
                    )
                )
            if isinstance(node.func, ast.Attribute):
                for registered in method_aliases.get(node.func.attr, ()):
                    findings.append(
                        _finding(
                            path=path,
                            line_number=line_number,
                            node=node,
                            registered=registered,
                        )
                    )
                for registered in classmethod_aliases.get(node.func.attr, ()):
                    if _owner_matches(node.func.value, registered.source_use.owner):
                        findings.append(
                            _finding(
                                path=path,
                                line_number=line_number,
                                node=node,
                                registered=registered,
                            )
                        )

    return findings


def _markdown_inline_deprecation_findings(
    path: Path,
    source: str,
) -> list[_DeprecationFinding]:
    findings: list[_DeprecationFinding] = []
    for line_number, snippet in _inline_code_spans(source):
        snippet_findings = _python_deprecation_findings(
            path=path,
            source=_parseable_inline_python(snippet),
            line_offset=line_number - 1,
        )
        findings.extend(snippet_findings)
        findings.extend(_literal_inline_source_findings(path, line_number, snippet))
    return _dedupe_findings(findings)


def _literal_inline_source_findings(
    path: Path,
    line_number: int,
    snippet: str,
) -> list[_DeprecationFinding]:
    stripped = snippet.strip()
    normalized_call = stripped[:-2] if stripped.endswith("()") else stripped
    findings: list[_DeprecationFinding] = []
    for registered in _checked_source_uses():
        source_use = registered.source_use
        if source_use.kind in {"exact-string-value", "import-route"}:
            continue
        if stripped == source_use.token or normalized_call == source_use.token:
            findings.append(
                _finding(
                    path=path,
                    line_number=line_number,
                    node=None,
                    registered=registered,
                )
            )
    return findings


def _parseable_inline_python(snippet: str) -> str:
    try:
        ast.parse(snippet)
        return snippet
    except SyntaxError:
        return f"__phospy_inline__({snippet})"


def _python_code_blocks(source: str) -> tuple[tuple[int, str], ...]:
    blocks: list[tuple[int, str]] = []
    for match in re.finditer(
        r"```python[^\n]*\n(?P<code>.*?)\n```",
        source,
        flags=re.DOTALL,
    ):
        line_number = source[: match.start("code")].count("\n") + 1
        blocks.append((line_number, match.group("code")))
    return tuple(blocks)


def _inline_code_spans(source: str) -> tuple[tuple[int, str], ...]:
    spans: list[tuple[int, str]] = []
    fenced_ranges = [
        (match.start(), match.end())
        for match in re.finditer(r"```.*?```", source, flags=re.DOTALL)
    ]
    for match in re.finditer(r"(?<!`)`([^`\n]+)`(?!`)", source):
        if any(start <= match.start() < end for start, end in fenced_ranges):
            continue
        line_number = source[: match.start(1)].count("\n") + 1
        spans.append((line_number, match.group(1)))
    return tuple(spans)


def _scan_inline_code(path: Path) -> bool:
    try:
        relative = path.relative_to(ROOT)
    except ValueError:
        relative = path
    parts = relative.parts
    return (
        relative == Path("README.md")
        or (len(parts) >= 1 and parts[0] == "examples")
        or (len(parts) >= 2 and parts[0] == "docs" and parts[1] == "api")
        or relative
        in {
            Path("docs/reference_bundles.md"),
            Path("docs/scientific-coverage.md"),
            Path("docs/workflow_contracts.md"),
        }
    )


def _python_exemptions(source: str) -> tuple[_Exemption, ...]:
    exemptions: list[_Exemption] = []
    try:
        tokens = tokenize.generate_tokens(io.StringIO(source).readline)
        for token in tokens:
            if token.type != tokenize.COMMENT:
                continue
            exemptions.extend(
                _exemptions_from_line(
                    token.string,
                    line_number=token.start[0],
                )
            )
    except tokenize.TokenError:
        return ()
    return tuple(exemptions)


def _text_exemptions(source: str) -> tuple[_Exemption, ...]:
    exemptions: list[_Exemption] = []
    for line_number, line in enumerate(source.splitlines(), start=1):
        if "<!--" not in line:
            continue
        exemptions.extend(_exemptions_from_line(line, line_number=line_number))
    return tuple(exemptions)


def _exemptions_from_line(line: str, *, line_number: int) -> tuple[_Exemption, ...]:
    if EXEMPTION_MARKER not in line:
        return ()
    raw = line.split(EXEMPTION_MARKER, 1)[1].replace("-->", "")
    return tuple(
        _Exemption(identifier=identifier.strip(), line_number=line_number)
        for identifier in re.split(r"[,\s]+", raw.strip())
        if identifier.strip()
    )


def _matching_exemption(
    finding: _DeprecationFinding,
    exemptions: tuple[_Exemption, ...],
) -> _Exemption | None:
    for exemption in exemptions:
        if exemption.identifier != finding.identifier:
            continue
        if exemption.line_number in {finding.line_number, finding.line_number - 1}:
            return exemption
    return None


def _inside_inspected_phospy_deprecation_warning_capture(
    node: ast.AST,
    parent_map: dict[ast.AST, ast.AST],
) -> bool:
    current = node
    while current in parent_map:
        current = parent_map[current]
        if not isinstance(current, ast.With):
            continue
        for item in current.items:
            if _is_inspected_phospy_deprecation_warning_capture(item):
                return True
    return False


def _is_inspected_phospy_deprecation_warning_capture(item: ast.withitem) -> bool:
    context = item.context_expr
    if not isinstance(context, ast.Call):
        return False
    if not _is_pytest_warns_call(context):
        return False
    if not _warns_phospy_deprecation_warning(context):
        return False
    return item.optional_vars is not None or any(
        keyword.arg == "match" and not _is_none_literal(keyword.value)
        for keyword in context.keywords
    )


def _is_pytest_warns_call(node: ast.Call) -> bool:
    if isinstance(node.func, ast.Attribute):
        return (
            isinstance(node.func.value, ast.Name)
            and node.func.value.id == "pytest"
            and node.func.attr == "warns"
        )
    return isinstance(node.func, ast.Name) and node.func.id == "warns"


def _warns_phospy_deprecation_warning(node: ast.Call) -> bool:
    candidates = list(node.args[:1])
    candidates.extend(
        keyword.value
        for keyword in node.keywords
        if keyword.arg in {"expected_warning", "warning"}
    )
    return any(
        _is_phospy_deprecation_warning_name(candidate) for candidate in candidates
    )


def _is_phospy_deprecation_warning_name(node: ast.AST) -> bool:
    if isinstance(node, ast.Name):
        return node.id == "PhosPyDeprecationWarning"
    if isinstance(node, ast.Attribute):
        return node.attr == "PhosPyDeprecationWarning"
    return False


def _is_none_literal(node: ast.AST) -> bool:
    return isinstance(node, ast.Constant) and node.value is None


def _parent_map(tree: ast.AST) -> dict[ast.AST, ast.AST]:
    parents: dict[ast.AST, ast.AST] = {}
    for parent in ast.walk(tree):
        for child in ast.iter_child_nodes(parent):
            parents[child] = parent
    return parents


def _call_name(node: ast.AST) -> str | None:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return node.attr
    return None


def _owner_matches(node: ast.AST, owner: str | None) -> bool:
    if owner is None:
        return True
    if isinstance(node, ast.Name):
        return node.id == owner
    if isinstance(node, ast.Attribute):
        return node.attr == owner
    return False


def _module_aliases(tree: ast.AST) -> dict[str, str]:
    aliases: dict[str, str] = {}
    for node in ast.walk(tree):
        if not isinstance(node, ast.Import):
            continue
        for alias in node.names:
            if alias.asname is not None:
                aliases[alias.asname] = alias.name
    return aliases


def _source_uses_by_import_route() -> dict[tuple[str, str], _RegisteredSourceUse]:
    routes: dict[tuple[str, str], _RegisteredSourceUse] = {}
    for registered in _checked_source_uses():
        source_use = registered.source_use
        if source_use.kind != "import-route":
            continue
        if source_use.module is None:
            raise AssertionError(
                f"{registered.record.identifier} import-route source use lacks module"
            )
        routes[(source_use.module, source_use.token)] = registered
    return routes


def _source_uses_by_token(
    source_use_kind: str,
) -> dict[str, tuple[_RegisteredSourceUse, ...]]:
    grouped: dict[str, list[_RegisteredSourceUse]] = {}
    for registered in _checked_source_uses():
        if registered.source_use.kind != source_use_kind:
            continue
        grouped.setdefault(registered.source_use.token, []).append(registered)
    return {token: tuple(items) for token, items in grouped.items()}


def _checked_source_uses() -> tuple[_RegisteredSourceUse, ...]:
    return tuple(
        _RegisteredSourceUse(record=record, source_use=source_use)
        for record in retained_deprecations()
        for source_use in record.source_uses
        if source_use.check_consumer_sources
    )


def _finding(
    *,
    path: Path,
    line_number: int,
    node: ast.AST | None,
    registered: _RegisteredSourceUse,
) -> _DeprecationFinding:
    record = registered.record
    source_use = registered.source_use
    module_prefix = f"{source_use.module}." if source_use.module else ""
    owner_prefix = f"{source_use.owner}." if source_use.owner else ""
    return _DeprecationFinding(
        identifier=record.identifier,
        line_number=line_number,
        node=node,
        message=(
            f"{path}:{line_number}: {record.identifier} uses deprecated "
            f"{source_use.kind} {module_prefix}{owner_prefix}{source_use.token!r}; "
            f"use {record.replacement}"
        ),
    )


def _dedupe_findings(
    findings: list[_DeprecationFinding],
) -> list[_DeprecationFinding]:
    seen: set[tuple[str, int, str]] = set()
    deduped: list[_DeprecationFinding] = []
    for finding in findings:
        key = (finding.identifier, finding.line_number, finding.message)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(finding)
    return deduped
