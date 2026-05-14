from __future__ import annotations

import ast
import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class TestFunction:
    file_path: Path
    relpath: str
    name: str
    start_line: int
    end_line: int
    lines: tuple[str, ...]
    assert_count: int


@dataclass(frozen=True)
class Candidate:
    relpath: str
    test_name: str
    line_range: str
    score: int
    assert_count: int
    payload_assert_count: int
    key_check_count: int
    dataframe_numeric_count: int
    large_comparison_count: int
    reason: str
    stability_assessment: str
    stability_note: str


REPO_ROOT = Path(__file__).resolve().parents[2]
TESTS_ROOT = REPO_ROOT / "tests"
OUTPUT_PATH = REPO_ROOT / "docs" / "testing" / "diagnostic_payload_test_report.md"
PRIMARY_TARGET = TESTS_ROOT / "unit" / "test_signalome_workflow_diagnostics.py"

PAYLOAD_ASSERT_PATTERN = re.compile(
    r"assert\s+.*\b(?:payload|diagnostic|diagnostics|report|summary|provenance)\s*\[",
    re.IGNORECASE,
)
KEY_CHECK_PATTERN = re.compile(
    r"assert\s+['\"][^'\"]+['\"]\s+in\s+[A-Za-z_][A-Za-z0-9_]*",
    re.IGNORECASE,
)
DATAFRAME_NUMERIC_PATTERN = re.compile(
    r"(pd\.testing|pytest\.approx|np\.(isclose|allclose)|finite numeric|DataFrame|to_dict|equals\()",
    re.IGNORECASE,
)
LARGE_COMPARISON_PATTERN = re.compile(
    r"assert\s+.*==\s*[\{\[]|assert\s+.*==\s*expected_[A-Za-z0-9_]*",
    re.IGNORECASE,
)

INTERNAL_MARKERS: tuple[re.Pattern[str], ...] = (
    re.compile(r"\._(?!_)[A-Za-z0-9_]+\b"),
    re.compile(r"\bmonkeypatch\b"),
    re.compile(r"\bpatch(?:\.object|\.dict|\.multiple)?\s*\("),
    re.compile(r"\bMock\s*\(|\bMagicMock\s*\("),
    re.compile(r"_borrow_[A-Za-z0-9_]+"),
    re.compile(r"from\s+[A-Za-z0-9_.]+\._[A-Za-z0-9_.]+\s+import"),
    re.compile(r"\binternal\b", re.IGNORECASE),
    re.compile(r"\bboundary[_ ]error\b", re.IGNORECASE),
    re.compile(r"\bseam\b", re.IGNORECASE),
    re.compile(r"\bexecutor\b", re.IGNORECASE),
    re.compile(r"\bworkflow_diagnostics\b", re.IGNORECASE),
)

STABLE_MARKERS: tuple[re.Pattern[str], ...] = (
    re.compile(r"\bpublic\b", re.IGNORECASE),
    re.compile(r"\bscientific\b", re.IGNORECASE),
    re.compile(r"\bprovenance\b", re.IGNORECASE),
    re.compile(r"\bcontract\b", re.IGNORECASE),
    re.compile(r"\bschema\b", re.IGNORECASE),
    re.compile(r"\boutput\b", re.IGNORECASE),
)


def _discover_primary_and_diagnostic_files() -> tuple[
    Path | None, list[Path], list[Path]
]:
    test_files = sorted(TESTS_ROOT.rglob("test*.py"), key=lambda path: path.as_posix())
    diagnostic_files = [
        path
        for path in test_files
        if "diagnostic" in path.name.lower() or "diagnostics" in path.name.lower()
    ]

    primary_path: Path | None = None
    if PRIMARY_TARGET.exists():
        primary_path = PRIMARY_TARGET
    else:
        alternatives = [
            path
            for path in test_files
            if "signalome" in path.name.lower() and "diagnostic" in path.name.lower()
        ]
        if alternatives:
            primary_path = alternatives[0]
            if primary_path not in diagnostic_files:
                diagnostic_files.append(primary_path)

    diagnostic_files = sorted(set(diagnostic_files), key=lambda path: path.as_posix())
    return primary_path, diagnostic_files, test_files


def _extract_test_functions(path: Path) -> list[TestFunction]:
    source = path.read_text(encoding="utf-8-sig", errors="replace")
    lines = source.splitlines()
    module = ast.parse(source)
    relpath = path.relative_to(REPO_ROOT).as_posix()

    functions: list[TestFunction] = []
    for node in module.body:
        if not isinstance(node, ast.FunctionDef):
            continue
        if not node.name.startswith("test_"):
            continue
        if node.end_lineno is None:
            continue
        start = node.lineno
        end = node.end_lineno
        fn_lines = tuple(lines[start - 1 : end])
        assert_count = sum(isinstance(inner, ast.Assert) for inner in ast.walk(node))
        functions.append(
            TestFunction(
                file_path=path,
                relpath=relpath,
                name=node.name,
                start_line=start,
                end_line=end,
                lines=fn_lines,
                assert_count=assert_count,
            )
        )
    return functions


def _count_matches(lines: tuple[str, ...], pattern: re.Pattern[str]) -> int:
    return sum(1 for line in lines if pattern.search(line))


def _count_marker_hits(text: str, patterns: tuple[re.Pattern[str], ...]) -> int:
    return sum(1 for pattern in patterns if pattern.search(text))


def _build_reason(
    *,
    assert_count: int,
    payload_assert_count: int,
    key_check_count: int,
    dataframe_numeric_count: int,
    large_comparison_count: int,
) -> str:
    reasons: list[str] = []
    if payload_assert_count >= 4:
        reasons.append("repeated `assert payload[...]`-style indexing checks")
    if key_check_count >= 3:
        reasons.append("repeated key-presence assertions")
    if dataframe_numeric_count >= 3:
        reasons.append("many dataframe/numeric diagnostic field assertions")
    if large_comparison_count >= 1:
        reasons.append("large expected dict/list comparison pattern")
    if assert_count >= 15:
        reasons.append("high assertion volume in one test")
    if not reasons:
        reasons.append("multiple diagnostic assertions suggest high maintenance load")
    return "; ".join(reasons)


def _assess_stability(test_function: TestFunction) -> tuple[str, str]:
    text = "\n".join(test_function.lines)
    stable_hits = _count_marker_hits(
        test_function.name, STABLE_MARKERS
    ) + _count_marker_hits(text, STABLE_MARKERS)
    internal_hits = _count_marker_hits(
        test_function.relpath.lower(), INTERNAL_MARKERS
    ) + _count_marker_hits(text, INTERNAL_MARKERS)

    if internal_hits >= 2:
        return (
            "internal/unstable",
            "Touches private/internal seams or heavy test doubles; snapshot locking may be brittle.",
        )
    if stable_hits >= 2 and internal_hits == 0:
        return (
            "public/stable",
            "Appears to assert public/scientific/contract-facing payload semantics.",
        )
    return (
        "mixed/needs review",
        "Contains both potentially stable contract checks and implementation-coupled details.",
    )


def _score_cluster(
    *,
    assert_count: int,
    payload_assert_count: int,
    key_check_count: int,
    dataframe_numeric_count: int,
    large_comparison_count: int,
) -> int:
    score = 0
    score += payload_assert_count * 3
    score += key_check_count * 2
    score += dataframe_numeric_count * 2
    score += large_comparison_count * 3
    score += max(assert_count - 8, 0) // 2
    return score


def _is_candidate(
    *,
    score: int,
    assert_count: int,
    payload_assert_count: int,
    key_check_count: int,
    dataframe_numeric_count: int,
    large_comparison_count: int,
) -> bool:
    cluster_signals = (
        payload_assert_count
        + key_check_count
        + dataframe_numeric_count
        + large_comparison_count
    )
    if assert_count < 8:
        return False
    if score < 14:
        return False
    return cluster_signals >= 4


def _collect_candidates(diagnostic_files: list[Path]) -> list[Candidate]:
    candidates: list[Candidate] = []
    for path in diagnostic_files:
        for test_function in _extract_test_functions(path):
            payload_assert_count = _count_matches(
                test_function.lines, PAYLOAD_ASSERT_PATTERN
            )
            key_check_count = _count_matches(test_function.lines, KEY_CHECK_PATTERN)
            dataframe_numeric_count = _count_matches(
                test_function.lines, DATAFRAME_NUMERIC_PATTERN
            )
            large_comparison_count = _count_matches(
                test_function.lines, LARGE_COMPARISON_PATTERN
            )

            score = _score_cluster(
                assert_count=test_function.assert_count,
                payload_assert_count=payload_assert_count,
                key_check_count=key_check_count,
                dataframe_numeric_count=dataframe_numeric_count,
                large_comparison_count=large_comparison_count,
            )
            if not _is_candidate(
                score=score,
                assert_count=test_function.assert_count,
                payload_assert_count=payload_assert_count,
                key_check_count=key_check_count,
                dataframe_numeric_count=dataframe_numeric_count,
                large_comparison_count=large_comparison_count,
            ):
                continue

            stability_assessment, stability_note = _assess_stability(test_function)
            reason = _build_reason(
                assert_count=test_function.assert_count,
                payload_assert_count=payload_assert_count,
                key_check_count=key_check_count,
                dataframe_numeric_count=dataframe_numeric_count,
                large_comparison_count=large_comparison_count,
            )
            candidates.append(
                Candidate(
                    relpath=test_function.relpath,
                    test_name=test_function.name,
                    line_range=f"{test_function.start_line}-{test_function.end_line}",
                    score=score,
                    assert_count=test_function.assert_count,
                    payload_assert_count=payload_assert_count,
                    key_check_count=key_check_count,
                    dataframe_numeric_count=dataframe_numeric_count,
                    large_comparison_count=large_comparison_count,
                    reason=reason,
                    stability_assessment=stability_assessment,
                    stability_note=stability_note,
                )
            )

    candidates.sort(
        key=lambda item: (
            {"public/stable": 0, "mixed/needs review": 1, "internal/unstable": 2}.get(
                item.stability_assessment, 3
            ),
            -item.score,
            item.relpath,
            item.line_range,
        )
    )
    return candidates


def _render_report(
    *,
    primary_path: Path | None,
    diagnostic_files: list[Path],
    scanned_test_files: list[Path],
    candidates: list[Candidate],
) -> str:
    primary_display = (
        primary_path.relative_to(REPO_ROOT).as_posix()
        if primary_path is not None
        else "(not found)"
    )

    lines: list[str] = []
    lines.append("# Diagnostic Payload Test Report")
    lines.append("")
    lines.append(
        "Generated by `python tools/testing/find_diagnostic_assertion_clusters.py`."
    )
    lines.append("")
    lines.append(
        "> Candidate report only. Tests with dense assertions are not automatically wrong."
    )
    lines.append("")
    lines.append(
        "> Caution: snapshot/fixture-style assertions should be reserved for stable public/scientific payload contracts, not volatile internal diagnostics."
    )
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append(f"- Total `test*.py` files scanned: {len(scanned_test_files)}")
    lines.append(f"- Diagnostic-focused files scanned: {len(diagnostic_files)}")
    lines.append(f"- Primary target path used: `{primary_display}`")
    lines.append("")
    lines.append("### Diagnostic Files Used")
    lines.append("")
    for path in diagnostic_files:
        lines.append(f"- `{path.relative_to(REPO_ROOT).as_posix()}`")
    lines.append("")

    if not candidates:
        lines.append("## Candidate Clusters")
        lines.append("")
        lines.append(
            "No high-density diagnostic assertion clusters were detected with current heuristics."
        )
        lines.append("")
        return "\n".join(lines)

    stable = [c for c in candidates if c.stability_assessment == "public/stable"]
    mixed = [c for c in candidates if c.stability_assessment == "mixed/needs review"]
    unstable = [c for c in candidates if c.stability_assessment == "internal/unstable"]

    lines.append("## Candidate Clusters")
    lines.append("")
    lines.append("| File | Test | Approx Lines | Why Candidate | Stability |")
    lines.append("| --- | --- | ---: | --- | --- |")
    for candidate in candidates:
        lines.append(
            "| "
            f"{candidate.relpath} | "
            f"{candidate.test_name} | "
            f"{candidate.line_range} | "
            f"{candidate.reason} | "
            f"{candidate.stability_assessment} |"
        )
    lines.append("")

    def _section(title: str, items: list[Candidate]) -> None:
        lines.append(f"## {title}")
        lines.append("")
        if not items:
            lines.append("- None identified by current heuristics.")
            lines.append("")
            return
        for item in items:
            lines.append(
                f"- `{item.relpath}` :: `{item.test_name}` (`lines {item.line_range}`): {item.stability_note}"
            )
        lines.append("")

    _section("Likely Stable-Contract Candidates", stable)
    _section("Mixed Candidates", mixed)
    _section("Likely Internal/Unstable Candidates", unstable)

    return "\n".join(lines)


def main() -> None:
    primary_path, diagnostic_files, scanned_test_files = (
        _discover_primary_and_diagnostic_files()
    )
    candidates = _collect_candidates(diagnostic_files)
    report = _render_report(
        primary_path=primary_path,
        diagnostic_files=diagnostic_files,
        scanned_test_files=scanned_test_files,
        candidates=candidates,
    )
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
