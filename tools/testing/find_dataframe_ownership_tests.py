from __future__ import annotations

import ast
import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class TestFunction:
    relpath: str
    name: str
    start_line: int
    end_line: int
    text: str


@dataclass(frozen=True)
class PatternGroup:
    key: str
    label: str
    detect: tuple[re.Pattern[str], ...]
    description: str


REPO_ROOT = Path(__file__).resolve().parents[2]
TESTS_ROOT = REPO_ROOT / "tests"
PRIMARY_TARGET = TESTS_ROOT / "unit" / "test_frame_ownership_policy.py"
OUTPUT_PATH = REPO_ROOT / "docs" / "testing" / "dataframe_ownership_test_report.md"

OWNERSHIP_FILE_KEYWORDS: tuple[re.Pattern[str], ...] = (
    re.compile(r"\b_from_owned\s*\("),
    re.compile(r"\b_borrow_[A-Za-z0-9_]+\s*\("),
    re.compile(r"\bto_dataframe\s*\("),
    re.compile(r"\bto_pandas\s*\("),
    re.compile(r"\bfingerprint_table\s*\("),
    re.compile(r"defensive_snapshot", re.IGNORECASE),
    re.compile(r"isolated_from", re.IGNORECASE),
    re.compile(r"copy=False"),
    re.compile(r"dataframe_deep"),
)

PATTERN_GROUPS: tuple[PatternGroup, ...] = (
    PatternGroup(
        key="export_isolation",
        label="Export isolation",
        detect=(
            re.compile(r"isolated_from_mutation", re.IGNORECASE),
            re.compile(r"public_export", re.IGNORECASE),
            re.compile(r"defensive_snapshot", re.IGNORECASE),
            re.compile(r"\bto_dataframe\s*\("),
            re.compile(r"\bto_pandas\s*\("),
        ),
        description=(
            "Public export mutations should not leak back into owned internal frames."
        ),
    ),
    PatternGroup(
        key="mutation_safety",
        label="Mutation safety",
        detect=(
            re.compile(r"caller_mutation", re.IGNORECASE),
            re.compile(r"after_build", re.IGNORECASE),
            re.compile(r"can_be_mutated_after_owned_transfer", re.IGNORECASE),
            re.compile(r"phospho\.iloc|site_metadata\.iloc"),
        ),
        description=(
            "Mutating caller-owned inputs after construction/build should not corrupt model state."
        ),
    ),
    PatternGroup(
        key="accessor_copy_behavior",
        label="Accessor copy behavior",
        detect=(
            re.compile(r"defensive_snapshot", re.IGNORECASE),
            re.compile(r"copy keyword", re.IGNORECASE),
            re.compile(r"unexpected keyword argument"),
            re.compile(r"\bcopy=False\b"),
            re.compile(r"\bnot\s+re(read|ad)"),
        ),
        description=(
            "Accessors should expose safe snapshots and reject legacy copy-control API."
        ),
    ),
    PatternGroup(
        key="constructor_copy_behavior",
        label="Constructor copy behavior",
        detect=(
            re.compile(r"\b_from_owned\s*\("),
            re.compile(r"boundary_copy_and_owned_transfer_modes", re.IGNORECASE),
            re.compile(r"\bis\s+not\b.*\bpred_mat\b|\bis\b.*\bpred_mat\b"),
        ),
        description=(
            "Public constructors copy inputs; internal owned constructors transfer by alias."
        ),
    ),
    PatternGroup(
        key="representative_result_type_matrix",
        label="Representative result type matrix",
        detect=(
            re.compile(r"result_table_properties", re.IGNORECASE),
            re.compile(r"DatasetPreprocessingReport"),
            re.compile(r"Kinase(Activity|Prediction|Scoring)Result"),
            re.compile(
                r"SignalomeWorkflowResult|SignalomeAssignments|SignalomeModules|KinaseNetwork"
            ),
        ),
        description=(
            "Cross-result contract matrix validating consistent ownership semantics across types."
        ),
    ),
)

HIGH_RISK_MARKERS: tuple[re.Pattern[str], ...] = (
    re.compile(r"\b_count_dataframe_deep_copies\b"),
    re.compile(r"\bdataframe_deep\s*==\s*0\b"),
    re.compile(r"\bdataframe_copy_churn_regression_budget\b"),
    re.compile(r"\bfingerprint_table\b"),
    re.compile(r"\b_borrow_[A-Za-z0-9_]+\b"),
    re.compile(r"\bvalidator\b.*\bdoes_not_mutate\b", re.IGNORECASE),
)

RESULT_TYPES: tuple[str, ...] = (
    "AnalysisReadyPhosphoDataset",
    "DatasetPreprocessingReport",
    "KinasePredictionResult",
    "KinaseScoringResult",
    "KinaseActivityResult",
    "KinaseWorkflowResult",
    "SignalomeWorkflowResult",
    "SignalomeAssignments",
    "SignalomeModules",
    "KinaseNetwork",
    "ReferenceBundle",
    "TableSchema",
    "PhosphoIntensityMatrix",
)


def _discover_test_files() -> list[Path]:
    return sorted(TESTS_ROOT.rglob("test*.py"), key=lambda path: path.as_posix())


def _load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8-sig", errors="replace")


def _ownership_score(content: str) -> int:
    return sum(len(pattern.findall(content)) for pattern in OWNERSHIP_FILE_KEYWORDS)


def _discover_ownership_files(test_files: list[Path]) -> tuple[Path | None, list[Path]]:
    primary: Path | None = PRIMARY_TARGET if PRIMARY_TARGET.exists() else None
    discovered: list[Path] = []

    for path in test_files:
        rel = path.relative_to(REPO_ROOT).as_posix().lower()
        score = _ownership_score(_load_text(path))
        if "ownership" in path.name.lower() or "frame_ownership" in rel:
            discovered.append(path)
            continue
        if score >= 6:
            discovered.append(path)

    if primary is None:
        fallback = [
            path
            for path in discovered
            if "frame" in path.name.lower() and "ownership" in path.name.lower()
        ]
        if fallback:
            primary = sorted(fallback, key=lambda p: p.as_posix())[0]

    discovered = sorted(set(discovered), key=lambda path: path.as_posix())
    if primary is not None and primary not in discovered:
        discovered.append(primary)
        discovered.sort(key=lambda path: path.as_posix())
    return primary, discovered


def _extract_tests(path: Path) -> list[TestFunction]:
    source = _load_text(path)
    lines = source.splitlines()
    module = ast.parse(source)
    relpath = path.relative_to(REPO_ROOT).as_posix()

    tests: list[TestFunction] = []
    for node in module.body:
        if not isinstance(node, ast.FunctionDef):
            continue
        if not node.name.startswith("test_"):
            continue
        if node.end_lineno is None:
            continue
        start = node.lineno
        end = node.end_lineno
        block = "\n".join(lines[start - 1 : end])
        tests.append(
            TestFunction(
                relpath=relpath,
                name=node.name,
                start_line=start,
                end_line=end,
                text=block,
            )
        )
    return tests


def _detect_groups(test: TestFunction) -> set[str]:
    detected: set[str] = set()
    text = f"{test.name}\n{test.text}"
    for group in PATTERN_GROUPS:
        if any(pattern.search(text) for pattern in group.detect):
            detected.add(group.key)
    return detected


def _detect_result_types(test: TestFunction) -> set[str]:
    found: set[str] = set()
    text = test.text
    for result_type in RESULT_TYPES:
        if re.search(rf"\b{re.escape(result_type)}\b", text):
            found.add(result_type)
    return found


def _is_high_risk(test: TestFunction) -> bool:
    text = f"{test.name}\n{test.text}"
    return any(pattern.search(text) for pattern in HIGH_RISK_MARKERS)


def _render_report(
    *,
    primary: Path | None,
    scanned_files: list[Path],
    ownership_files: list[Path],
    tests: list[TestFunction],
    groups_by_test: dict[tuple[str, str, int], set[str]],
    types_by_test: dict[tuple[str, str, int], set[str]],
) -> str:
    primary_display = (
        primary.relative_to(REPO_ROOT).as_posix()
        if primary is not None
        else "(not found)"
    )

    lines: list[str] = []
    lines.append("# DataFrame Ownership Test Report")
    lines.append("")
    lines.append(
        "Generated by `python tools/testing/find_dataframe_ownership_tests.py`."
    )
    lines.append("")
    lines.append(
        "> Candidate consolidation report only. No test behavior changes are implied."
    )
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append(f"- Total `test*.py` files scanned: {len(scanned_files)}")
    lines.append(f"- Ownership/copy-policy files discovered: {len(ownership_files)}")
    lines.append(f"- Primary target path used: `{primary_display}`")
    lines.append("")
    lines.append("### Discovered Files")
    lines.append("")
    for path in ownership_files:
        lines.append(f"- `{path.relative_to(REPO_ROOT).as_posix()}`")
    lines.append("")

    group_rows: list[tuple[PatternGroup, list[TestFunction]]] = []
    for group in PATTERN_GROUPS:
        matched_tests = [
            test
            for test in tests
            if group.key in groups_by_test[(test.relpath, test.name, test.start_line)]
        ]
        if matched_tests:
            group_rows.append((group, matched_tests))

    lines.append("## Repeated Pattern Groups")
    lines.append("")
    lines.append("| Group | Test Count | Distinct Files | Notes |")
    lines.append("| --- | ---: | ---: | --- |")
    for group, matched in group_rows:
        files = {test.relpath for test in matched}
        lines.append(
            f"| {group.label} | {len(matched)} | {len(files)} | {group.description} |"
        )
    if not group_rows:
        lines.append("| _No repeated groups detected_ | 0 | 0 | n/a |")
    lines.append("")

    lines.append("## Example Clusters")
    lines.append("")
    for group, matched in group_rows:
        lines.append(f"### {group.label}")
        lines.append("")
        for test in sorted(
            matched,
            key=lambda t: (t.relpath, t.start_line, t.name),
        )[:8]:
            lines.append(
                f"- `{test.relpath}` :: `{test.name}` (`lines {test.start_line}-{test.end_line}`)"
            )
        if len(matched) > 8:
            lines.append(f"- ... and {len(matched) - 8} more tests")
        lines.append("")

    type_counts: dict[str, int] = {result_type: 0 for result_type in RESULT_TYPES}
    type_files: dict[str, set[str]] = {
        result_type: set() for result_type in RESULT_TYPES
    }
    for test in tests:
        key = (test.relpath, test.name, test.start_line)
        for result_type in types_by_test[key]:
            type_counts[result_type] += 1
            type_files[result_type].add(test.relpath)

    lines.append("## Result Types Covered")
    lines.append("")
    lines.append("| Result/Data Type | Tests Mentioning Type | Files |")
    lines.append("| --- | ---: | ---: |")
    for result_type in RESULT_TYPES:
        if type_counts[result_type] == 0:
            continue
        lines.append(
            f"| {result_type} | {type_counts[result_type]} | {len(type_files[result_type])} |"
        )
    lines.append("")

    high_risk_tests = [test for test in tests if _is_high_risk(test)]
    high_risk_tests.sort(key=lambda t: (t.relpath, t.start_line, t.name))

    lines.append("## High-Risk Paths To Keep Targeted")
    lines.append("")
    lines.append(
        "These paths protect copy-budget, borrow-alias, and non-mutation guarantees and should not be removed casually:"
    )
    lines.append("")
    for test in high_risk_tests[:15]:
        lines.append(
            f"- `{test.relpath}` :: `{test.name}` (`lines {test.start_line}-{test.end_line}`)"
        )
    if len(high_risk_tests) > 15:
        lines.append(
            f"- ... and {len(high_risk_tests) - 15} additional high-risk tests"
        )
    lines.append("")

    lines.append("## Candidate Generic Contract Matrix")
    lines.append("")
    lines.append(
        "| Contract Row | Minimal Representative Types | Keep Targeted Extras |"
    )
    lines.append("| --- | --- | --- |")
    lines.append(
        "| Constructor input isolation (public ctor copies) | AnalysisReadyPhosphoDataset, KinasePredictionResult, SignalomeWorkflowResult | Keep one per workflow boundary where builder/interpreter handoff semantics differ |"
    )
    lines.append(
        "| Owned transfer aliasing (`_from_owned`) | KinasePredictionResult, KinaseScoringResult, KinaseActivityResult | Keep borrow-path alias checks that guard internal zero-copy guarantees |"
    )
    lines.append(
        "| Public accessor defensive snapshots | AnalysisReadyPhosphoDataset, DatasetPreprocessingReport, KinaseActivityResult, SignalomeWorkflowResult | Keep table-specific cases with nullable/optional outputs (`nodes`, `candidate_correlations`) |"
    )
    lines.append(
        "| Legacy `copy=` keyword rejection | Dataset, prediction/scoring/activity result exports, signalome/table exports | Keep at least one integration-level check for each public export family |"
    )
    lines.append(
        "| Read-path mutation safety via fingerprints | Signalome validator/interpreter borrowed-read paths | Keep full targeted checks (high-risk) due regression blast radius |"
    )
    lines.append("")

    lines.append("## Follow-Up Consolidation Plan")
    lines.append("")
    lines.append(
        "1. Build shared parametrized helpers for export-isolation and accessor-snapshot assertions across result types."
    )
    lines.append(
        "2. Reduce repeated per-field mutation checks to representative table subsets per result family."
    )
    lines.append(
        "3. Preserve explicit high-risk tests for borrow aliasing, copy-churn budgets, and fingerprint non-mutation guarantees."
    )
    lines.append(
        "4. Track consolidation in dedicated ticket(s) and require before/after coverage mapping for removed duplicates."
    )
    lines.append("")

    return "\n".join(lines)


def main() -> None:
    scanned_files = _discover_test_files()
    primary, ownership_files = _discover_ownership_files(scanned_files)

    tests: list[TestFunction] = []
    for path in ownership_files:
        tests.extend(_extract_tests(path))
    tests.sort(key=lambda t: (t.relpath, t.start_line, t.name))

    groups_by_test: dict[tuple[str, str, int], set[str]] = {}
    types_by_test: dict[tuple[str, str, int], set[str]] = {}
    for test in tests:
        key = (test.relpath, test.name, test.start_line)
        groups_by_test[key] = _detect_groups(test)
        types_by_test[key] = _detect_result_types(test)

    report = _render_report(
        primary=primary,
        scanned_files=scanned_files,
        ownership_files=ownership_files,
        tests=tests,
        groups_by_test=groups_by_test,
        types_by_test=types_by_test,
    )
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
