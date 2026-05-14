from __future__ import annotations

import ast
import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class TestFunction:
    file_path: Path
    file_relpath: str
    name: str
    start_line: int
    end_line: int
    lines: tuple[str, ...]


@dataclass(frozen=True)
class Theme:
    key: str
    label: str
    kind: str  # "primitive" or "domain"
    patterns: tuple[re.Pattern[str], ...]
    recommendation: str
    ticket_hint: str


@dataclass(frozen=True)
class ThemeHit:
    theme_key: str
    file_relpath: str
    test_name: str
    test_start_line: int
    match_line: int
    snippet: str


REPO_ROOT = Path(__file__).resolve().parents[2]
TESTS_ROOT = REPO_ROOT / "tests"
OUTPUT_PATH = REPO_ROOT / "docs" / "testing" / "validator_test_pattern_report.md"
PRIMARY_VALIDATOR_PATH = TESTS_ROOT / "unit" / "test_validator_boundaries.py"


PRIMITIVE_THEMES: tuple[Theme, ...] = (
    Theme(
        key="optional_positive_integer",
        label="Optional Positive Integer",
        kind="primitive",
        patterns=(
            re.compile(r"\b(min_substrates|top_n_substrates|module_count)\b"),
            re.compile(r"must be greater than or equal to"),
            re.compile(r"require_optional_int_at_least"),
            re.compile(r"(minimum_observed_values|min_observed_values)"),
        ),
        recommendation=(
            "Consolidate into shared parameterized boundary-case matrices for "
            "None, wrong type, zero/negative, and valid positive values."
        ),
        ticket_hint="TST-FOLLOWUP-VAL-PRIMITIVE-INT",
    ),
    Theme(
        key="non_empty_string",
        label="Non-Empty String",
        kind="primitive",
        patterns=(
            re.compile(r"non-empty string"),
            re.compile(r"must contain non-empty string values"),
            re.compile(r"must be non-empty"),
            re.compile(r"empty string"),
        ),
        recommendation=(
            "Group repeated string-empty/whitespace/type checks into one helper-driven "
            "test matrix per validator boundary."
        ),
        ticket_hint="TST-FOLLOWUP-VAL-PRIMITIVE-STRING",
    ),
    Theme(
        key="bounded_numeric_range",
        label="Bounded Numeric Range",
        kind="primitive",
        patterns=(
            re.compile(r"\b(threshold|cutoff|pseudocount)\b"),
            re.compile(r"must be finite"),
            re.compile(r"require_optional_real_between"),
            re.compile(r"must be greater than or equal to"),
            re.compile(r"must be .*?between"),
        ),
        recommendation=(
            "Centralize lower/upper-bound range checks with one parameterized "
            "table per field family (probabilities, thresholds, pseudocounts)."
        ),
        ticket_hint="TST-FOLLOWUP-VAL-PRIMITIVE-RANGE",
    ),
    Theme(
        key="enum_validation",
        label="Enum / Supported Literal Validation",
        kind="primitive",
        patterns=(
            re.compile(r"must be one of"),
            re.compile(r"\bunsupported\b"),
            re.compile(r"\binvalid\b"),
            re.compile(r"coerce_policy_enum"),
            re.compile(r"require_supported_literal"),
        ),
        recommendation=(
            "Use shared enum-contract test helpers to validate unknown values and "
            "error-message quality consistently."
        ),
        ticket_hint="TST-FOLLOWUP-VAL-PRIMITIVE-ENUM",
    ),
    Theme(
        key="nullable_collection",
        label="Nullable Collection / Optional Field Coupling",
        kind="primitive",
        patterns=(
            re.compile(r"must be None"),
            re.compile(r"requires .*? when"),
            re.compile(r"\bpairs\b"),
            re.compile(r"\boptional\b"),
        ),
        recommendation=(
            "Consolidate policy-coupled optional fields into parameterized cases "
            "that explicitly encode allowed/forbidden None combinations."
        ),
        ticket_hint="TST-FOLLOWUP-VAL-PRIMITIVE-NULLABLE",
    ),
    Theme(
        key="dataframe_shape",
        label="DataFrame Shape / Schema / Alignment",
        kind="primitive",
        patterns=(
            re.compile(r"\bDataFrame\b"),
            re.compile(r"missing required columns"),
            re.compile(r"shape must align"),
            re.compile(r"index must be unique"),
            re.compile(r"columns must be unique"),
            re.compile(r"non-numeric columns"),
            re.compile(r"shared_count="),
        ),
        recommendation=(
            "Collapse repeated dataframe validation checks into reusable schema/shape "
            "case sets while keeping domain-specific columns in dedicated tests."
        ),
        ticket_hint="TST-FOLLOWUP-VAL-PRIMITIVE-DATAFRAME",
    ),
    Theme(
        key="path_like_value",
        label="Path-Like Value",
        kind="primitive",
        patterns=(
            re.compile(r"file path"),
            re.compile(r"filesystem path"),
            re.compile(r"path-like"),
            re.compile(r"require_local_filesystem_path"),
        ),
        recommendation=(
            "Create a small shared suite for local-path acceptance vs URL/remote "
            "rejection semantics."
        ),
        ticket_hint="TST-FOLLOWUP-VAL-PRIMITIVE-PATH",
    ),
)


DOMAIN_THEMES: tuple[Theme, ...] = (
    Theme(
        key="mixed_total_protein_quantitative_meaning",
        label="Mixed Total-Protein Quantitative Meaning",
        kind="domain",
        patterns=(
            re.compile(r"mixed total-protein quantitative meaning"),
            re.compile(r"allow_mixed_total_protein_quantitative_meaning"),
            re.compile(r"unmatched_policy='allow_uncorrected'"),
        ),
        recommendation=(
            "Keep domain-specific assertions explicit but reduce duplicate setup by "
            "sharing mixed-dataset fixture builders across workflow validators."
        ),
        ticket_hint="TST-FOLLOWUP-VAL-DOMAIN-MIXED-TOTAL",
    ),
    Theme(
        key="signalome_config_boundary_contracts",
        label="Signalome Config Boundary Contracts",
        kind="domain",
        patterns=(
            re.compile(r"signalome workflow request config\."),
            re.compile(r"Signalome(Scientific|Output|Validation|Clustering)Config"),
            re.compile(r"module_selection_"),
        ),
        recommendation=(
            "Consolidate repeated invalid-config boundary checks with a signalome "
            "config invalid-case table, preserving clear field-level expectations."
        ),
        ticket_hint="TST-FOLLOWUP-VAL-DOMAIN-SIGNALOME-CONFIG",
    ),
    Theme(
        key="reference_compatibility_and_bundle_contracts",
        label="Reference Compatibility / Bundle Contracts",
        kind="domain",
        patterns=(
            re.compile(r"ReferenceCompatibilityError"),
            re.compile(r"ReferencePreset"),
            re.compile(r"ReferenceBundle"),
        ),
        recommendation=(
            "Keep this domain logic separate from primitive validators; consider one "
            "shared compatibility scenario matrix for request-validator boundaries."
        ),
        ticket_hint="TST-FOLLOWUP-VAL-DOMAIN-REFERENCE",
    ),
)


ALL_THEMES: tuple[Theme, ...] = PRIMITIVE_THEMES + DOMAIN_THEMES


def _discover_test_files() -> list[Path]:
    return sorted(TESTS_ROOT.rglob("test*.py"), key=lambda path: path.as_posix())


def _resolve_primary_validator_file() -> tuple[Path | None, str]:
    if PRIMARY_VALIDATOR_PATH.exists():
        return PRIMARY_VALIDATOR_PATH, PRIMARY_VALIDATOR_PATH.relative_to(
            REPO_ROOT
        ).as_posix()

    candidates = sorted(
        TESTS_ROOT.rglob("test*validator*boundar*.py"),
        key=lambda path: path.as_posix(),
    )
    if candidates:
        replacement = candidates[0]
        return replacement, replacement.relative_to(REPO_ROOT).as_posix()
    return None, "(not found)"


def _is_related_test_file(path: Path, primary_file: Path | None) -> bool:
    name = path.name.lower()
    relpath = path.relative_to(REPO_ROOT).as_posix().lower()
    if primary_file is not None and path == primary_file:
        return True
    if "validator" in name or "validation" in name:
        return True
    if "validator" in relpath or "validation" in relpath:
        return True
    return False


def _extract_test_functions(path: Path) -> list[TestFunction]:
    source = path.read_text(encoding="utf-8-sig", errors="replace")
    lines = source.splitlines()
    module = ast.parse(source)
    relpath = path.relative_to(REPO_ROOT).as_posix()

    test_functions: list[TestFunction] = []
    for node in module.body:
        if not isinstance(node, ast.FunctionDef):
            continue
        if not node.name.startswith("test_"):
            continue
        if node.end_lineno is None:
            continue
        start = node.lineno
        end = node.end_lineno
        snippet = tuple(lines[start - 1 : end])
        test_functions.append(
            TestFunction(
                file_path=path,
                file_relpath=relpath,
                name=node.name,
                start_line=start,
                end_line=end,
                lines=snippet,
            )
        )
    return test_functions


def _clean_snippet(line: str) -> str:
    compact = " ".join(line.strip().split())
    compact = compact.replace("|", "\\|")
    return compact[:160]


def _find_theme_hits(test_function: TestFunction, theme: Theme) -> list[ThemeHit]:
    hits: list[ThemeHit] = []
    for offset, line in enumerate(test_function.lines):
        for pattern in theme.patterns:
            if pattern.search(line):
                hits.append(
                    ThemeHit(
                        theme_key=theme.key,
                        file_relpath=test_function.file_relpath,
                        test_name=test_function.name,
                        test_start_line=test_function.start_line,
                        match_line=test_function.start_line + offset,
                        snippet=_clean_snippet(line),
                    )
                )
                break
    return hits


def _collect_hits(test_functions: list[TestFunction]) -> dict[str, list[ThemeHit]]:
    hits_by_theme: dict[str, list[ThemeHit]] = {theme.key: [] for theme in ALL_THEMES}
    for test_function in test_functions:
        for theme in ALL_THEMES:
            hits = _find_theme_hits(test_function, theme)
            if not hits:
                continue
            first_hit = min(hits, key=lambda hit: hit.match_line)
            hits_by_theme[theme.key].append(first_hit)

    for key in hits_by_theme:
        hits_by_theme[key].sort(
            key=lambda hit: (hit.file_relpath, hit.match_line, hit.test_name)
        )
    return hits_by_theme


def _count_hits_for_file(hits: list[ThemeHit], file_relpath: str) -> int:
    return sum(1 for hit in hits if hit.file_relpath == file_relpath)


def _render_theme_section(
    themes: tuple[Theme, ...],
    hits_by_theme: dict[str, list[ThemeHit]],
    primary_relpath: str,
) -> list[str]:
    lines: list[str] = []
    lines.append("| Theme | Total Hits | In Primary File | Distinct Files |")
    lines.append("| --- | ---: | ---: | ---: |")
    for theme in themes:
        hits = hits_by_theme[theme.key]
        if len(hits) < 2:
            continue
        distinct_files = len({hit.file_relpath for hit in hits})
        primary_hits = _count_hits_for_file(hits, primary_relpath)
        lines.append(
            f"| {theme.label} | {len(hits)} | {primary_hits} | {distinct_files} |"
        )
    if lines[-1].startswith("| ---"):
        lines.append("| _No repeated themes detected_ | 0 | 0 | 0 |")
    return lines


def _render_theme_details(theme: Theme, hits: list[ThemeHit]) -> list[str]:
    lines: list[str] = []
    lines.append(f"### {theme.label}")
    lines.append("")
    lines.append(f"- Recommendation: {theme.recommendation}")
    lines.append(f"- Follow-up ticket candidate: `{theme.ticket_hint}`")
    lines.append("")
    lines.append("| File | Test | Line | Evidence |")
    lines.append("| --- | --- | ---: | --- |")
    for hit in hits:
        lines.append(
            f"| {hit.file_relpath} | {hit.test_name} | {hit.match_line} | `{hit.snippet}` |"
        )
    lines.append("")
    return lines


def _build_report(
    *,
    primary_file: Path | None,
    primary_relpath: str,
    scanned_files: list[Path],
    related_files: list[Path],
    test_functions: list[TestFunction],
    hits_by_theme: dict[str, list[ThemeHit]],
) -> str:
    lines: list[str] = []
    lines.append("# Validator Test Pattern Report")
    lines.append("")
    lines.append("Generated by `python tools/testing/find_validator_test_patterns.py`.")
    lines.append("")
    lines.append(
        "> Candidate consolidation report only. This does not change tests and does not "
        "claim any individual test is incorrect."
    )
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append(f"- Total `test*.py` files scanned: {len(scanned_files)}")
    lines.append(f"- Related validator/validation files analyzed: {len(related_files)}")
    lines.append(f"- Test functions analyzed: {len(test_functions)}")
    if primary_file is not None:
        lines.append(f"- Primary boundary file used: `{primary_relpath}`")
    else:
        lines.append(
            "- Primary boundary file used: not found; fallback discovery did not locate an equivalent."
        )
    lines.append("")
    lines.append("## Shared Primitive Validation Patterns")
    lines.append("")
    lines.extend(
        _render_theme_section(PRIMITIVE_THEMES, hits_by_theme, primary_relpath)
    )
    lines.append("")
    lines.append("## Domain-Specific Validation Patterns")
    lines.append("")
    lines.extend(_render_theme_section(DOMAIN_THEMES, hits_by_theme, primary_relpath))
    lines.append("")

    repeated_primitive = [
        theme for theme in PRIMITIVE_THEMES if len(hits_by_theme[theme.key]) >= 2
    ]
    repeated_domain = [
        theme for theme in DOMAIN_THEMES if len(hits_by_theme[theme.key]) >= 2
    ]

    lines.append("## Consolidation Candidates")
    lines.append("")
    if not repeated_primitive and not repeated_domain:
        lines.append("No repeated patterns met the reporting threshold (>=2 hits).")
        lines.append("")
    else:
        lines.append(
            "These are candidate consolidations for follow-up tickets, not immediate rewrites."
        )
        lines.append("")
        for theme in repeated_primitive:
            lines.extend(_render_theme_details(theme, hits_by_theme[theme.key]))
        for theme in repeated_domain:
            lines.extend(_render_theme_details(theme, hits_by_theme[theme.key]))

    lines.append("## Suggested Follow-Up Tickets")
    lines.append("")
    if not repeated_primitive and not repeated_domain:
        lines.append("- No follow-up tickets proposed based on current threshold.")
    else:
        for theme in repeated_primitive + repeated_domain:
            lines.append(
                f"- `{theme.ticket_hint}`: Consolidation pass for {theme.label.lower()} tests."
            )
    lines.append("")

    return "\n".join(lines)


def main() -> None:
    scanned_files = _discover_test_files()
    primary_file, primary_relpath = _resolve_primary_validator_file()
    related_files = [
        path for path in scanned_files if _is_related_test_file(path, primary_file)
    ]

    test_functions: list[TestFunction] = []
    for path in related_files:
        test_functions.extend(_extract_test_functions(path))
    test_functions.sort(key=lambda fn: (fn.file_relpath, fn.start_line, fn.name))

    hits_by_theme = _collect_hits(test_functions)
    report = _build_report(
        primary_file=primary_file,
        primary_relpath=primary_relpath,
        scanned_files=scanned_files,
        related_files=related_files,
        test_functions=test_functions,
        hits_by_theme=hits_by_theme,
    )

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
