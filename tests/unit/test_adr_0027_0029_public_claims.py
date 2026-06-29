"""Repository-wide ADR-0027/ADR-0029 public terminology guardrails."""

from __future__ import annotations

import ast
import re
from dataclasses import dataclass
from pathlib import Path
from re import Pattern

ROOT = Path(__file__).resolve().parents[2]
README = ROOT / "README.md"
DOCS_ROOT = ROOT / "docs"
SRC_ROOT = ROOT / "src" / "phospy"
EXAMPLES_ROOT = ROOT / "examples"
TESTS_ROOT = ROOT / "tests"


@dataclass(frozen=True)
class TextBlock:
    path: Path
    line: int
    text: str


@dataclass(frozen=True)
class PhraseRule:
    phrase: str
    pattern: Pattern[str]


def _rx(pattern: str) -> Pattern[str]:
    return re.compile(pattern, re.IGNORECASE)


UNSUPPORTED_PUBLIC_CLAIMS: tuple[PhraseRule, ...] = (
    PhraseRule("PhosR-equivalent", _rx(r"\bphosr[-\s]+equivalent\b")),
    PhraseRule("PhosR equivalent", _rx(r"\bphosr\s+equivalent\b")),
    PhraseRule("SPS/RUV-III parity", _rx(r"\bsps[/\s-]+ruv[-\s]*iii\s+parity\b")),
    PhraseRule("SPS RUV III parity", _rx(r"\bsps\s+ruv\s+iii\s+parity\b")),
    PhraseRule("RUV-III parity", _rx(r"\bruv[-\s]*iii\s+parity\b")),
    PhraseRule("RUV III parity", _rx(r"\bruv\s+iii\s+parity\b")),
    PhraseRule("full RUV support", _rx(r"\bfull\s+ruv\s+support\b")),
    PhraseRule("general RUV support", _rx(r"\bgeneral\s+ruv\s+support\b")),
    PhraseRule(
        "automatic control-site selection",
        _rx(r"\bautomatic\s+control[-\s]+site\s+selection\b"),
    ),
    PhraseRule(
        "automatic control site selection",
        _rx(r"\bautomatic\s+control\s+site\s+selection\b"),
    ),
    PhraseRule("built-in control set", _rx(r"\bbuilt[-\s]+in\s+control\s+set\b")),
    PhraseRule("bundled control set", _rx(r"\bbundled\s+control\s+set\b")),
    PhraseRule("online control lookup", _rx(r"\bonline\s+control\s+lookup\b")),
)

NEGATIVE_SCOPE_QUALIFIERS: tuple[Pattern[str], ...] = (
    _rx(r"\bnot\b"),
    _rx(r"\bunsupported\b"),
    _rx(r"\bdoes\s+not\b"),
    _rx(r"\bno\b"),
    _rx(r"\bnon[-\s]+equivalent\b"),
    _rx(r"\bnon[-\s]+parity\b"),
    _rx(r"\bnot\s+phosr[-\s]+equivalent\b"),
    _rx(r"\bnot\s+ruv[-\s]*iii\b"),
    _rx(r"\bnot\s+replicate-aware\s+ruv[-\s]*iii\b"),
)

REPORT_ONLY_RUV_READINESS = (
    "report-only",
    "report-only ruv-readiness metadata",
    "metadata-only",
    "does not modify the matrix",
    "not executable correction",
    "readiness only",
)
FIXED_EFFECT_LINEAR_RESIDUALIZE = (
    "fixed-effect",
    "fixed effect",
    "residualisation",
    "residualization",
    "not ruv",
    "not sps/ruv-style",
)
NATIVE_SPS_RUV_CONFIG = (
    "native phospy sps/ruv-style",
    "native sps/ruv-style",
    "phospy-native sps/ruv-style",
    "not phosr-equivalent",
    "not ruv-iii",
)
AMBIGUOUS_EXECUTABLE_TEMPORARY_IMPUTATION = _rx(
    r"\bexecutable\s+temporary[-\s]+imputation\b"
)
PUBLIC_NATIVE_MISSINGNESS_REJECTION_RULES: tuple[Pattern[str], ...] = (
    _rx(r"\bpublic\s+native\b.{0,180}\bcomplete\s+correction-stage\s+matrix\b"),
    _rx(r"\brejects\s+actual\s+missing\s+values\s+\(nans?\)"),
)

RAW_PUBLIC_BANS: tuple[PhraseRule, ...] = (
    PhraseRule("use_ruv", _rx(r"\buse_ruv\b")),
    PhraseRule("RUV-compatible", _rx(r"\bruv[-\s]+compatible\b")),
    PhraseRule(
        "control_site_ruv_style",
        _rx(r"\bcontrol_site_ruv_style\b"),
    ),
)


def _public_doc_paths() -> tuple[Path, ...]:
    return (README, *sorted(DOCS_ROOT.rglob("*.md")))


def _public_python_paths() -> tuple[Path, ...]:
    roots = (
        SRC_ROOT / "contracts" / "configs",
        SRC_ROOT / "workflows",
        SRC_ROOT / "science",
    )
    return tuple(path for root in roots for path in sorted(root.rglob("*.py")))


def _public_raw_paths() -> tuple[Path, ...]:
    examples = (
        tuple(sorted(EXAMPLES_ROOT.rglob("*.py")))
        + tuple(sorted(EXAMPLES_ROOT.rglob("*.md")))
        if EXAMPLES_ROOT.exists()
        else ()
    )
    return (*_public_doc_paths(), *_public_python_paths(), *examples)


def _relative_path(path: Path) -> str:
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def _line_for_offset(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def _normalise(text: str) -> str:
    text = text.replace("`", "")
    text = text.replace("*", "")
    return " ".join(text.split()).casefold()


def _sentence_blocks_from_text(path: Path, source: str) -> tuple[TextBlock, ...]:
    blocks: list[TextBlock] = []
    for match in re.finditer(r"[^.!?\n]+(?:[.!?]+|\n|$)", source):
        sentence = match.group(0).strip()
        if sentence:
            blocks.append(
                TextBlock(
                    path=path,
                    line=_line_for_offset(source, match.start()),
                    text=_normalise(sentence),
                )
            )
    return tuple(blocks)


def _markdown_prose_blocks(path: Path) -> tuple[TextBlock, ...]:
    source = path.read_text(encoding="utf-8")
    blocks: list[TextBlock] = []
    current: list[str] = []
    start_line: int | None = None
    in_fence = False

    def flush() -> None:
        nonlocal current, start_line
        if current and start_line is not None:
            blocks.append(TextBlock(path, start_line, _normalise(" ".join(current))))
        current = []
        start_line = None

    for line_number, raw_line in enumerate(source.splitlines(), start=1):
        stripped = raw_line.strip()
        if stripped.startswith("```"):
            flush()
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        if not stripped:
            flush()
            continue
        if stripped.startswith("|") or stripped.startswith(("-", "*")):
            flush()
            blocks.append(TextBlock(path, line_number, _normalise(stripped)))
            continue
        if start_line is None:
            start_line = line_number
        current.append(stripped)

    flush()
    return tuple(blocks)


def _python_prose_blocks(path: Path) -> tuple[TextBlock, ...]:
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(path))
    blocks: list[TextBlock] = []

    def append_docstring(node: ast.AST) -> None:
        docstring = ast.get_docstring(node)
        if not docstring:
            return
        line = getattr(node, "lineno", 1)
        blocks.append(TextBlock(path, line, _normalise(docstring)))

    append_docstring(tree)
    for node in ast.walk(tree):
        if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            if getattr(node, "name", "").startswith("_"):
                continue
            append_docstring(node)
    return tuple(blocks)


def _public_prose_blocks() -> tuple[TextBlock, ...]:
    markdown = tuple(
        block for path in _public_doc_paths() for block in _markdown_prose_blocks(path)
    )
    python = tuple(
        block for path in _public_python_paths() for block in _python_prose_blocks(path)
    )
    return (*markdown, *python)


def _has_negative_scope(text: str) -> bool:
    return any(pattern.search(text) for pattern in NEGATIVE_SCOPE_QUALIFIERS)


def _offender(path: Path, line: int, phrase: str, text: str) -> str:
    return f"{_relative_path(path)}:{line}: {phrase}: {text}"


def _find_positive_unsupported_claims(blocks: tuple[TextBlock, ...]) -> tuple[str, ...]:
    failures: list[str] = []
    for block in blocks:
        for sentence in _sentence_blocks_from_text(block.path, block.text):
            for rule in UNSUPPORTED_PUBLIC_CLAIMS:
                if not rule.pattern.search(sentence.text):
                    continue
                if _has_negative_scope(sentence.text):
                    continue
                failures.append(
                    _offender(block.path, block.line, rule.phrase, sentence.text)
                )
    return tuple(failures)


def _find_raw_banned_terms(paths: tuple[Path, ...]) -> tuple[str, ...]:
    failures: list[str] = []
    for path in paths:
        source = path.read_text(encoding="utf-8")
        for rule in RAW_PUBLIC_BANS:
            for match in rule.pattern.finditer(source):
                failures.append(
                    _offender(
                        path,
                        _line_for_offset(source, match.start()),
                        rule.phrase,
                        match.group(0),
                    )
                )
    return tuple(failures)


def _has_nearby_qualifier(
    text: str,
    occurrence: re.Match[str],
    accepted_phrases: tuple[str, ...],
) -> bool:
    start = max(0, occurrence.start() - 220)
    end = min(len(text), occurrence.end() + 220)
    window = text[start:end]
    return any(phrase in window for phrase in accepted_phrases)


def _find_unqualified_mentions(
    blocks: tuple[TextBlock, ...],
    mention_pattern: Pattern[str],
    name: str,
    accepted_phrases: tuple[str, ...],
) -> tuple[str, ...]:
    failures: list[str] = []
    for block in blocks:
        for occurrence in mention_pattern.finditer(block.text):
            if _has_nearby_qualifier(block.text, occurrence, accepted_phrases):
                continue
            failures.append(_offender(block.path, block.line, name, block.text))
    return tuple(failures)


def _has_public_native_missingness_rejection_rule(text: str) -> bool:
    return all(
        pattern.search(text) for pattern in PUBLIC_NATIVE_MISSINGNESS_REJECTION_RULES
    )


def _find_ambiguous_executable_temporary_imputation_mentions(
    blocks: tuple[TextBlock, ...],
) -> tuple[str, ...]:
    failures: list[str] = []
    for block in blocks:
        if not AMBIGUOUS_EXECUTABLE_TEMPORARY_IMPUTATION.search(block.text):
            continue
        if _has_public_native_missingness_rejection_rule(block.text):
            continue
        failures.append(
            _offender(
                block.path,
                block.line,
                "executable temporary imputation",
                block.text,
            )
        )
    return tuple(failures)


def _test_function_for_line(path: Path, line_number: int) -> str:
    source_lines = path.read_text(encoding="utf-8").splitlines()
    for index in range(line_number - 1, -1, -1):
        match = re.match(r"def\s+(test_[a-zA-Z0-9_]+)\(", source_lines[index])
        if match:
            return match.group(1)
    return ""


def test_public_docs_and_source_do_not_make_unsupported_positive_claims() -> None:
    failures = _find_positive_unsupported_claims(_public_prose_blocks())

    assert not failures, "\n".join(failures)


def test_public_docs_examples_and_source_do_not_expose_removed_ruv_shortcuts() -> None:
    failures = _find_raw_banned_terms(_public_raw_paths())

    assert not failures, "\n".join(failures)


def test_public_mentions_keep_ruv_readiness_report_only() -> None:
    failures = _find_unqualified_mentions(
        _public_prose_blocks(),
        _rx(r"\bruv_readiness\b"),
        "ruv_readiness",
        REPORT_ONLY_RUV_READINESS,
    )

    assert not failures, "\n".join(failures)


def test_public_mentions_keep_linear_residualize_batch_fixed_effect() -> None:
    failures = _find_unqualified_mentions(
        _public_prose_blocks(),
        _rx(r"\blinear_residualize_batch\b"),
        "linear_residualize_batch",
        FIXED_EFFECT_LINEAR_RESIDUALIZE,
    )

    assert not failures, "\n".join(failures)


def test_public_mentions_keep_sps_ruv_config_native_or_non_parity_scoped() -> None:
    failures = _find_unqualified_mentions(
        _public_prose_blocks(),
        _rx(r"\bspsruvbatchcorrectionconfig\b"),
        "SpsRuvBatchCorrectionConfig",
        NATIVE_SPS_RUV_CONFIG,
    )

    assert not failures, "\n".join(failures)


def test_public_mentions_do_not_use_ambiguous_temporary_imputation_wording() -> None:
    failures = _find_ambiguous_executable_temporary_imputation_mentions(
        _public_prose_blocks()
    )

    assert not failures, "\n".join(failures)


def test_removed_internal_method_alias_is_absent_outside_rejection_tests() -> None:
    alias = "control_site" + "_ruv_style"
    alias_pattern = re.compile(rf"\b{re.escape(alias)}\b")
    search_paths = (
        *_public_doc_paths(),
        *_public_python_paths(),
        *sorted(TESTS_ROOT.rglob("*.py")),
    )
    failures: list[str] = []
    for path in search_paths:
        if path == Path(__file__).resolve():
            continue
        source = path.read_text(encoding="utf-8")
        for match in alias_pattern.finditer(source):
            line = _line_for_offset(source, match.start())
            if path.is_relative_to(TESTS_ROOT):
                test_name = _test_function_for_line(path, line)
                if "reject" in test_name or "rejection" in test_name:
                    continue
            failures.append(_offender(path, line, alias, match.group(0)))

    assert not failures, "\n".join(failures)


def test_adr_0027_0029_guard_rejects_synthetic_violations() -> None:
    blocks = tuple(
        TextBlock(Path("synthetic.md"), index, _normalise(text))
        for index, text in enumerate(
            (
                "PhosPy provides PhosR-equivalent SPS/RUV-III parity.",
                "ruv_readiness applies executable correction to the matrix.",
                "ruv_readiness produces RUV-compatible preprocessing.",
                "linear_residualize_batch is SPS/RUV-style correction.",
                "SpsRuvBatchCorrectionConfig runs correction.",
                "Executable temporary imputation labels include row_median_temporary.",
            ),
            start=1,
        )
    )

    assert _find_positive_unsupported_claims(blocks)
    assert _find_unqualified_mentions(
        blocks, _rx(r"\bruv_readiness\b"), "ruv_readiness", REPORT_ONLY_RUV_READINESS
    )
    assert _find_unqualified_mentions(
        blocks,
        _rx(r"\blinear_residualize_batch\b"),
        "linear_residualize_batch",
        FIXED_EFFECT_LINEAR_RESIDUALIZE,
    )
    assert _find_unqualified_mentions(
        blocks,
        _rx(r"\bspsruvbatchcorrectionconfig\b"),
        "SpsRuvBatchCorrectionConfig",
        NATIVE_SPS_RUV_CONFIG,
    )
    assert _find_ambiguous_executable_temporary_imputation_mentions(blocks)


def test_adr_0027_0029_guard_allows_compliant_wording() -> None:
    blocks = tuple(
        TextBlock(Path("synthetic.md"), index, _normalise(text))
        for index, text in enumerate(
            (
                "PhosR-equivalent SPS/RUV-III parity is not supported.",
                "ruv_readiness is report-only RUV-readiness metadata and does not modify the matrix.",
                "linear_residualize_batch is fixed-effect residualisation, not SPS/RUV-style correction.",
                "SpsRuvBatchCorrectionConfig provides native PhosPy SPS/RUV-style correction, not RUV-III.",
                "Executable temporary imputation labels include row_median_temporary, but the public native workflow requires a complete correction-stage matrix and rejects actual missing values (NaNs) before executor invocation.",
            ),
            start=1,
        )
    )

    assert not _find_positive_unsupported_claims(blocks)
    assert not _find_unqualified_mentions(
        blocks, _rx(r"\bruv_readiness\b"), "ruv_readiness", REPORT_ONLY_RUV_READINESS
    )
    assert not _find_unqualified_mentions(
        blocks,
        _rx(r"\blinear_residualize_batch\b"),
        "linear_residualize_batch",
        FIXED_EFFECT_LINEAR_RESIDUALIZE,
    )
    assert not _find_unqualified_mentions(
        blocks,
        _rx(r"\bspsruvbatchcorrectionconfig\b"),
        "SpsRuvBatchCorrectionConfig",
        NATIVE_SPS_RUV_CONFIG,
    )
    assert not _find_ambiguous_executable_temporary_imputation_mentions(blocks)
