"""Guard public docs against unsupported scientific parity claims."""

from __future__ import annotations

import ast
import re
from dataclasses import dataclass
from pathlib import Path
from re import Pattern

ROOT = Path(__file__).resolve().parents[2]
DOCS_ROOT = ROOT / "docs"
README = ROOT / "README.md"
SRC_ROOT = ROOT / "src" / "phospy"


@dataclass(frozen=True)
class ClaimBlock:
    path: Path
    start_line: int
    text: str


@dataclass(frozen=True)
class SupportClaimRule:
    name: str
    unsupported_claim: Pattern[str]
    allowed_contexts: tuple[Pattern[str], ...]
    update_when_supported: str


def _rx(pattern: str) -> Pattern[str]:
    return re.compile(pattern, re.IGNORECASE)


NEGATED_OR_LIMITED = (
    _rx(r"\b(?:not|no|never|without|unsupported|outside|out of scope)\b"),
    _rx(r"\bdoes\s+not\b"),
    _rx(r"\bdo\s+not\b"),
    _rx(r"\bmust\s+not\b"),
    _rx(r"\bnot\s+(?:claimed|provided|implemented|supported|executable)\b"),
    _rx(r"\bnot\s+a\s+(?:claim|current feature-support claim)\b"),
    _rx(r"\bopen\s+gap\b"),
    _rx(r"\blimitation(?:s)?\b"),
    _rx(r"\broadmap\b"),
)

ROADMAP_OR_UNSUPPORTED = NEGATED_OR_LIMITED + (
    _rx(r"\bfuture\b"),
    _rx(r"\bunsafe\b"),
    _rx(r"\bnon-supported\b"),
    _rx(r"\bnot\s+equivalent\b"),
)


CLAIM_WORDS = (
    r"support(?:s|ed|ing)?",
    r"implement(?:s|ed|ation)?",
    r"provid(?:e|es|ed|ing)",
    r"bundl(?:e|es|ed|ing)",
    r"parity",
    r"equivalence",
    r"equivalent",
    r"compatib(?:le|ility)",
    r"replacement",
    r"clone",
    r"fitted",
    r"run(?:s|ning)?",
    r"execut(?:e|es|ed|able|ion)",
)
CLAIM_WORD = rf"(?:{'|'.join(CLAIM_WORDS)})"


SUPPORT_CLAIM_RULES: tuple[SupportClaimRule, ...] = (
    SupportClaimRule(
        name="global PhosR parity",
        unsupported_claim=_rx(
            r"\bphosr\b.{0,90}\bparity\b.{0,90}\b(?:support\w*|claim\w*|"
            r"provid\w*|full|global|complete|compatib\w*)\b|"
            r"\b(?:support\w*|claim\w*|provid\w*|full|global|complete)\b"
            r".{0,90}\bphosr\b.{0,90}\bparity\b|"
            r"\bphosr\b.{0,90}\b(?:equivalence|equivalent|compatib\w*|"
            r"replacement|clone)\b|"
            r"\b(?:equivalence|equivalent|compatib\w*|replacement|clone)\b"
            r".{0,90}\bphosr\b"
        ),
        allowed_contexts=(
            *NEGATED_OR_LIMITED,
            _rx(r"\b(?:phosr-style|phosr-inspired)\b"),
            _rx(r"\bfeature-(?:scoped|specific)\b"),
            _rx(r"\bevidence-scoped\b"),
            _rx(r"\bfixture-backed\b"),
            _rx(r"\bcomparison evidence\b"),
            _rx(r"\bselected phosr-style\b"),
            _rx(r"\bphosr-derived lineage\b"),
            _rx(r"\bupstream phosr project\b"),
            _rx(r"\bnon-phosr-equivalence\b"),
            _rx(r"\btests/parity/"),
            _rx(r"\bparity\.md\b"),
        ),
        update_when_supported=(
            "Only add a positive PhosR parity context when the support category "
            "and fixture-backed evidence scope are named explicitly."
        ),
    ),
    SupportClaimRule(
        name="broad limma parity",
        unsupported_claim=_rx(
            rf"\blimma\b.{{0,90}}\b{CLAIM_WORD}\b|"
            rf"\b{CLAIM_WORD}\b.{{0,90}}\blimma\b"
        ),
        allowed_contexts=(
            *NEGATED_OR_LIMITED,
            _rx(r"\binstead of implying\b"),
            _rx(r"\blimma-(?:style|trend|backed|envelope)\b"),
            _rx(r"\blimma\s+`?(?:duplicatecorrelation|removebatcheffect)`?\b"),
            _rx(r"\bcontract difference vs limma\b"),
            _rx(r"\bdifferential limma parity fixtures\b"),
            _rx(r"\btests/parity/test_differential_limma_parity\.py\b"),
            _rx(r"\btests/fixtures/rewrite_parity/differential_limma_envelope/\b"),
        ),
        update_when_supported=(
            "If a broader limma-compatible lane is added, list its exact model "
            "scope and tests here instead of allowing generic limma parity text."
        ),
    ),
    SupportClaimRule(
        name="MSstatsPTM equivalence",
        unsupported_claim=_rx(
            rf"\bmsstatsptm\b.{{0,90}}\b{CLAIM_WORD}\b|"
            rf"\b{CLAIM_WORD}\b.{{0,90}}\bmsstatsptm\b"
        ),
        allowed_contexts=(
            *NEGATED_OR_LIMITED,
            _rx(r"\bpreparation-only\b"),
            _rx(r"\bmodel-input preparation\b"),
        ),
        update_when_supported=(
            "If joint PTM/protein modelling is implemented, replace the "
            "negative-only allowance with the exact supported MSstatsPTM-style "
            "scope and evidence."
        ),
    ),
    SupportClaimRule(
        name="PTM-SEA support or parity",
        unsupported_claim=_rx(
            rf"\bptm-sea\b.{{0,90}}\b{CLAIM_WORD}\b|"
            rf"\b{CLAIM_WORD}\b.{{0,90}}\bptm-sea\b"
        ),
        allowed_contexts=(
            *NEGATED_OR_LIMITED,
            _rx(r"\bora\s+is\s+not\b"),
            _rx(r"\boffline over-representation analysis\b"),
            _rx(r"\bcaller-supplied\b"),
            _rx(r"\bordinary local collections\b"),
            _rx(r"\bseparate from\b"),
        ),
        update_when_supported=(
            "If PTM-SEA support is added, allow only a context that names the "
            "implemented method, resource provenance, and validation tests."
        ),
    ),
    SupportClaimRule(
        name="official Kinase Library bundled support or parity",
        unsupported_claim=_rx(
            rf"\bofficial\b.{{0,60}}\bkinase library\b.{{0,90}}\b{CLAIM_WORD}\b|"
            rf"\bkinase library\b.{{0,90}}\b(?:official|bundl\w*|parity|"
            rf"compatib\w*)\b.{{0,90}}\b{CLAIM_WORD}\b|"
            rf"\b{CLAIM_WORD}\b.{{0,90}}\b(?:official|bundl\w*|parity|"
            rf"compatib\w*)\b.{{0,90}}\bkinase library\b"
        ),
        allowed_contexts=(
            *NEGATED_OR_LIMITED,
            _rx(r"\bkinase library-style\b"),
            _rx(r"\bcaller-supplied\b"),
            _rx(r"\blocal resources?\b"),
            _rx(r"\bcompatible local `?kinaselibraryresource`?\b"),
            _rx(r"\bpure science-layer\b"),
            _rx(r"\bvalidated phospy implementation\b"),
            _rx(r"\bnot an official kinase library implementation\b"),
            _rx(r"\bdoes not bundle official kinase library data\b"),
            _rx(r"\bno official kinase library\b"),
        ),
        update_when_supported=(
            "If official bundled Kinase Library support is added, require a "
            "license/provenance statement and official parity evidence before "
            "allowing positive bundled/official language."
        ),
    ),
    SupportClaimRule(
        name="ComBat support",
        unsupported_claim=_rx(
            rf"\bcombat\b.{{0,90}}\b{CLAIM_WORD}\b|"
            rf"\b{CLAIM_WORD}\b.{{0,90}}\bcombat\b"
        ),
        allowed_contexts=NEGATED_OR_LIMITED
        + (
            _rx(r"\blinear_residualize_batch\b"),
            _rx(r"\bbatch-effect methods remain outside\b"),
        ),
        update_when_supported=(
            "If ComBat support is added, allow only text naming the exact "
            "algorithm, provenance, validation, and supported workflow boundary."
        ),
    ),
    SupportClaimRule(
        name="RUV-III support or implementation",
        unsupported_claim=_rx(
            rf"\b(?:ruv-iii|ruv_iii_style)\b.{{0,90}}\b{CLAIM_WORD}\b|"
            rf"\b{CLAIM_WORD}\b.{{0,90}}\b(?:ruv-iii|ruv_iii_style)\b"
        ),
        allowed_contexts=ROADMAP_OR_UNSUPPORTED
        + (
            _rx(r"\breplicate-aware ruv-iii\b"),
            _rx(r"\bmethod label\b"),
        ),
        update_when_supported=(
            "If RUV-III support is added, allow only text naming the exact "
            "replicate-aware method semantics, validation scope, and evidence."
        ),
    ),
    SupportClaimRule(
        name="RUV support",
        unsupported_claim=_rx(
            rf"\bruv\b.{{0,90}}\b{CLAIM_WORD}\b|"
            rf"\b{CLAIM_WORD}\b.{{0,90}}\bruv\b"
        ),
        allowed_contexts=NEGATED_OR_LIMITED
        + (
            _rx(r"\bruv-compatible\b"),
            _rx(r"\bnative sps/ruv-style\b"),
            _rx(r"\bspsruvbatchcorrectionconfig\b"),
            _rx(r"\bmetadata readiness reporting\b"),
            _rx(r"\bbatch-effect methods remain outside\b"),
        ),
        update_when_supported=(
            "Generic RUV support language must be replaced with the exact "
            "implemented native method, validation boundary, and evidence."
        ),
    ),
    SupportClaimRule(
        name="limma duplicateCorrelation support",
        unsupported_claim=_rx(
            rf"\bduplicatecorrelation\b.{{0,90}}\b{CLAIM_WORD}\b|"
            rf"\b{CLAIM_WORD}\b.{{0,90}}\bduplicatecorrelation\b"
        ),
        allowed_contexts=NEGATED_OR_LIMITED
        + (
            _rx(r"\bstyle\b"),
            _rx(r"\bcorrelated-replicate\b"),
            _rx(r"\bunsupported-design rejection policy\b"),
            _rx(r"\blinear_residualize_batch\b"),
        ),
        update_when_supported=(
            "If duplicateCorrelation-style modelling is implemented, allow only "
            "text tied to its explicit design/result contract and parity tests."
        ),
    ),
    SupportClaimRule(
        name="mixed-effects support",
        unsupported_claim=_rx(
            rf"\bmixed[-\s]effects?\b.{{0,90}}\b{CLAIM_WORD}\b|"
            rf"\b{CLAIM_WORD}\b.{{0,90}}\bmixed[-\s]effects?\b"
        ),
        allowed_contexts=NEGATED_OR_LIMITED
        + (
            _rx(r"\brandom subject\b"),
            _rx(r"\borderinary fixed effects\b"),
            _rx(r"\blinear_residualize_batch\b"),
            _rx(r"\bbefore adding\b"),
        ),
        update_when_supported=(
            "If mixed-effects modelling is implemented, allow only text naming "
            "the model family, design contract, and validation/parity evidence."
        ),
    ),
    SupportClaimRule(
        name="limma removeBatchEffect support",
        unsupported_claim=_rx(
            rf"\bremovebatcheffect\b.{{0,90}}\b{CLAIM_WORD}\b|"
            rf"\b{CLAIM_WORD}\b.{{0,90}}\bremovebatcheffect\b"
        ),
        allowed_contexts=NEGATED_OR_LIMITED
        + (
            _rx(r"\blinear_residualize_batch\b"),
            _rx(r"\bbroader batch-effect modelling\b"),
        ),
        update_when_supported=(
            "If removeBatchEffect parity is implemented, allow only text tied "
            "to that exact method, boundary, and evidence."
        ),
    ),
)


DISTINCTION_RULES: tuple[SupportClaimRule, ...] = (
    SupportClaimRule(
        name="PhosR-equivalent SPS/RUV-III claim",
        unsupported_claim=_rx(r"\bphosr-equivalent\s+sps/ruv-iii\b"),
        allowed_contexts=ROADMAP_OR_UNSUPPORTED,
        update_when_supported=(
            "Do not use PhosR-equivalent SPS/RUV-III wording as a positive "
            "claim unless parity evidence and scope are explicit."
        ),
    ),
    SupportClaimRule(
        name="stale future-only native correction lane",
        unsupported_claim=_rx(
            r"\bfuture\s+correction\b|\bwithout\s+numerical\s+correction\b"
        ),
        allowed_contexts=(),
        update_when_supported=(
            "Active correction APIs must describe native SPS/RUV-style "
            "preprocessing correction and resolved numerical execution plans."
        ),
    ),
    SupportClaimRule(
        name="generic RUV support shortcut",
        unsupported_claim=_rx(r"\bruv\s+support\b"),
        allowed_contexts=ROADMAP_OR_UNSUPPORTED
        + (
            _rx(r"\bnative\s+phospy\s+sps/ruv-style\b"),
            _rx(r"\bnative\s+sps/ruv-style\b"),
            _rx(r"\bspsruvbatchcorrectionconfig\b"),
        ),
        update_when_supported=(
            "Say native PhosPy SPS/RUV-style preprocessing correction, not "
            "generic RUV support."
        ),
    ),
    SupportClaimRule(
        name="boolean RUV shortcut",
        unsupported_claim=_rx(r"\buse_ruv\b"),
        allowed_contexts=ROADMAP_OR_UNSUPPORTED,
        update_when_supported=(
            "Do not introduce a boolean use_ruv public method or example; use "
            "the structured native SPS/RUV-style config instead."
        ),
    ),
    SupportClaimRule(
        name="automatic control-site selection",
        unsupported_claim=_rx(
            r"\b(?:automatic(?:ally)?|hidden|implicit(?:ly)?)\b.{0,80}"
            r"\b(?:control[-\s]?site|controls?)\b.{0,80}"
            r"\b(?:select(?:s|ed|ion)?|lookup|fetch(?:es|ed|ing)?)\b|"
            r"\b(?:control[-\s]?site|controls?)\b.{0,80}"
            r"\b(?:select(?:s|ed|ion)?|lookup|fetch(?:es|ed|ing)?)\b.{0,80}"
            r"\b(?:automatic(?:ally)?|hidden|implicit(?:ly)?)\b"
        ),
        allowed_contexts=ROADMAP_OR_UNSUPPORTED
        + (
            _rx(r"\bcaller-supplied\b"),
            _rx(r"\bexplicit\b"),
        ),
        update_when_supported=(
            "Control-site selection must be described as caller-supplied or "
            "as a separately implemented explicit method with provenance."
        ),
    ),
    SupportClaimRule(
        name="ruv_readiness as correction",
        unsupported_claim=_rx(
            r"\bruv_readiness\b.{0,120}\b(?:correct(?:s|ed|ion)?|apply|"
            r"applies|applied|modify|modifies|modified|matrix|selects?)\b|"
            r"\b(?:correct(?:s|ed|ion)?|apply|applies|applied|modify|"
            r"modifies|modified|matrix|selects?)\b.{0,120}\bruv_readiness\b"
        ),
        allowed_contexts=ROADMAP_OR_UNSUPPORTED
        + (
            _rx(r"\breport-only\b"),
            _rx(r"\bdiagnostic(?:s)?\b"),
            _rx(r"\bmetadata readiness\b"),
        ),
        update_when_supported=(
            "`ruv_readiness` must remain diagnostic/report-only language, not "
            "correction language."
        ),
    ),
    SupportClaimRule(
        name="linear residualisation as SPS/RUV-style correction",
        unsupported_claim=_rx(
            r"\blinear_residualize_batch\b.{0,140}\b(?:sps/ruv|ruv-iii|"
            r"native\s+sps/ruv-style|unwanted-factor|control[-\s]?site)\b|"
            r"\b(?:sps/ruv|ruv-iii|native\s+sps/ruv-style|unwanted-factor|"
            r"control[-\s]?site)\b.{0,140}\blinear_residualize_batch\b"
        ),
        allowed_contexts=ROADMAP_OR_UNSUPPORTED
        + (
            _rx(r"\bfixed-effect residualisation\b"),
            _rx(r"\bdistinct\b"),
        ),
        update_when_supported=(
            "`linear_residualize_batch` must be described as fixed-effect "
            "residualisation, not SPS/RUV-style correction."
        ),
    ),
    SupportClaimRule(
        name="differential batch covariates as preprocessing correction",
        unsupported_claim=_rx(
            r"\b(?:differential\b.{0,80}\bbatch\s+covariates?|"
            r"batch\s+covariates?\b.{0,80}\bdifferential|"
            r"differential\b.{0,80}\bbatch\s+(?:fixed\s+)?effects?|"
            r"batch\s+(?:fixed\s+)?effects?\b.{0,80}\bdifferential)"
            r".{0,120}\b(?:preprocessing\s+correction|correct(?:s|ed)?\s+"
            r"(?:the\s+)?(?:input\s+)?(?:data|matrix)|batch-correction method)\b"
        ),
        allowed_contexts=ROADMAP_OR_UNSUPPORTED
        + (
            _rx(r"\bmodel terms?\b"),
            _rx(r"\bmodel covariates?\b"),
        ),
        update_when_supported=(
            "Differential batch covariates must be model terms, not documented "
            "as preprocessing correction."
        ),
    ),
)


def _public_docs_paths() -> tuple[Path, ...]:
    docs = [
        path
        for path in DOCS_ROOT.rglob("*.md")
        if not _is_under(path, DOCS_ROOT / "testing")
    ]
    return (README, *sorted(docs))


def _api_doc_paths() -> tuple[Path, ...]:
    return (
        SRC_ROOT
        / "contracts"
        / "configs"
        / "preprocessing"
        / "correction_missingness.py",
        SRC_ROOT / "contracts" / "configs" / "preprocessing" / "control_sites.py",
        SRC_ROOT / "contracts" / "configs" / "preprocessing" / "batch_correction.py",
        SRC_ROOT
        / "contracts"
        / "configs"
        / "preprocessing"
        / "internal_batch_correction.py",
        SRC_ROOT / "validation" / "configs" / "preprocessing.py",
        SRC_ROOT / "workflows" / "batch_correction" / "interpreter.py",
    )


def _is_under(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def _is_adr_path(path: Path) -> bool:
    return _is_under(path, DOCS_ROOT / "adr")


def _iter_claim_blocks(path: Path) -> tuple[ClaimBlock, ...]:
    blocks: list[ClaimBlock] = []
    current: list[str] = []
    start_line: int | None = None
    in_fence = False

    def flush() -> None:
        nonlocal current, start_line
        if current and start_line is not None:
            blocks.append(
                ClaimBlock(
                    path=path,
                    start_line=start_line,
                    text=_normalise_claim_text(" ".join(current)),
                )
            )
        current = []
        start_line = None

    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
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
        if stripped.startswith("|"):
            flush()
            blocks.append(
                ClaimBlock(
                    path=path,
                    start_line=line_number,
                    text=_normalise_claim_text(stripped),
                )
            )
            continue
        if start_line is None:
            start_line = line_number
        current.append(stripped)

    flush()
    return tuple(blocks)


def _iter_api_doc_blocks(path: Path) -> tuple[ClaimBlock, ...]:
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(path))
    blocks: list[ClaimBlock] = []

    def append_docstring(node: ast.AST) -> None:
        docstring = ast.get_docstring(node)
        if not docstring:
            return
        lineno = getattr(node, "lineno", 1)
        blocks.append(
            ClaimBlock(
                path=path,
                start_line=lineno,
                text=_normalise_claim_text(docstring),
            )
        )

    append_docstring(tree)
    for node in ast.walk(tree):
        if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name.startswith("_"):
                continue
            append_docstring(node)

    return tuple(blocks)


def _normalise_claim_text(text: str) -> str:
    text = text.replace("`", "")
    text = text.replace("*", "")
    return " ".join(text.split()).casefold()


def _is_allowed(block: ClaimBlock, rule: SupportClaimRule) -> bool:
    if any(pattern.search(block.text) for pattern in rule.allowed_contexts):
        return True
    if not _is_adr_path(block.path):
        return False
    return any(pattern.search(block.text) for pattern in ROADMAP_OR_UNSUPPORTED)


def _find_unsupported_claims(
    blocks: tuple[ClaimBlock, ...],
) -> tuple[str, ...]:
    failures: list[str] = []
    for block in blocks:
        for rule in SUPPORT_CLAIM_RULES + DISTINCTION_RULES:
            if not rule.unsupported_claim.search(block.text):
                continue
            if _is_allowed(block, rule):
                continue
            try:
                rel_path = block.path.relative_to(ROOT).as_posix()
            except ValueError:
                rel_path = block.path.as_posix()
            failures.append(
                f"{rel_path}:{block.start_line}: {rule.name}: {block.text}\n"
                f"Update rule: {rule.update_when_supported}"
            )
    return tuple(failures)


def _assert_each_block_rejected(blocks: tuple[ClaimBlock, ...]) -> None:
    failures = _find_unsupported_claims(blocks)
    missing: list[str] = []
    for block in blocks:
        try:
            rel_path = block.path.relative_to(ROOT).as_posix()
        except ValueError:
            rel_path = block.path.as_posix()
        prefix = f"{rel_path}:{block.start_line}:"
        if not any(failure.startswith(prefix) for failure in failures):
            missing.append(f"{prefix} {block.text}")

    assert not missing, "\n".join(missing)


def test_public_docs_do_not_make_unsupported_scientific_parity_claims() -> None:
    blocks = tuple(
        block for path in _public_docs_paths() for block in _iter_claim_blocks(path)
    ) + tuple(
        block for path in _api_doc_paths() for block in _iter_api_doc_blocks(path)
    )

    failures = _find_unsupported_claims(blocks)

    assert not failures, "\n\n".join(failures)


def test_public_docs_describe_ruv_iii_style_only_as_non_executable() -> None:
    allowed = re.compile(
        r"\b(?:not|no)\s+executable\b|\bnot\s+currently\s+supported\b",
        re.IGNORECASE,
    )
    offenders: list[str] = []
    for path in _public_docs_paths():
        for line_number, line in enumerate(
            path.read_text(encoding="utf-8").splitlines(),
            start=1,
        ):
            if "ruv_iii_style" not in line:
                continue
            if allowed.search(line):
                continue
            offenders.append(f"{path.relative_to(ROOT).as_posix()}:{line_number}")

    assert not offenders


def test_public_docs_do_not_promote_control_site_ruv_style_method_label() -> None:
    offenders: list[str] = []
    for path in _public_docs_paths():
        for line_number, line in enumerate(
            path.read_text(encoding="utf-8").splitlines(),
            start=1,
        ):
            if "control_site_ruv_style" in line:
                offenders.append(f"{path.relative_to(ROOT).as_posix()}:{line_number}")

    assert not offenders


def test_claim_guard_rejects_unsupported_positive_claims() -> None:
    bad_blocks = tuple(
        ClaimBlock(Path("synthetic.md"), index, text)
        for index, text in enumerate(
            (
                "PhosPy supports full PhosR parity.",
                "The differential workflow provides broad limma parity.",
                "Protein-aware preparation provides MSstatsPTM equivalence.",
                "EnrichmentWorkflow supports PTM-SEA parity.",
                "PhosPy bundles official Kinase Library data.",
                "Dataset preprocessing supports ComBat batch correction.",
                "Dataset preprocessing implements RUV correction.",
                "Dataset preprocessing implements RUV-III correction.",
                "Differential analysis runs duplicateCorrelation.",
                "Differential analysis supports mixed-effects modelling.",
                "Batch correction supports removeBatchEffect parity.",
            ),
            start=1,
        )
    )

    _assert_each_block_rejected(
        tuple(
            ClaimBlock(
                path=block.path,
                start_line=block.start_line,
                text=_normalise_claim_text(block.text),
            )
            for block in bad_blocks
        )
    )


def test_claim_guard_rejects_batch_correction_shortcuts_and_blurred_boundaries() -> (
    None
):
    bad_blocks = tuple(
        ClaimBlock(Path("synthetic.md"), index, _normalise_claim_text(text))
        for index, text in enumerate(
            (
                "Enable PhosR-equivalent SPS/RUV-III batch correction.",
                "Set use_ruv=True to run batch correction.",
                "Missingness and temporary-imputation contracts for future correction.",
                "Execution plan for SPS/RUV-style correction without numerical correction.",
                "PhosPy provides automatic control-site selection for RUV.",
                "ruv_readiness applies correction to the matrix.",
                "linear_residualize_batch is SPS/RUV-style correction.",
                "linear_residualize_batch is native SPS/RUV-style correction.",
                "Differential batch covariates are preprocessing correction.",
                "Differential batch fixed effects correct the input matrix.",
            ),
            start=1,
        )
    )

    _assert_each_block_rejected(bad_blocks)


def test_claim_guard_allows_current_limitation_language() -> None:
    allowed_blocks = tuple(
        ClaimBlock(Path("synthetic.md"), index, _normalise_claim_text(text))
        for index, text in enumerate(
            (
                "PhosPy does not claim global PhosR parity.",
                "This is not ComBat, not RUV, not limma removeBatchEffect parity.",
                "No executable RUV, ComBat, or limma removeBatchEffect parity lane is supported.",
                "Native PhosPy SPS/RUV-style preprocessing correction is available through SpsRuvBatchCorrectionConfig.",
                "The interpreter returns a resolved numerical execution plan for native SPS/RUV-style preprocessing correction.",
                "Protein-aware preparation is preparation-only and does not claim MSstatsPTM equivalence.",
                "ORA is not GSEA, ssGSEA, or PTM-SEA support.",
                "No official Kinase Library compatibility or parity claim is made.",
                "Fixed-block terms are not limma duplicateCorrelation or mixed-effects modelling.",
                "limma-style moderated variance is supported for the scoped differential lane.",
                "Do not interpret ruv_readiness as RUV support.",
                "ruv_readiness remains diagnostic/report-only metadata readiness reporting.",
                "linear_residualize_batch is limited fixed-effect residualisation, not SPS/RUV-style correction.",
                "Differential batch covariates are model terms, not preprocessing correction.",
            ),
            start=1,
        )
    )

    assert not _find_unsupported_claims(allowed_blocks)


def test_adr_context_allows_unsupported_method_discussion_only_when_qualified() -> None:
    allowed_adr_blocks = tuple(
        ClaimBlock(Path("docs/adr/synthetic.md"), index, _normalise_claim_text(text))
        for index, text in enumerate(
            (
                "A future use_ruv shortcut is unsafe and unsupported.",
                "This ADR discusses RUV-III support as a non-supported roadmap item.",
                "PhosR-equivalent SPS/RUV-III correction is not supported.",
            ),
            start=1,
        )
    )
    bad_adr_blocks = tuple(
        ClaimBlock(Path("docs/adr/synthetic.md"), index, _normalise_claim_text(text))
        for index, text in enumerate(
            (
                "Add use_ruv to simplify examples.",
                "This ADR records RUV support.",
                "PhosR-equivalent SPS/RUV-III correction is available.",
            ),
            start=1,
        )
    )

    assert not _find_unsupported_claims(allowed_adr_blocks)
    _assert_each_block_rejected(bad_adr_blocks)
