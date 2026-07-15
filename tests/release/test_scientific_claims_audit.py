from __future__ import annotations

import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.release_gate

ROOT = Path(__file__).resolve().parents[2]

README = ROOT / "README.md"
SCIENTIFIC_COVERAGE = ROOT / "docs" / "scientific-coverage.md"
WORKFLOW_CONTRACTS = ROOT / "docs" / "workflow_contracts.md"
DIFFERENTIAL_API = ROOT / "docs" / "api" / "differential-analysis.md"
KINASE_API = ROOT / "docs" / "api" / "kinase.md"
ENRICHMENT_API = ROOT / "docs" / "api" / "enrichment.md"
DATASET_BUILD_API = ROOT / "docs" / "api" / "dataset-build-workflow.md"

DIFFERENTIAL_CAVEATS = (
    ROOT / "src" / "phospy" / "workflows" / "differential" / "caveats.py"
)
KINASE_CAVEATS = ROOT / "src" / "phospy" / "workflows" / "kinase" / "caveats.py"
ENRICHMENT_CAVEATS = ROOT / "src" / "phospy" / "workflows" / "enrichment" / "caveats.py"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _normalized(path: Path) -> str:
    text = _read(path).replace("`", "")
    return re.sub(r"\s+", " ", text.casefold())


def _assert_contains_all(path: Path, required: tuple[str, ...]) -> None:
    normalized = _normalized(path)
    missing = [phrase for phrase in required if phrase.casefold() not in normalized]
    assert not missing, f"{path.relative_to(ROOT).as_posix()} missing: {missing}"


def _assert_any(path: Path, alternatives: tuple[str, ...]) -> None:
    normalized = _normalized(path)
    assert any(phrase.casefold() in normalized for phrase in alternatives), (
        f"{path.relative_to(ROOT).as_posix()} missing one of: {alternatives}"
    )


def test_release_docs_keep_differential_scope_limited_to_tested_envelopes() -> None:
    for path in (README, SCIENTIFIC_COVERAGE, WORKFLOW_CONTRACTS, DIFFERENTIAL_API):
        _assert_contains_all(path, ("tested design and contrast envelopes",))
        _assert_any(
            path,
            (
                "not full limma or phosr parity",
                "not full phosr or limma parity",
            ),
        )


def test_release_docs_keep_kinase_scores_relative_not_causal_proof() -> None:
    _assert_contains_all(
        README,
        (
            "relative support within a run",
            "not direct proof of kinase activation or causal pathway activity",
        ),
    )
    for path in (SCIENTIFIC_COVERAGE, WORKFLOW_CONTRACTS, KINASE_API):
        _assert_contains_all(path, ("relative support",))
        _assert_any(
            path,
            (
                "not calibrated causal inference",
                "not calibrated probabilities or proof of causal regulation",
                "not causal proof",
            ),
        )


def test_release_docs_keep_enrichment_ora_only() -> None:
    for path in (README, SCIENTIFIC_COVERAGE, WORKFLOW_CONTRACTS, ENRICHMENT_API):
        _assert_contains_all(path, ("ora",))
        _assert_any(
            path,
            (
                "offline over-representation analysis",
                "offline ora",
            ),
        )
        _assert_any(
            path,
            (
                "does not imply gsea or ptm-sea support",
                "ora does not imply gsea, ssgsea, or ptm-sea support",
                "ora is not gsea, ssgsea, or ptm-sea support",
                "does not implement ranked-list enrichment, gsea, ssgsea, or ptm-sea",
            ),
        )


def test_release_docs_keep_bundled_references_rat_only() -> None:
    for path in (README, SCIENTIFIC_COVERAGE, WORKFLOW_CONTRACTS, KINASE_API):
        _assert_contains_all(path, ("rat-only", "referencebundle"))

    combined = " ".join(
        _normalized(path) for path in (README, SCIENTIFIC_COVERAGE, KINASE_API)
    )
    assert "rat-first" not in combined
    assert "bundled human" not in combined
    assert "bundled mouse" not in combined


def test_release_docs_keep_batch_correction_limits_explicit() -> None:
    for path in (README, SCIENTIFIC_COVERAGE, WORKFLOW_CONTRACTS, DATASET_BUILD_API):
        _assert_contains_all(
            path,
            (
                "not combat",
                "not ruv",
                "not limma removebatcheffect parity",
                "not mixed-effects",
            ),
        )

    _assert_contains_all(
        SCIENTIFIC_COVERAGE,
        (
            "not phosr-equivalent sps/ruv-iii parity",
            "replicate-aware ruv-iii semantics",
            "requires a complete correction-stage matrix",
            "rejects actual missing values (nans) before executor invocation",
        ),
    )
    _assert_contains_all(
        WORKFLOW_CONTRACTS,
        (
            "replicate-aware ruv-iii correction semantics",
            "not used for numerical unwanted-factor estimation",
            "rejects actual missing values (nans)",
        ),
    )


def test_release_caveat_builders_expose_scientific_interpretation_limits() -> None:
    _assert_contains_all(
        DIFFERENTIAL_CAVEATS,
        (
            "not full",
            "limma or phosr parity",
            "tested_design_and_contrast_envelope",
            "full_limma_or_phosr_parity_claimed",
        ),
    )
    _assert_contains_all(
        KINASE_CAVEATS,
        (
            "relative support values",
            "not causal kinase activity proof",
            "relative_support_within_run",
            "not_causal_activity_proof",
        ),
    )
    _assert_contains_all(
        ENRICHMENT_CAVEATS,
        (
            "offline over-representation analysis",
            "ora overlap counts and p-values only",
            "gsea, ssgsea, ptm-sea",
            "rank_based_enrichment_supported",
        ),
    )
