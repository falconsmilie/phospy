from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

from tests.support.fixture_byte_policy import (
    assert_canonical_text_bytes,
    assert_lf_gitattributes_coverage,
    assert_text_fixture_matches_sha256,
    iter_importer_fixture_index_references,
    iter_manifest_governed_text_fixture_references,
)

ROOT = Path(__file__).resolve().parents[2]
FIXTURE_ROOT = ROOT / "tests" / "fixtures"
RELEASE_VALIDATION_ROOT = FIXTURE_ROOT / "release_validation_regression"
LARGE_LIMMA_TREND_ROOT = (
    FIXTURE_ROOT / "rewrite_parity" / "differential_limma_trend_large"
)
CANONICAL_BYTE_POLICY = "utf-8 LF with final newline"

MANIFEST_GOVERNED_FIXTURE_DIRS = (
    RELEASE_VALIDATION_ROOT / "evidence_resolution",
    RELEASE_VALIDATION_ROOT / "kinase_sparse_support",
    RELEASE_VALIDATION_ROOT / "signalome_safety",
    RELEASE_VALIDATION_ROOT / "sps_ruv_planted_unwanted_factor",
    RELEASE_VALIDATION_ROOT / "peptide_site_bias_regimes",
    RELEASE_VALIDATION_ROOT / "kinase_activity_known_membership",
    RELEASE_VALIDATION_ROOT / "signalome_planted_modules",
    RELEASE_VALIDATION_ROOT / "importer_edge_cases",
    LARGE_LIMMA_TREND_ROOT,
)

pytestmark = [pytest.mark.release_gate, pytest.mark.reproducibility]


def _read_manifest(fixture_dir: Path) -> dict[str, Any]:
    return json.loads((fixture_dir / "MANIFEST.json").read_text(encoding="utf-8"))


def _tree_files(root: Path) -> tuple[Path, ...]:
    return tuple(
        sorted(path.relative_to(root) for path in root.rglob("*") if path.is_file())
    )


def _first_differing_byte(left: bytes, right: bytes) -> int | None:
    for index, (left_byte, right_byte) in enumerate(zip(left, right, strict=False)):
        if left_byte != right_byte:
            return index
    if len(left) != len(right):
        return min(len(left), len(right))
    return None


def _byte_at(data: bytes, index: int) -> str:
    if index >= len(data):
        return "<EOF>"
    return f"0x{data[index]:02x}"


def _display_path(path: Path) -> Path:
    try:
        return path.relative_to(ROOT)
    except ValueError:
        return path


def _validate_manifest_hashes(fixture_dir: Path) -> None:
    manifest_path = fixture_dir / "MANIFEST.json"
    manifest_bytes = manifest_path.read_bytes()
    manifest = _read_manifest(fixture_dir)
    byte_policy = str(manifest.get("byte_policy"))

    assert byte_policy == CANONICAL_BYTE_POLICY
    assert_canonical_text_bytes(
        _display_path(manifest_path),
        manifest_bytes,
        byte_policy=byte_policy,
    )

    declared_paths = {
        Path(str(file_entry["relative_path"])) for file_entry in manifest["files"]
    }
    expected_paths = tuple(sorted(declared_paths | {Path("MANIFEST.json")}))
    assert _tree_files(fixture_dir) == expected_paths

    for file_entry in manifest["files"]:
        relative_path = Path(str(file_entry["relative_path"]))
        path = fixture_dir / relative_path
        payload = path.read_bytes()
        assert_canonical_text_bytes(
            _display_path(path),
            payload,
            byte_policy=byte_policy,
        )
        actual = hashlib.sha256(payload).hexdigest()
        expected = str(file_entry["sha256"])
        assert actual == expected, (
            "manifest fixture digest mismatch: "
            f"file={path.relative_to(ROOT).as_posix()} "
            f"byte_policy={byte_policy!r} expected={expected} actual={actual}"
        )


def _assert_generated_tree_matches_checked_in(
    *,
    generated_root: Path,
    checked_in_root: Path,
    byte_policy: str,
) -> None:
    generated_files = _tree_files(generated_root)
    checked_in_files = _tree_files(checked_in_root)
    assert generated_files == checked_in_files

    for relative_path in generated_files:
        generated_path = generated_root / relative_path
        checked_in_path = checked_in_root / relative_path
        generated = generated_path.read_bytes()
        checked_in = checked_in_path.read_bytes()
        assert_canonical_text_bytes(
            generated_path,
            generated,
            byte_policy=byte_policy,
        )
        if generated != checked_in:
            index = _first_differing_byte(generated, checked_in)
            assert index is not None
            pytest.fail(
                "generated manifest-governed fixture bytes differ from "
                "checked-in bytes: "
                f"file={checked_in_path.relative_to(ROOT).as_posix()} "
                f"byte_policy={byte_policy!r} first_differing_byte={index} "
                f"generated={_byte_at(generated, index)} "
                f"checked_in={_byte_at(checked_in, index)}"
            )


def _run_checked(command: list[str]) -> None:
    result = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, (
        f"command failed: {subprocess.list2cmdline(command)}\n"
        f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"
    )


def test_manifest_governed_fixtures_use_canonical_lf_bytes_and_valid_hashes() -> None:
    for fixture_dir in MANIFEST_GOVERNED_FIXTURE_DIRS:
        _validate_manifest_hashes(fixture_dir)
    for reference in iter_importer_fixture_index_references(ROOT):
        assert reference.path.is_file()
        assert_text_fixture_matches_sha256(
            reference.path,
            expected_sha256=reference.expected_sha256,
            repo_root=ROOT,
        )


def test_manifest_governed_text_fixtures_have_lf_gitattributes_coverage() -> None:
    references = iter_manifest_governed_text_fixture_references(ROOT)
    assert references, "expected manifest-governed text fixture references"
    assert_lf_gitattributes_coverage(
        ROOT,
        tuple(reference.path for reference in references),
    )


def test_release_validation_generator_reproduces_checked_in_bytes(
    tmp_path: Path,
) -> None:
    generated_root = tmp_path / "release_validation_regression"
    _run_checked(
        [
            sys.executable,
            "scripts/active/generate_release_validation_regression_fixtures.py",
            "--outdir",
            str(generated_root),
            "--manifest-outdir-label",
            "tests/fixtures/release_validation_regression",
            "--timestamp",
            "2026-07-24T00:00:00Z",
            "--seed",
            "20260724",
        ]
    )

    _assert_generated_tree_matches_checked_in(
        generated_root=generated_root,
        checked_in_root=RELEASE_VALIDATION_ROOT,
        byte_policy=CANONICAL_BYTE_POLICY,
    )
    for fixture_dir in sorted(
        path for path in generated_root.iterdir() if path.is_dir()
    ):
        _validate_manifest_hashes(fixture_dir)


def _matching_rscript_or_skip() -> str:
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is not available for large limma fixture regeneration")

    expected_manifest = _read_manifest(LARGE_LIMMA_TREND_ROOT)
    expected_r_version = expected_manifest["external_implementation"]["r_version"]
    expected_limma_version = expected_manifest["external_implementation"][
        "limma_version"
    ]
    result = subprocess.run(
        [
            rscript,
            "-e",
            (
                "if (!requireNamespace('limma', quietly = TRUE)) quit(status = 42); "
                "cat(R.version$version.string, '\\n', "
                "as.character(packageVersion('limma')), '\\n', sep = '')"
            ),
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode == 42:
        pytest.skip("R package 'limma' is not available for fixture regeneration")
    assert result.returncode == 0, result.stderr
    observed = result.stdout.splitlines()
    assert len(observed) >= 2, result.stdout
    if observed[0] != expected_r_version or observed[1] != expected_limma_version:
        pytest.skip(
            "large limma fixture exact-byte regeneration requires matching "
            f"R/limma versions: expected {expected_r_version!r}/"
            f"{expected_limma_version!r}, observed {observed[0]!r}/{observed[1]!r}"
        )
    return rscript


def test_large_limma_trend_generator_reproduces_checked_in_bytes(
    tmp_path: Path,
) -> None:
    rscript = _matching_rscript_or_skip()
    generated_root = tmp_path / "differential_limma_trend_large"
    _run_checked(
        [
            rscript,
            "scripts/active/generate_large_differential_limma_trend_fixture.R",
            "--outdir",
            str(generated_root),
            "--manifest-outdir-label",
            "tests/fixtures/rewrite_parity/differential_limma_trend_large",
            "--seed",
            "20260724",
            "--timestamp",
            "2026-07-24T00:00:00Z",
            "--n_features",
            "1600",
        ]
    )

    _assert_generated_tree_matches_checked_in(
        generated_root=generated_root,
        checked_in_root=LARGE_LIMMA_TREND_ROOT,
        byte_policy=CANONICAL_BYTE_POLICY,
    )
    _validate_manifest_hashes(generated_root)
