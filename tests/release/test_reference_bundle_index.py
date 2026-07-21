from __future__ import annotations

import json
import subprocess
from hashlib import sha256
from pathlib import Path

import pytest

from scripts.validate_reference_bundle_index import (
    ReferenceBundleIndexError,
    validate_reference_bundle_index,
)

pytestmark = pytest.mark.release_gate

ROOT = Path(__file__).resolve().parents[2]
RAT_MANIFEST_PATH = "src/phospy/data/reference_bundles/rat/l6_native/manifest.json"
EXPECTED_RAT_DIGESTS = {
    "substrate_map.csv": (
        "dc5ad357a6aaf2ef29d1fa931f117b06522ff2abe51146f0cfbc84aa4f36c32d"
    ),
    "site_sequences.csv": (
        "11afdb3aefd82d58f6d3c04af481afeb0ea0c3a3627679cfaa5e648371d8af3f"
    ),
    "motif_scores.csv": (
        "66905f60e65294b20dc7300faaf461c8d58943def63b488515096e79f65e790c"
    ),
    "motif_sizes.csv": (
        "9ab23117b6d489f2ab8a131fa590a8d26cce549162b479238b992274bde8fda0"
    ),
    "ATTRIBUTION.md": (
        "341f6dcf517bda8b90f69ec09fa80c49c50db9c2728f3569a7dff57bb2a9c01e"
    ),
}


def test_real_rat_staged_reference_bundle_files_match_manifest() -> None:
    validated_files = validate_reference_bundle_index(repo_root=ROOT)

    rat_digests = {
        item.relative_path: item.actual_sha256
        for item in validated_files
        if item.manifest_path == RAT_MANIFEST_PATH
    }

    assert rat_digests == EXPECTED_RAT_DIGESTS


def test_staged_index_validation_rejects_crlf_payload_mismatch(
    tmp_path: Path,
) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    _run_git(repo, "init")
    _run_git(repo, "config", "core.autocrlf", "false")
    bundle_root = repo / "src" / "phospy" / "data" / "reference_bundles" / "rat"
    bundle_root = bundle_root / "l6_native"
    bundle_root.mkdir(parents=True)

    intended_bytes = b"Unit attribution\n"
    converted_bytes = b"Unit attribution\r\n"
    attribution_path = bundle_root / "ATTRIBUTION.md"
    attribution_path.write_bytes(converted_bytes)
    manifest_path = bundle_root / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "reference_id": "l6_native",
                "files": [
                    {
                        "relative_path": "ATTRIBUTION.md",
                        "sha256": sha256(intended_bytes).hexdigest(),
                    }
                ],
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    _run_git(
        repo,
        "add",
        "src/phospy/data/reference_bundles/rat/l6_native/manifest.json",
        "src/phospy/data/reference_bundles/rat/l6_native/ATTRIBUTION.md",
    )

    with pytest.raises(ReferenceBundleIndexError) as exc_info:
        validate_reference_bundle_index(repo_root=repo)

    message = str(exc_info.value)
    assert "reference bundle staged blob digest mismatch" in message
    assert "affected file=ATTRIBUTION.md" in message
    assert f"expected digest={sha256(intended_bytes).hexdigest()}" in message
    assert f"actual digest={sha256(converted_bytes).hexdigest()}" in message


def test_staged_index_validation_rejects_extra_reference_file(
    tmp_path: Path,
) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    _run_git(repo, "init")
    bundle_root = repo / "src" / "phospy" / "data" / "reference_bundles" / "rat"
    bundle_root = bundle_root / "l6_native"
    bundle_root.mkdir(parents=True)
    attribution_path = bundle_root / "ATTRIBUTION.md"
    attribution_path.write_text("Unit attribution\n", encoding="utf-8")
    extra_path = bundle_root / "extra.csv"
    extra_path.write_text("extra\n", encoding="utf-8")
    manifest_path = bundle_root / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "reference_id": "l6_native",
                "files": [
                    {
                        "relative_path": "ATTRIBUTION.md",
                        "sha256": sha256(attribution_path.read_bytes()).hexdigest(),
                    }
                ],
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    _run_git(repo, "add", "src/phospy/data/reference_bundles")

    with pytest.raises(ReferenceBundleIndexError) as exc_info:
        validate_reference_bundle_index(repo_root=repo)

    message = str(exc_info.value)
    assert "Git index contains undeclared reference-bundle file" in message
    assert (
        "affected file=src/phospy/data/reference_bundles/rat/l6_native/extra.csv"
        in (message)
    )


def _run_git(repo: Path, *args: str) -> None:
    subprocess.run(
        ("git", *args),
        cwd=repo,
        capture_output=True,
        check=True,
    )
