from __future__ import annotations

import io
import json
import subprocess
import tarfile
from hashlib import sha256
from pathlib import Path, PurePosixPath
from zipfile import ZIP_DEFLATED, ZipFile

import pytest

from phospy.science.references.validation import validate_bundled_reference_manifests
from scripts.validate_reference_bundle_distribution import (
    ReferenceBundleDistributionError,
    validate_reference_bundle_archive,
    validate_reference_bundle_wheel,
)

pytestmark = pytest.mark.release_gate

_BUNDLE_ROOT = PurePosixPath("phospy/data/reference_bundles/rat/l6_native")
_SOURCE_BUNDLE_ROOT = (
    Path(__file__).resolve().parents[2]
    / "src"
    / "phospy"
    / "data"
    / "reference_bundles"
    / "rat"
    / "l6_native"
)
_SOURCE_REFERENCE_BUNDLES_ROOT = _SOURCE_BUNDLE_ROOT.parents[1]
_BASE_FILES = {
    "substrate_map.csv": b"kinase,site_id\nAKT1,MAPK1;S123;\n",
    "site_sequences.csv": b"site_id,centralized_sequence\nMAPK1;S123;,ACD\n",
    "ATTRIBUTION.md": b"Unit attribution\n",
}


def test_real_rat_source_bundle_validates_exact_file_hashes() -> None:
    manifests = validate_bundled_reference_manifests(_SOURCE_REFERENCE_BUNDLES_ROOT)

    rat_manifest = next(
        manifest for manifest in manifests if manifest.reference_id == "l6_native"
    )
    attribution = next(
        item for item in rat_manifest.files if item.relative_path == "ATTRIBUTION.md"
    )

    assert attribution.sha256 == _sha256(
        (_SOURCE_BUNDLE_ROOT / "ATTRIBUTION.md").read_bytes()
    )


def test_valid_temporary_wheel_archive_validates(tmp_path: Path) -> None:
    wheel_path = _write_reference_wheel(tmp_path=tmp_path)

    validate_reference_bundle_wheel(wheel_path)


def test_valid_temporary_sdist_archive_validates(tmp_path: Path) -> None:
    sdist_path = _write_reference_sdist(tmp_path=tmp_path)

    validate_reference_bundle_archive(sdist_path)


def test_distribution_validation_compares_archives_to_git_index(
    tmp_path: Path,
) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    _run_git(repo, "init")
    _write_git_source_bundle(repo)
    _run_git(repo, "add", "src/phospy/data/reference_bundles")
    changed_bytes = b"Changed attribution\n"
    wheel_path = _write_reference_wheel(
        tmp_path=tmp_path,
        packaged_file_overrides={"ATTRIBUTION.md": changed_bytes},
    )

    with pytest.raises(ReferenceBundleDistributionError) as exc_info:
        validate_reference_bundle_archive(
            wheel_path,
            compare_git_index=True,
            repo_root=repo,
        )

    message = str(exc_info.value)
    assert (
        "distribution reference-bundle file does not reproduce committed Git index blob"
    ) in message
    assert (
        "affected file=src/phospy/data/reference_bundles/rat/l6_native/ATTRIBUTION.md"
    ) in message
    assert f"expected digest={_sha256(_BASE_FILES['ATTRIBUTION.md'])}" in message
    assert f"actual digest={_sha256(changed_bytes)}" in message


def test_distribution_validation_fails_for_changed_attribution_file(
    tmp_path: Path,
) -> None:
    changed_bytes = b"Changed attribution\n"
    wheel_path = _write_reference_wheel(
        tmp_path=tmp_path,
        packaged_file_overrides={"ATTRIBUTION.md": changed_bytes},
    )

    message = _distribution_error(wheel_path)

    assert "reference bundle file digest mismatch" in message
    assert "affected file=ATTRIBUTION.md" in message
    assert f"expected digest={_sha256(_BASE_FILES['ATTRIBUTION.md'])}" in message
    assert f"actual digest={_sha256(changed_bytes)}" in message


def test_distribution_validation_fails_for_missing_attribution_file(
    tmp_path: Path,
) -> None:
    wheel_path = _write_reference_wheel(
        tmp_path=tmp_path,
        omitted_files={"ATTRIBUTION.md"},
    )

    message = _distribution_error(wheel_path)

    assert "required bundle-local attribution file is missing from archive" in message
    assert "affected file=ATTRIBUTION.md" in message
    assert f"expected digest={_sha256(_BASE_FILES['ATTRIBUTION.md'])}" in message
    assert "actual digest=missing" in message


def test_distribution_validation_fails_for_missing_csv(tmp_path: Path) -> None:
    wheel_path = _write_reference_wheel(
        tmp_path=tmp_path,
        omitted_files={"substrate_map.csv"},
    )

    message = _distribution_error(wheel_path)

    assert "manifest-listed file is missing from archive" in message
    assert "affected file=substrate_map.csv" in message
    assert f"expected digest={_sha256(_BASE_FILES['substrate_map.csv'])}" in message
    assert "actual digest=missing" in message


def test_distribution_validation_fails_for_digest_mismatch(tmp_path: Path) -> None:
    changed_bytes = b"site_id,centralized_sequence\nMAPK1;S123;,XYZ\n"
    wheel_path = _write_reference_wheel(
        tmp_path=tmp_path,
        packaged_file_overrides={"site_sequences.csv": changed_bytes},
    )

    message = _distribution_error(wheel_path)

    assert "reference bundle file digest mismatch" in message
    assert "affected file=site_sequences.csv" in message
    assert f"expected digest={_sha256(_BASE_FILES['site_sequences.csv'])}" in message
    assert f"actual digest={_sha256(changed_bytes)}" in message


def test_distribution_validation_error_message_identifies_wheel_bundle_and_reference(
    tmp_path: Path,
) -> None:
    changed_bytes = b"kinase,site_id\nAKT1,MAPK1;S999;\n"
    wheel_path = _write_reference_wheel(
        tmp_path=tmp_path,
        packaged_file_overrides={"substrate_map.csv": changed_bytes},
    )

    message = _distribution_error(wheel_path)

    assert f"archive path={wheel_path}" in message
    assert "bundle path=src/phospy/data/reference_bundles/rat/l6_native" in message
    assert "reference ID=l6_native" in message
    assert "affected file=substrate_map.csv" in message
    assert f"expected digest={_sha256(_BASE_FILES['substrate_map.csv'])}" in message
    assert f"actual digest={_sha256(changed_bytes)}" in message


@pytest.mark.parametrize(
    (
        "packaged_file_overrides",
        "omitted_files",
        "expected_file",
        "expected_digest",
        "actual_digest",
        "expected_reason",
    ),
    [
        pytest.param(
            None,
            {"site_sequences.csv"},
            "site_sequences.csv",
            sha256(_BASE_FILES["site_sequences.csv"]).hexdigest(),
            "missing",
            "manifest-listed file is missing from archive",
            id="wheel-missing-file",
        ),
        pytest.param(
            {"site_sequences.csv": b"site_id,centralized_sequence\nMAPK1;S123;,XYZ\n"},
            None,
            "site_sequences.csv",
            sha256(_BASE_FILES["site_sequences.csv"]).hexdigest(),
            sha256(b"site_id,centralized_sequence\nMAPK1;S123;,XYZ\n").hexdigest(),
            "reference bundle file digest mismatch",
            id="wheel-hash-mismatch",
        ),
    ],
)
def test_wheel_file_errors_keep_actionable_digest_contract(
    tmp_path: Path,
    packaged_file_overrides: dict[str, bytes] | None,
    omitted_files: set[str] | None,
    expected_file: str,
    expected_digest: str,
    actual_digest: str,
    expected_reason: str,
) -> None:
    wheel_path = _write_reference_wheel(
        tmp_path=tmp_path,
        packaged_file_overrides=packaged_file_overrides,
        omitted_files=omitted_files,
    )

    message = _distribution_error(wheel_path)

    assert expected_reason in message
    assert f"archive path={wheel_path}" in message
    assert "bundle path=src/phospy/data/reference_bundles/rat/l6_native" in message
    assert "reference ID=l6_native" in message
    assert f"affected file={expected_file}" in message
    assert f"expected digest={expected_digest}" in message
    assert f"actual digest={actual_digest}" in message


def _distribution_error(wheel_path: Path) -> str:
    with pytest.raises(ReferenceBundleDistributionError) as exc_info:
        validate_reference_bundle_wheel(wheel_path)
    return str(exc_info.value)


def _write_reference_wheel(
    *,
    tmp_path: Path,
    packaged_file_overrides: dict[str, bytes] | None = None,
    omitted_files: set[str] | None = None,
) -> Path:
    packaged_files = dict(_BASE_FILES)
    if packaged_file_overrides is not None:
        packaged_files.update(packaged_file_overrides)
    omitted = omitted_files or set()
    wheel_path = tmp_path / "phospy-1.6.0-py3-none-any.whl"
    with ZipFile(wheel_path, "w", compression=ZIP_DEFLATED) as archive:
        archive.writestr(
            (_BUNDLE_ROOT / "manifest.json").as_posix(),
            json.dumps(_manifest_payload(), indent=2) + "\n",
        )
        for relative_path, data in packaged_files.items():
            if relative_path in omitted:
                continue
            archive.writestr((_BUNDLE_ROOT / relative_path).as_posix(), data)
        archive.writestr("phospy/__init__.py", b"")
        archive.writestr(
            "phospy-1.6.0.dist-info/WHEEL",
            "Wheel-Version: 1.0\nGenerator: unit-test\nRoot-Is-Purelib: true\n",
        )
    return wheel_path


def _write_reference_sdist(*, tmp_path: Path) -> Path:
    sdist_path = tmp_path / "phospy-1.6.0.tar.gz"
    with tarfile.open(sdist_path, "w:gz") as archive:
        root = PurePosixPath("phospy-1.6.0") / "src" / _BUNDLE_ROOT
        _add_tar_bytes(
            archive,
            (root / "manifest.json").as_posix(),
            (json.dumps(_manifest_payload(), indent=2) + "\n").encode("utf-8"),
        )
        for relative_path, data in _BASE_FILES.items():
            _add_tar_bytes(archive, (root / relative_path).as_posix(), data)
    return sdist_path


def _add_tar_bytes(
    archive: tarfile.TarFile,
    name: str,
    data: bytes,
) -> None:
    info = tarfile.TarInfo(name)
    info.size = len(data)
    archive.addfile(info, io.BytesIO(data))


def _write_git_source_bundle(repo: Path) -> None:
    bundle_root = repo / "src" / "phospy" / "data" / "reference_bundles"
    bundle_root = bundle_root / "rat" / "l6_native"
    bundle_root.mkdir(parents=True)
    (bundle_root / "manifest.json").write_text(
        json.dumps(_manifest_payload(), indent=2) + "\n",
        encoding="utf-8",
    )
    for relative_path, data in _BASE_FILES.items():
        (bundle_root / relative_path).write_bytes(data)


def _run_git(repo: Path, *args: str) -> None:
    subprocess.run(
        ("git", *args),
        cwd=repo,
        capture_output=True,
        check=True,
    )


def _manifest_payload() -> dict[str, object]:
    return {
        "reference_id": "l6_native",
        "reference_version": "bundled-snapshot-2026-04-16",
        "redistribution_evidence": {
            "attribution": {
                "bundle_attribution_path": "ATTRIBUTION.md",
            },
        },
        "files": [
            _file_payload(
                "substrate_map.csv",
                role="kinase_substrate",
                file_format="csv",
            ),
            _file_payload(
                "site_sequences.csv",
                role="site_sequences",
                file_format="csv",
            ),
            _file_payload(
                "ATTRIBUTION.md",
                role="attribution",
                file_format="markdown",
            ),
        ],
    }


def _file_payload(
    relative_path: str,
    *,
    role: str,
    file_format: str,
) -> dict[str, object]:
    return {
        "relative_path": relative_path,
        "role": role,
        "format": file_format,
        "sha256": _sha256(_BASE_FILES[relative_path]),
        "row_count": None,
        "column_names": None,
    }


def _sha256(data: bytes) -> str:
    return sha256(data).hexdigest()
