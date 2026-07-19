from __future__ import annotations

import importlib.util
import json
import sys
from hashlib import sha256
from pathlib import Path
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parents[2]
VERIFIER_PATH = ROOT / "scripts" / "verify_distribution_artifact.py"

pytestmark = pytest.mark.release_gate


class _FakeDistribution:
    files = (Path("phospy") / "__init__.py",)

    def __init__(self, distribution_root: Path) -> None:
        self._distribution_root = distribution_root

    def locate_file(self, file_path: Path) -> Path:
        return self._distribution_root / file_path


def _load_verifier() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "verify_distribution_artifact_under_test",
        VERIFIER_PATH,
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _clear_phospy_modules(monkeypatch: pytest.MonkeyPatch) -> None:
    for module_name in list(sys.modules):
        if module_name == "phospy" or module_name.startswith("phospy."):
            monkeypatch.delitem(sys.modules, module_name, raising=False)


def _fake_phospy_module(package_init: Path, *, version: str) -> ModuleType:
    package = ModuleType("phospy")
    package.__file__ = str(package_init)
    package.__version__ = version
    package.__path__ = [str(package_init.parent)]  # type: ignore[attr-defined]
    return package


def _identity_context(verifier: ModuleType, tmp_path: Path):
    return verifier.VerificationContext(
        artifact_kind="wheel",
        artifact_path=tmp_path / "phospy-1.6.0-py3-none-any.whl",
        build_manifest_path=tmp_path / "build-manifest.json",
        repository_root=tmp_path / "repo",
        report_json=tmp_path / "report.json",
    )


def _write_artifact_attestation_inputs(
    tmp_path: Path,
    *,
    artifact_kind: str = "wheel",
    artifact_name: str = "phospy-1.6.0-py3-none-any.whl",
) -> tuple[Path, Path, str]:
    artifact_path = tmp_path / artifact_name
    artifact_path.write_bytes(b"distribution artifact bytes")
    artifact_sha256 = sha256(artifact_path.read_bytes()).hexdigest()
    manifest_path = tmp_path / "build-manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "schema": "phospy.build-manifest/v1",
                "source_identity_digest": "sha256:" + ("1" * 64),
                "package_version": "1.6.0",
                "artifacts": [
                    {
                        "kind": artifact_kind,
                        "filename": artifact_name,
                        "sha256": artifact_sha256,
                    }
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return artifact_path, manifest_path, artifact_sha256


def test_verifier_rejects_package_imported_from_checkout(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    verifier = _load_verifier()
    context = _identity_context(verifier, tmp_path)
    package_init = context.repository_root / "src" / "phospy" / "__init__.py"
    package_init.parent.mkdir(parents=True)
    package_init.write_text("", encoding="utf-8")

    _clear_phospy_modules(monkeypatch)
    monkeypatch.setitem(
        sys.modules,
        "phospy",
        _fake_phospy_module(package_init, version="1.6.0"),
    )
    monkeypatch.setattr(verifier, "_require_isolated_interpreter", lambda: None)
    monkeypatch.setattr(
        verifier,
        "_require_external_working_directory",
        lambda repository_root: None,
    )
    monkeypatch.setattr(
        verifier.metadata,
        "distribution",
        lambda name: _FakeDistribution(package_init.parents[1]),
    )
    monkeypatch.setattr(verifier.metadata, "version", lambda name: "1.6.0")

    with pytest.raises(verifier.VerificationError, match="source checkout"):
        verifier._check_installed_package_identity(context)


def test_verifier_rejects_metadata_runtime_version_disagreement(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    verifier = _load_verifier()
    context = _identity_context(verifier, tmp_path)
    package_init = tmp_path / "site-packages" / "phospy" / "__init__.py"
    package_init.parent.mkdir(parents=True)
    package_init.write_text("", encoding="utf-8")

    _clear_phospy_modules(monkeypatch)
    monkeypatch.setitem(
        sys.modules,
        "phospy",
        _fake_phospy_module(package_init, version="0.0.0"),
    )
    monkeypatch.setattr(verifier, "_require_isolated_interpreter", lambda: None)
    monkeypatch.setattr(
        verifier,
        "_require_external_working_directory",
        lambda repository_root: None,
    )
    monkeypatch.setattr(
        verifier.metadata,
        "distribution",
        lambda name: _FakeDistribution(package_init.parents[1]),
    )
    monkeypatch.setattr(verifier.metadata, "version", lambda name: "1.6.0")

    with pytest.raises(verifier.VerificationError, match="runtime version"):
        verifier._check_installed_package_identity(context)


def test_verifier_rejects_digest_invalid_scientific_resource(
    tmp_path: Path,
) -> None:
    verifier = _load_verifier()
    bundle_root = tmp_path / "bundle"
    bundle_root.mkdir()
    resource = bundle_root / "resource.csv"
    resource.write_text("changed\n", encoding="utf-8")
    manifest_payload = {
        "files": [
            {
                "relative_path": "resource.csv",
                "sha256": "0" * 64,
            }
        ]
    }

    with pytest.raises(verifier.VerificationError, match="digest mismatch"):
        verifier._validate_manifest_resource_digests(bundle_root, manifest_payload)


def test_verifier_report_exposes_stable_artifact_attestation_contract(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    verifier = _load_verifier()
    artifact_path, manifest_path, artifact_sha256 = _write_artifact_attestation_inputs(
        tmp_path
    )
    report_path = tmp_path / "artifact-report.json"

    def _identity_check(context: object) -> dict[str, object]:
        context.package_root = tmp_path / "site-packages" / "phospy"
        context.distribution_version = "1.6.0"
        return {
            "metadata_version": "1.6.0",
            "runtime_version": "1.6.0",
        }

    def _pass_check(context: object) -> dict[str, object]:
        return {"ok": True}

    monkeypatch.setattr(
        verifier,
        "_assert_loaded_phospy_modules_not_from_repository",
        lambda context: None,
    )
    monkeypatch.setattr(
        verifier,
        "_CHECKS",
        (
            ("artifact-manifest-binding", verifier._check_artifact_manifest_binding),
            ("installed-package-identity", _identity_check),
            ("packaged-scientific-resources", _pass_check),
            ("corrected-construction-and-provenance-path", _pass_check),
            ("corrected-derived-and-ownership-path", _pass_check),
            ("corrected-differential-path", _pass_check),
        ),
    )

    exit_status = verifier.main(
        [
            "--artifact-kind",
            "wheel",
            "--artifact-path",
            str(artifact_path),
            "--build-manifest",
            str(manifest_path),
            "--repository-root",
            str(tmp_path / "repo"),
            "--report-json",
            str(report_path),
        ]
    )

    payload = json.loads(report_path.read_text(encoding="utf-8"))
    assert exit_status == 0
    assert payload["schema"] == "phospy.artifact-verification/v1"
    assert payload["status"] == "success"
    assert payload["source_identity_digest"] == "sha256:" + ("1" * 64)
    assert payload["artifact"] == {
        "kind": "wheel",
        "filename": artifact_path.name,
        "sha256": artifact_sha256,
    }
    assert payload["package"] == {"name": "phospy", "version": "1.6.0"}
    assert payload["environment"]["python"] == sys.version.split()[0]
    assert payload["environment"]["dependency_snapshot_sha256"].startswith("sha256:")
    assert payload["checks"] == {
        "artifact_manifest_binding": "pass",
        "installed_import_origin": "pass",
        "package_metadata": "pass",
        "scientific_resources": "pass",
        "trusted_construction": "pass",
        "provenance_immutability": "pass",
        "derived_lineage": "pass",
        "dataframe_ownership": "pass",
        "differential_execution": "pass",
    }
    assert payload["check_details"][0]["name"] == "artifact-manifest-binding"


def test_verifier_returns_nonzero_when_corrected_path_invariant_fails(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    verifier = _load_verifier()
    artifact_path, manifest_path, _artifact_sha256 = _write_artifact_attestation_inputs(
        tmp_path
    )
    report_path = tmp_path / "artifact-report.json"

    def _pass_check(context: object) -> dict[str, object]:
        return {"ok": True}

    def _fail_check(context: object) -> dict[str, object]:
        raise verifier.VerificationError("trusted assertion evidence missing")

    monkeypatch.setattr(
        verifier,
        "_assert_loaded_phospy_modules_not_from_repository",
        lambda context: None,
    )
    monkeypatch.setattr(
        verifier,
        "_CHECKS",
        (
            ("installed-package-identity", _pass_check),
            ("corrected-construction-and-provenance-path", _fail_check),
        ),
    )

    exit_status = verifier.main(
        [
            "--artifact-kind",
            "wheel",
            "--artifact-path",
            str(artifact_path),
            "--build-manifest",
            str(manifest_path),
            "--repository-root",
            str(tmp_path / "repo"),
            "--report-json",
            str(report_path),
        ]
    )

    payload = json.loads(report_path.read_text(encoding="utf-8"))
    checks = {check["name"]: check for check in payload["check_details"]}
    assert exit_status == 1
    assert payload["status"] == "failure"
    assert payload["checks"]["trusted_construction"] == "fail"
    assert checks["corrected-construction-and-provenance-path"]["status"] == "fail"
    assert (
        "trusted assertion evidence missing"
        in (checks["corrected-construction-and-provenance-path"]["message"])
    )
