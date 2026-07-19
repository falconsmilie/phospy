from __future__ import annotations

import importlib.util
import json
import sys
from hashlib import sha256
from pathlib import Path
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parents[2]
MANIFEST_WRITER_PATH = ROOT / "scripts" / "write_build_manifest.py"

pytestmark = pytest.mark.release_gate


def _load_manifest_writer() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "write_build_manifest_under_test",
        MANIFEST_WRITER_PATH,
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_build_manifest_binds_source_identity_to_wheel_and_sdist(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    writer = _load_manifest_writer()
    dist_dir = tmp_path / "dist"
    dist_dir.mkdir()
    wheel = dist_dir / "phospy-1.6.0-py3-none-any.whl"
    sdist = dist_dir / "phospy-1.6.0.tar.gz"
    wheel.write_bytes(b"wheel bytes")
    sdist.write_bytes(b"sdist bytes")
    output_path = tmp_path / "build" / "reports" / "build-manifest.json"

    monkeypatch.setattr(
        writer,
        "_source_identity_digest",
        lambda repository_root: "sha256:" + ("a" * 64),
    )

    written_path = writer.write_build_manifest(
        dist_dir=dist_dir,
        output_path=output_path,
        repository_root=tmp_path,
        package_version="1.6.0",
    )

    payload = json.loads(written_path.read_text(encoding="utf-8"))
    artifacts = {artifact["kind"]: artifact for artifact in payload["artifacts"]}
    assert written_path == output_path
    assert payload["schema"] == "phospy.build-manifest/v1"
    assert payload["source_identity_digest"] == "sha256:" + ("a" * 64)
    assert payload["package_version"] == "1.6.0"
    assert artifacts == {
        "wheel": {
            "kind": "wheel",
            "filename": wheel.name,
            "sha256": sha256(wheel.read_bytes()).hexdigest(),
        },
        "sdist": {
            "kind": "sdist",
            "filename": sdist.name,
            "sha256": sha256(sdist.read_bytes()).hexdigest(),
        },
    }
