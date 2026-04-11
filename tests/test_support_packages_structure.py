from __future__ import annotations

from pathlib import Path


def test_support_packages_exist() -> None:
    package_root = Path(__file__).resolve().parents[1] / "src" / "phospy"

    assert (package_root / "io" / "readers.py").exists()
    assert (package_root / "io" / "mappings.py").exists()
    assert (package_root / "io" / "publishing.py").exists()
    assert (package_root / "io" / "writers.py").exists()
    assert (package_root / "errors" / "base.py").exists()
    assert (package_root / "errors" / "validation.py").exists()
    assert (package_root / "internal" / "constants.py").exists()
    assert (package_root / "internal" / "types.py").exists()


def test_legacy_support_modules_have_been_removed() -> None:
    package_root = Path(__file__).resolve().parents[1] / "src" / "phospy"

    assert not (package_root / "io.py").exists()
    assert not (package_root / "publishing.py").exists()
    assert not (package_root / "writers.py").exists()
    assert not (package_root / "constants.py").exists()
    assert not (package_root / "types.py").exists()
    assert not (package_root / "validation" / "errors.py").exists()
