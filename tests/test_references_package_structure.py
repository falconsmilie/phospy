from __future__ import annotations

from pathlib import Path


def test_reference_domain_modules_exist() -> None:
    package_dir = Path(__file__).resolve().parents[1] / "src" / "phospy" / "references"
    assert (package_dir / "models.py").exists()
    assert (package_dir / "resources.py").exists()
    assert (package_dir / "resolution.py").exists()


def test_reference_bundle_provider_is_no_longer_owned_by_motifs_module() -> None:
    motifs_module = Path(__file__).resolve().parents[1] / "src" / "phospy" / "motifs.py"
    source = motifs_module.read_text(encoding="utf-8")
    assert "class ReferenceBundle(" not in source
    assert "class BundledReferenceProvider" not in source
