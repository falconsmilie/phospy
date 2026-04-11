from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "src" / "phospy"


def test_signalome_domain_modules_live_under_signalomes_package() -> None:
    for relative_path in (
        "signalomes/analysis.py",
        "signalomes/assignments.py",
        "signalomes/clustering.py",
        "signalomes/maps.py",
        "signalomes/networks.py",
        "signalomes/results.py",
        "signalomes/site_ids.py",
    ):
        assert (ROOT / relative_path).exists()


def test_legacy_signalome_root_modules_have_been_removed() -> None:
    for relative_path in (
        "signalome_assignments.py",
        "signalome_clustering.py",
        "signalome_construction.py",
        "signalome_maps.py",
        "signalome_models.py",
        "signalome_networks.py",
        "signalome_site_ids.py",
        "signalomes.py",
    ):
        assert not (ROOT / relative_path).exists()
