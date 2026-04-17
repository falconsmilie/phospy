from __future__ import annotations

from pathlib import Path


def test_preprocessing_domain_contains_moved_core_modules() -> None:
    root = Path(__file__).resolve().parents[1] / "src" / "phospy"
    expected = [
        root / "preprocessing" / "core.py",
        root / "preprocessing" / "dataset.py",
        root / "preprocessing" / "services.py",
        root / "preprocessing" / "site_matrix.py",
        root / "preprocessing" / "steps.py",
        root / "preprocessing" / "modes.py",
        root / "preprocessing" / "primitives.py",
        root / "preprocessing" / "protein_correction.py",
    ]

    assert all(path.exists() for path in expected)


def test_legacy_flat_preprocessing_modules_have_been_removed() -> None:
    root = Path(__file__).resolve().parents[1] / "src" / "phospy"
    legacy = [
        root / "_preprocessing_primitives.py",
        root / "_protein_correction.py",
        root / "core_processing.py",
        root / "dataset_preprocessing.py",
        root / "preprocessing.py",
        root / "preprocessing_services.py",
        root / "site_matrix_builder.py",
    ]

    assert all(not path.exists() for path in legacy)
