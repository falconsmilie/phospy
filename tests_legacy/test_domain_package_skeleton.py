from __future__ import annotations

from importlib import import_module
from pathlib import Path

PACKAGE_NAMES = [
    "phospy.api",
    "phospy.datasets",
    "phospy.preprocessing",
    "phospy.prediction",
    "phospy.activities",
    "phospy.signalomes",
    "phospy.references",
    "phospy.io",
    "phospy.validation",
    "phospy.errors",
    "phospy.internal",
]


def test_domain_packages_exist_with_package_docstrings() -> None:
    for package_name in PACKAGE_NAMES:
        module = import_module(package_name)
        assert module.__doc__
        assert module.__doc__.strip()


def test_root_package_migration_map_exists() -> None:
    migration_map = Path("docs/architecture/root-package-migration-map.md")
    assert migration_map.exists()
    contents = migration_map.read_text(encoding="utf-8")
    assert "src/phospy/workflow.py" in contents
    assert "src/phospy/motifs.py" in contents
