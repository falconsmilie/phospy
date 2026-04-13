from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_package_layout_doc_describes_domain_packages() -> None:
    contents = (REPO_ROOT / "docs" / "architecture" / "package-layout.md").read_text(
        encoding="utf-8"
    )

    expected_packages = (
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
    )
    for package_name in expected_packages:
        assert package_name in contents


def test_contributing_links_to_package_layout_guidance() -> None:
    contents = (REPO_ROOT / ".github" / "CONTRIBUTING.md").read_text(encoding="utf-8")

    assert "docs/architecture/package-layout.md" in contents
    assert "phospy.api" in contents
    assert "phospy.internal" in contents


def test_root_package_docstring_marks_root_surface_as_intentional() -> None:
    contents = (REPO_ROOT / "src" / "phospy" / "__init__.py").read_text(
        encoding="utf-8"
    )

    assert "deliberately small convenience surface" in contents
    assert "Intentionally retained convenience exports" in contents


def test_contributing_documents_scientific_policy_expectations() -> None:
    contents = (REPO_ROOT / ".github" / "CONTRIBUTING.md").read_text(encoding="utf-8")

    assert "duplicate_site_strategy" in contents
    assert "missing_value_strategy" in contents
    assert "module_selection_strategy" in contents
