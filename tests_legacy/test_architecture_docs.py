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


def test_root_package_docstring_marks_no_reexport_policy_as_intentional() -> None:
    contents = (REPO_ROOT / "src" / "phospy" / "__init__.py").read_text(
        encoding="utf-8"
    )

    assert "does not re-export domain APIs" in contents
    assert "phospy.api" in contents


def test_contributing_documents_scientific_policy_expectations() -> None:
    contents = (REPO_ROOT / ".github" / "CONTRIBUTING.md").read_text(encoding="utf-8")

    assert "duplicate_site_strategy" in contents
    assert "missing_data_policy" in contents
    assert "missing_value_strategy" in contents
    assert "module_selection_strategy" in contents


def test_api_docs_describe_supported_workflow_result_contracts() -> None:
    contents = (REPO_ROOT / "docs" / "api.md").read_text(encoding="utf-8")

    assert "Supported public lane" in contents
    assert "Internal lane (advanced contributors only)" not in contents
    assert "SimpleKinaseWorkflowResult" in contents
    assert "`KinaseWorkflowResult`" not in contents
    assert "You do not need a separate predMat workflow." in contents
    assert "prediction_result.pred_mat_result" in contents
    assert "prediction_result" in contents
    assert "scoring_result" in contents
    assert "kinase_activity_result" in contents
    assert "SimpleKinaseExecutionGraph" in contents
    assert "EnsemblePredictorContract" in contents
    assert "KinasePredictionDebugTrace" in contents


def test_readme_describes_supported_result_access_patterns() -> None:
    contents = (REPO_ROOT / "README.md").read_text(encoding="utf-8")

    assert (
        "SimpleKinaseWorkflow.run(...)` returns `SimpleKinaseWorkflowResult`"
        in contents
    )
    assert (
        "Use `result.pred_mat_result` for the canonical predMat table contract."
        in contents
    )
    assert (
        "Use `result.prediction_result` when you need full prediction payload details"
        in contents
    )
    assert (
        "Use `result.scoring_result` for `profile_scores`, `combined_scores`, and `weights`."
        in contents
    )
    assert "Use `result.kinase_activity_result` for activity summaries" in contents
    assert "supported public lane is `phospy.api.SimpleKinaseWorkflow`" in contents
