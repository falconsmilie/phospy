"""Documentation contracts for differential handling of imputed datasets."""

from __future__ import annotations

import inspect
from pathlib import Path

from phospy.differential.public import DifferentialAnalysisWorkflow
from phospy.preprocessing.contracts.configuration.missing_data import (
    DatasetMissingDataConfig,
)

_ROOT = Path(__file__).resolve().parents[2]
_DOCS = _ROOT / "docs"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _section(text: str, start: str, end: str) -> str:
    return text.split(start, maxsplit=1)[1].split(end, maxsplit=1)[0]


def _normalise(text: str) -> str:
    return " ".join(text.split()).casefold()


def test_public_workflow_docstring_exposes_imputed_dataset_boundary() -> None:
    docstring = inspect.getdoc(DifferentialAnalysisWorkflow)

    assert docstring is not None
    assert 'imputed_value_policy="reject"' in docstring
    assert '"withhold_imputed_features"' in docstring
    assert "by default" in docstring.casefold()


def test_missing_data_config_docstring_conditions_imputed_dataset_state() -> None:
    docstring = inspect.getdoc(DatasetMissingDataConfig)

    assert docstring is not None
    normalised = _normalise(docstring)
    assert "when it actually fills cells" in normalised
    assert "datasets carrying those imputed cells" in normalised
    assert 'imputed_value_policy="withhold_imputed_features"' in docstring
    assert "can retain imputed cells in tested features" in normalised
    assert "not observed-only fitting" in normalised


def test_differential_guide_warns_before_detailed_policy_section() -> None:
    guide = _read(_DOCS / "api" / "differential-analysis.md")
    headings = [line for line in guide.splitlines() if line.startswith("## ")]
    input_requirements = _section(
        guide,
        "## Input Requirements and Policies",
        "## Empirical Bayes Variance Moderation",
    )
    detailed_policy_position = guide.index(
        "### Missing, Imputed, and Authoritative Matrix Policy"
    )

    assert "## Before You Begin" not in headings
    assert "## Input Requirements and Policies" in headings
    assert len(headings) == len(set(headings))
    assert 'imputed_value_policy="reject"' in input_requirements
    assert guide.index('imputed_value_policy="reject"') < detailed_policy_position
    assert '"withhold_imputed_features"' in input_requirements
    assert "imputed_value_max_fraction" in input_requirements
    assert "minimum_condition_replicates" in input_requirements
    assert "not observed-only fitting" in input_requirements
    assert "retained tested features can still contain imputed cells" in (
        input_requirements.casefold()
    )


def test_dataset_missing_data_section_connects_every_imputation_policy() -> None:
    guide = _read(_DOCS / "api" / "dataset-build-workflow.md")
    missing_data = _section(guide, "### Missing Data", "#### Group-Aware KNN + MinProb")

    for policy in (
        '"impute_row_median"',
        '"impute_minprob"',
        '"impute_knn"',
        '"impute_group_aware"',
    ):
        assert policy in missing_data
    assert "Successful dataset construction does not itself authorise" in missing_data
    assert "DifferentialAnalysisWorkflow" in missing_data
    assert (
        "differential-analysis.md#missing-imputed-and-authoritative-matrix-policy"
        in (missing_data)
    )


def test_quickstart_directs_imputation_users_to_differential_policy() -> None:
    quickstart = _read(_DOCS / "quickstart.md")
    normalised = _normalise(quickstart)

    assert "This quickstart uses non-imputed data" in quickstart
    assert "if you select any dataset imputation policy" in normalised
    assert (
        "api/differential-analysis.md#missing-imputed-and-authoritative-matrix-policy"
        in (quickstart)
    )
