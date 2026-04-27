from __future__ import annotations

import importlib.util
import inspect
from dataclasses import fields
from pathlib import Path
from typing import get_type_hints

import pandas as pd
import pytest

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11
    import tomli as tomllib

import phospy
import phospy.api as public_api
import phospy.io as phospy_io
from phospy import (
    AnalysisReadyDatasetBuilder,
)
from phospy.api import (
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    Organism,
    PhosPyInputError,
)
from phospy.api.configs import (
    DATASET_SITE_MATRIX_MISSING_DATA_POLICIES,
    DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING,
    DatasetMissingDataConfig,
)
from phospy.datasets.models import AnalysisReadyPhosphoDataset

ROOT = Path(__file__).resolve().parents[2]


def _public_methods(cls: type[object]) -> set[str]:
    return {
        name
        for name, value in inspect.getmembers(cls)
        if callable(value) and not name.startswith("_")
    }


def test_packaging_contract_excludes_archived_package_namespace() -> None:
    pyproject = ROOT / "pyproject.toml"
    config = tomllib.loads(pyproject.read_text(encoding="utf-8"))
    package_find = config["tool"]["setuptools"]["packages"]["find"]
    assert set(package_find["include"]) == {"phospy", "phospy.*"}
    assert {"phospy_legacy", "phospy_legacy.*"}.issubset(set(package_find["exclude"]))


def test_source_tree_does_not_expose_archived_namespace() -> None:
    assert importlib.util.find_spec("phospy_legacy") is None


def test_public_dataset_ingestion_story_is_builder_only() -> None:
    dataset_exports = {
        name for name in phospy.__all__ if "Dataset" in name or name.endswith("Builder")
    }
    dataset_exports = {name for name in dataset_exports if not name.endswith("Error")}
    assert dataset_exports == {
        "AnalysisReadyDatasetBuilder",
        "AnalysisReadyPhosphoDataset",
    }
    assert not hasattr(phospy, "build_dataset_from_files")
    assert not hasattr(phospy_io, "build_dataset_from_files")


def test_builder_exposes_only_run_request_contract() -> None:
    assert _public_methods(AnalysisReadyDatasetBuilder) == {"run"}
    assert not hasattr(AnalysisReadyDatasetBuilder, "execute")
    hints = get_type_hints(AnalysisReadyDatasetBuilder.run)
    assert hints["request"] is DatasetBuildRequest
    assert hints["return"] is AnalysisReadyPhosphoDataset


def test_dataset_build_request_excludes_user_declared_intensity_scale_state() -> None:
    request_fields = {field.name for field in fields(DatasetBuildRequest)}
    assert "intensity_scale_state" not in request_fields


def test_dataset_build_request_rejects_user_declared_intensity_scale_state() -> None:
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        DatasetBuildRequest(
            phospho=object(),
            site_metadata=object(),
            intensity_scale_state=object(),  # type: ignore[call-arg]
        )


def test_dataset_preprocessing_config_is_top_level_public_type() -> None:
    assert "DatasetPreprocessingConfig" not in phospy.__all__
    assert "DatasetPreprocessingConfig" in public_api.__all__
    assert DatasetPreprocessingConfig().missing_data.policy == "forbid"


def test_site_matrix_missing_data_public_contract_is_complete_case_only() -> None:
    assert DATASET_SITE_MATRIX_MISSING_DATA_POLICIES == frozenset(
        {DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING}
    )
    config = DatasetSiteMatrixConfig(policy="build_from_metadata")
    assert (
        config.missing_data_policy
        == DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING
    )
    assert config.minimum_observed_values is None


def test_public_dataset_configs_reject_unsupported_local_state_on_construction() -> (
    None
):
    with pytest.raises(PhosPyInputError, match="minimum_observed_values"):
        DatasetSiteMatrixConfig(
            policy="build_from_metadata",
            minimum_observed_values=1,
        )
    with pytest.raises(PhosPyInputError, match="missing_data.min_observed_values"):
        DatasetMissingDataConfig(policy="forbid", min_observed_values=2)


def test_site_matrix_config_docstring_foregrounds_strict_public_lane() -> None:
    docstring = DatasetSiteMatrixConfig.__doc__
    assert docstring is not None
    assert "missing-value-free `AnalysisReadyPhosphoDataset`" in docstring
    assert "internal-only state" in docstring


def test_site_matrix_build_contract_is_row_wise_for_mixed_sequence_support() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [1.5, 2.5, 3.5],
        },
        index=pd.Index(["MAPK14;Y182;", "FAKE1;S1;", "GSK3B;S9;"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "FAKE1", "GSK3B"],
            "site": ["Y182", "S1", "S9"],
        },
        index=phospho.index.copy(),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
            ),
        )
    )

    assert built.phospho.index.tolist() == ["GSK3B;S9;", "MAPK14;Y182;"]
    assert built.preprocessing_report is not None
    assert built.preprocessing_report.row_audit is not None
    dropped = built.preprocessing_report.row_audit.loc[
        (built.preprocessing_report.row_audit.loc[:, "stage"] == "site_matrix")
        & (built.preprocessing_report.row_audit.loc[:, "action"] == "dropped")
    ]
    assert set(dropped.loc[:, "source_row_id"].astype(str)) == {"FAKE1;S1;"}


def test_site_matrix_build_contract_retains_all_rows_when_fully_resolvable() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=pd.Index(["row_a", "row_b"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
        },
        index=phospho.index.copy(),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
            ),
        )
    )

    assert built.phospho.index.tolist() == ["GSK3B;S9;", "MAPK14;Y182;"]
    assert built.preprocessing_report is not None
    assert built.preprocessing_report.row_audit is not None
    dropped = built.preprocessing_report.row_audit.loc[
        (built.preprocessing_report.row_audit.loc[:, "stage"] == "site_matrix")
        & (built.preprocessing_report.row_audit.loc[:, "action"] == "dropped")
    ]
    assert dropped.empty


def test_site_matrix_build_contract_reports_empty_when_fully_unresolvable() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=pd.Index(["row_a", "row_b"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["FAKE1", "FAKE2"],
            "site": ["S1", "T2"],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "site-matrix construction produced no retained rows after filtering; "
            "input_rows=2, dropped_missing_sequence=2"
        ),
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                preprocessing_config=DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
                ),
            )
        )
