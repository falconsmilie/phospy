from __future__ import annotations

import importlib
import json
from dataclasses import asdict
from typing import get_type_hints

import numpy as np
import pandas as pd
import pytest

import phospy
import phospy.advanced as advanced_api
import phospy.api as public_api
import phospy.api.configs as stable_config_api
from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy._api_inventory import (
    ADVANCED_API_STABILITY_JUSTIFICATIONS,
    ADVANCED_CONFIG_API,
    ADVANCED_PUBLIC_API,
    STABLE_CONFIG_API,
    STABLE_PUBLIC_API,
)
from phospy.advanced import (
    DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1,
    DifferentialAnalysisConfig,
    DifferentialProteinAwareModelConfig,
)
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.errors import ContractValidationError, PhosPyInputError
from phospy.errors.validation import WorkflowValidationError
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.differential.models import EmpiricalBayesPriorDiagnostics
from phospy.science.differential.models.protein_aware import (
    ProteinAwareDifferentialComputationRequest,
    ProteinAwareDifferentialComputationResult,
)
from phospy.science.differential.models.tables import (
    DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_CONTRAST_NON_ESTIMABLE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_NON_FINITE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_MODEL_FIT_QUANTITIES_NON_FINITE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_MODEL_RESIDUAL_VARIANCE_NON_FINITE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_MODEL_RESIDUAL_VARIANCE_NON_POSITIVE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_EXCLUDED,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK,
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_MODEL_FIT_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
    validate_result_table_contract,
)
from phospy.workflows.differential.interpreter import _resolve_execution_config
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

METHOD_ID = DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1
ADVANCED_PROTEIN_AWARE_SYMBOLS = {
    "DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1",
    "DifferentialProteinAwareModelConfig",
}
UNSUPPORTED_PROTEIN_AWARE_CONFIG_SYMBOLS = {
    "DifferentialProteinAwareModelMethod",
    "SUPPORTED_DIFFERENTIAL_PROTEIN_AWARE_MODEL_METHODS",
}
SAMPLE_ORDER = ("A_1", "A_2", "B_1", "B_2")


def test_protein_aware_config_defaults_and_json_compatible_payload() -> None:
    ordinary = DifferentialAnalysisConfig()
    protein_aware = DifferentialProteinAwareModelConfig()
    selected = DifferentialAnalysisConfig(protein_aware_model=protein_aware)

    assert ordinary.protein_aware_model is None
    assert protein_aware.method == METHOD_ID
    assert selected.protein_aware_model is protein_aware

    payload = asdict(protein_aware)
    assert payload == {"method": METHOD_ID}
    assert json.loads(json.dumps(payload)) == payload


def test_protein_aware_config_rejects_unsupported_method() -> None:
    with pytest.raises(
        ContractValidationError,
        match="differential.protein_aware_model.method must be one of",
    ):
        DifferentialProteinAwareModelConfig(method="unsupported")  # type: ignore[arg-type]


def test_differential_config_rejects_wrong_protein_aware_config_type() -> None:
    with pytest.raises(
        ContractValidationError,
        match="differential.protein_aware_model must be DifferentialProteinAwareModelConfig",
    ):
        DifferentialAnalysisConfig(protein_aware_model=object())  # type: ignore[arg-type]


def test_protein_aware_config_is_advanced_only_api_surface() -> None:
    assert advanced_api.DifferentialProteinAwareModelConfig is (
        DifferentialProteinAwareModelConfig
    )
    assert (
        advanced_api.DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1
        == METHOD_ID
    )
    assert ADVANCED_PROTEIN_AWARE_SYMBOLS <= set(ADVANCED_CONFIG_API)
    assert ADVANCED_PROTEIN_AWARE_SYMBOLS <= set(ADVANCED_PUBLIC_API)
    assert ADVANCED_PROTEIN_AWARE_SYMBOLS <= set(advanced_api.__all__)
    assert ADVANCED_PROTEIN_AWARE_SYMBOLS <= set(ADVANCED_API_STABILITY_JUSTIFICATIONS)

    assert ADVANCED_PROTEIN_AWARE_SYMBOLS.isdisjoint(STABLE_CONFIG_API)
    assert ADVANCED_PROTEIN_AWARE_SYMBOLS.isdisjoint(STABLE_PUBLIC_API)
    assert ADVANCED_PROTEIN_AWARE_SYMBOLS.isdisjoint(public_api.__all__)
    assert ADVANCED_PROTEIN_AWARE_SYMBOLS.isdisjoint(stable_config_api.__all__)
    for symbol_name in ADVANCED_PROTEIN_AWARE_SYMBOLS:
        assert not hasattr(phospy, symbol_name)
        with pytest.raises(ImportError):
            exec(f"from phospy import {symbol_name}", {})


def test_protein_aware_config_has_no_stable_compatibility_import_routes() -> None:
    blocked_symbols = (
        ADVANCED_PROTEIN_AWARE_SYMBOLS | UNSUPPORTED_PROTEIN_AWARE_CONFIG_SYMBOLS
    )
    for module_name in (
        "phospy.api",
        "phospy.api.configs",
        "phospy.api.configs.differential",
    ):
        module = importlib.import_module(module_name)
        for symbol_name in blocked_symbols:
            assert not hasattr(module, symbol_name)
            with pytest.raises(ImportError):
                exec(f"from {module_name} import {symbol_name}", {})


def test_differential_request_fields_remain_unchanged() -> None:
    hints = get_type_hints(DifferentialAnalysisRequest)

    assert set(hints) == {"dataset", "design", "contrasts", "config"}
    assert hints["config"] is DifferentialAnalysisConfig


def test_resolved_execution_config_carries_normalized_protein_aware_method() -> None:
    ordinary = _resolve_execution_config(DifferentialAnalysisConfig())
    selected = _resolve_execution_config(
        DifferentialAnalysisConfig(
            protein_aware_model=DifferentialProteinAwareModelConfig()
        )
    )

    assert ordinary.protein_aware_method is None
    assert selected.protein_aware_method == METHOD_ID
    assert not isinstance(
        selected.protein_aware_method,
        DifferentialProteinAwareModelConfig,
    )


def test_selected_protein_aware_workflow_lane_requires_dataset_sidecar() -> None:
    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.protein_aware_inputs.sidecar_missing",
    ):
        DifferentialAnalysisWorkflow().run(
            DifferentialAnalysisRequest(
                dataset=_public_dataset(),
                design=_public_design(),
                contrasts=_public_contrasts(),
                config=DifferentialAnalysisConfig(
                    protein_aware_model=DifferentialProteinAwareModelConfig()
                ),
            )
        )


@pytest.mark.parametrize(
    ("mutation", "pattern"),
    [
        (
            "wrong_type",
            "protein_aware_model must be DifferentialProteinAwareModelConfig",
        ),
        ("unsupported_method", "protein_aware_model.method must be one of"),
    ],
)
def test_differential_validator_rejects_corrupted_protein_aware_config(
    mutation: str,
    pattern: str,
) -> None:
    if mutation == "wrong_type":
        config = DifferentialAnalysisConfig()
        object.__setattr__(config, "protein_aware_model", object())
    else:
        protein_aware = DifferentialProteinAwareModelConfig()
        object.__setattr__(protein_aware, "method", "unsupported")
        config = DifferentialAnalysisConfig(protein_aware_model=protein_aware)

    with pytest.raises(WorkflowValidationError, match=pattern):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=_public_dataset(),
                design=_public_design(),
                contrasts=_public_contrasts(),
                config=config,
            )
        )


def test_private_protein_aware_request_owns_and_validates_alignment() -> None:
    kwargs = _request_kwargs()
    request = ProteinAwareDifferentialComputationRequest(**kwargs)  # type: ignore[arg-type]

    assert request.sample_order == SAMPLE_ORDER
    assert request.method_id == METHOD_ID
    assert request.phosphosite_matrix.index.tolist() == ["site_a", "site_b"]
    assert request.matched_pairs["site_key"].tolist() == ["site_a", "site_b"]

    phosphosite_matrix = kwargs["phosphosite_matrix"]
    matched_pairs = kwargs["matched_pairs"]
    assert isinstance(phosphosite_matrix, pd.DataFrame)
    assert isinstance(matched_pairs, pd.DataFrame)
    phosphosite_matrix.iloc[0, 0] = 999.0
    matched_pairs.loc[0, "total_protein_row_key"] = "changed"

    assert float(request.phosphosite_matrix.iloc[0, 0]) == 10.0
    assert request.matched_pairs.loc[0, "total_protein_row_key"] == "protein_a"


@pytest.mark.parametrize(
    ("mutation", "pattern"),
    [
        ("reordered_samples", "phosphosite_matrix.columns must exactly match"),
        ("reordered_matched_pairs", "matched_pairs.site_key must exactly match"),
        ("duplicate_matched_pair_columns", "matched_pairs.columns must be unique"),
        (
            "missing_protein_covariate",
            "must contain each matched total_protein_row_key",
        ),
        ("unsupported_method", "method_id must be one of"),
    ],
)
def test_private_protein_aware_request_rejects_invalid_contracts(
    mutation: str,
    pattern: str,
) -> None:
    kwargs = _request_kwargs()
    if mutation == "reordered_samples":
        matrix = kwargs["phosphosite_matrix"]
        assert isinstance(matrix, pd.DataFrame)
        kwargs["phosphosite_matrix"] = matrix.loc[:, list(reversed(SAMPLE_ORDER))]
    elif mutation == "reordered_matched_pairs":
        pairs = kwargs["matched_pairs"]
        assert isinstance(pairs, pd.DataFrame)
        kwargs["matched_pairs"] = pairs.iloc[[1, 0], :].reset_index(drop=True)
    elif mutation == "duplicate_matched_pair_columns":
        pairs = kwargs["matched_pairs"]
        assert isinstance(pairs, pd.DataFrame)
        duplicate_pairs = pairs.copy(deep=True)
        duplicate_pairs.insert(
            1,
            "site_key",
            duplicate_pairs.loc[:, "site_key"],
            allow_duplicates=True,
        )
        kwargs["matched_pairs"] = duplicate_pairs
    elif mutation == "missing_protein_covariate":
        covariates = kwargs["resolved_protein_covariates"]
        assert isinstance(covariates, pd.DataFrame)
        kwargs["resolved_protein_covariates"] = covariates.drop(index="protein_b")
    elif mutation == "unsupported_method":
        kwargs["method_id"] = "unsupported"

    with pytest.raises(PhosPyInputError, match=pattern):
        ProteinAwareDifferentialComputationRequest(**kwargs)  # type: ignore[arg-type]


def test_private_protein_aware_request_from_owned_still_validates_alignment() -> None:
    kwargs = _request_kwargs()
    matrix = kwargs["phosphosite_matrix"]
    assert isinstance(matrix, pd.DataFrame)
    kwargs["phosphosite_matrix"] = matrix.loc[:, list(reversed(SAMPLE_ORDER))]

    with pytest.raises(PhosPyInputError, match="phosphosite_matrix.columns"):
        ProteinAwareDifferentialComputationRequest._from_owned(**kwargs)  # type: ignore[arg-type]


def test_private_protein_aware_result_owns_alignment_without_global_decomposition() -> (
    None
):
    kwargs = _result_kwargs()
    result = ProteinAwareDifferentialComputationResult(**kwargs)  # type: ignore[arg-type]

    assert result.tested_site_ids == ("site_a", "site_b")
    assert result.method_id == METHOD_ID
    assert not hasattr(result, "design_decomposition")
    assert result.table_for("B_vs_A").index.tolist() == ["site_a", "site_b"]
    assert result.protein_coefficient_series().tolist() == [0.25, -0.5]

    exported_table = result.table_for("B_vs_A")
    exported_site_diagnostics = result.site_diagnostics_dataframe()
    exported_table.iloc[0, 0] = 999.0
    exported_site_diagnostics.loc["site_a", "protein_coefficient"] = 999.0

    assert float(result.table_for("B_vs_A").iloc[0, 0]) == 1.0
    assert (
        float(
            result.site_diagnostics_dataframe().loc[
                "site_a",
                "protein_coefficient",
            ]
        )
        == 0.25
    )


@pytest.mark.parametrize(
    ("mutation", "pattern"),
    [
        ("contrast_index", "contrast_tables\\['B_vs_A'\\].index"),
        ("site_protein_coefficient", "site_diagnostics.protein_coefficient"),
        (
            "missing_group_diagnostic",
            "must contain each per-site total_protein_row_key",
        ),
    ],
)
def test_private_protein_aware_result_rejects_invalid_alignment(
    mutation: str,
    pattern: str,
) -> None:
    kwargs = _result_kwargs()
    if mutation == "contrast_index":
        tables = kwargs["contrast_tables"]
        assert isinstance(tables, dict)
        tables["B_vs_A"] = tables["B_vs_A"].iloc[[1, 0], :]
    elif mutation == "site_protein_coefficient":
        site_diagnostics = kwargs["site_diagnostics"]
        assert isinstance(site_diagnostics, pd.DataFrame)
        site_diagnostics.loc["site_a", "protein_coefficient"] = 999.0
    elif mutation == "missing_group_diagnostic":
        group_diagnostics = kwargs["augmented_design_diagnostics"]
        assert isinstance(group_diagnostics, pd.DataFrame)
        kwargs["augmented_design_diagnostics"] = group_diagnostics.drop(
            index="protein_b"
        )

    with pytest.raises(PhosPyInputError, match=pattern):
        ProteinAwareDifferentialComputationResult(**kwargs)  # type: ignore[arg-type]


def test_private_protein_aware_result_assume_owned_still_validates_alignment() -> None:
    kwargs = _result_kwargs()
    contrast_tables = kwargs["contrast_tables"]
    assert isinstance(contrast_tables, dict)
    contrast_tables["B_vs_A"] = contrast_tables["B_vs_A"].iloc[[1, 0], :]

    with pytest.raises(PhosPyInputError, match="contrast_tables\\['B_vs_A'\\].index"):
        ProteinAwareDifferentialComputationResult(
            **kwargs,
            _assume_owned=True,
        )  # type: ignore[arg-type]


def test_protein_aware_withheld_statuses_pass_result_table_contract() -> None:
    table = _protein_aware_status_table()

    validate_result_table_contract(
        table,
        field_name="protein_aware_differential_result_table",
    )


def test_protein_aware_status_requires_matching_stable_reason_code() -> None:
    table = _protein_aware_status_table()
    table.loc[table.index[1], DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] = (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_NON_FINITE
    )

    with pytest.raises(
        PhosPyInputError,
        match="unsupported protein-aware reason codes",
    ):
        validate_result_table_contract(
            table,
            field_name="protein_aware_differential_result_table",
        )


def test_protein_aware_reason_code_requires_protein_aware_status() -> None:
    table = _protein_aware_status_table()
    table.loc[table.index[0], DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] = (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK
    )

    with pytest.raises(
        PhosPyInputError,
        match="unsupported protein-aware reason codes",
    ):
        validate_result_table_contract(
            table,
            field_name="protein_aware_differential_result_table",
        )


def test_protein_aware_withheld_status_requires_missing_statistics() -> None:
    table = _protein_aware_status_table()
    table.loc[table.index[1], "P.Value"] = 0.5

    with pytest.raises(
        PhosPyInputError,
        match="withheld rows must contain missing values",
    ):
        validate_result_table_contract(
            table,
            field_name="protein_aware_differential_result_table",
        )


def _request_kwargs() -> dict[str, object]:
    phosphosite_matrix = pd.DataFrame(
        {
            "A_1": [10.0, 8.0],
            "A_2": [10.2, 8.1],
            "B_1": [11.0, 8.8],
            "B_2": [11.1, 8.7],
        },
        index=pd.Index(["site_a", "site_b"], name="site_key"),
    )
    design = pd.DataFrame(
        {
            "condition_A": [1.0, 1.0, 0.0, 0.0],
            "condition_B": [0.0, 0.0, 1.0, 1.0],
        },
        index=pd.Index(SAMPLE_ORDER, name="sample_id"),
    )
    contrasts = pd.DataFrame(
        {"B_vs_A": [-1.0, 1.0]},
        index=pd.Index(["condition_A", "condition_B"], name="coefficient"),
    )
    matched_pairs = pd.DataFrame(
        {
            "site_key": ["site_a", "site_b"],
            "protein_identifier": ["P_A", "P_B"],
            "total_protein_row_key": ["protein_a", "protein_b"],
        }
    )
    resolved_protein_covariates = pd.DataFrame(
        {
            "A_1": [9.0, 7.0],
            "A_2": [9.1, 7.2],
            "B_1": [9.4, 7.4],
            "B_2": [9.5, 7.5],
        },
        index=pd.Index(["protein_a", "protein_b"], name="total_protein_row_key"),
    )
    return {
        "phosphosite_matrix": phosphosite_matrix,
        "base_design": design,
        "base_contrasts": contrasts,
        "sample_order": SAMPLE_ORDER,
        "matched_pairs": matched_pairs,
        "resolved_protein_covariates": resolved_protein_covariates,
        "method_id": METHOD_ID,
    }


def _result_kwargs() -> dict[str, object]:
    index = pd.Index(["site_a", "site_b"], name="site_key")
    residual_variance = pd.Series([0.2, 0.4], index=index.copy())
    prior_variance = pd.Series([0.3, 0.3], index=index.copy())
    prior_dof = pd.Series([8.0, 8.0], index=index.copy())
    protein_coefficient = pd.Series(
        [0.25, -0.5],
        index=index.copy(),
        name="protein_coefficient",
    )
    return {
        "residual_variance": residual_variance,
        "posterior_residual_variance": pd.Series([0.25, 0.35], index=index.copy()),
        "prior_residual_variance": prior_variance,
        "prior_degrees_of_freedom_series_value": prior_dof,
        "prior_variance": 0.3,
        "prior_degrees_of_freedom": 8.0,
        "residual_degrees_of_freedom": 1.0,
        "empirical_bayes_method": "standard",
        "empirical_bayes_robust": False,
        "empirical_bayes_trend": False,
        "prior_diagnostics": EmpiricalBayesPriorDiagnostics(
            method="standard",
            robust=False,
            trend=False,
            winsor_tail_p=(0.05, 0.1),
            base_prior_variance=0.3,
            base_prior_degrees_of_freedom=8.0,
            robust_outlier_count=0,
            robust_outlier_fraction=0.0,
            winsorized_low_count=0,
            winsorized_high_count=0,
            prior_variance=prior_variance,
            prior_degrees_of_freedom=prior_dof,
        ),
        "mean_variance_trend_diagnostics": None,
        "contrast_tables": {
            "B_vs_A": pd.DataFrame(
                {
                    "logFC": [1.0, -0.2],
                    "t": [3.0, -1.5],
                    "P.Value": [0.02, 0.2],
                    "adj.P.Val": [0.04, 0.2],
                },
                index=index.copy(),
            )
        },
        "protein_coefficient": protein_coefficient,
        "site_diagnostics": pd.DataFrame(
            {
                "site_key": ["site_a", "site_b"],
                "total_protein_row_key": ["protein_a", "protein_b"],
                "protein_coefficient": [0.25, -0.5],
            },
            index=index.copy(),
        ),
        "augmented_design_diagnostics": pd.DataFrame(
            {
                "total_protein_row_key": ["protein_a", "protein_b"],
                "sample_count": [4, 4],
                "coefficient_count": [3, 3],
                "rank": [3, 3],
                "residual_degrees_of_freedom": [1.0, 1.0],
                "condition_number": [5.0, 6.0],
                "max_condition_number": [1.0e10, 1.0e10],
            },
            index=pd.Index(
                ["protein_a", "protein_b"],
                name="total_protein_row_key",
            ),
        ),
        "tested_site_ids": ("site_a", "site_b"),
        "method_id": METHOD_ID,
    }


def _protein_aware_status_table() -> pd.DataFrame:
    site_index = protein_site_key_index(
        protein_identifiers=(
            "MAPK14",
            "AKT1",
            "GSK3B",
            "PRKACA",
            "RPS6",
            "EIF4EBP1",
            "RAF1",
            "MTOR",
            "BAD",
        ),
        sites=("Y182", "T308", "S9", "S339", "S235", "T37", "S259", "S2448", "S136"),
    )
    statuses = [
        DIFFERENTIAL_RESULT_STATUS_TESTED,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_MODEL_FIT_INVALID,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_MODEL_FIT_INVALID,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_MODEL_FIT_INVALID,
    ]
    reasons = [
        "",
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_EXCLUDED,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_CONTRAST_NON_ESTIMABLE,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_MODEL_RESIDUAL_VARIANCE_NON_POSITIVE,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_MODEL_RESIDUAL_VARIANCE_NON_FINITE,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_MODEL_FIT_QUANTITIES_NON_FINITE,
    ]
    table = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": [
                "MAPK14;Y182;",
                "AKT1;T308;",
                "GSK3B;S9;",
                "PRKACA;S339;",
                "RPS6;S235;",
                "EIF4EBP1;T37;",
                "RAF1;S259;",
                "MTOR;S2448;",
                "BAD;S136;",
            ],
            **site_key_context_columns(site_index),
            "gene_symbol": [
                "MAPK14",
                "AKT1",
                "GSK3B",
                "PRKACA",
                "RPS6",
                "EIF4EBP1",
                "RAF1",
                "MTOR",
                "BAD",
            ],
            "site": [
                "Y182",
                "T308",
                "S9",
                "S339",
                "S235",
                "T37",
                "S259",
                "S2448",
                "S136",
            ],
            "logFC": [
                1.0,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            ],
            "t": [2.0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            "P.Value": [
                0.04,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            ],
            "adj.P.Val": [
                0.04,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            ],
            DIFFERENTIAL_RESULT_STATUS_COLUMN: statuses,
            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN: reasons,
        },
        index=site_index.copy(),
    )
    return table


def _public_dataset():
    phospho = pd.DataFrame(
        {
            "A_1": [10.0, 8.0],
            "A_2": [10.2, 8.1],
            "B_1": [11.0, 8.8],
            "B_2": [11.1, 8.7],
        },
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "protein_id": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "localisation_confidence": [0.95, 0.95],
            "site_sequence": [
                ("A" * 15) + "Y" + ("A" * 15),
                ("A" * 15) + "S" + ("A" * 15),
            ],
        },
        index=phospho.index.copy(),
    )
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="log2",
        )
    )


def _public_design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_r2",
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_r2",
            ),
        )
    )


def _public_contrasts() -> tuple[Contrast, ...]:
    return (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
    )
