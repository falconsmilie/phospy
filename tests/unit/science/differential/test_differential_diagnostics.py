from __future__ import annotations

import json

import pandas as pd
import pytest

from phospy.errors import PhosPyInputError
from phospy.science.differential.models import (
    DifferentialAnalysisRequest,
    DifferentialContrastDefinition,
    DifferentialModelDiagnostics,
)
from phospy.science.differential.validator import DifferentialAnalysisRequestValidator


def test_differential_model_diagnostics_payload_records_required_fields() -> None:
    contrast = DifferentialContrastDefinition(
        name="B_vs_A",
        numerator_condition="B",
        denominator_condition="A",
        coefficients=(("A", -1.0), ("B", 1.0)),
    )
    diagnostics = DifferentialModelDiagnostics(
        model_type="moderated_ols_fixed_effect",
        design_columns=("A", "B"),
        contrast_definitions=(contrast,),
        rank=2,
        n_samples=4,
        n_sites=3,
        residual_degrees_of_freedom=2.0,
        variance_method="ordinary_least_squares_residual_variance",
        moderation_method="empirical_bayes_standard",
        multiple_testing_method="benjamini_hochberg",
        imputation_policy="reject",
        missing_value_policy="reject_missing_values_before_differential_execution",
        intensity_scale="log2",
        normalisation_state="none",
        batch_or_covariate_terms=(),
        unsupported_assumptions=(
            "mixed-effects differential modelling is unsupported",
        ),
        warnings=("scope is limited to validated fixed-effect designs",),
    )

    payload = diagnostics.to_payload()

    assert set(payload) == {
        "model_type",
        "design_columns",
        "contrast_definitions",
        "rank",
        "n_samples",
        "n_sites",
        "residual_degrees_of_freedom",
        "decomposition_method",
        "solver",
        "column_scale_method",
        "rank_tolerance_policy",
        "rank_tolerance",
        "condition_number",
        "max_condition_number",
        "singular_values",
        "variance_method",
        "moderation_method",
        "multiple_testing_method",
        "imputation_policy",
        "missing_value_policy",
        "intensity_scale",
        "normalisation_state",
        "batch_or_covariate_terms",
        "unsupported_assumptions",
        "warnings",
    }
    assert payload["contrast_definitions"] == [
        {
            "name": "B_vs_A",
            "numerator_condition": "B",
            "denominator_condition": "A",
            "coefficients": [
                {"coefficient": "A", "weight": -1.0},
                {"coefficient": "B", "weight": 1.0},
            ],
            "description": "",
        }
    ]
    json.dumps(payload)


def test_science_validator_rejects_contrast_with_unknown_design_column() -> None:
    matrix = pd.DataFrame(
        {
            "A_1": [1.0, 2.0],
            "A_2": [1.2, 2.1],
            "B_1": [2.0, 2.2],
            "B_2": [2.1, 2.4],
        },
        index=pd.Index(["site_1", "site_2"], name="site_id"),
    )
    design = pd.DataFrame(
        {
            "A": [1.0, 1.0, 0.0, 0.0],
            "B": [0.0, 0.0, 1.0, 1.0],
        },
        index=pd.Index(matrix.columns, name="sample"),
    )
    contrasts = pd.DataFrame(
        {"B_vs_A": [-1.0, 1.0]},
        index=pd.Index(["A", "unknown_column"], name="coefficient"),
    )
    request = DifferentialAnalysisRequest(
        matrix=matrix,
        design=design,
        contrasts=contrasts,
    )

    with pytest.raises(PhosPyInputError, match="contrasts.index.*design.columns"):
        DifferentialAnalysisRequestValidator().run(request)
