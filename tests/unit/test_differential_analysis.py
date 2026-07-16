from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api import (
    CategoricalCovariate,
    ContinuousCovariate,
    Contrast,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    EmpiricalBayesConfig,
    ExperimentalDesign,
    MultipleTestingConfig,
    MultipleTestingMethod,
    Organism,
    SampleDesignRecord,
)
from phospy.api.configs import SUPPORTED_MULTIPLE_TESTING_METHODS
from phospy.errors import (
    ContractValidationError,
    PhosPyInputError,
    WorkflowValidationError,
)
from phospy.science.datasets.builders.preprocessing import (
    build_dataset_processing_state,
)
from phospy.science.datasets.builders.transformation_resolver import (
    DatasetIntensityScaleResolver,
)
from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT,
)
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    encode_site_key,
)
from phospy.science.statistics.multiple_testing import adjust_p_values
from phospy.science.transformations.models import (
    IntensityScaleEstablishmentMode,
    IntensityScaleEstablishmentSource,
    IntensityScaleState,
    MatrixIntensityScaleState,
)
from phospy.science.transformations.transformers import (
    IdentityTransformer,
    Log2Transformer,
)
from phospy.workflows.differential.interpreter import DifferentialAnalysisInterpreter
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
)
from tests.support.site_keys import site_key_context_columns


def _matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "A_1": [1.0, 2.0, 1.0, 1.2, 2.2],
            "A_2": [1.1, 2.1, 1.1, 1.1, 2.0],
            "B_1": [2.1, 2.0, 1.0, 1.3, 2.1],
            "B_2": [2.0, 2.2, 0.9, 1.2, 2.3],
            "C_1": [1.0, 2.0, 3.0, 1.5, 2.4],
            "C_2": [1.1, 2.1, 3.1, 1.4, 2.2],
        },
        index=pd.Index(
            [
                "MAPK14;Y182;",
                "GSK3B;S9;",
                "AKT1;T308;",
                "RPS6KB1;T389;",
                "MTOR;S2448;",
            ],
            name="site_id",
        ),
    )


def _site_key_for_display_id(
    display_id: str,
    *,
    protein_id: str | None = None,
) -> str:
    parts = [token.strip() for token in display_id.split(";") if token.strip()]
    gene_symbol = parts[0]
    site = parts[1]
    key = build_protein_scoped_site_key(
        organism="rat",
        protein_namespace="protein_id",
        protein_identifier=protein_id or gene_symbol,
        residue=site.upper()[0],
        position=int(site[1:]),
        field_name="tests.unit.test_differential_analysis.site_key",
        error_type=ValueError,
    )
    return encode_site_key(key)


def _dataset(
    matrix: pd.DataFrame | None = None,
    *,
    intensity_scale_state: IntensityScaleState | None = None,
):
    phospho = _matrix() if matrix is None else matrix
    display_ids = phospho.index.astype(str).tolist()
    gene_site = [site_id.split(";") for site_id in display_ids]
    protein_ids = [parts[0] for parts in gene_site]
    site_keys = [
        _site_key_for_display_id(display_id, protein_id=protein_id)
        for display_id, protein_id in zip(display_ids, protein_ids, strict=True)
    ]
    phospho = phospho.copy(deep=True)
    phospho.index = pd.Index(site_keys, name="site_key")
    site_metadata = pd.DataFrame(
        {
            "site_key": site_keys,
            "display_id": display_ids,
            **site_key_context_columns(site_keys),
            "gene_symbol": [parts[0] for parts in gene_site],
            "site": [parts[1] for parts in gene_site],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in [parts[1] for parts in gene_site]
            ],
            "protein_id": protein_ids,
        },
        index=phospho.index.copy(),
    )
    return supported_dataset(
        phospho=phospho,
        site_metadata=site_metadata,
        intensity_scale_state=intensity_scale_state,
    )


def supported_dataset(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    intensity_scale_state: IntensityScaleState | None = None,
):
    from phospy import AnalysisReadyPhosphoDataset

    if intensity_scale_state is None:
        intensity_scale_state = supported_log2_intensity_scale_state(
            has_total_matrix=False
        )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=intensity_scale_state,
        processing_state=build_dataset_processing_state(
            plan=PreprocessingPlan.default(),
            intensity_scale_state=intensity_scale_state,
        ),
    )


def _declared_log2_intensity_scale_state(
    phospho: pd.DataFrame,
) -> IntensityScaleState:
    declared_state = IntensityScaleState(
        phospho=MatrixIntensityScaleState.log2(established_by="test.declaration"),
        total=None,
    )
    return (
        DatasetIntensityScaleResolver(transformer=IdentityTransformer())
        .run(
            phospho=phospho,
            total=None,
            declared_input_scale_state=declared_state,
            declared_input_establishment_mode=IntensityScaleEstablishmentMode.DECLARED,
            input_declaration_source="tests.unit.test_differential_analysis",
        )
        .intensity_scale_state
    )


def _phospy_transformed_log2_dataset(phospho: pd.DataFrame):
    resolved = DatasetIntensityScaleResolver(
        transformer=Log2Transformer(pseudocount=1.0)
    ).run(
        phospho=phospho,
        total=None,
    )
    return _dataset(
        resolved.phospho,
        intensity_scale_state=resolved.intensity_scale_state,
    )


def _design_from_conditions(
    entries: tuple[tuple[str, str], ...],
) -> ExperimentalDesign:
    replicate_counts: defaultdict[str, int] = defaultdict(int)
    records: list[SampleDesignRecord] = []
    for sample_id, condition in entries:
        replicate_counts[condition] += 1
        records.append(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=condition,
                biological_replicate_id=(f"{condition}_r{replicate_counts[condition]}"),
            )
        )
    return ExperimentalDesign(samples=tuple(records))


def _design() -> ExperimentalDesign:
    return _design_from_conditions(
        (
            ("A_1", "A"),
            ("A_2", "A"),
            ("B_1", "B"),
            ("B_2", "B"),
            ("C_1", "C"),
            ("C_2", "C"),
        )
    )


def _contrasts() -> tuple[Contrast, ...]:
    return (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
        Contrast(
            name="C_vs_A",
            numerator_condition="C",
            denominator_condition="A",
        ),
    )


def _request(
    *,
    dataset=None,
    design: ExperimentalDesign | None = None,
    contrasts: tuple[Contrast, ...] | None = None,
    empirical_bayes: EmpiricalBayesConfig | None = None,
    multiple_testing: MultipleTestingConfig | None = None,
    minimum_condition_replicates: int = 2,
    allow_suspicious_declared_input_scale: bool = False,
) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=_dataset() if dataset is None else dataset,
        design=_design() if design is None else design,
        contrasts=_contrasts() if contrasts is None else contrasts,
        config=DifferentialAnalysisConfig(
            minimum_condition_replicates=minimum_condition_replicates,
            allow_suspicious_declared_input_scale=(
                allow_suspicious_declared_input_scale
            ),
            empirical_bayes=(
                EmpiricalBayesConfig(method="standard")
                if empirical_bayes is None
                else empirical_bayes
            ),
            multiple_testing=(
                MultipleTestingConfig()
                if multiple_testing is None
                else multiple_testing
            ),
        ),
    )


def test_differential_analysis_returns_per_contrast_moderated_tables() -> None:
    result = DifferentialAnalysisWorkflow().run(_request())

    assert set(result.contrast_tables) == {"B_vs_A", "C_vs_A"}
    assert result.empirical_bayes_method == "standard"
    assert result.empirical_bayes_robust is False
    assert result.empirical_bayes_trend is False
    assert result.mean_variance_trend_diagnostics is None
    assert result.policy_provenance is not None
    assert result.policy_provenance.design.formula == "~0 + condition"
    assert result.policy_provenance.design.rank == 3
    assert result.policy_provenance.design.residual_degrees_of_freedom == pytest.approx(
        3.0
    )
    assert result.policy_provenance.statistical_testing.p_value_method == (
        "two_sided_t_distribution_survival_function"
    )
    assert result.policy_provenance.statistical_testing.adjusted_p_value_method == (
        "benjamini_hochberg"
    )
    assert result.policy_provenance.statistical_testing.input_intensity_scale == "log2"
    assert result.policy_provenance.statistical_testing.logfc_interpretation == (
        "fitted condition contrast on the established log2 phosphosite intensity scale"
    )
    assert result.policy_provenance.missing_values.policy == (
        "reject_missing_values_before_differential_execution"
    )
    assert (
        "correlated repeated-measure differential modelling beyond explicit fixed blocks"
        in result.policy_provenance.unsupported_design.intentionally_rejected_features
    )
    assert result.policy_provenance.replicates.condition_replicate_counts == (
        ("A", 2),
        ("B", 2),
        ("C", 2),
    )
    assert [item.name for item in result.policy_provenance.contrasts] == [
        "B_vs_A",
        "C_vs_A",
    ]
    for contrast_name in ("B_vs_A", "C_vs_A"):
        table = result.table_for(contrast_name)
        assert list(table.columns) == [
            "site_key",
            "display_id",
            "organism",
            "protein_namespace",
            "protein_identifier",
            "gene_symbol",
            "site",
            "protein_id",
            "logFC",
            "t",
            "P.Value",
            "adj.P.Val",
        ]
        assert table.shape[0] == 5
        assert table.index.tolist() == _dataset().phospho.index.tolist()
        assert (table.loc[:, "P.Value"] >= 0.0).all()
        assert (table.loc[:, "P.Value"] <= 1.0).all()
        assert (table.loc[:, "adj.P.Val"] >= 0.0).all()
        assert (table.loc[:, "adj.P.Val"] <= 1.0).all()


def test_differential_rejects_suspicious_declared_log2_by_default() -> None:
    suspicious_matrix = _matrix() * 10000.0
    state = _declared_log2_intensity_scale_state(suspicious_matrix)
    provenance = state.establishment_provenance
    assert provenance is not None
    assert provenance.mode is IntensityScaleEstablishmentMode.DECLARED
    assert provenance.diagnostic_warnings

    with pytest.raises(WorkflowValidationError) as exc_info:
        DifferentialAnalysisWorkflow().run(
            _request(
                dataset=_dataset(
                    suspicious_matrix,
                    intensity_scale_state=state,
                )
            )
        )

    message = str(exc_info.value)
    assert (
        "differential analysis rejects suspicious declared log2 intensity scale "
        "by default"
    ) in message
    assert provenance.diagnostic_warnings[0] in message
    assert "rebuild dataset with correct input scale" in message
    assert "apply supported log2 transformation" in message
    assert "explicitly set differential override" in message


def test_differential_accepts_declared_log2_without_warnings() -> None:
    state = _declared_log2_intensity_scale_state(_matrix())
    provenance = state.establishment_provenance
    assert provenance is not None
    assert provenance.mode is IntensityScaleEstablishmentMode.DECLARED
    assert provenance.diagnostic_warnings == ()

    result = DifferentialAnalysisWorkflow().run(
        _request(dataset=_dataset(_matrix(), intensity_scale_state=state))
    )

    assert result.policy_provenance is not None
    testing_policy = result.policy_provenance.statistical_testing
    assert testing_policy.allow_suspicious_declared_input_scale is False


def test_differential_accepts_phospy_transformed_log2_state() -> None:
    dataset = _phospy_transformed_log2_dataset(_matrix() + 10.0)
    provenance = dataset.intensity_scale_state.establishment_provenance
    assert provenance is not None
    assert provenance.source is IntensityScaleEstablishmentSource.TRANSFORMED_BY_PHOSPY

    result = DifferentialAnalysisWorkflow().run(_request(dataset=dataset))

    assert result.policy_provenance is not None
    assert result.policy_provenance.statistical_testing.input_intensity_scale == "log2"


def test_differential_accepts_suspicious_declared_log2_with_explicit_override() -> None:
    suspicious_matrix = _matrix() * 10000.0
    state = _declared_log2_intensity_scale_state(suspicious_matrix)

    result = DifferentialAnalysisWorkflow().run(
        _request(
            dataset=_dataset(
                suspicious_matrix,
                intensity_scale_state=state,
            ),
            allow_suspicious_declared_input_scale=True,
        )
    )

    assert result.policy_provenance is not None
    testing_policy = result.policy_provenance.statistical_testing
    assert testing_policy.allow_suspicious_declared_input_scale is True


def test_differential_default_multiple_testing_correction_is_unchanged() -> None:
    default = DifferentialAnalysisWorkflow().run(_request())
    explicit = DifferentialAnalysisWorkflow().run(
        _request(
            multiple_testing=MultipleTestingConfig(method="benjamini_hochberg"),
        )
    )

    assert default.policy_provenance is not None
    assert (
        default.policy_provenance.statistical_testing.adjusted_p_value_method
        == "benjamini_hochberg"
    )
    for contrast_name in ("B_vs_A", "C_vs_A"):
        pdt.assert_frame_equal(
            default.table_for(contrast_name),
            explicit.table_for(contrast_name),
            check_exact=False,
            rtol=1e-12,
            atol=0.0,
        )


@pytest.mark.parametrize("method", SUPPORTED_MULTIPLE_TESTING_METHODS)
def test_differential_analysis_supports_configured_multiple_testing_methods(
    method: MultipleTestingMethod,
) -> None:
    result = DifferentialAnalysisWorkflow().run(
        _request(multiple_testing=MultipleTestingConfig(method=method))
    )

    assert result.policy_provenance is not None
    assert (
        result.policy_provenance.statistical_testing.adjusted_p_value_method == method
    )
    assert result.policy_provenance.missing_values.adjusted_p_value_scope == (
        "adjustment_over_tested_features_only_per_contrast"
    )
    for contrast_name in ("B_vs_A", "C_vs_A"):
        table = result.table_for(contrast_name)
        expected = adjust_p_values(
            table.loc[:, "P.Value"].to_numpy(dtype=float),
            method=method,
        )
        np.testing.assert_allclose(
            table.loc[:, "adj.P.Val"].to_numpy(dtype=float),
            expected,
            rtol=1e-12,
            atol=0.0,
        )


def test_differential_correction_is_applied_per_contrast() -> None:
    result = DifferentialAnalysisWorkflow().run(
        _request(multiple_testing=MultipleTestingConfig(method="bonferroni"))
    )

    b_vs_a = result.table_for("B_vs_A")
    c_vs_a = result.table_for("C_vs_A")
    pooled_adjusted = adjust_p_values(
        np.concatenate(
            [
                b_vs_a.loc[:, "P.Value"].to_numpy(dtype=float),
                c_vs_a.loc[:, "P.Value"].to_numpy(dtype=float),
            ]
        ),
        method="bonferroni",
    )

    np.testing.assert_allclose(
        b_vs_a.loc[:, "adj.P.Val"].to_numpy(dtype=float),
        adjust_p_values(
            b_vs_a.loc[:, "P.Value"].to_numpy(dtype=float),
            method="bonferroni",
        ),
        rtol=1e-12,
        atol=0.0,
    )
    np.testing.assert_allclose(
        c_vs_a.loc[:, "adj.P.Val"].to_numpy(dtype=float),
        adjust_p_values(
            c_vs_a.loc[:, "P.Value"].to_numpy(dtype=float),
            method="bonferroni",
        ),
        rtol=1e-12,
        atol=0.0,
    )
    assert not np.allclose(
        b_vs_a.loc[:, "adj.P.Val"].to_numpy(dtype=float),
        pooled_adjusted[: b_vs_a.shape[0]],
        rtol=1e-12,
        atol=0.0,
    )


def test_differential_validator_rejects_invalid_multiple_testing_method() -> None:
    config = MultipleTestingConfig()
    object.__setattr__(config, "method", "storey")

    with pytest.raises(
        WorkflowValidationError,
        match="multiple_testing.method must be one of",
    ):
        DifferentialAnalysisValidator().run(_request(multiple_testing=config))


def test_differential_interpreter_condition_only_design_inputs_remain_unchanged() -> (
    None
):
    interpreted = DifferentialAnalysisInterpreter().run(
        DifferentialAnalysisValidator().run(_request())
    )

    expected_design = pd.DataFrame(
        {
            "A": [1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            "B": [0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
            "C": [0.0, 0.0, 0.0, 0.0, 1.0, 1.0],
        },
        index=pd.Index(["A_1", "A_2", "B_1", "B_2", "C_1", "C_2"], name="sample"),
    )
    expected_design.columns = pd.Index(expected_design.columns, name="coefficient")
    expected_contrasts = pd.DataFrame(
        {
            "B_vs_A": [-1.0, 1.0, 0.0],
            "C_vs_A": [-1.0, 0.0, 1.0],
        },
        index=pd.Index(["A", "B", "C"], name="coefficient"),
    )
    expected_contrasts.columns = pd.Index(expected_contrasts.columns, name="contrast")

    pdt.assert_frame_equal(
        interpreted.computation_request.design.to_dataframe(),
        expected_design,
    )
    pdt.assert_frame_equal(
        interpreted.computation_request.contrasts.to_dataframe(),
        expected_contrasts,
    )
    execution_design = interpreted.execution_design
    assert execution_design is not None
    pdt.assert_frame_equal(
        execution_design.design_matrix.to_dataframe(),
        expected_design,
    )
    pdt.assert_frame_equal(
        execution_design.contrast_matrix.to_dataframe(),
        expected_contrasts,
    )
    assert execution_design.condition_contrast_vectors[0].coefficients == (
        ("A", -1.0),
        ("B", 1.0),
        ("C", 0.0),
    )
    assert execution_design.covariate_columns == ()
    assert execution_design.formula == "~0 + condition"
    assert execution_design.description == "condition-only fixed-effect design"
    assert execution_design.sample_order == ("A_1", "A_2", "B_1", "B_2", "C_1", "C_2")
    assert execution_design.paired_design_policy == "reject"
    assert execution_design.block_column_metadata is None


def test_differential_interpreter_builds_categorical_covariate_inputs() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_r2",
                covariates={"sex": "M"},
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_r2",
                covariates={"sex": "M"},
            ),
            SampleDesignRecord(
                sample_id="C_1",
                condition="C",
                biological_replicate_id="C_r1",
                covariates={"sex": "F"},
            ),
            SampleDesignRecord(
                sample_id="C_2",
                condition="C",
                biological_replicate_id="C_r2",
                covariates={"sex": "M"},
            ),
        ),
        fixed_effects=(CategoricalCovariate("sex"),),
    )

    interpreted = DifferentialAnalysisInterpreter().run(
        DifferentialAnalysisValidator().run(_request(design=design))
    )

    expected_design = pd.DataFrame(
        {
            "A": [1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            "B": [0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
            "C": [0.0, 0.0, 0.0, 0.0, 1.0, 1.0],
            "sex[M]": [0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
        },
        index=pd.Index(["A_1", "A_2", "B_1", "B_2", "C_1", "C_2"], name="sample"),
    )
    expected_design.columns = pd.Index(expected_design.columns, name="coefficient")
    expected_contrasts = pd.DataFrame(
        {
            "B_vs_A": [-1.0, 1.0, 0.0, 0.0],
            "C_vs_A": [-1.0, 0.0, 1.0, 0.0],
        },
        index=pd.Index(["A", "B", "C", "sex[M]"], name="coefficient"),
    )
    expected_contrasts.columns = pd.Index(expected_contrasts.columns, name="contrast")
    execution_design = interpreted.execution_design
    assert execution_design is not None

    pdt.assert_frame_equal(
        interpreted.computation_request.design.to_dataframe(),
        expected_design,
    )
    pdt.assert_frame_equal(
        execution_design.design_matrix.to_dataframe(),
        expected_design,
    )
    pdt.assert_frame_equal(
        execution_design.contrast_matrix.to_dataframe(),
        expected_contrasts,
    )
    assert execution_design.condition_contrast_vectors[0].coefficients == (
        ("A", -1.0),
        ("B", 1.0),
        ("C", 0.0),
        ("sex[M]", 0.0),
    )
    assert execution_design.condition_contrast_vectors[1].coefficients == (
        ("A", -1.0),
        ("B", 0.0),
        ("C", 1.0),
        ("sex[M]", 0.0),
    )
    assert len(execution_design.covariate_columns) == 1
    covariate = execution_design.covariate_columns[0]
    assert covariate.name == "sex"
    assert covariate.kind == "categorical"
    assert covariate.columns == ("sex[M]",)
    assert covariate.levels == ("F", "M")
    assert covariate.reference_level == "F"
    assert covariate.unused_levels == ()
    assert execution_design.formula == "~0 + condition + sex"
    assert execution_design.description == "fixed-effect design: ~0 + condition + sex"


def test_differential_interpreter_builds_continuous_covariate_inputs() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
                covariates={"dose": 0.0},
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_r2",
                covariates={"dose": 1.0},
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
                covariates={"dose": 0.0},
            ),
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_r2",
                covariates={"dose": 1.0},
            ),
            SampleDesignRecord(
                sample_id="C_1",
                condition="C",
                biological_replicate_id="C_r1",
                covariates={"dose": 0.0},
            ),
            SampleDesignRecord(
                sample_id="C_2",
                condition="C",
                biological_replicate_id="C_r2",
                covariates={"dose": 1.0},
            ),
        ),
        fixed_effects=(ContinuousCovariate("dose"),),
    )

    interpreted = DifferentialAnalysisInterpreter().run(
        DifferentialAnalysisValidator().run(_request(design=design))
    )

    expected_design = pd.DataFrame(
        {
            "A": [1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            "B": [0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
            "C": [0.0, 0.0, 0.0, 0.0, 1.0, 1.0],
            "dose": [0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
        },
        index=pd.Index(["A_1", "A_2", "B_1", "B_2", "C_1", "C_2"], name="sample"),
    )
    expected_design.columns = pd.Index(expected_design.columns, name="coefficient")
    execution_design = interpreted.execution_design
    assert execution_design is not None

    pdt.assert_frame_equal(
        interpreted.computation_request.design.to_dataframe(),
        expected_design,
    )
    pdt.assert_frame_equal(
        execution_design.design_matrix.to_dataframe(),
        expected_design,
    )
    assert execution_design.condition_contrast_vectors[0].coefficients == (
        ("A", -1.0),
        ("B", 1.0),
        ("C", 0.0),
        ("dose", 0.0),
    )
    assert len(execution_design.covariate_columns) == 1
    covariate = execution_design.covariate_columns[0]
    assert covariate.name == "dose"
    assert covariate.kind == "continuous"
    assert covariate.columns == ("dose",)
    assert covariate.levels == ()
    assert covariate.reference_level is None
    assert execution_design.formula == "~0 + condition + dose"


def test_differential_interpreter_preserves_execution_sample_order() -> None:
    sample_order = ("B_2", "A_1", "C_2", "B_1", "A_2", "C_1")
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="B_2",
                condition="B",
                biological_replicate_id="B_r2",
                covariates={"dose": 1.0},
            ),
            SampleDesignRecord(
                sample_id="A_1",
                condition="A",
                biological_replicate_id="A_r1",
                covariates={"dose": 0.0},
            ),
            SampleDesignRecord(
                sample_id="C_2",
                condition="C",
                biological_replicate_id="C_r2",
                covariates={"dose": 1.0},
            ),
            SampleDesignRecord(
                sample_id="B_1",
                condition="B",
                biological_replicate_id="B_r1",
                covariates={"dose": 0.0},
            ),
            SampleDesignRecord(
                sample_id="A_2",
                condition="A",
                biological_replicate_id="A_r2",
                covariates={"dose": 1.0},
            ),
            SampleDesignRecord(
                sample_id="C_1",
                condition="C",
                biological_replicate_id="C_r1",
                covariates={"dose": 0.0},
            ),
        ),
        fixed_effects=(ContinuousCovariate("dose"),),
    )
    matrix = _matrix().loc[:, ["C_1", "B_1", "A_2", "C_2", "B_2", "A_1"]]

    interpreted = DifferentialAnalysisInterpreter().run(
        DifferentialAnalysisValidator().run(
            _request(dataset=_dataset(matrix), design=design)
        )
    )

    execution_design = interpreted.execution_design
    assert execution_design is not None
    assert execution_design.sample_order == sample_order
    assert interpreted.computation_request.matrix.columns.tolist() == list(sample_order)
    assert execution_design.design_matrix.to_dataframe().index.tolist() == list(
        sample_order
    )


def test_empirical_bayes_config_rejects_invalid_winsor_tail_values() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="empirical_bayes.winsor_tail_p values must sum to less than 1.0",
    ):
        EmpiricalBayesConfig(method="robust", winsor_tail_p=(0.5, 0.5))


def test_robust_mode_downweights_variance_outlier() -> None:
    matrix = _matrix().copy()
    matrix.loc["MAPK14;Y182;", "C_1"] = 9.5
    matrix.loc["MAPK14;Y182;", "C_2"] = -6.0
    standard = DifferentialAnalysisWorkflow().run(
        _request(
            dataset=_dataset(matrix),
            empirical_bayes=EmpiricalBayesConfig(method="standard"),
        )
    )
    robust = DifferentialAnalysisWorkflow().run(
        _request(
            dataset=_dataset(matrix),
            empirical_bayes=EmpiricalBayesConfig(
                method="robust",
                winsor_tail_p=(0.05, 0.1),
            ),
        )
    )

    assert robust.empirical_bayes_robust is True
    assert robust.policy_provenance is not None
    assert robust.policy_provenance.empirical_bayes.method == "robust"
    assert robust.policy_provenance.empirical_bayes.robust is True
    assert robust.policy_provenance.empirical_bayes.winsor_tail_p == (0.05, 0.1)
    assert robust.prior_diagnostics.robust_outlier_count >= 1
    outlier_site = _site_key_for_display_id("MAPK14;Y182;")
    assert (
        robust.prior_degrees_of_freedom_series().loc[outlier_site]
        <= standard.prior_degrees_of_freedom_series().loc[outlier_site]
    )


def test_trend_mode_stores_mean_variance_diagnostics() -> None:
    matrix = _matrix().copy()
    matrix.loc["MAPK14;Y182;"] = matrix.loc["MAPK14;Y182;"] * 0.1
    matrix.loc["MTOR;S2448;"] = matrix.loc["MTOR;S2448;"] * 4.0
    result = DifferentialAnalysisWorkflow().run(
        _request(
            dataset=_dataset(matrix),
            empirical_bayes=EmpiricalBayesConfig(
                method="standard",
                trend=True,
            ),
        )
    )
    assert result.empirical_bayes_trend is True
    assert result.mean_variance_trend_diagnostics is not None
    diagnostics = result.mean_variance_trend_diagnostics
    expected_index = _dataset(matrix).phospho.index.tolist()
    assert diagnostics.mean_intensity.index.tolist() == expected_index
    assert diagnostics.fitted_log_prior_variance.index.tolist() == expected_index
    assert not diagnostics.fitted_log_prior_variance.equals(
        diagnostics.log_residual_variance
    )


def test_low_replicate_mode_remains_stable_with_robust_and_trend() -> None:
    matrix = pd.DataFrame(
        {
            "A_1": [1.0, 2.0, 0.8, 0.5],
            "A_2": [1.1, 2.3, 1.0, 0.4],
            "B_1": [1.8, 2.2, 0.7, 2.1],
            "B_2": [1.9, 2.5, 0.9, 2.2],
        },
        index=pd.Index(
            ["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;", "RPS6KB1;T389;"],
            name="site_id",
        ),
    )
    design = _design_from_conditions(
        (
            ("A_1", "A"),
            ("A_2", "A"),
            ("B_1", "B"),
            ("B_2", "B"),
        )
    )
    contrasts = (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
    )

    result = DifferentialAnalysisWorkflow().run(
        _request(
            dataset=_dataset(matrix),
            design=design,
            contrasts=contrasts,
            empirical_bayes=EmpiricalBayesConfig(
                method="robust",
                trend=True,
                winsor_tail_p=(0.05, 0.1),
            ),
        )
    )
    table = result.table_for("B_vs_A")
    assert np.isfinite(table.loc[:, "t"]).all()
    assert np.isfinite(table.loc[:, "P.Value"]).all()
    assert (table.loc[:, "P.Value"] >= 0.0).all()
    assert (table.loc[:, "P.Value"] <= 1.0).all()


def test_differential_analysis_fails_on_sample_design_mismatch() -> None:
    design = _design_from_conditions(
        (
            ("X1", "A"),
            ("X2", "A"),
            ("X3", "B"),
            ("X4", "B"),
            ("X5", "C"),
            ("X6", "C"),
        )
    )
    with pytest.raises(
        WorkflowValidationError,
        match="samples not present in dataset",
    ):
        DifferentialAnalysisWorkflow().run(_request(design=design))


def test_differential_analysis_fails_on_contrast_design_term_mismatch() -> None:
    contrasts = (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A_wrong",
        ),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="unknown denominator condition",
    ):
        DifferentialAnalysisWorkflow().run(_request(contrasts=contrasts))


def test_differential_analysis_fails_when_residual_dof_is_non_positive() -> None:
    matrix = _matrix().loc[:, ["A_1", "B_1", "C_1"]].copy(deep=True)
    matrix.iloc[1, 1] += 0.05
    design = _design_from_conditions(
        (
            ("A_1", "A"),
            ("B_1", "B"),
            ("C_1", "C"),
        )
    )
    with pytest.raises(
        WorkflowValidationError,
        match="residual degrees of freedom must be positive",
    ):
        DifferentialAnalysisWorkflow().run(
            _request(
                dataset=_dataset(matrix),
                design=design,
                minimum_condition_replicates=1,
            )
        )


def test_differential_analysis_withholds_all_constant_site_intensities() -> None:
    matrix = pd.DataFrame(
        {
            "A_1": [5.0, 1.0, 2.0],
            "A_2": [5.0, 1.2, 2.1],
            "B_1": [5.0, 2.1, 2.3],
            "B_2": [5.0, 2.2, 2.4],
        },
        index=pd.Index(
            ["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"],
            name="site_id",
        ),
    )
    design = _design_from_conditions(
        (
            ("A_1", "A"),
            ("A_2", "A"),
            ("B_1", "B"),
            ("B_2", "B"),
        )
    )
    contrasts = (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
    )

    result = DifferentialAnalysisWorkflow().run(
        _request(
            dataset=_dataset(matrix),
            design=design,
            contrasts=contrasts,
        )
    )
    table = result.table_for("B_vs_A")

    assert table.iloc[0][DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT
    )
    assert table.iloc[0][["logFC", "t", "P.Value", "adj.P.Val"]].isna().all()
    assert table.iloc[1:][DIFFERENTIAL_RESULT_STATUS_COLUMN].tolist() == [
        DIFFERENTIAL_RESULT_STATUS_TESTED,
        DIFFERENTIAL_RESULT_STATUS_TESTED,
    ]
    assert (
        np.isfinite(table.iloc[1:][["logFC", "t", "P.Value", "adj.P.Val"]]).all().all()
    )


def test_differential_analysis_rejects_empty_condition_labels() -> None:
    with pytest.raises(
        ContractValidationError,
        match="condition",
    ):
        _request(
            design=ExperimentalDesign(
                samples=(
                    SampleDesignRecord(sample_id="A_1", condition="A"),
                    SampleDesignRecord(sample_id="A_2", condition=""),
                    SampleDesignRecord(sample_id="B_1", condition="B"),
                    SampleDesignRecord(sample_id="B_2", condition="B"),
                )
            )
        )


def test_differential_analysis_sample_order_mismatch_is_resolved_by_label() -> None:
    matrix = pd.DataFrame(
        {
            "A_1": [1.0, 2.0],
            "A_2": [1.1, 2.1],
            "B_1": [2.0, 1.9],
            "B_2": [2.2, 2.2],
        },
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    design = _design_from_conditions(
        (
            ("A_1", "A"),
            ("A_2", "A"),
            ("B_1", "B"),
            ("B_2", "B"),
        )
    )
    contrasts = (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
    )
    reordered_samples = ["B_2", "A_1", "B_1", "A_2"]
    reordered_design = _design_from_conditions(
        (
            ("B_2", "B"),
            ("A_1", "A"),
            ("B_1", "B"),
            ("A_2", "A"),
        )
    )

    aligned = (
        DifferentialAnalysisWorkflow()
        .run(
            _request(
                dataset=_dataset(matrix),
                design=design,
                contrasts=contrasts,
            )
        )
        .table_for("B_vs_A")
    )
    reordered = (
        DifferentialAnalysisWorkflow()
        .run(
            _request(
                dataset=_dataset(matrix.loc[:, reordered_samples]),
                design=reordered_design,
                contrasts=contrasts,
            )
        )
        .table_for("B_vs_A")
    )

    pdt.assert_frame_equal(
        aligned,
        reordered,
        check_exact=False,
        rtol=1e-12,
        atol=1e-12,
    )
