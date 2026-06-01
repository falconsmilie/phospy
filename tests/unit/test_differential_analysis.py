from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api import (
    Contrast,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    EmpiricalBayesConfig,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.errors import (
    PhosPyInputError,
    WorkflowBoundaryError,
    WorkflowValidationError,
)
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    encode_site_key,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)


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


def _dataset(matrix: pd.DataFrame | None = None):
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
    return supported_dataset(phospho=phospho, site_metadata=site_metadata)


def supported_dataset(*, phospho: pd.DataFrame, site_metadata: pd.DataFrame):
    from phospy import AnalysisReadyPhosphoDataset

    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
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
    minimum_condition_replicates: int = 2,
) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=_dataset() if dataset is None else dataset,
        design=_design() if design is None else design,
        contrasts=_contrasts() if contrasts is None else contrasts,
        config=DifferentialAnalysisConfig(
            minimum_condition_replicates=minimum_condition_replicates,
            empirical_bayes=(
                EmpiricalBayesConfig(method="standard")
                if empirical_bayes is None
                else empirical_bayes
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
    assert result.policy_provenance.missing_values.policy == (
        "reject_missing_values_before_differential_execution"
    )
    assert "batch-aware differential modelling" in (
        result.policy_provenance.unsupported_design.intentionally_rejected_features
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
    matrix = _matrix().loc[:, ["A_1", "B_1", "C_1"]]
    design = _design_from_conditions(
        (
            ("A_1", "A"),
            ("B_1", "B"),
            ("C_1", "C"),
        )
    )
    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.interpreter.residual_dof",
    ):
        DifferentialAnalysisWorkflow().run(
            _request(
                dataset=_dataset(matrix),
                design=design,
                minimum_condition_replicates=1,
            )
        )


def test_differential_analysis_handles_zero_variance_features() -> None:
    matrix = pd.DataFrame(
        {
            "A_1": [5.0, 1.0],
            "A_2": [5.0, 1.2],
            "B_1": [5.0, 2.1],
            "B_2": [5.0, 2.2],
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

    result = DifferentialAnalysisWorkflow().run(
        _request(
            dataset=_dataset(matrix),
            design=design,
            contrasts=contrasts,
        )
    )
    table = result.table_for("B_vs_A")

    site_key = _site_key_for_display_id("MAPK14;Y182;")
    assert table.at[site_key, "logFC"] == pytest.approx(0.0)
    assert table.at[site_key, "t"] == pytest.approx(0.0)
    assert table.at[site_key, "P.Value"] == pytest.approx(1.0)
    assert np.isfinite(table.loc[:, "t"]).all()
    assert np.isfinite(table.loc[:, "P.Value"]).all()


def test_differential_analysis_rejects_empty_condition_labels() -> None:
    with pytest.raises(
        WorkflowValidationError,
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

    pdt.assert_frame_equal(aligned, reordered, check_exact=False, rtol=1e-12, atol=0.0)
