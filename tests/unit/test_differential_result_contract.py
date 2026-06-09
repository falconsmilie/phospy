from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
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
from phospy.errors import DatasetValidationError, PhosPyInputError
from phospy.science.differential.executor import (
    DifferentialAnalysisExecutor as DifferentialComputationExecutor,
)
from phospy.science.differential.models import (
    DifferentialAnalysisRequest as DifferentialComputationRequest,
)
from phospy.science.differential.models import (
    DifferentialAnalysisResult,
    DifferentialComputationResult,
    EmpiricalBayesPriorDiagnostics,
)
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    encode_site_key,
)
from phospy.workflows.differential.executor import (
    DifferentialAnalysisExecutor as DifferentialWorkflowExecutor,
)
from phospy.workflows.differential.interpreter import DifferentialAnalysisInterpreter
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import site_key_context_columns

FIXTURE_DIR = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "rewrite_parity"
    / "differential_limma_envelope"
)
NEGATIVE_FIXTURE_DIR = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "rewrite_parity"
    / "differential_contract_negative_cases"
)
IDENTITY_COLUMNS = [
    "site_key",
    "display_id",
    "gene_symbol",
    "site",
    "organism",
    "protein_namespace",
    "protein_identifier",
    "protein_id",
]
STATISTIC_COLUMNS = ["logFC", "t", "P.Value", "adj.P.Val"]


def _canonical_site_id(raw_site_id: str, *, ordinal: int) -> str:
    tokens = [token.strip() for token in raw_site_id.split(";") if token.strip()]
    if len(tokens) >= 2:
        return f"{tokens[0]};{tokens[1]};"
    return f"SITE{ordinal};S{ordinal};"


def _site_key(
    *,
    gene_symbol: str,
    site: str,
    protein_identifier: str | None = None,
    organism: str = "rat",
    protein_namespace: str = "protein_id",
) -> str:
    residue = site.strip().upper()[0]
    position = int(site.strip()[1:])
    key = build_protein_scoped_site_key(
        organism=organism,
        protein_namespace=protein_namespace,
        protein_identifier=protein_identifier or gene_symbol,
        residue=residue,
        position=position,
        field_name="tests.unit.test_differential_result_contract.site_key",
        error_type=ValueError,
    )
    return encode_site_key(key)


def _build_dataset(matrix: pd.DataFrame):
    from phospy import AnalysisReadyPhosphoDataset

    display_ids = matrix.index.astype(str).tolist()
    parsed = [site_id.split(";") for site_id in display_ids]
    site_keys = [_site_key(gene_symbol=parts[0], site=parts[1]) for parts in parsed]
    phospho = matrix.copy(deep=True)
    phospho.index = pd.Index(site_keys, name="site_key")
    site_metadata = pd.DataFrame(
        {
            "site_key": site_keys,
            "display_id": display_ids,
            **site_key_context_columns(site_keys),
            "gene_symbol": [parts[0] for parts in parsed],
            "site": [parts[1] for parts in parsed],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in [parts[1] for parts in parsed]
            ],
            "protein_id": [parts[0] for parts in parsed],
        },
        index=pd.Index(site_keys, name="site_key"),
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


def _load_matrix() -> pd.DataFrame:
    frame = pd.read_csv(FIXTURE_DIR / "matrix.csv")
    matrix = frame.set_index("site_id")
    raw_site_ids = matrix.index.astype(str).tolist()
    canonical_ids = [
        _canonical_site_id(site_id, ordinal=idx)
        for idx, site_id in enumerate(raw_site_ids, start=1)
    ]
    matrix.index = pd.Index(canonical_ids, name="display_id")
    return matrix


def _request_for_reverse_contrasts(
    matrix: pd.DataFrame,
) -> DifferentialAnalysisRequest:
    design = ExperimentalDesign(
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
    contrasts = (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
        Contrast(
            name="A_vs_B",
            numerator_condition="A",
            denominator_condition="B",
        ),
    )
    return DifferentialAnalysisRequest(
        dataset=_build_dataset(matrix),
        design=design,
        contrasts=contrasts,
        config=DifferentialAnalysisConfig(
            empirical_bayes=EmpiricalBayesConfig(method="standard")
        ),
    )


def _prior_diagnostics(index: pd.Index) -> EmpiricalBayesPriorDiagnostics:
    return EmpiricalBayesPriorDiagnostics(
        method="standard",
        robust=False,
        trend=False,
        winsor_tail_p=(0.05, 0.1),
        base_prior_variance=1.0,
        base_prior_degrees_of_freedom=10.0,
        robust_outlier_count=0,
        robust_outlier_fraction=0.0,
        winsorized_low_count=0,
        winsorized_high_count=0,
        prior_variance=pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="prior_residual_variance",
        ),
        prior_degrees_of_freedom=pd.Series(
            np.full(index.size, 10.0),
            index=index.copy(),
            name="prior_degrees_of_freedom",
        ),
    )


def _manual_result_with_table(table: pd.DataFrame) -> DifferentialAnalysisResult:
    index = table.index.copy()
    payload = {
        "residual_variance": pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="residual_variance",
        ),
        "posterior_residual_variance": pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="posterior_residual_variance",
        ),
        "prior_residual_variance": pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="prior_residual_variance",
        ),
        "prior_degrees_of_freedom_series_value": pd.Series(
            np.full(index.size, 10.0),
            index=index.copy(),
            name="prior_degrees_of_freedom",
        ),
        "prior_variance": 1.0,
        "prior_degrees_of_freedom": 10.0,
        "residual_degrees_of_freedom": 4.0,
        "empirical_bayes_method": "standard",
        "empirical_bayes_robust": False,
        "empirical_bayes_trend": False,
        "prior_diagnostics": _prior_diagnostics(index),
        "mean_variance_trend_diagnostics": None,
        "contrast_tables": {"B_vs_A": table},
    }
    return DifferentialAnalysisResult(**payload)


def _manual_computation_result_with_table(
    table: pd.DataFrame,
) -> DifferentialComputationResult:
    index = table.index.copy()
    return DifferentialComputationResult(
        residual_variance=pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="residual_variance",
        ),
        posterior_residual_variance=pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="posterior_residual_variance",
        ),
        prior_residual_variance=pd.Series(
            np.full(index.size, 1.0),
            index=index.copy(),
            name="prior_residual_variance",
        ),
        prior_degrees_of_freedom_series_value=pd.Series(
            np.full(index.size, 10.0),
            index=index.copy(),
            name="prior_degrees_of_freedom",
        ),
        prior_variance=1.0,
        prior_degrees_of_freedom=10.0,
        residual_degrees_of_freedom=4.0,
        empirical_bayes_method="standard",
        empirical_bayes_robust=False,
        empirical_bayes_trend=False,
        prior_diagnostics=_prior_diagnostics(index),
        mean_variance_trend_diagnostics=None,
        contrast_tables={"B_vs_A": table},
    )


def _strict_result_table() -> pd.DataFrame:
    site_keys = [
        _site_key(gene_symbol="MAPK14", site="Y182"),
        _site_key(gene_symbol="GSK3B", site="S9"),
    ]
    context = site_key_context_columns(site_keys)
    return pd.DataFrame(
        {
            "site_key": site_keys,
            "display_id": ["MAPK14;Y182;", "GSK3B;S9;"],
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "organism": context["organism"],
            "protein_namespace": context["protein_namespace"],
            "protein_identifier": context["protein_identifier"],
            "protein_id": ["MAPK14", "GSK3B"],
            "logFC": [1.0, -1.0],
            "t": [2.0, -2.0],
            "P.Value": [0.05, 0.10],
            "adj.P.Val": [0.10, 0.10],
        },
        index=pd.Index(site_keys, name="site_key"),
    )


def _akt1_t309_result_table() -> pd.DataFrame:
    site_keys = [
        _site_key(
            gene_symbol="AKT1",
            site="T309",
            protein_identifier="P31749",
            organism="human",
            protein_namespace="uniprot",
        )
    ]
    return pd.DataFrame(
        {
            "site_key": site_keys,
            "display_id": ["AKT1;T309;"],
            "gene_symbol": ["AKT1"],
            "site": ["T309"],
            "organism": ["human"],
            "protein_namespace": ["uniprot"],
            "protein_identifier": ["P31749"],
            "protein_id": ["P31749"],
            "logFC": [1.0],
            "t": [2.0],
            "P.Value": [0.05],
            "adj.P.Val": [0.10],
        },
        index=pd.Index(site_keys, name="site_key"),
    )


def _stat_only_display_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "logFC": [1.0, -1.0],
            "t": [2.0, -2.0],
            "P.Value": [0.05, 0.10],
            "adj.P.Val": [0.10, 0.10],
        },
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="display_id"),
    )


def _stat_only_site_key_table() -> pd.DataFrame:
    strict_table = _strict_result_table()
    return strict_table.loc[:, STATISTIC_COLUMNS].copy()


def test_direct_result_construction_accepts_site_key_identity_table() -> None:
    result = _manual_result_with_table(_strict_result_table())
    table = result.table_for("B_vs_A")

    assert list(table.columns) == IDENTITY_COLUMNS + STATISTIC_COLUMNS
    assert table.index.name == "site_key"
    assert table.loc[:, "site_key"].tolist() == table.index.tolist()


def test_direct_result_construction_rejects_display_indexed_stat_only_table() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="missing required columns: site_key, display_id, gene_symbol, site",
    ):
        _manual_result_with_table(_stat_only_display_table())


def test_direct_result_construction_rejects_site_key_indexed_stat_only_table() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="missing required columns: site_key, display_id, gene_symbol, site",
    ):
        _manual_result_with_table(_stat_only_site_key_table())


def test_direct_result_construction_rejects_missing_site_key_column() -> None:
    table = _strict_result_table().drop(columns=["site_key"])

    with pytest.raises(PhosPyInputError, match="missing required columns: site_key"):
        _manual_result_with_table(table)


def test_direct_result_construction_rejects_missing_display_id_column() -> None:
    table = _strict_result_table().drop(columns=["display_id"])

    with pytest.raises(PhosPyInputError, match="missing required columns: display_id"):
        _manual_result_with_table(table)


def test_direct_result_construction_rejects_missing_gene_symbol_column() -> None:
    table = _strict_result_table().drop(columns=["gene_symbol"])

    with pytest.raises(PhosPyInputError, match="missing required columns: gene_symbol"):
        _manual_result_with_table(table)


def test_direct_result_construction_rejects_missing_site_column() -> None:
    table = _strict_result_table().drop(columns=["site"])

    with pytest.raises(PhosPyInputError, match="missing required columns: site"):
        _manual_result_with_table(table)


def test_direct_result_construction_allows_duplicate_display_ids() -> None:
    site_keys = [
        _site_key(
            gene_symbol="MAPK14",
            site="Y182",
            protein_identifier="MAPK14_CANONICAL",
        ),
        _site_key(
            gene_symbol="MAPK14",
            site="Y182",
            protein_identifier="MAPK14_ISOFORM_2",
        ),
    ]
    context = site_key_context_columns(site_keys)
    table = pd.DataFrame(
        {
            "site_key": site_keys,
            "display_id": ["MAPK14;Y182;", "MAPK14;Y182;"],
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "organism": context["organism"],
            "protein_namespace": context["protein_namespace"],
            "protein_identifier": context["protein_identifier"],
            "protein_id": ["MAPK14_A", "MAPK14_B"],
            "logFC": [1.0, -1.0],
            "t": [2.0, -2.0],
            "P.Value": [0.05, 0.10],
            "adj.P.Val": [0.10, 0.10],
        },
        index=pd.Index(site_keys, name="site_key"),
    )

    result = _manual_result_with_table(table)
    exported = result.table_for("B_vs_A")

    assert exported.loc[:, "display_id"].tolist() == [
        "MAPK14;Y182;",
        "MAPK14;Y182;",
    ]
    assert exported.loc[:, "site_key"].nunique() == exported.shape[0]
    assert exported.index.tolist() == table.index.tolist()


def test_direct_result_construction_rejects_site_mismatch() -> None:
    table = _akt1_t309_result_table()
    table.loc[:, "site"] = "T308"
    table.loc[:, "display_id"] = "AKT1;T308;"

    with pytest.raises(
        PhosPyInputError,
        match="site_key encodes T309 but row metadata site is 'T308'",
    ):
        _manual_result_with_table(table)


def test_direct_result_construction_rejects_display_id_site_mismatch() -> None:
    table = _akt1_t309_result_table()
    table.loc[:, "display_id"] = "AKT1;T308;"

    with pytest.raises(
        PhosPyInputError,
        match=r"display_id does not match gene_symbol \+ site",
    ):
        _manual_result_with_table(table)


def test_direct_result_construction_rejects_organism_mismatched_to_site_key() -> None:
    table = _akt1_t309_result_table()
    table.loc[:, "organism"] = "mouse"

    with pytest.raises(PhosPyInputError, match="organism is incoherent"):
        _manual_result_with_table(table)


def test_direct_result_construction_rejects_protein_namespace_mismatch() -> None:
    table = _akt1_t309_result_table()
    table.loc[:, "protein_namespace"] = "refseq"

    with pytest.raises(PhosPyInputError, match="protein_namespace is incoherent"):
        _manual_result_with_table(table)


def test_direct_result_construction_rejects_protein_identifier_mismatch() -> None:
    table = _akt1_t309_result_table()
    table.loc[:, "protein_identifier"] = "P99999"

    with pytest.raises(PhosPyInputError, match="protein_identifier is incoherent"):
        _manual_result_with_table(table)


def test_direct_result_construction_rejects_display_id_gene_site_mismatch() -> None:
    table = _akt1_t309_result_table()
    table.loc[:, "display_id"] = "AKT2;T309;"

    with pytest.raises(
        PhosPyInputError,
        match=r"display_id does not match gene_symbol \+ site",
    ):
        _manual_result_with_table(table)


def test_direct_result_construction_rejects_display_labels_as_site_key_column() -> None:
    table = _strict_result_table()
    table.loc[:, "site_key"] = table.loc[:, "display_id"].astype(str).tolist()

    with pytest.raises(PhosPyInputError, match="site_key must exactly match"):
        _manual_result_with_table(table)


def test_direct_result_construction_rejects_arbitrary_site_key_values() -> None:
    table = _strict_result_table()
    arbitrary_index = pd.Index(["row-1", "row-2"], name="site_key")
    table.index = arbitrary_index
    table.loc[:, "site_key"] = arbitrary_index.tolist()

    with pytest.raises(PhosPyInputError, match="must start with 'phospy:v1'"):
        _manual_result_with_table(table)


def test_direct_result_construction_rejects_duplicate_site_key_rows() -> None:
    table = _strict_result_table().iloc[[0, 0], :].copy()

    with pytest.raises(PhosPyInputError, match="must be unique"):
        _manual_result_with_table(table)


def test_direct_result_construction_rejects_site_key_index_column_mismatch() -> None:
    table = _strict_result_table()
    site_keys = table.index.astype(str).tolist()
    table.loc[site_keys[0], "site_key"] = site_keys[1]

    with pytest.raises(PhosPyInputError, match="site_key must exactly match"):
        _manual_result_with_table(table)


def test_internal_stat_only_result_type_is_explicit() -> None:
    result = _manual_computation_result_with_table(_stat_only_display_table())
    table = result.table_for("B_vs_A")

    assert isinstance(result, DifferentialComputationResult)
    assert not isinstance(result, DifferentialAnalysisResult)
    assert list(table.columns) == STATISTIC_COLUMNS
    assert table.index.name == "display_id"


def test_low_level_differential_executor_returns_computation_result() -> None:
    interpreted = DifferentialAnalysisInterpreter().run(
        DifferentialAnalysisValidator().run(
            _request_for_reverse_contrasts(_load_matrix())
        )
    )
    result = DifferentialComputationExecutor().run(interpreted.computation_request)

    assert isinstance(result, DifferentialComputationResult)
    assert not isinstance(result, DifferentialAnalysisResult)
    for contrast_name in ("B_vs_A", "A_vs_B"):
        table = result.table_for(contrast_name)
        assert list(table.columns) == STATISTIC_COLUMNS


def test_result_tables_follow_public_differential_contract() -> None:
    result = DifferentialAnalysisWorkflow().run(
        _request_for_reverse_contrasts(_load_matrix())
    )

    assert isinstance(result, DifferentialAnalysisResult)
    assert list(result.contrast_tables) == ["B_vs_A", "A_vs_B"]
    for contrast_name in ("B_vs_A", "A_vs_B"):
        table = result.table_for(contrast_name)
        assert list(table.columns) == IDENTITY_COLUMNS + STATISTIC_COLUMNS
        assert table.index.name == "site_key"
        assert table.loc[:, "site_key"].tolist() == table.index.tolist()
        assert np.isfinite(table.loc[:, "logFC"]).all()
        assert np.isfinite(table.loc[:, "t"]).all()
        assert (table.loc[:, "P.Value"] >= 0.0).all()
        assert (table.loc[:, "P.Value"] <= 1.0).all()
        assert (table.loc[:, "adj.P.Val"] >= 0.0).all()
        assert (table.loc[:, "adj.P.Val"] <= 1.0).all()


def test_reverse_contrasts_are_directionally_consistent() -> None:
    result = DifferentialAnalysisWorkflow().run(
        _request_for_reverse_contrasts(_load_matrix())
    )
    b_vs_a = result.table_for("B_vs_A")
    a_vs_b = result.table_for("A_vs_B")

    np.testing.assert_allclose(
        b_vs_a.loc[:, "logFC"].to_numpy(dtype=float),
        -a_vs_b.loc[:, "logFC"].to_numpy(dtype=float),
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        b_vs_a.loc[:, "t"].to_numpy(dtype=float),
        -a_vs_b.loc[:, "t"].to_numpy(dtype=float),
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        b_vs_a.loc[:, "P.Value"].to_numpy(dtype=float),
        a_vs_b.loc[:, "P.Value"].to_numpy(dtype=float),
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        b_vs_a.loc[:, "adj.P.Val"].to_numpy(dtype=float),
        a_vs_b.loc[:, "adj.P.Val"].to_numpy(dtype=float),
        rtol=1e-12,
        atol=1e-12,
    )


def test_workflow_rejects_incoherent_duplicate_display_ids() -> None:
    request = _request_for_reverse_contrasts(_load_matrix())
    interpreted = DifferentialAnalysisInterpreter().run(
        DifferentialAnalysisValidator().run(request)
    )
    duplicate_display_id = interpreted.result_identity_metadata.loc[
        interpreted.result_identity_metadata.index[0], "display_id"
    ]
    mutated_identity_metadata = interpreted.result_identity_metadata.copy(deep=True)
    mutated_identity_metadata.loc[:, "display_id"] = duplicate_display_id
    interpreted_with_duplicate_display_ids = replace(
        interpreted,
        result_identity_metadata=mutated_identity_metadata,
    )

    with pytest.raises(
        PhosPyInputError,
        match=r"display_id does not match gene_symbol \+ site",
    ):
        DifferentialWorkflowExecutor().run(interpreted_with_duplicate_display_ids)


def test_workflow_keeps_duplicate_display_ids_with_distinct_site_keys() -> None:
    from phospy import AnalysisReadyPhosphoDataset

    site_keys = [
        _site_key(
            gene_symbol="MAPK14",
            site="Y182",
            protein_identifier="MAPK14_CANONICAL",
        ),
        _site_key(
            gene_symbol="MAPK14",
            site="Y182",
            protein_identifier="MAPK14_ISOFORM_2",
        ),
        _site_key(gene_symbol="GSK3B", site="S9"),
        _site_key(gene_symbol="AKT1", site="T308"),
    ]
    display_ids = ["MAPK14;Y182;", "MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"]
    gene_symbols = ["MAPK14", "MAPK14", "GSK3B", "AKT1"]
    sites = ["Y182", "Y182", "S9", "T308"]
    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 1.2, 2.0, 1.0],
            "A_2": [1.1, 1.3, 2.1, 1.1],
            "B_1": [2.1, 2.2, 2.0, 0.9],
            "B_2": [2.0, 2.1, 2.2, 1.0],
        },
        index=pd.Index(site_keys, name="site_key"),
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_keys,
            "display_id": display_ids,
            **site_key_context_columns(site_keys),
            "gene_symbol": gene_symbols,
            "site": sites,
            "site_sequence": [("A" * 15) + site[0] + ("A" * 15) for site in sites],
            "protein_id": ["MAPK14_A", "MAPK14_B", "GSK3B", "AKT1"],
        },
        index=phospho.index.copy(),
    )
    dataset = AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )

    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=ExperimentalDesign(
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
            ),
            contrasts=(
                Contrast(
                    name="B_vs_A",
                    numerator_condition="B",
                    denominator_condition="A",
                ),
            ),
        )
    )
    table = result.table_for("B_vs_A")

    duplicate_rows = table.loc[site_keys[:2], :]
    assert duplicate_rows.loc[:, "display_id"].tolist() == [
        "MAPK14;Y182;",
        "MAPK14;Y182;",
    ]
    assert duplicate_rows.loc[:, "site_key"].tolist() == site_keys[:2]
    assert table.loc[:, "site_key"].nunique() == table.shape[0]
    assert table.index.tolist() == site_keys
    assert (
        table.loc[:, "organism"].tolist() == site_metadata.loc[:, "organism"].tolist()
    )
    assert (
        table.loc[:, "protein_namespace"].tolist()
        == site_metadata.loc[:, "protein_namespace"].tolist()
    )
    assert (
        table.loc[:, "protein_identifier"].tolist()
        == site_metadata.loc[:, "protein_identifier"].tolist()
    )
    assert (
        table.loc[:, "protein_id"].tolist()
        == site_metadata.loc[:, "protein_id"].tolist()
    )


def test_identity_required_result_rejects_display_indexed_table() -> None:
    display_index = pd.Index(["AKT1;T308;", "GSK3B;S9;"], name="site_key")
    site_keys = [
        _site_key(gene_symbol="AKT1", site="T308"),
        _site_key(gene_symbol="GSK3B", site="S9"),
    ]
    table = pd.DataFrame(
        {
            "site_key": site_keys,
            "display_id": display_index.astype(str).tolist(),
            "gene_symbol": ["AKT1", "GSK3B"],
            "site": ["T308", "S9"],
            "logFC": [1.0, -1.0],
            "t": [2.0, -2.0],
            "P.Value": [0.05, 0.10],
            "adj.P.Val": [0.10, 0.10],
        },
        index=display_index,
    )

    with pytest.raises(PhosPyInputError, match="must start with 'phospy:v1'"):
        _manual_result_with_table(table)


def test_display_unique_statistics_are_unchanged_by_site_key_result_identity() -> None:
    matrix = _load_matrix()
    request = _request_for_reverse_contrasts(matrix)
    interpreted = DifferentialAnalysisInterpreter().run(
        DifferentialAnalysisValidator().run(request)
    )
    site_key_result = DifferentialWorkflowExecutor().run(interpreted)
    display_index_request = DifferentialComputationRequest(
        matrix=matrix.copy(deep=True),
        design=interpreted.computation_request.design,
        contrasts=interpreted.computation_request.contrasts,
        empirical_bayes=interpreted.computation_request.empirical_bayes,
    )
    display_index_result = DifferentialComputationExecutor().run(display_index_request)

    assert isinstance(display_index_result, DifferentialComputationResult)
    for contrast_name in ("B_vs_A", "A_vs_B"):
        site_key_table = site_key_result.table_for(contrast_name)
        display_index_table = display_index_result.table_for(contrast_name)
        assert (
            site_key_table.loc[:, "display_id"].astype(str).tolist()
            == display_index_table.index.astype(str).tolist()
        )
        np.testing.assert_allclose(
            site_key_table.loc[:, STATISTIC_COLUMNS].to_numpy(dtype=float),
            display_index_table.loc[:, STATISTIC_COLUMNS].to_numpy(dtype=float),
            rtol=1e-12,
            atol=1e-12,
        )


def test_missing_values_are_rejected_before_differential_execution() -> None:
    matrix = pd.read_csv(NEGATIVE_FIXTURE_DIR / "matrix_with_missing.csv").set_index(
        "site_id"
    )
    raw_site_ids = matrix.index.astype(str).tolist()
    canonical_ids = [
        _canonical_site_id(site_id, ordinal=idx)
        for idx, site_id in enumerate(raw_site_ids, start=1)
    ]
    matrix.index = pd.Index(canonical_ids, name=matrix.index.name)

    with pytest.raises(DatasetValidationError, match="dataset.phospho"):
        _build_dataset(matrix)
