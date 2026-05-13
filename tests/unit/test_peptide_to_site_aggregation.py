from __future__ import annotations

import pandas as pd

from phospy.differential.aggregation import (
    PEPTIDE_TO_SITE_STRATEGY_COMPAT_BEST_P_VALUE,
    PEPTIDE_TO_SITE_STRATEGY_FIXED_EFFECT_META,
    PeptideToSiteAggregationConfig,
    PeptideToSiteAggregator,
)
from phospy.differential.aggregation.scientific_policies import (
    build_peptide_to_site_aggregation_policy,
)
from phospy.evidence import PeptideEvidenceTable


def _evidence_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "peptide_row_id": ["pep_1", "pep_2", "pep_3"],
            "site_id": ["MAPK1;S10;", "MAPK1;S10;", "MAPK1;T12;"],
            "unique_feature_id": ["feat_1", "feat_2", "feat_3"],
            "gene_symbol": ["MAPK1", "MAPK1", "MAPK1"],
            "protein_accession": ["P28482", "P28482", "P28482"],
            "site_string": ["S10", "S10", "T12"],
            "sample_a": [10.0, 8.0, 12.0],
            "sample_b": [11.0, 9.0, 13.0],
            "peptide_sequence": ["AAAAA", "BBBBB", "CCCCC"],
            "modified_peptide_sequence": ["AA[+80]AAA", "BB[+80]BBB", "CC[+80]CCC"],
            "multi_site": [False, False, False],
            "provenance_source": ["maxquant", "maxquant", "maxquant"],
        }
    )


def _differential_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "logFC": [1.0, -0.8, 0.5],
            "t": [5.0, -4.0, 3.0],
            "P.Value": [1e-4, 2e-4, 5e-3],
            "adj.P.Val": [3e-4, 4e-4, 7e-3],
        },
        index=pd.Index(["feat_1", "feat_2", "feat_3"], name="feature_id"),
    )


def test_default_strategy_is_not_minimum_p_value() -> None:
    assert (
        PeptideToSiteAggregationConfig().strategy
        == PEPTIDE_TO_SITE_STRATEGY_FIXED_EFFECT_META
    )
    assert (
        PeptideToSiteAggregationConfig().strategy
        != PEPTIDE_TO_SITE_STRATEGY_COMPAT_BEST_P_VALUE
    )


def test_peptide_to_site_scientific_policy_ownership_is_explicit() -> None:
    assert (
        build_peptide_to_site_aggregation_policy.__module__
        == "phospy.differential.aggregation.scientific_policies"
    )


def test_aggregation_handles_one_and_multiple_peptides_per_site() -> None:
    evidence = PeptideEvidenceTable(
        frame=_evidence_frame(),
        sample_intensity_columns=("sample_a", "sample_b"),
    )
    aggregated = PeptideToSiteAggregator().run_table(
        peptide_differential_table=_differential_table(),
        evidence=evidence,
        config=PeptideToSiteAggregationConfig(
            strategy=PEPTIDE_TO_SITE_STRATEGY_FIXED_EFFECT_META
        ),
        contrast_name="B_vs_A",
    )
    table = aggregated.to_dataframe()

    assert set(table.index.tolist()) == {"MAPK1;S10;", "MAPK1;T12;"}
    assert int(table.loc["MAPK1;S10;", "n_peptide_observations"]) == 2
    assert int(table.loc["MAPK1;T12;", "n_peptide_observations"]) == 1
    assert float(table.loc["MAPK1;T12;", "logFC"]) == 0.5


def test_conflicting_peptide_effects_reduce_site_level_significance() -> None:
    evidence = PeptideEvidenceTable(
        frame=_evidence_frame(),
        sample_intensity_columns=("sample_a", "sample_b"),
    )
    aggregated = PeptideToSiteAggregator().run_table(
        peptide_differential_table=_differential_table(),
        evidence=evidence,
        config=PeptideToSiteAggregationConfig(
            strategy=PEPTIDE_TO_SITE_STRATEGY_FIXED_EFFECT_META
        ),
        contrast_name="B_vs_A",
    )
    table = aggregated.to_dataframe()
    site_row = table.loc["MAPK1;S10;"]

    assert abs(float(site_row["logFC"])) < 0.2
    assert float(site_row["P.Value"]) > 0.05


def test_compatibility_mode_reproduces_minimum_p_value_selection() -> None:
    evidence = PeptideEvidenceTable(
        frame=_evidence_frame(),
        sample_intensity_columns=("sample_a", "sample_b"),
    )
    aggregated = PeptideToSiteAggregator().run_table(
        peptide_differential_table=_differential_table(),
        evidence=evidence,
        config=PeptideToSiteAggregationConfig(
            strategy=PEPTIDE_TO_SITE_STRATEGY_COMPAT_BEST_P_VALUE
        ),
        contrast_name="B_vs_A",
    )
    table = aggregated.to_dataframe()
    site_row = table.loc["MAPK1;S10;"]

    assert float(site_row["logFC"]) == 1.0
    assert float(site_row["P.Value"]) == 1e-4
    assert int(site_row["n_peptides_used"]) == 1


def test_missing_variance_rows_are_excluded_from_variance_weighted_strategies() -> None:
    evidence = PeptideEvidenceTable(
        frame=_evidence_frame(),
        sample_intensity_columns=("sample_a", "sample_b"),
    )
    table = _differential_table()
    table.loc["feat_2", "t"] = 0.0

    aggregated = PeptideToSiteAggregator().run_table(
        peptide_differential_table=table,
        evidence=evidence,
        config=PeptideToSiteAggregationConfig(
            strategy=PEPTIDE_TO_SITE_STRATEGY_FIXED_EFFECT_META
        ),
        contrast_name="B_vs_A",
    )
    site_row = aggregated.to_dataframe().loc["MAPK1;S10;"]
    assert int(site_row["n_peptide_observations"]) == 2
    assert int(site_row["n_peptides_used"]) == 1


def test_scientific_policy_metadata_warns_for_compatibility_min_p_mode() -> None:
    evidence = PeptideEvidenceTable(
        frame=_evidence_frame(),
        sample_intensity_columns=("sample_a", "sample_b"),
    )
    result = PeptideToSiteAggregator().run_table(
        peptide_differential_table=_differential_table(),
        evidence=evidence,
        config=PeptideToSiteAggregationConfig(
            strategy=PEPTIDE_TO_SITE_STRATEGY_COMPAT_BEST_P_VALUE
        ),
        contrast_name="B_vs_A",
    )
    policy = result.scientific_policies[0]
    assert result.warnings
    assert policy.parameters["compatibility_mode_warning"] is True
    assert policy.parameters["strategy"] == PEPTIDE_TO_SITE_STRATEGY_COMPAT_BEST_P_VALUE
    assert "multi_site_handling" in result.provenance
    handling = result.provenance["multi_site_handling"]
    assert isinstance(handling, dict)
    assert handling["statistical_modeling_policy"] == "keep_joint"


def test_aggregation_preserves_joint_multi_site_ids_by_default() -> None:
    multi_frame = pd.DataFrame(
        {
            "peptide_row_id": ["pep_joint"],
            "site_id": ["MAPK1;S10;"],
            "unique_feature_id": ["feat_joint"],
            "gene_symbol": ["MAPK1"],
            "protein_accession": ["P28482"],
            "site_string": ["S10,T12"],
            "sample_a": [10.0],
            "sample_b": [11.0],
            "peptide_sequence": ["AAAAA"],
            "modified_peptide_sequence": ["AA[+80]AAA"],
            "multi_site": [True],
            "provenance_source": ["maxquant"],
        }
    )
    evidence = PeptideEvidenceTable(
        frame=multi_frame,
        sample_intensity_columns=("sample_a", "sample_b"),
    )
    differential = pd.DataFrame(
        {
            "logFC": [0.75],
            "t": [4.0],
            "P.Value": [0.01],
            "adj.P.Val": [0.01],
        },
        index=pd.Index(["feat_joint"], name="feature_id"),
    )

    result = PeptideToSiteAggregator().run_table(
        peptide_differential_table=differential,
        evidence=evidence,
        config=PeptideToSiteAggregationConfig(
            strategy=PEPTIDE_TO_SITE_STRATEGY_FIXED_EFFECT_META
        ),
        contrast_name="B_vs_A",
    )
    table = result.to_dataframe()
    assert list(table.index.tolist()) == ["MAPK1;S10,T12;"]
