from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import PhosPyInputError
from phospy.science.evidence import (
    MULTI_SITE_POLICY_ERROR,
    MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING,
    MULTI_SITE_POLICY_FIRST_SITE_COMPATIBILITY,
    MULTI_SITE_POLICY_KEEP_JOINT,
    MULTI_SITE_POLICY_SPLIT_EQUAL_WEIGHT,
    MULTI_SITE_POLICY_SPLIT_WITH_UNCERTAINTY,
    MultiSiteHandlingConfig,
    PeptideEvidenceTable,
)
from phospy.science.evidence.multi_site import parse_phospho_site_tokens


def _base_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "peptide_row_id": ["row_1", "row_2"],
            "site_id": ["MAPK1;S1246;", "MAPK1;S10;"],
            "unique_feature_id": ["feat_1", "feat_2"],
            "gene_symbol": ["MAPK1", "MAPK1"],
            "protein_accession": ["P28482", "P28482"],
            "site_string": ["S1246,T1247", "S10"],
            "sample_a": [10.0, 11.0],
            "sample_b": [12.0, 13.0],
            "peptide_sequence": ["AAAAA", "BBBBB"],
            "modified_peptide_sequence": ["AA[+80]AAA", "BB[+80]BBB"],
            "multi_site": [True, False],
            "provenance_source": ["maxquant", "maxquant"],
        }
    )


def test_parse_site_string_supports_explicit_multi_site_tokens() -> None:
    parsed = parse_phospho_site_tokens(
        "S1246,T1247",
        field_name="site_string",
    )
    assert [item.token for item in parsed] == ["S1246", "T1247"]


def test_parse_site_string_supports_three_site_and_mixed_residue_tokens() -> None:
    parsed = parse_phospho_site_tokens(
        "S12;T13;Y14",
        field_name="site_string",
    )
    assert [item.token for item in parsed] == ["S12", "T13", "Y14"]


@pytest.mark.parametrize(
    "bad_site_string",
    [
        "",
        "S0,T2",
        "S12,,T13",
        "A12,T13",
        "S12|T13",
    ],
)
def test_parse_site_string_rejects_malformed_values(bad_site_string: str) -> None:
    with pytest.raises(PhosPyInputError):
        parse_phospho_site_tokens(
            bad_site_string,
            field_name="site_string",
        )


def test_default_statistical_policy_keeps_joint_site_id_for_multi_site_rows() -> None:
    evidence = PeptideEvidenceTable(
        frame=_base_frame(),
        sample_intensity_columns=("sample_a", "sample_b"),
    )
    mapping = evidence.site_mapping.to_dataframe()
    assert "MAPK1;S1246,T1247;" in set(mapping.loc[:, "site_id"].tolist())
    assert "MAPK1;S10;" in set(mapping.loc[:, "site_id"].tolist())


def test_first_site_compatibility_is_only_enabled_explicitly() -> None:
    evidence = PeptideEvidenceTable(
        frame=_base_frame(),
        sample_intensity_columns=("sample_a", "sample_b"),
        multi_site_handling_config=MultiSiteHandlingConfig(
            statistical_modeling_policy=MULTI_SITE_POLICY_FIRST_SITE_COMPATIBILITY,
            kinase_sequence_scoring_policy=MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING,
        ),
    )
    mapping = evidence.site_mapping.to_dataframe()
    assert "MAPK1;S1246;" in set(mapping.loc[:, "site_id"].tolist())
    assert "MAPK1;S1246,T1247;" not in set(mapping.loc[:, "site_id"].tolist())


def test_split_policies_emit_row_level_sites_without_row_explosion() -> None:
    evidence = PeptideEvidenceTable(
        frame=_base_frame(),
        sample_intensity_columns=("sample_a", "sample_b"),
        multi_site_handling_config=MultiSiteHandlingConfig(
            statistical_modeling_policy=MULTI_SITE_POLICY_SPLIT_EQUAL_WEIGHT,
            kinase_sequence_scoring_policy=MULTI_SITE_POLICY_SPLIT_WITH_UNCERTAINTY,
        ),
    )
    assert evidence.to_dataframe().shape[0] == 2
    mapping = evidence.site_mapping.to_dataframe()
    assert {"MAPK1;S1246;", "MAPK1;T1247;", "MAPK1;S10;"} == set(
        mapping.loc[:, "site_id"].tolist()
    )
    split_weights = mapping.loc[
        mapping.loc[:, "site_id"].isin({"MAPK1;S1246;", "MAPK1;T1247;"}),
        "mapping_weight",
    ]
    assert split_weights.tolist() == [0.5, 0.5]
    kinase_mapping = evidence.kinase_sequence_site_mapping()
    uncertain = kinase_mapping.loc[
        kinase_mapping.loc[:, "site_id"].isin({"MAPK1;S1246;", "MAPK1;T1247;"}),
        "mapping_uncertainty",
    ]
    assert uncertain.tolist() == [True, True]


def test_kinase_sequence_policy_can_exclude_or_keep_joint_multi_site_rows() -> None:
    excluded = PeptideEvidenceTable(
        frame=_base_frame(),
        sample_intensity_columns=("sample_a", "sample_b"),
        multi_site_handling_config=MultiSiteHandlingConfig(
            statistical_modeling_policy=MULTI_SITE_POLICY_KEEP_JOINT,
            kinase_sequence_scoring_policy=MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING,
        ),
    )
    excluded_mapping = excluded.kinase_sequence_site_mapping()
    assert set(excluded_mapping.loc[:, "site_id"].tolist()) == {"MAPK1;S10;"}

    joint = PeptideEvidenceTable(
        frame=_base_frame(),
        sample_intensity_columns=("sample_a", "sample_b"),
        multi_site_handling_config=MultiSiteHandlingConfig(
            statistical_modeling_policy=MULTI_SITE_POLICY_KEEP_JOINT,
            kinase_sequence_scoring_policy=MULTI_SITE_POLICY_KEEP_JOINT,
        ),
    )
    joint_mapping = joint.kinase_sequence_site_mapping()
    assert "MAPK1;S1246,T1247;" in set(joint_mapping.loc[:, "site_id"].tolist())


def test_multi_site_policy_error_rejects_multi_site_rows() -> None:
    with pytest.raises(PhosPyInputError, match="policy='error'"):
        PeptideEvidenceTable(
            frame=_base_frame(),
            sample_intensity_columns=("sample_a", "sample_b"),
            multi_site_handling_config=MultiSiteHandlingConfig(
                statistical_modeling_policy=MULTI_SITE_POLICY_ERROR,
                kinase_sequence_scoring_policy=MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING,
            ),
        )


def test_multi_site_provenance_records_policy_selection() -> None:
    evidence = PeptideEvidenceTable(
        frame=_base_frame(),
        sample_intensity_columns=("sample_a", "sample_b"),
        multi_site_handling_config=MultiSiteHandlingConfig(
            statistical_modeling_policy=MULTI_SITE_POLICY_KEEP_JOINT,
            kinase_sequence_scoring_policy=MULTI_SITE_POLICY_FIRST_SITE_COMPATIBILITY,
        ),
    )
    provenance = evidence.multi_site_policy_provenance()
    assert provenance["statistical_modeling_policy"] == MULTI_SITE_POLICY_KEEP_JOINT
    assert (
        provenance["kinase_sequence_scoring_policy"]
        == MULTI_SITE_POLICY_FIRST_SITE_COMPATIBILITY
    )
    assert int(provenance["multi_site_rows"]) == 1
    assert int(provenance["single_site_rows"]) == 1
