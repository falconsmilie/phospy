from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import DatasetBuildRequest, PhosPyInputError
from phospy.api.requests import (
    DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING,
    DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
    DATASET_MULTI_SITE_POLICY_REJECT,
    DATASET_MULTI_SITE_POLICY_SPLIT,
    DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
    DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED,
)
from phospy.science.datasets.builders.validator import DatasetBuildRequestValidator


def _site_level_phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )


def _site_level_site_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            "protein_id": ["MAPK14"],
        },
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )


def _peptide_evidence_frame(*, include_single_site: bool = True) -> pd.DataFrame:
    rows = [
        {
            "peptide_row_id": "pep_joint",
            "site_id": "MAPK1;S10;",
            "unique_feature_id": "feat_joint",
            "gene_symbol": "MAPK1",
            "protein_accession": "P28482",
            "site_string": "S10,T12",
            "sample_a": 10.0,
            "sample_b": 12.0,
            "peptide_sequence": "AAAAA",
            "modified_peptide_sequence": "AA[+80]AAA",
            "multi_site": True,
            "provenance_source": "maxquant",
            "site_sequence": "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
            "localisation_confidence": 0.95,
        }
    ]
    if include_single_site:
        rows.append(
            {
                "peptide_row_id": "pep_single",
                "site_id": "AKT1;S473;",
                "unique_feature_id": "feat_single",
                "gene_symbol": "AKT1",
                "protein_accession": "P31749",
                "site_string": "S473",
                "sample_a": 7.0,
                "sample_b": 9.0,
                "peptide_sequence": "BBBBB",
                "modified_peptide_sequence": "BB[+80]BBB",
                "multi_site": False,
                "provenance_source": "maxquant",
                "site_sequence": "BBBBBBBBBBBBBBBSBBBBBBBBBBBBBBB",
                "localisation_confidence": 0.9,
            }
        )
    return pd.DataFrame(rows)


def test_site_level_input_works_with_safe_default_or_explicit_declaration() -> None:
    safe_default_request = DatasetBuildRequest(
        phospho=_site_level_phospho(),
        site_metadata=_site_level_site_metadata(),
        input_intensity_scale="linear",
    )
    explicit_request = DatasetBuildRequest(
        phospho=_site_level_phospho(),
        site_metadata=_site_level_site_metadata(),
        input_intensity_scale="linear",
        site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED,
    )
    assert (
        DatasetBuildRequestValidator().run(safe_default_request) is safe_default_request
    )
    assert DatasetBuildRequestValidator().run(explicit_request) is explicit_request


def test_peptide_evidence_requires_multi_site_policy() -> None:
    request = DatasetBuildRequest(
        site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
        peptide_evidence=_peptide_evidence_frame(),
        peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
        input_intensity_scale="linear",
    )
    with pytest.raises(
        PhosPyInputError,
        match="peptide_evidence input requires multi_site_policy",
    ):
        DatasetBuildRequestValidator().run(request)


def test_reject_policy_fails_on_ambiguous_peptide_evidence() -> None:
    request = DatasetBuildRequest(
        site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
        peptide_evidence=_peptide_evidence_frame(include_single_site=False),
        peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
        multi_site_policy=DATASET_MULTI_SITE_POLICY_REJECT,
        input_intensity_scale="linear",
    )
    with pytest.raises(
        PhosPyInputError,
        match="multi-site observation cannot be resolved under policy='error'",
    ):
        AnalysisReadyDatasetBuilder().run(request)


def test_exclude_policy_records_exclusions_in_report_and_provenance() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            peptide_evidence=_peptide_evidence_frame(),
            peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING,
            input_intensity_scale="linear",
        )
    )
    assert list(built.phospho.index.tolist()) == ["AKT1;S473;"]
    assert built.preprocessing_report is not None
    resolution_rows = built.preprocessing_report.operations.loc[
        built.preprocessing_report.operations.loc[:, "stage"]
        == "peptide_evidence_resolution"
    ]
    assert int(resolution_rows.shape[0]) == 1
    parameters = resolution_rows.iloc[0]["parameters"]
    assert isinstance(parameters, dict)
    assert int(parameters["peptide_observations_received"]) == 2
    assert int(parameters["ambiguous_observations"]) == 1
    assert int(parameters["excluded_observations"]) == 1
    assert int(parameters["split_observations"]) == 0
    assert built.provenance is not None
    payload = built.provenance.workflow_parameters["peptide_evidence_resolution"]
    assert isinstance(payload, dict)
    assert (
        payload["multi_site_policy"]
        == DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING
    )


def test_keep_joint_policy_preserves_joint_ambiguous_site_representation() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            peptide_evidence=_peptide_evidence_frame(include_single_site=False),
            peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_KEEP_JOINT,
            input_intensity_scale="linear",
        )
    )
    assert list(built.phospho.index.tolist()) == ["MAPK1;S10,T12;"]
    site_value = built.site_metadata.loc["MAPK1;S10,T12;", "site"]
    assert str(site_value) == "S10,T12"


def test_split_policy_applies_deterministic_equal_split() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            site_resolution_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
            peptide_evidence=_peptide_evidence_frame(include_single_site=False),
            peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
            multi_site_policy=DATASET_MULTI_SITE_POLICY_SPLIT,
            input_intensity_scale="linear",
        )
    )
    assert set(built.phospho.index.tolist()) == {"MAPK1;S10;", "MAPK1;T12;"}
    assert float(built.phospho.loc["MAPK1;S10;", "sample_a"]) == pytest.approx(5.0)
    assert float(built.phospho.loc["MAPK1;T12;", "sample_a"]) == pytest.approx(5.0)
    assert built.provenance is not None
    payload = built.provenance.workflow_parameters["peptide_evidence_resolution"]
    assert isinstance(payload, dict)
    assert int(payload["split_observations"]) == 1
