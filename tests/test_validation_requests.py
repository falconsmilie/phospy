from __future__ import annotations

import pandas as pd
import pytest

from phospy.validation.errors import RequestValidationError
from phospy.validation.requests import CorePipelineRequest, KinaseWorkflowRequest


def test_core_pipeline_request_requires_existing_paths(tmp_path) -> None:
    total_path = tmp_path / "missing_total.tsv"
    phospho_path = tmp_path / "missing_phospho.tsv"

    with pytest.raises(RequestValidationError, match="Path does not exist"):
        CorePipelineRequest.validate_request(
            total_path=total_path,
            phospho_path=phospho_path,
        )


def test_core_pipeline_request_rejects_duplicate_comparisons(tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    total_path.write_text("genes\tgroup1\nPRKACA\t1\n")
    phospho_path.write_text("uid\tgene_names\n1\tPRKACA\n")

    with pytest.raises(RequestValidationError, match="Duplicate comparison pair"):
        CorePipelineRequest.validate_request(
            total_path=total_path,
            phospho_path=phospho_path,
            comparisons=(("group1", "group4"), ("group1", "group4")),
        )


def test_core_pipeline_request_rejects_unknown_comparison_group(tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    total_path.write_text("genes\tgroup1\nPRKACA\t1\n")
    phospho_path.write_text("uid\tgene_names\n1\tPRKACA\n")

    with pytest.raises(RequestValidationError, match="Unknown comparison group"):
        CorePipelineRequest.validate_request(
            total_path=total_path,
            phospho_path=phospho_path,
            comparisons=(("group1", "group9"),),
        )


def test_kinase_workflow_request_rejects_invalid_threshold() -> None:
    phospho_matrix = pd.DataFrame({"sample_1": [1.0]}, index=["SITE_1"])

    with pytest.raises(RequestValidationError, match="score_threshold"):
        KinaseWorkflowRequest.validate_request(
            phospho_matrix=phospho_matrix,
            substrate_map={"KINASE_A": ["SITE_1"]},
            motif_sequences=None,
            allow_profile_only_fallback=True,
            score_threshold=1.5,
        )


def test_kinase_workflow_request_requires_site_sequences_with_motifs() -> None:
    phospho_matrix = pd.DataFrame({"sample_1": [1.0]}, index=["SITE_1"])

    with pytest.raises(RequestValidationError, match="site_sequences are required"):
        KinaseWorkflowRequest.validate_request(
            phospho_matrix=phospho_matrix,
            substrate_map={"KINASE_A": ["SITE_1"]},
            motif_sequences={"KINASE_A": ["QQAAAAAYY"]},
        )


def test_core_pipeline_request_rejects_invalid_max_unmatched_fraction(tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    total_path.write_text("genes\tgroup1\nPRKACA\t1\n")
    phospho_path.write_text("uid\tgene_names\n1\tPRKACA\n")

    with pytest.raises(RequestValidationError, match="max_unmatched_fraction"):
        CorePipelineRequest.validate_request(
            total_path=total_path,
            phospho_path=phospho_path,
            max_unmatched_fraction=1.5,
        )


def test_core_pipeline_request_rejects_directory_paths(tmp_path) -> None:
    total_dir = tmp_path / "total_dir"
    phospho_path = tmp_path / "phospho.tsv"
    total_dir.mkdir()
    phospho_path.write_text(
        "uid	gene_names	gene_p_site	localization_prob	centralized_sequence	p_group1	p_group2	p_group3	p_group4	p_group5	p_group6\n1	PRKACA	PRKACA_S339	0.95	AAAAAA	1	1	1	1	1	1\n"
    )

    with pytest.raises(RequestValidationError, match="Path is not a file"):
        CorePipelineRequest.validate_request(
            total_path=total_dir,
            phospho_path=phospho_path,
        )


def test_kinase_workflow_request_rejects_plain_site_sequence_lists() -> None:
    phospho_matrix = pd.DataFrame({"sample_1": [1.0]}, index=["SITE_1"])

    with pytest.raises(RequestValidationError, match="site_sequences must be provided"):
        KinaseWorkflowRequest.validate_request(
            phospho_matrix=phospho_matrix,
            substrate_map={"KINASE_A": ["SITE_1"]},
            site_sequences=["QQAAAAAYY"],
            motif_sequences={"KINASE_A": ["QQAAAAAYY"]},
        )
