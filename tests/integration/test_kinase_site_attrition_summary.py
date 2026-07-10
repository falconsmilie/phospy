from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
from phospy.api import (
    DatasetBuildRequest,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferenceContextCompatibilityPolicy,
)
from phospy.io.publishers.workflows import publish_kinase_workflow

pytestmark = pytest.mark.integration


def _centred_sequence(residue: str) -> str:
    return ("A" * 15) + residue + ("A" * 15)


def _build_dataset_with_attrition_mix() -> object:
    phospho = pd.DataFrame(
        {
            "sample_a": [8.0, 7.0, 5.0, 4.0, 3.0, 2.0, 9.0, 6.0],
            "sample_b": [9.0, 6.0, 5.5, 4.5, 3.5, 2.5, np.nan, 6.5],
        },
        index=pd.Index(
            [
                "source_1",
                "source_2",
                "source_3",
                "source_4",
                "source_5",
                "source_6",
                "source_7",
                "source_8",
            ],
            name="source_row_id",
        ),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": [
                "MAPK14",
                "MAPK14",
                "GSK3B",
                "AKT1",
                "BADID",
                "NOREF",
                "MISSQ",
                "REF2",
            ],
            "protein_id": [
                "MAPK14",
                "MAPK14",
                "GSK3B",
                "AKT1",
                "BADID",
                "NOREF",
                "MISSQ",
                "REF2",
            ],
            "site": ["Y182", "Y183", "S9", "T308", "S123", "S1", "S1", "T1"],
            "site_sequence": [
                _centred_sequence("Y"),
                _centred_sequence("Y"),
                _centred_sequence("S"),
                pd.NA,
                _centred_sequence("S"),
                _centred_sequence("S"),
                _centred_sequence("S"),
                _centred_sequence("T"),
            ],
            "localisation_confidence": [0.95] * phospho.shape[0],
        },
        index=phospho.index.copy(),
    )
    preprocessing = DatasetPreprocessingConfig(
        missing_data=DatasetMissingDataConfig(
            policy="impute_row_median",
            min_observed_values=1,
        ),
        site_matrix=DatasetSiteMatrixConfig(
            policy="build_from_metadata",
            duplicate_site_policy="first",
        ),
    )
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
            allow_opaque_site_values=True,
            preprocessing_config=preprocessing,
        )
    )


def _build_references() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1", "K1", "K1", "K2", "K2", "K2"],
                "substrate_site": [
                    "MAPK14;Y182;",
                    "GSK3B;S9;",
                    "REF2;T1;",
                    "MAPK14;Y182;",
                    "GSK3B;S9;",
                    "REF2;T1;",
                ],
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    _centred_sequence("Y"),
                    _centred_sequence("S"),
                    _centred_sequence("T"),
                    _centred_sequence("S"),
                    _centred_sequence("S"),
                ]
            },
            index=pd.Index(
                [
                    "MAPK14;Y182;",
                    "GSK3B;S9;",
                    "REF2;T1;",
                    "BADID;123;",
                    "NOREF;S1;",
                ],
                name="site_id",
            ),
        ),
    )


def test_kinase_workflow_exposes_compact_site_attrition_summary() -> None:
    dataset = _build_dataset_with_attrition_mix()
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=_build_references(),
            scoring_config=KinaseScoringConfig(
                min_substrates=2,
                reference_context_compatibility_policy=(
                    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                ),
            ),
            prediction_config=KinasePredictionConfig(
                top_k=3,
                deterministic_max_selected_kinases=3,
                adaptive_ensemble_runs=3,
            ),
            activity_config=KinaseActivityConfig(
                enabled=True,
                threshold=0.0,
                min_substrates=1,
                top_n_substrates=3,
            ),
        )
    )

    assert result.site_attrition_summary is not None
    assert result.eligibility_report is not None
    preprocessing = result.site_attrition_summary.preprocessing
    scoring = result.site_attrition_summary.scoring
    eligibility = result.eligibility_report

    assert preprocessing.input_rows == 8
    assert preprocessing.rows_removed_during_preprocessing == 1
    assert preprocessing.rows_removed_invalid_or_missing_site_identifiers == 0
    assert preprocessing.duplicate_sites_merged_or_resolved == 0
    assert preprocessing.output_rows == 7
    assert preprocessing.sequence_complete_sites == 7

    assert eligibility.total_dataset_sites == 7
    assert eligibility.sequence_complete_sites == 7
    assert eligibility.localisation_eligible_sites == 7
    assert eligibility.reference_overlap_sites == 3
    assert eligibility.excluded_no_reference_match == 4
    assert eligibility.excluded_low_localisation == 0
    assert eligibility.eligible_kinases == 2
    assert eligibility.excluded_kinases_below_min_substrates == 0

    assert scoring.rows_removed_invalid_or_missing_site_identifiers == 0
    assert scoring.final_quantitative_sites_entering_scoring == 7
    assert scoring.sites_with_valid_site_sequence == 7
    assert scoring.sites_without_usable_site_sequence == 0
    assert scoring.sites_eligible_for_motif_scoring == 7
    assert scoring.sites_with_kinase_substrate_reference_profile_evidence == 3
    assert scoring.sites_contributing_to_final_fused_prediction_scoring_output == int(
        result.prediction_result.pred_mat.notna().any(axis=1).sum()
    )
    assert scoring.sites_contributing_to_activity_scoring is not None
    assert scoring.sites_contributing_to_activity_scoring == int(
        result.prediction_result.pred_mat.loc[
            result.prediction_result.pred_mat.index.intersection(
                result.dataset.phospho.index
            )
        ]
        .notna()
        .any(axis=1)
        .sum()
    )

    assert result.dataset.preprocessing_report is not None
    preprocessing_report_summary = (
        result.dataset.preprocessing_report.site_attrition_summary()
    )
    assert preprocessing_report_summary.input_rows == preprocessing.input_rows
    assert (
        preprocessing_report_summary.rows_removed_during_preprocessing
        == preprocessing.rows_removed_during_preprocessing
    )
    assert (
        preprocessing_report_summary.rows_removed_invalid_or_missing_site_identifiers
        == preprocessing.rows_removed_invalid_or_missing_site_identifiers
    )
    assert (
        preprocessing_report_summary.duplicate_sites_merged_or_resolved
        == preprocessing.duplicate_sites_merged_or_resolved
    )


def test_published_kinase_manifest_includes_site_attrition_summary(
    tmp_path: Path,
) -> None:
    dataset = _build_dataset_with_attrition_mix()
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=_build_references(),
            scoring_config=KinaseScoringConfig(
                min_substrates=2,
                reference_context_compatibility_policy=(
                    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                ),
            ),
            prediction_config=KinasePredictionConfig(
                top_k=3,
                deterministic_max_selected_kinases=3,
                adaptive_ensemble_runs=3,
            ),
            activity_config=None,
        )
    )

    written = publish_kinase_workflow(
        result, tmp_path / "published", output_format="csv"
    )
    manifest = json.loads(written["kinase.manifest"].read_text(encoding="utf-8"))
    payload = manifest.get("site_attrition_summary")
    eligibility_payload = manifest.get("eligibility_report")

    assert payload is not None
    assert eligibility_payload is not None
    assert payload["preprocessing"]["input_rows"] == 8
    assert payload["preprocessing"]["rows_removed_during_preprocessing"] == 1
    assert payload["preprocessing"]["duplicate_sites_merged_or_resolved"] == 0
    assert payload["preprocessing"]["sequence_complete_sites"] == 7
    assert payload["scoring"]["final_quantitative_sites_entering_scoring"] == 7
    assert eligibility_payload["total_dataset_sites"] == 7
    assert eligibility_payload["sequence_complete_sites"] == 7
    assert eligibility_payload["localisation_eligible_sites"] == 7
    assert eligibility_payload["reference_overlap_sites"] == 3
    assert eligibility_payload["excluded_no_reference_match"] == 4
    assert eligibility_payload["excluded_low_localisation"] == 0
    assert eligibility_payload["eligible_kinases"] == 2
    assert eligibility_payload["excluded_kinases_below_min_substrates"] == 0
    assert json.loads(json.dumps(eligibility_payload)) == eligibility_payload
