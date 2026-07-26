from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from phospy.api import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api.configs import ProfileSelfInclusionPolicy
from phospy.errors.validation import WorkflowValidationError
from phospy.provenance.models import ReferenceProvenance
from phospy.science.activities.methods import (
    SSGSEA_SIGNIFICANCE_STATUS_AVAILABLE,
    SSGSEA_SIGNIFICANCE_STATUS_UNAVAILABLE_NO_PERMUTATIONS,
    SsgseaSubstrateEnrichmentActivityMethod,
)
from phospy.science.references.models import ReferenceContext
from phospy.tables.kinase import (
    KINASE_PROFILE_SCORE_DIAGNOSTIC_REASON_INSUFFICIENT_SUBSTRATES_AFTER_LEAVE_ONE_OUT,
    KINASE_PROFILE_SCORE_DIAGNOSTIC_STATUS_UNSCORED,
)
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)

pytestmark = pytest.mark.release_gate

FIXTURE_DIR = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "release_validation_regression"
    / "kinase_sparse_support"
)


def _manifest() -> dict[str, Any]:
    return json.loads((FIXTURE_DIR / "MANIFEST.json").read_text(encoding="utf-8"))


def _contracts() -> dict[str, Any]:
    return json.loads(
        (FIXTURE_DIR / "expected_contracts.json").read_text(encoding="utf-8")
    )


def _context() -> ReferenceContext:
    return ReferenceContext(
        organism="rat",
        protein_namespace="protein_id",
        source_name="release-validation-regression",
        source_version="2026-07-24",
        proteome_version=None,
        reference_table_sha256="b" * 64,
    )


def _dataset(*, with_reference_context: bool = False):
    phospho = pd.read_csv(FIXTURE_DIR / "phospho.csv").set_index("site_id")
    display_ids = phospho.index.astype(str).tolist()
    site_index = site_key_index_from_display_ids(display_ids)
    phospho.index = site_index.copy()
    metadata = pd.read_csv(FIXTURE_DIR / "site_metadata.csv").set_index("site_id")
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": display_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": metadata.loc[display_ids, "gene_symbol"]
            .astype(str)
            .tolist(),
            "protein_id": metadata.loc[display_ids, "protein_id"].astype(str).tolist(),
            "site": metadata.loc[display_ids, "site"].astype(str).tolist(),
            "site_sequence": metadata.loc[display_ids, "site_sequence"]
            .astype(str)
            .tolist(),
            "localisation_confidence": metadata.loc[
                display_ids, "localisation_confidence"
            ]
            .astype(float)
            .tolist(),
        },
        index=site_index.copy(),
    )
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )
    if not with_reference_context:
        return dataset
    if dataset.provenance is None:
        raise AssertionError("test dataset must carry construction provenance")
    return trusted_analysis_ready_dataset_from_tables(
        phospho=dataset.phospho,
        site_metadata=dataset.site_metadata,
        organism=dataset.organism,
        intensity_scale_state=dataset.intensity_scale_state,
        processing_state=dataset.processing_state,
        provenance=replace(dataset.provenance, reference_context=_context()),
    )


def _references(*, with_reference_context: bool = False) -> ReferenceBundle:
    sequences = pd.read_csv(FIXTURE_DIR / "site_sequences.csv").set_index("site_id")
    provenance = None
    if with_reference_context:
        context = _context()
        provenance = ReferenceProvenance(
            source_type="explicit",
            organism=Organism.RAT.value,
            bundle_id=None,
            source_name=context.source_name,
            source_version=context.source_version,
            identifier_namespace=context.protein_namespace,
            table_fingerprints=(),
            reference_context=context,
        )
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.read_csv(FIXTURE_DIR / "substrate_map.csv"),
        site_sequences=sequences,
        provenance=provenance,
    )


def _request(
    *,
    scoring_config: KinaseScoringConfig | None = None,
    with_reference_context: bool = False,
) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=_dataset(with_reference_context=with_reference_context),
        references=_references(with_reference_context=with_reference_context),
        scoring_config=scoring_config
        or KinaseScoringConfig(
            min_substrates=2,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=4,
            deterministic_max_selected_kinases=4,
            adaptive_ensemble_runs=2,
        ),
        activity_config=None,
    )


def _ssgsea_effect_matrix() -> pd.DataFrame:
    frame = pd.read_csv(FIXTURE_DIR / "ssgsea_effect_matrix.csv").set_index("site_id")
    display_ids = frame.index.astype(str).tolist()
    frame.index = site_key_index_from_display_ids(display_ids)
    return frame


def _ssgsea_membership() -> pd.DataFrame:
    substrate_map = pd.read_csv(FIXTURE_DIR / "substrate_map.csv")
    rows = substrate_map.loc[
        substrate_map.loc[:, "kinase"].isin(["K_SSGSEA_HIGH", "K_SSGSEA_LOW"]),
        :,
    ].copy(deep=True)
    rows.loc[:, "substrate_site"] = (
        site_key_index_from_display_ids(
            rows.loc[:, "substrate_site"].astype(str).tolist()
        )
        .astype(str)
        .tolist()
    )
    return rows


def test_kinase_sparse_support_fixture_manifest_hashes_match_files() -> None:
    manifest = _manifest()

    assert manifest["classification"] == "regression"
    assert manifest["fixture_family"] == "kinase_sparse_support"
    assert "not external parity" in manifest["source_policy"]
    assert manifest["seed"] == 20260724
    assert _contracts()["external_reference"] is None

    for file_entry in manifest["files"]:
        path = FIXTURE_DIR / str(file_entry["relative_path"])
        assert path.is_file()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == file_entry["sha256"]


def test_kinase_sparse_support_fixture_exercises_substrate_floor_and_overlap() -> None:
    result = KinaseWorkflow().run(_request())
    score_columns = set(result.scoring_result.profile_scores.columns.astype(str))
    contracts = _contracts()["substrate_support_classes"]

    assert set(contracts["below_minimum"]).isdisjoint(score_columns)
    assert set(contracts["at_minimum"]) <= score_columns
    assert set(contracts["above_minimum"]) <= score_columns
    assert "K_SPARSE" not in score_columns
    assert result.site_attrition_summary is not None
    assert (
        result.site_attrition_summary.scoring.sites_with_kinase_substrate_reference_profile_evidence
        >= 4
    )


def test_kinase_leave_one_out_changes_scores_and_reports_minimum_support_loss() -> None:
    allow = KinaseWorkflow().run(_request())
    leave_one_out = KinaseWorkflow().run(
        _request(
            scoring_config=KinaseScoringConfig(
                min_substrates=2,
                profile_self_inclusion_policy=ProfileSelfInclusionPolicy.LEAVE_ONE_OUT,
                reference_context_compatibility_policy=(
                    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                ),
            )
        )
    )
    scored_site = allow.scoring_result.profile_scores.index[0]

    assert (
        allow.scoring_result.profile_scores.at[scored_site, "K_AT"]
        != leave_one_out.scoring_result.profile_scores.at[scored_site, "K_AT"]
    )
    diagnostics = leave_one_out.scoring_result.profile_score_diagnostics
    assert diagnostics is not None
    at_rows = diagnostics.loc[diagnostics.loc[:, "kinase"] == "K_AT", :]
    assert set(at_rows.loc[:, "status"]) == {
        KINASE_PROFILE_SCORE_DIAGNOSTIC_STATUS_UNSCORED
    }
    assert set(at_rows.loc[:, "reason"]) == {
        KINASE_PROFILE_SCORE_DIAGNOSTIC_REASON_INSUFFICIENT_SUBSTRATES_AFTER_LEAVE_ONE_OUT
    }


def test_kinase_production_profile_threshold_failure_is_explicit() -> None:
    with pytest.raises(WorkflowValidationError) as exc_info:
        KinaseWorkflow().run(
            _request(
                scoring_config=KinaseScoringConfig.production(
                    minimum_reference_overlap_fraction=0.90,
                    minimum_sequence_supported_fraction=0.90,
                    minimum_scored_fraction=0.90,
                ),
                with_reference_context=True,
            )
        )

    message = str(exc_info.value)
    assert "policy=require_threshold" in message
    assert "localisation_confidence must be >= 0.750" in message
    assert "affected_rows=1" in message


def test_kinase_ssgsea_fixture_reports_optional_permutation_significance() -> None:
    effect_matrix = _ssgsea_effect_matrix()
    membership = _ssgsea_membership()

    no_permutation = SsgseaSubstrateEnrichmentActivityMethod(
        min_substrates=2,
        permutation_count=0,
        random_seed=19,
    ).run(
        effect_matrix=effect_matrix,
        kinase_substrate_membership=membership,
    )
    with_permutation = SsgseaSubstrateEnrichmentActivityMethod(
        min_substrates=2,
        permutation_count=12,
        random_seed=19,
    ).run(
        effect_matrix=effect_matrix,
        kinase_substrate_membership=membership,
    )

    assert no_permutation.p_value_matrix is None
    assert no_permutation.q_value_matrix is None
    assert set(no_permutation.statistics_table["significance_status"]) == {
        SSGSEA_SIGNIFICANCE_STATUS_UNAVAILABLE_NO_PERMUTATIONS
    }
    assert with_permutation.p_value_matrix is not None
    assert with_permutation.q_value_matrix is not None
    assert set(with_permutation.statistics_table["significance_status"]) == {
        SSGSEA_SIGNIFICANCE_STATUS_AVAILABLE
    }
    expected_support = _contracts()["ssgsea_expected_support"]
    for kinase, count in expected_support.items():
        assert int(with_permutation.substrate_count_matrix.loc[kinase].max()) == int(
            count
        )
