from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phospy.api import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api.configs import (
    KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
    KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
    KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT,
    KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
    KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
    KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
    KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
)
from phospy.errors import WorkflowBoundaryError, WorkflowValidationError
from phospy.science.activities.method_contracts import (
    all_kinase_activity_method_quantitative_contracts,
    kinase_activity_method_quantitative_input_contract,
)
from phospy.science.activities.methods import (
    SimplifiedWeightedSubstrateActivityMethod,
    SsgseaSubstrateEnrichmentActivityMethod,
)
from phospy.science.activities.models import (
    KinaseActivityInputs,
    PredMatOverlapSummary,
)
from phospy.science.activities.semantics import (
    ActivityInputMatrix,
    ActivityProfileAxis,
    ActivityQuantitativeSemantics,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.quantitative_method_contracts import (
    MethodQuantitativeInputContract,
    method_contracts_to_markdown_table,
)
from phospy.science.transformations.models import (
    IntensityScaleKind,
    QuantitativeMeaning,
)
from phospy.validation.workflows.method_quantitative import (
    MethodQuantitativeInputValidator,
)
from phospy.workflows.kinase.science import (
    build_kinase_profiles,
    score_profile_correlations,
)
from phospy.workflows.kinase.scoring_mode_contracts import (
    all_kinase_scoring_method_quantitative_contracts,
    kinase_scoring_method_quantitative_input_contract,
)
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
    supported_log2_intensity_scale_state,
    supported_log2_intensity_scale_state_with_meaning,
    supported_log2_processing_state,
    supported_log2_processing_state_with_meaning,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)

ROOT = Path(__file__).resolve().parents[4]
KINASE_DOC = ROOT / "docs" / "api" / "kinase.md"

_PROFILE_SCORING_MODES = (
    KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
    KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
    KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
)
_ALL_SCORING_MODES = (
    *_PROFILE_SCORING_MODES,
    KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
)
_ACTIVITY_METHODS = (
    KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
    KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
    KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT,
)


@pytest.mark.parametrize("scoring_mode", _ALL_SCORING_MODES)
def test_scoring_methods_declare_quantitative_contracts(scoring_mode: str) -> None:
    contract = kinase_scoring_method_quantitative_input_contract(scoring_mode)
    payload = contract.to_payload()

    assert payload["method_id"] == f"kinase_scoring.{scoring_mode}"
    assert payload["accepted_scales"]
    assert payload["accepted_meanings"]
    assert payload["required_centring"]
    assert payload["required_standardisation"]
    assert payload["missing_value_treatment"]
    assert payload["profile_axis_requirements"]
    assert payload["statistical_interpretation"]
    assert payload["no_implicit_transformation"] is True


@pytest.mark.parametrize("scoring_mode", _PROFILE_SCORING_MODES)
def test_profile_scoring_contracts_accept_abundance_and_reject_effect_meanings(
    scoring_mode: str,
) -> None:
    contract = kinase_scoring_method_quantitative_input_contract(scoring_mode)

    assert IntensityScaleKind.LINEAR in contract.accepted_scales
    assert IntensityScaleKind.LOG2 in contract.accepted_scales
    assert QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE in contract.accepted_meanings
    assert QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE in contract.accepted_meanings
    assert QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO in contract.accepted_meanings
    assert (
        QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE not in contract.accepted_meanings
    )
    assert (
        QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE not in contract.accepted_meanings
    )
    assert QuantitativeMeaning.UNKNOWN not in contract.accepted_meanings
    assert contract.quantitative_input_required is True


def test_motif_only_scoring_contract_records_scale_without_consuming_values() -> None:
    contract = kinase_scoring_method_quantitative_input_contract(
        KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY
    )

    assert contract.quantitative_input_required is False
    assert QuantitativeMeaning.UNKNOWN in contract.accepted_meanings
    assert "not consumed" in contract.scale_sensitivity


@pytest.mark.parametrize("method", _ACTIVITY_METHODS)
def test_activity_methods_declare_quantitative_contracts(method: str) -> None:
    contract = kinase_activity_method_quantitative_input_contract(method)
    payload = contract.to_payload()

    assert payload["method_id"].endswith("_v1")
    assert payload["accepted_scales"]
    assert payload["accepted_meanings"]
    assert payload["accepted_activity_profile_axes"]
    assert payload["accepted_activity_quantitative_semantics"]
    assert payload["required_centring"]
    assert payload["required_standardisation"]
    assert payload["missing_value_treatment"]
    assert payload["profile_axis_requirements"]
    assert payload["statistical_interpretation"]
    assert payload["no_implicit_transformation"] is True


def test_weighted_activity_contract_accepts_abundance_but_not_effect_semantics() -> (
    None
):
    contract = kinase_activity_method_quantitative_input_contract(
        KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY
    )

    assert contract.accepted_scales == (
        IntensityScaleKind.LINEAR,
        IntensityScaleKind.LOG2,
    )
    assert QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE in contract.accepted_meanings
    assert (
        QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE not in contract.accepted_meanings
    )
    assert contract.accepted_activity_profile_axes == (
        ActivityProfileAxis.SAMPLE,
        ActivityProfileAxis.CONDITION_SUMMARY,
    )
    assert ActivityQuantitativeSemantics.SAMPLE_LEVEL_ABUNDANCE in (
        contract.accepted_activity_quantitative_semantics
    )
    assert ActivityQuantitativeSemantics.STANDARDISED_EFFECT not in (
        contract.accepted_activity_quantitative_semantics
    )


def test_ksea_activity_contract_accepts_log2_abundance_and_effect_inputs() -> None:
    contract = kinase_activity_method_quantitative_input_contract(
        KINASE_ACTIVITY_METHOD_KSEA_ZSCORE
    )

    assert contract.accepted_scales == (IntensityScaleKind.LOG2,)
    assert contract.accepted_meanings == (
        QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO,
        QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE,
        QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE,
    )
    assert contract.accepted_activity_profile_axes == (
        ActivityProfileAxis.SAMPLE,
        ActivityProfileAxis.CONTRAST,
        ActivityProfileAxis.EFFECT,
    )
    assert contract.accepted_activity_quantitative_semantics == (
        ActivityQuantitativeSemantics.SAMPLE_LEVEL_ABUNDANCE,
        ActivityQuantitativeSemantics.CONTRAST_LOG_FOLD_CHANGE,
        ActivityQuantitativeSemantics.STANDARDISED_EFFECT,
    )


def test_ssgsea_activity_contract_requires_log2_effect_inputs() -> None:
    contract = kinase_activity_method_quantitative_input_contract(
        KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT
    )

    assert contract.accepted_scales == (IntensityScaleKind.LOG2,)
    assert contract.accepted_meanings == (
        QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE,
        QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE,
    )
    assert contract.accepted_activity_profile_axes == (
        ActivityProfileAxis.CONTRAST,
        ActivityProfileAxis.EFFECT,
    )
    assert contract.accepted_activity_quantitative_semantics == (
        ActivityQuantitativeSemantics.CONTRAST_LOG_FOLD_CHANGE,
        ActivityQuantitativeSemantics.STANDARDISED_EFFECT,
    )


@pytest.mark.parametrize(
    ("dataset", "method"),
    [
        (
            lambda: _linear_dataset(),
            KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
        ),
        (
            lambda: _log2_abundance_dataset(),
            KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
        ),
        (lambda: _log2_abundance_dataset(), KINASE_ACTIVITY_METHOD_KSEA_ZSCORE),
        (lambda: _effect_dataset(), KINASE_ACTIVITY_METHOD_KSEA_ZSCORE),
        (
            lambda: _effect_dataset(),
            KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT,
        ),
    ],
    ids=[
        "weighted-linear-abundance",
        "weighted-log2-abundance",
        "ksea-log2-abundance",
        "ksea-log2-effect",
        "ssgsea-log2-effect",
    ],
)
def test_activity_method_contract_validator_accepts_declared_inputs(
    dataset: Callable[[], AnalysisReadyPhosphoDataset],
    method: str,
) -> None:
    resolved = MethodQuantitativeInputValidator().run(
        dataset=dataset(),
        contract=kinase_activity_method_quantitative_input_contract(method),
        context="test activity contract",
    )

    assert resolved.to_payload()["method_id"].endswith("_v1")


@pytest.mark.parametrize(
    ("dataset", "method", "message"),
    [
        (
            lambda: _linear_dataset(),
            KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
            "requires quantitative meaning",
        ),
        (
            lambda: _linear_dataset(),
            KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT,
            "requires quantitative meaning",
        ),
        (
            lambda: _effect_dataset(),
            KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
            "requires quantitative meaning",
        ),
    ],
    ids=[
        "ksea-rejects-linear-abundance",
        "ssgsea-rejects-linear-abundance",
        "weighted-rejects-effect",
    ],
)
def test_activity_method_contract_validator_rejects_invalid_scale_or_meaning(
    dataset: Callable[[], AnalysisReadyPhosphoDataset],
    method: str,
    message: str,
) -> None:
    with pytest.raises(WorkflowValidationError, match=message):
        MethodQuantitativeInputValidator().run(
            dataset=dataset(),
            contract=kinase_activity_method_quantitative_input_contract(method),
            context="test activity contract",
        )


def test_invalid_activity_scale_meaning_fails_at_workflow_validator_before_execution() -> (
    None
):
    request = _request(
        dataset=_linear_dataset(),
        activity_config=KinaseActivityConfig(
            method=KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
            ksea_min_substrates=2,
        ),
    )

    with pytest.raises(WorkflowValidationError) as exc_info:
        KinaseWorkflowValidator().run(request)

    message = str(exc_info.value)
    assert "kinase activity request dataset" in message
    assert "ksea_zscore_v1" in message
    assert "differential_effect_size" in message
    assert "phosphosite_abundance" in message


@pytest.mark.parametrize(
    "contract",
    (
        *all_kinase_scoring_method_quantitative_contracts(),
        *all_kinase_activity_method_quantitative_contracts(),
    ),
    ids=lambda contract: contract.method_id,
)
def test_method_contracts_declare_missingness_policy_without_imputation(
    contract: MethodQuantitativeInputContract,
) -> None:
    payload = contract.to_payload()
    policy = str(payload["missing_value_treatment"]).lower()

    assert "missing" in policy
    assert "imputation" in policy
    assert "no_implicit_transformation" in payload


def test_activity_p_value_interpretation_metadata_is_method_specific() -> None:
    weighted = kinase_activity_method_quantitative_input_contract(
        KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY
    ).to_payload()
    ksea = kinase_activity_method_quantitative_input_contract(
        KINASE_ACTIVITY_METHOD_KSEA_ZSCORE
    ).to_payload()
    ssgsea = kinase_activity_method_quantitative_input_contract(
        KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT
    ).to_payload()

    assert weighted["p_value_interpretation"] is None
    assert "normal-approximation p-values" in str(ksea["p_value_interpretation"])
    assert "Benjamini-Hochberg" in str(ksea["p_value_interpretation"])
    assert "No p-values are produced unless seeded permutations are requested" in str(
        ssgsea["p_value_interpretation"]
    )
    assert "permutation p-values" in str(ssgsea["p_value_interpretation"])


def test_activity_methods_reject_interchangeable_typed_input_semantics() -> None:
    inputs = _activity_inputs(
        ActivityInputMatrix.standardised_effect(_activity_matrix(), _assume_owned=True)
    )

    with pytest.raises(WorkflowBoundaryError, match="no implicit transformation"):
        SimplifiedWeightedSubstrateActivityMethod(
            threshold=0.5,
            min_substrates=2,
            top_n_substrates=2,
        ).run(inputs)

    abundance_inputs = _activity_inputs(
        ActivityInputMatrix.sample_level_abundance(
            _activity_matrix(),
            _assume_owned=True,
        )
    )
    with pytest.raises(WorkflowBoundaryError, match="no implicit transformation"):
        SsgseaSubstrateEnrichmentActivityMethod(
            min_substrates=2,
            permutation_count=0,
        ).run(
            activity_input=abundance_inputs.activity_input,
            kinase_substrate_membership=_activity_membership(),
        )


def test_linear_and_log2_profile_scores_are_not_treated_as_interchangeable() -> None:
    sites = pd.Index(["s1", "s2", "s3"], name="site_key")
    linear = pd.DataFrame(
        {
            "profile_a": [1.0, 2.0, 4.0],
            "profile_b": [8.0, 4.0, 2.0],
            "profile_c": [3.0, 6.0, 12.0],
        },
        index=sites,
    )
    log2 = np.log2(linear)
    kinase_substrate_map = pd.DataFrame(
        {
            "kinase": ["K1", "K1", "K2", "K2"],
            "substrate_site": ["s1", "s2", "s2", "s3"],
        }
    )

    linear_profiles = build_kinase_profiles(
        phospho=linear,
        kinase_substrate_map=kinase_substrate_map,
        min_substrates=2,
    ).profile_matrix
    log2_profiles = build_kinase_profiles(
        phospho=log2,
        kinase_substrate_map=kinase_substrate_map,
        min_substrates=2,
    ).profile_matrix
    linear_scores = score_profile_correlations(
        phospho=linear,
        profile_matrix=linear_profiles,
    )
    log2_scores = score_profile_correlations(
        phospho=log2,
        profile_matrix=log2_profiles,
    )

    assert not np.allclose(
        linear_scores.to_numpy(),
        log2_scores.to_numpy(),
        equal_nan=True,
    )


@pytest.mark.parametrize(
    ("dataset", "expected_scale", "expected_meaning"),
    [
        (lambda: _linear_dataset(), "linear", "phosphosite_abundance"),
        (lambda: _log2_abundance_dataset(), "log2", "phosphosite_log_abundance"),
    ],
    ids=["linear-provenance", "log2-provenance"],
)
def test_result_provenance_records_resolved_scoring_method_contract(
    dataset: Callable[[], AnalysisReadyPhosphoDataset],
    expected_scale: str,
    expected_meaning: str,
) -> None:
    result = KinaseWorkflow().run(_request(dataset=dataset()))

    payload = result.provenance.workflow_parameters["scoring_config"][
        "method_input_contract"
    ]
    assert payload["method_id"] == "kinase_scoring.phosr_rank_weighted"
    assert payload["resolved_scale"] == expected_scale
    assert payload["resolved_meaning"] == expected_meaning
    assert payload["no_implicit_transformation"] is True
    assert "not numerically interchangeable" in payload["statistical_interpretation"]


def test_result_provenance_records_resolved_activity_method_contract() -> None:
    result = KinaseWorkflow().run(
        _request(
            dataset=_log2_abundance_dataset(),
            activity_config=KinaseActivityConfig(
                method=KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
                ksea_min_substrates=2,
            ),
        )
    )

    payload = result.provenance.workflow_parameters["activity_config"][
        "method_input_contract"
    ]
    assert payload["method_id"] == "ksea_zscore_v1"
    assert payload["resolved_scale"] == "log2"
    assert payload["resolved_meaning"] == "phosphosite_log_abundance"
    assert payload["resolved_activity_profile_axis"] == "sample"
    assert (
        payload["resolved_activity_quantitative_semantics"] == "sample_level_abundance"
    )
    assert "normal-approximation p-values" in payload["p_value_interpretation"]


def test_kinase_docs_quantitative_contract_table_matches_method_declarations() -> None:
    expected_table = method_contracts_to_markdown_table(
        (
            *all_kinase_scoring_method_quantitative_contracts(),
            *all_kinase_activity_method_quantitative_contracts(),
        )
    )
    docs = KINASE_DOC.read_text(encoding="utf-8")

    assert expected_table in docs


def _linear_dataset() -> AnalysisReadyPhosphoDataset:
    return _dataset_from_phospho(
        pd.DataFrame(
            {
                "profile_a": [1.0, 2.0],
                "profile_b": [1.5, 2.5],
            },
            index=_site_index(),
        ),
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _log2_abundance_dataset() -> AnalysisReadyPhosphoDataset:
    linear = pd.DataFrame(
        {
            "profile_a": [1.0, 2.0],
            "profile_b": [1.5, 2.5],
        },
        index=_site_index(),
    )
    return _dataset_from_phospho(
        np.log2(linear),
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


def _effect_dataset() -> AnalysisReadyPhosphoDataset:
    return _dataset_from_phospho(
        pd.DataFrame(
            {
                "contrast_a": [1.0, -0.5],
                "contrast_b": [0.25, 1.25],
            },
            index=_site_index(),
        ),
        intensity_scale_state=supported_log2_intensity_scale_state_with_meaning(
            has_total_matrix=False,
            meaning=QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE,
        ),
        processing_state=supported_log2_processing_state_with_meaning(
            has_total_matrix=False,
            meaning=QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE,
        ),
    )


def _dataset_from_phospho(
    phospho: pd.DataFrame,
    *,
    intensity_scale_state: object,
    processing_state: object,
) -> AnalysisReadyPhosphoDataset:
    display_ids = ["GENE1;S10;", "GENE2;T20;"]
    site_index = _site_index()
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": display_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": ["GENE1", "GENE2"],
            "site": ["S10", "T20"],
            "protein_id": ["GENE1", "GENE2"],
            "site_sequence": [_window("S"), _window("T")],
        },
        index=site_index.copy(),
    )
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=intensity_scale_state,
        processing_state=processing_state,
    )


def _site_index() -> pd.Index:
    return site_key_index_from_display_ids(("GENE1;S10;", "GENE2;T20;"))


def _references() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["KPROFILE", "KPROFILE"],
                "substrate_site": ["GENE1;S10;", "GENE2;T20;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": [_window("S"), _window("T")]},
            index=pd.Index(["GENE1;S10;", "GENE2;T20;"], name="site_id"),
        ),
    )


def _request(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    activity_config: KinaseActivityConfig | None = None,
) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=dataset,
        references=_references(),
        scoring_config=KinaseScoringConfig(
            reliability_profile="custom",
            min_substrates=2,
            scoring_mode=KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
            mode="deterministic_ranking",
        ),
        activity_config=activity_config,
    )


def _activity_inputs(activity_input: ActivityInputMatrix) -> KinaseActivityInputs:
    frame = activity_input.frame
    pred_mat = pd.DataFrame(
        {"K1": [0.8, 0.7]},
        index=frame.index.copy(),
    )
    return KinaseActivityInputs(
        pred_mat=pred_mat,
        phospho_matrix=frame.copy(deep=True),
        threshold=0.5,
        min_substrates=2,
        top_n_substrates=2,
        overlap_summary=PredMatOverlapSummary(
            overlap_count=int(frame.shape[0]),
            pred_mat_rows=int(frame.shape[0]),
            phospho_rows=int(frame.shape[0]),
        ),
        activity_input=activity_input,
    )


def _activity_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "profile_a": [1.0, np.nan],
            "profile_b": [2.0, 3.0],
        },
        index=_site_index(),
    )


def _activity_membership() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "kinase": ["K1", "K1"],
            "substrate_site": _site_index().astype(str).tolist(),
        }
    )


def _window(residue: str) -> str:
    return ("A" * 7) + residue + ("A" * 7)
