from __future__ import annotations

from dataclasses import FrozenInstanceError, asdict

import pandas as pd
import pytest

from phospy.api import (
    Contrast,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    EnrichmentConfig,
    ExperimentalDesign,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferenceContextCompatibilityPolicy,
    SampleDesignRecord,
)
from phospy.api.configs import ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
from phospy.api.results import (
    EnrichmentWorkflowResult,
    KinaseWorkflowResult,
    ResultCaveat,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.enrichment.models import GeneSetCollection
from phospy.workflows.differential.caveats import (
    DIFFERENTIAL_DIRECT_TRUSTED_DATASET_CAVEAT_CODE,
)
from phospy.workflows.kinase.caveats import (
    KINASE_DIRECT_TRUSTED_DATASET_CAVEAT_CODE,
    KINASE_SCORING_LIMITATION_CAVEAT_CODE,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

_GENES = ["MAPK14", "AKT1", "GSK3B", "RPS6"]
_SITES = ["Y182", "T308", "S9", "S235"]
_SAMPLES = ("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")


def _trusted_site_index() -> pd.Index:
    return protein_site_key_index(protein_identifiers=_GENES, sites=_SITES)


def _trusted_display_ids() -> list[str]:
    return [f"{gene};{site};" for gene, site in zip(_GENES, _SITES, strict=True)]


def _trusted_phospho(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "A_1": [10.0, 20.0, 5.0, 30.0],
            "A_2": [10.1, 20.2, 5.1, 30.2],
            "A_3": [9.9, 20.1, 5.2, 30.1],
            "B_1": [15.0, 24.0, 6.0, 29.0],
            "B_2": [15.1, 24.1, 6.2, 29.2],
            "B_3": [14.9, 24.2, 6.1, 29.1],
        },
        index=index.copy(),
    )


def _trusted_site_metadata(index: pd.Index) -> pd.DataFrame:
    display_ids = _trusted_display_ids()
    return pd.DataFrame(
        {
            "site_key": index.astype(str).tolist(),
            "display_id": display_ids,
            **site_key_context_columns(index),
            "gene_symbol": _GENES,
            "protein_id": _GENES,
            "site": _SITES,
            "site_sequence": [("A" * 15) + site[0] + ("A" * 15) for site in _SITES],
            "localisation_confidence": [0.95, 0.9, 0.92, 0.93],
        },
        index=index.copy(),
    )


def _trusted_dataset_without_assertion_metadata() -> AnalysisReadyPhosphoDataset:
    index = _trusted_site_index()
    return AnalysisReadyPhosphoDataset(
        phospho=_trusted_phospho(index),
        site_metadata=_trusted_site_metadata(index),
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


def _differential_design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=sample_id.split("_", maxsplit=1)[0],
                biological_replicate_id=sample_id,
            )
            for sample_id in _SAMPLES
        )
    )


def _kinase_references() -> ReferenceBundle:
    display_ids = _trusted_display_ids()
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["KINASE_A", "KINASE_A", "KINASE_B", "KINASE_B"],
                "substrate_site": display_ids,
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": [("A" * 15) + site[0] + ("A" * 15) for site in _SITES]},
            index=pd.Index(display_ids, name="site_id"),
        ),
    )


def _caveat_by_code(
    caveats: tuple[ResultCaveat, ...],
    code: str,
) -> ResultCaveat:
    matches = [caveat for caveat in caveats if caveat.code == code]
    assert len(matches) == 1
    return matches[0]


def test_result_caveat_is_immutable() -> None:
    caveat = ResultCaveat(
        code="low_scored_fraction",
        severity="warning",
        message="Only half of the input sites contributed to scoring.",
        details={"observed_fraction": 0.5},
    )

    with pytest.raises(FrozenInstanceError):
        caveat.code = "mutated"  # type: ignore[misc]

    with pytest.raises(TypeError):
        caveat.details["observed_fraction"] = 0.9  # type: ignore[index]


def test_workflow_result_defaults_to_empty_caveats() -> None:
    result = EnrichmentWorkflowResult(
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=GeneSetCollection(
            sets={"mapk_pathway": ("AKT1", "MAPK1")},
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        ),
        config=EnrichmentConfig(),
    )

    assert result.caveats == ()


def test_result_caveat_contains_code_severity_message_and_details() -> None:
    caveat = ResultCaveat(
        code="minimum_scored_fraction_not_met",
        severity="warning",
        message="The scored site fraction is below the configured threshold.",
        details={
            "threshold_name": "minimum_scored_fraction",
            "configured_threshold": 0.75,
            "observed_value": 0.5,
        },
    )

    assert caveat.code == "minimum_scored_fraction_not_met"
    assert caveat.severity == "warning"
    assert caveat.message == (
        "The scored site fraction is below the configured threshold."
    )
    assert caveat.details["threshold_name"] == "minimum_scored_fraction"
    assert caveat.to_payload() == {
        "code": "minimum_scored_fraction_not_met",
        "severity": "warning",
        "message": "The scored site fraction is below the configured threshold.",
        "details": {
            "threshold_name": "minimum_scored_fraction",
            "configured_threshold": 0.75,
            "observed_value": 0.5,
        },
    }
    assert asdict(caveat) == caveat.to_payload()


def test_differential_result_caveat_warns_for_missing_trusted_assertions() -> None:
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_trusted_dataset_without_assertion_metadata(),
            design=_differential_design(),
            contrasts=(
                Contrast(
                    name="B_vs_A",
                    numerator_condition="B",
                    denominator_condition="A",
                ),
            ),
        )
    )

    caveat = _caveat_by_code(
        result.caveats,
        DIFFERENTIAL_DIRECT_TRUSTED_DATASET_CAVEAT_CODE,
    )

    assert caveat.severity == "warning"
    assert caveat.details["trusted_assertion_metadata_provided"] is False
    assert caveat.details["missing_trusted_assertions"] == [
        "identity_user_asserted",
        "quantitative_meaning_user_asserted",
        "localisation_user_asserted",
        "sequence_user_asserted",
        "reference_context_user_asserted",
    ]


def test_kinase_result_caveat_warns_for_missing_trusted_assertions() -> None:
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_trusted_dataset_without_assertion_metadata(),
            references=_kinase_references(),
            scoring_config=KinaseScoringConfig(
                min_substrates=2,
                reference_context_compatibility_policy=(
                    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                ),
            ),
            prediction_config=KinasePredictionConfig(
                top_k=2,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
            ),
            activity_config=KinaseActivityConfig(enabled=False),
        )
    )
    assert isinstance(result, KinaseWorkflowResult)

    caveat = _caveat_by_code(result.caveats, KINASE_DIRECT_TRUSTED_DATASET_CAVEAT_CODE)

    assert caveat.severity == "warning"
    assert caveat.details["trusted_assertion_metadata_provided"] is False
    assert caveat.details["missing_trusted_assertions"] == [
        "identity_user_asserted",
        "quantitative_meaning_user_asserted",
        "localisation_user_asserted",
        "sequence_user_asserted",
        "reference_context_user_asserted",
    ]

    scoring_caveat = _caveat_by_code(
        result.caveats,
        KINASE_SCORING_LIMITATION_CAVEAT_CODE,
    )
    assert scoring_caveat.severity == "info"
    assert "relative support values" in scoring_caveat.message
    assert "not causal kinase activity proof" in scoring_caveat.message
    assert scoring_caveat.details["score_interpretation"] == (
        "relative_support_within_run"
    )
    assert scoring_caveat.details["not_causal_activity_proof"] is True
    assert scoring_caveat.details["not_calibrated_probability"] is True
