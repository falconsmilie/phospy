from __future__ import annotations

from dataclasses import FrozenInstanceError, asdict

import pytest

from phospy.api import EnrichmentConfig
from phospy.api.configs import ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
from phospy.api.results import EnrichmentWorkflowResult, ResultCaveat
from phospy.science.enrichment.models import GeneSetCollection


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
