from __future__ import annotations

import pandas as pd

from phospy.science.prediction.models import KinaseScoringResult
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from phospy.workflows.kinase.row_attrition import (
    build_kinase_row_attrition_provenance,
)
from tests.unit.workflows.kinase.test_kinase_attrition_policy import (
    _request,
    _strict_scored_fraction_policy,
    _validated_request,
)


def _resolved_request():
    return KinaseWorkflowInterpreter().run(
        _validated_request(
            _request(_strict_scored_fraction_policy(on_violation="warn"))
        )
    )


def _scoring_result(scores: pd.DataFrame) -> KinaseScoringResult:
    return KinaseScoringResult(profile_scores=scores)


def test_kinase_scoring_retention_record_uses_actual_output_index() -> None:
    request = _resolved_request()
    scoring_site_index = request.scoring_site_index
    retained_index = scoring_site_index.delete(1)
    scores = pd.DataFrame(
        {"KINASE_A": [1.0, 0.5, 0.25]},
        index=retained_index,
    )

    provenance = build_kinase_row_attrition_provenance(
        request=request,
        scoring_result=_scoring_result(scores),
    )

    assert provenance.row_attrition is not None
    records = provenance.row_attrition.records
    assert len(records) == 1
    assert records[0].stage == "kinase_scoring_retention"
    assert records[0].reason == "sites_removed_by_scoring_retention"
    assert records[0].input_rows == int(scoring_site_index.size)
    assert records[0].output_rows == int(retained_index.size)
    assert records[0].examples == (str(scoring_site_index[1]),)
    assert provenance.row_attrition.final_rows == int(scores.shape[0])
    assert all(record.removed_rows > 0 for record in records)


def test_kinase_pair_attrition_is_not_encoded_as_site_row_attrition() -> None:
    request = _resolved_request()
    scores = pd.DataFrame(
        {
            "KINASE_A": [float("nan"), 0.7, 0.3, 0.1],
            "KINASE_B": [0.2, float("nan"), 0.4, 0.8],
        },
        index=request.scoring_site_index,
    )

    provenance = build_kinase_row_attrition_provenance(
        request=request,
        scoring_result=_scoring_result(scores),
    )

    assert provenance.row_attrition is None
    assert provenance.metrics["site_kinase_pairs_considered"] == 8
    assert provenance.metrics["site_kinase_pairs_scored"] == 6
    assert (
        provenance.metrics["site_kinase_pairs_unscored_due_to_insufficient_evidence"]
        == 2
    )
