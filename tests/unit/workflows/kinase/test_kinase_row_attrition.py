from __future__ import annotations

from types import SimpleNamespace

import pandas as pd

from phospy.contracts.configs import (
    KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
    KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
    LocalisationRequirement,
)
from phospy.workflows.kinase.row_attrition import (
    build_kinase_row_attrition_provenance,
)


def _window(residue: str = "S") -> str:
    return ("A" * 15) + residue + ("A" * 15)


def _request(
    *,
    site_ids: tuple[str, ...] = ("S1", "S2", "S3"),
    scoring_site_ids: tuple[str, ...] | None = None,
    reference_site_ids: tuple[str, ...] | None = None,
    sequences: dict[str, object] | None = None,
    localisation: dict[str, object] | None = None,
    localisation_requirement: LocalisationRequirement | None = None,
    scoring_mode: str = KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
) -> SimpleNamespace:
    scoring_site_ids = site_ids if scoring_site_ids is None else scoring_site_ids
    reference_site_ids = (
        scoring_site_ids if reference_site_ids is None else reference_site_ids
    )
    sequence_values = {site_id: _window("S") for site_id in site_ids} | (
        sequences or {}
    )
    localisation_values = {site_id: 0.95 for site_id in site_ids} | (localisation or {})
    index = pd.Index(site_ids, name="site_key")
    return SimpleNamespace(
        dataset=SimpleNamespace(
            phospho=pd.DataFrame({"sample": range(len(site_ids))}, index=index),
            site_metadata=pd.DataFrame(
                {
                    "site": ["S1" for _ in site_ids],
                    "site_sequence": [sequence_values[site_id] for site_id in site_ids],
                    "localisation_probability": [
                        localisation_values[site_id] for site_id in site_ids
                    ],
                },
                index=index,
            ),
        ),
        scoring_site_index=pd.Index(scoring_site_ids, name="site_key"),
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1" for _ in reference_site_ids],
                "substrate_site": list(reference_site_ids),
            }
        ),
        execution_config=SimpleNamespace(
            scoring_mode=scoring_mode,
            localisation_requirement=(
                LocalisationRequirement()
                if localisation_requirement is None
                else localisation_requirement
            ),
        ),
        attrition_metrics=None,
    )


def _scoring_result(
    scores: pd.DataFrame,
) -> SimpleNamespace:
    return SimpleNamespace(
        authoritative_scores=scores,
        profile_score_diagnostics=None,
    )


def _records(provenance) -> tuple[dict[str, object], ...]:
    assert provenance.row_attrition is not None
    return tuple(record.to_payload() for record in provenance.row_attrition.records)


def test_kinase_row_attrition_records_sequence_context_drops() -> None:
    request = _request(
        scoring_site_ids=("S1", "S3"),
        reference_site_ids=("S1", "S3"),
        sequences={"S2": ""},
    )
    result = build_kinase_row_attrition_provenance(
        request=request,
        scoring_result=_scoring_result(
            pd.DataFrame({"K1": [1.0, 0.5]}, index=["S1", "S3"])
        ),
    )

    records = _records(result)

    assert records[0]["stage"] == "kinase_sequence_context"
    assert records[0]["reason"] == "sites_missing_valid_centered_sequence"
    assert records[0]["removed_rows"] == 1
    assert records[0]["examples"] == ["S2"]


def test_kinase_row_attrition_records_localisation_threshold_drops_when_filtering() -> (
    None
):
    request = _request(
        scoring_site_ids=("S1", "S3"),
        reference_site_ids=("S1", "S3"),
        localisation={"S2": 0.2},
        localisation_requirement=LocalisationRequirement(minimum_probability=0.75),
    )
    result = build_kinase_row_attrition_provenance(
        request=request,
        scoring_result=_scoring_result(
            pd.DataFrame({"K1": [1.0, 0.5]}, index=["S1", "S3"])
        ),
    )

    records = _records(result)

    assert records[0]["stage"] == "kinase_localisation_filter"
    assert records[0]["reason"] == "sites_below_localisation_threshold"
    assert records[0]["removed_rows"] == 1
    assert records[0]["examples"] == ["S2"]


def test_kinase_motif_only_does_not_record_reference_overlap_drops() -> None:
    request = _request(
        reference_site_ids=("S1",),
        scoring_mode=KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
    )
    result = build_kinase_row_attrition_provenance(
        request=request,
        scoring_result=_scoring_result(
            pd.DataFrame({"K1": [1.0, 0.5, 0.25]}, index=["S1", "S2", "S3"])
        ),
    )

    assert result.row_attrition is None
    assert result.metrics["sites_not_present_in_reference_resource"] == 0


def test_kinase_row_attrition_records_sites_without_any_authoritative_score() -> None:
    request = _request(reference_site_ids=("S1", "S2", "S3"))
    result = build_kinase_row_attrition_provenance(
        request=request,
        scoring_result=_scoring_result(
            pd.DataFrame(
                {"K1": [1.0, float("nan"), 0.5], "K2": [0.2, float("nan"), 0.1]},
                index=["S1", "S2", "S3"],
            )
        ),
    )

    records = _records(result)

    assert records[0]["stage"] == "kinase_scoring_retention"
    assert records[0]["reason"] == "sites_without_any_authoritative_score"
    assert records[0]["examples"] == ["S2"]


def test_kinase_pair_attrition_is_not_encoded_as_site_row_attrition() -> None:
    request = _request(site_ids=("S1", "S2"), reference_site_ids=("S1", "S2"))
    result = build_kinase_row_attrition_provenance(
        request=request,
        scoring_result=_scoring_result(
            pd.DataFrame(
                {"K1": [float("nan"), 0.5], "K2": [0.2, float("nan")]},
                index=["S1", "S2"],
            )
        ),
    )

    assert result.row_attrition is None
    assert result.metrics["site_kinase_pairs_considered"] == 4
    assert result.metrics["site_kinase_pairs_scored"] == 2
    assert (
        result.metrics["site_kinase_pairs_unscored_due_to_insufficient_evidence"] == 2
    )


def test_kinase_row_attrition_report_counts_are_continuous() -> None:
    request = _request(
        site_ids=("S1", "S2", "S3", "S4"),
        scoring_site_ids=("S1", "S3", "S4"),
        reference_site_ids=("S1", "S4"),
        sequences={"S2": ""},
    )
    result = build_kinase_row_attrition_provenance(
        request=request,
        scoring_result=_scoring_result(
            pd.DataFrame(
                {"K1": [1.0, 0.5, float("nan")]},
                index=["S1", "S3", "S4"],
            )
        ),
    )

    records = result.row_attrition.records

    assert [(record.input_rows, record.output_rows) for record in records] == [
        (4, 3),
        (3, 2),
        (2, 1),
    ]
    assert result.row_attrition.input_rows == 4
    assert result.row_attrition.final_rows == 1
