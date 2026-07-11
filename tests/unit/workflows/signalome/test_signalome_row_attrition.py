from __future__ import annotations

from types import SimpleNamespace

import pandas as pd

from phospy.contracts.configs import (
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    LocalisationRequirement,
)
from phospy.science.signalomes.models import SignalomeScorePreconditioningDiagnostics
from phospy.workflows.signalome.row_attrition import (
    build_signalome_row_attrition_provenance,
)


def _request(
    *,
    site_ids: tuple[str, ...] = ("S1", "S2", "S3"),
    retained_site_ids: tuple[str, ...] = ("S1", "S3"),
    prediction_site_ids: tuple[str, ...] | None = None,
    score_site_ids: tuple[str, ...] | None = None,
    sequences: dict[str, object] | None = None,
    localisation: dict[str, object] | None = None,
    proteins: dict[str, object] | None = None,
    score_values: dict[str, tuple[float, ...]] | None = None,
    preconditioning_drop_count: int = 0,
    localisation_requirement: LocalisationRequirement | None = None,
) -> SimpleNamespace:
    prediction_site_ids = (
        site_ids if prediction_site_ids is None else prediction_site_ids
    )
    score_site_ids = site_ids if score_site_ids is None else score_site_ids
    sequence_values = {
        site_id: ("A" * 15) + "S" + ("A" * 15) for site_id in site_ids
    } | (sequences or {})
    localisation_values = {site_id: 0.95 for site_id in site_ids} | (localisation or {})
    protein_values = {
        site_id: f"P{index}" for index, site_id in enumerate(site_ids)
    } | (proteins or {})
    score_values = score_values or {
        site_id: (1.0 + index,) for index, site_id in enumerate(score_site_ids)
    }
    raw_scores = pd.DataFrame.from_dict(
        {site_id: score_values[site_id] for site_id in score_site_ids},
        orient="index",
        columns=["K1"],
    )
    raw_scores.index.name = "site_key"
    prediction_matrix = pd.DataFrame(
        {"K1": [0.5 for _ in prediction_site_ids]},
        index=pd.Index(prediction_site_ids, name="site_key"),
    )
    retained_scores = raw_scores.reindex(index=list(retained_site_ids))
    input_row_count = len(set(prediction_site_ids).intersection(set(score_site_ids)))
    return SimpleNamespace(
        dataset=SimpleNamespace(
            phospho=pd.DataFrame(
                {"sample": range(len(site_ids))},
                index=pd.Index(site_ids, name="site_key"),
            ),
            site_metadata=pd.DataFrame(
                {
                    "site_sequence": [sequence_values[site_id] for site_id in site_ids],
                    "localisation_probability": [
                        localisation_values[site_id] for site_id in site_ids
                    ],
                    "protein_id": [protein_values[site_id] for site_id in site_ids],
                },
                index=pd.Index(site_ids, name="site_key"),
            ),
        ),
        kinase_result=SimpleNamespace(
            prediction_result=SimpleNamespace(pred_mat=prediction_matrix),
            scoring_result=SimpleNamespace(authoritative_scores=raw_scores),
        ),
        downstream_score_matrix=retained_scores,
        execution_config=SimpleNamespace(
            localisation_requirement=(
                LocalisationRequirement()
                if localisation_requirement is None
                else localisation_requirement
            )
        ),
        score_preconditioning_diagnostics=SignalomeScorePreconditioningDiagnostics(
            input_row_count=input_row_count,
            dropped_all_missing_row_count=preconditioning_drop_count,
            retained_row_count=len(retained_site_ids),
            policy=SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
        ),
    )


def _records(provenance) -> tuple[dict[str, object], ...]:
    assert provenance.row_attrition is not None
    return tuple(record.to_payload() for record in provenance.row_attrition.records)


def test_signalome_row_attrition_records_sequence_context_drops() -> None:
    request = _request(sequences={"S2": ""})
    result = build_signalome_row_attrition_provenance(request)

    records = _records(result)

    assert records[0]["stage"] == "signalome_sequence_context"
    assert records[0]["reason"] == "sites_missing_sequence_context"
    assert records[0]["examples"] == ["S2"]


def test_signalome_row_attrition_records_localisation_drops_when_applicable() -> None:
    request = _request(
        localisation={"S2": 0.2},
        localisation_requirement=LocalisationRequirement(minimum_probability=0.75),
    )
    result = build_signalome_row_attrition_provenance(request)

    records = _records(result)

    assert records[0]["stage"] == "signalome_localisation_metadata"
    assert records[0]["reason"] == "sites_below_localisation_threshold"
    assert records[0]["examples"] == ["S2"]


def test_signalome_row_attrition_records_missing_protein_grouping() -> None:
    request = _request(proteins={"S2": ""})
    result = build_signalome_row_attrition_provenance(request)

    records = _records(result)

    assert records[0]["stage"] == "signalome_protein_grouping"
    assert records[0]["reason"] == "sites_missing_protein_grouping_metadata"
    assert records[0]["examples"] == ["S2"]


def test_signalome_row_attrition_records_score_preconditioning_drops() -> None:
    request = _request(
        score_values={"S1": (1.0,), "S2": (float("nan"),), "S3": (0.5,)},
        preconditioning_drop_count=1,
    )
    result = build_signalome_row_attrition_provenance(request)

    records = _records(result)

    assert records[0]["stage"] == "signalome_score_preconditioning"
    assert records[0]["reason"] == "sites_removed_by_score_preconditioning"
    assert records[0]["examples"] == ["S2"]


def test_signalome_row_attrition_does_not_double_count_sites() -> None:
    request = _request(
        sequences={"S2": ""},
        score_values={"S1": (1.0,), "S2": (float("nan"),), "S3": (0.5,)},
        preconditioning_drop_count=1,
    )
    result = build_signalome_row_attrition_provenance(request)

    records = _records(result)

    assert [record["stage"] for record in records] == [
        "signalome_score_preconditioning"
    ]
    assert sum(int(record["removed_rows"]) for record in records) == 1
    assert result.metrics["sites_removed_by_score_preconditioning"] == 1


def test_signalome_row_attrition_report_counts_are_continuous() -> None:
    request = _request(
        site_ids=("S1", "S2", "S3", "S4", "S5"),
        retained_site_ids=("S1",),
        prediction_site_ids=("S1", "S4", "S5"),
        score_site_ids=("S1", "S3", "S4", "S5"),
        sequences={"S2": ""},
        score_values={
            "S1": (1.0,),
            "S3": (0.5,),
            "S4": (float("nan"),),
            "S5": (0.2,),
        },
        preconditioning_drop_count=1,
    )
    result = build_signalome_row_attrition_provenance(request)

    records = result.row_attrition.records

    assert [(record.input_rows, record.output_rows) for record in records] == [
        (5, 3),
        (3, 2),
        (2, 1),
    ]
    assert [record.stage for record in records] == [
        "signalome_site_alignment",
        "signalome_score_preconditioning",
        "signalome_scoring_clustering_retention",
    ]
    assert result.row_attrition.final_rows == 1
