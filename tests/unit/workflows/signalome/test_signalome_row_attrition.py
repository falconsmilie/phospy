from __future__ import annotations

import pandas as pd

from phospy.api import (
    Organism,
    ReferenceBundle,
    SignalomeWorkflowRequest,
)
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.workflows.signalome.interpreter import SignalomeWorkflowInterpreter
from phospy.workflows.signalome.row_attrition import (
    build_signalome_row_attrition_provenance,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.signalome_config import build_signalome_config
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)


def _window(display_id: str) -> str:
    residue = display_id.split(";")[1][0].upper()
    return ("A" * 15) + residue + ("A" * 15)


def _dataset(display_ids: list[str]) -> AnalysisReadyPhosphoDataset:
    site_index = site_key_index_from_display_ids(display_ids)
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {
                "sample_a": [1.0 + index for index, _ in enumerate(display_ids)],
                "sample_b": [2.0 + index for index, _ in enumerate(display_ids)],
            },
            index=site_index.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_index.astype(str).tolist(),
                "display_id": display_ids,
                **site_key_context_columns(site_index),
                "gene_symbol": [
                    display_id.split(";", 1)[0] for display_id in display_ids
                ],
                "protein_id": [
                    display_id.split(";", 1)[0] for display_id in display_ids
                ],
                "site": [display_id.split(";")[1] for display_id in display_ids],
                "site_sequence": [_window(display_id) for display_id in display_ids],
            },
            index=site_index.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _matrix(
    *,
    index: pd.Index,
    values: list[list[float]],
) -> pd.DataFrame:
    return pd.DataFrame(values, index=index.copy(), columns=["KINASE_A", "KINASE_B"])


def _kinase_result(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    prediction_matrix: pd.DataFrame,
    score_matrix: pd.DataFrame,
) -> KinaseWorkflowResult:
    display_ids = dataset.site_metadata.loc[:, "display_id"].astype(str).tolist()
    return KinaseWorkflowResult(
        dataset=dataset,
        references=ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {
                    "kinase": ["KINASE_A" for _ in display_ids],
                    "substrate_site": display_ids,
                }
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": [_window(display_id) for display_id in display_ids]},
                index=pd.Index(display_ids, name="site_id"),
            ),
        ),
        scoring_result=KinaseScoringResult(
            profile_scores=score_matrix,
            rank_weighted_fusion_scores=score_matrix,
        ),
        prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
        activity_result=None,
    )


def _interpreted_request(
    *,
    display_ids: list[str],
    prediction_site_index: pd.Index | None = None,
    score_site_index: pd.Index | None = None,
    score_values: list[list[float]] | None = None,
):
    dataset = _dataset(display_ids)
    site_index = dataset.phospho.index
    prediction_index = (
        site_index if prediction_site_index is None else prediction_site_index
    )
    score_index = site_index if score_site_index is None else score_site_index
    prediction_matrix = _matrix(
        index=prediction_index,
        values=[[0.8, 0.2] for _ in range(int(prediction_index.size))],
    )
    scores = _matrix(
        index=score_index,
        values=(
            score_values
            if score_values is not None
            else [[1.0, 0.5] for _ in range(int(score_index.size))]
        ),
    )
    return SignalomeWorkflowInterpreter().run(
        SignalomeWorkflowRequest(
            kinase_result=_kinase_result(
                dataset=dataset,
                prediction_matrix=prediction_matrix,
                score_matrix=scores,
            ),
            config=build_signalome_config(
                substrate_support_cutoff=0.5,
                score_preconditioning_policy="allow_and_report",
            ),
        )
    )


def test_signalome_score_preconditioning_record_uses_actual_stage_indexes() -> None:
    interpreted = _interpreted_request(
        display_ids=["P1;S1;", "P2;S2;", "P3;S3;"],
        score_values=[
            [float("nan"), float("nan")],
            [2.0, 3.0],
            [1.0, 2.0],
        ],
    )

    provenance = build_signalome_row_attrition_provenance(
        interpreted,
        final_site_ids=interpreted.downstream_score_matrix.index,
    )

    assert provenance.row_attrition is not None
    records = provenance.row_attrition.records
    assert len(records) == 1
    assert records[0].stage == "signalome_score_preconditioning"
    assert records[0].reason == "sites_removed_by_score_preconditioning"
    assert records[0].input_rows == 3
    assert records[0].output_rows == 2
    assert records[0].examples == (str(interpreted.dataset.phospho.index[0]),)
    assert provenance.metrics["sites_removed_by_score_preconditioning"] == 1
    assert provenance.row_attrition.final_rows == int(
        interpreted.downstream_score_matrix.shape[0]
    )


def test_signalome_internal_sequential_records_are_continuous_without_double_count() -> (
    None
):
    display_ids = ["P1;S1;", "P2;S2;", "P3;S3;", "P4;S4;"]
    dataset = _dataset(display_ids)
    site_index = dataset.phospho.index
    # The public validator rejects prediction/score index mismatches. This direct
    # interpreter test covers the private alignment filter without weakening that
    # public validation boundary.
    interpreted = _interpreted_request(
        display_ids=display_ids,
        prediction_site_index=site_index.delete(3),
        score_site_index=site_index,
        score_values=[
            [1.0, 0.5],
            [float("nan"), float("nan")],
            [2.0, 3.0],
            [4.0, 5.0],
        ],
    )

    provenance = build_signalome_row_attrition_provenance(
        interpreted,
        final_site_ids=interpreted.downstream_score_matrix.index,
    )

    assert provenance.row_attrition is not None
    records = provenance.row_attrition.records
    assert [record.stage for record in records] == [
        "signalome_site_alignment",
        "signalome_score_preconditioning",
    ]
    assert [(record.input_rows, record.output_rows) for record in records] == [
        (4, 3),
        (3, 2),
    ]
    assert sum(record.removed_rows for record in records) == 2
    assert all(record.removed_rows > 0 for record in records)


def test_signalome_no_zero_removal_record_when_no_rows_removed() -> None:
    interpreted = _interpreted_request(
        display_ids=["P1;S1;", "P2;S2;", "P3;S3;"],
    )

    provenance = build_signalome_row_attrition_provenance(
        interpreted,
        final_site_ids=interpreted.downstream_score_matrix.index,
    )

    assert provenance.row_attrition is None
    assert provenance.metrics["sites_removed_by_score_preconditioning"] == 0
