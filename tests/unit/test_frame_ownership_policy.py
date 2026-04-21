from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
)
from phospy.activities.models import KinaseActivityInputs, PredMatOverlapSummary
from phospy.api import (
    DatasetBuildRequest,
    Organism,
)
from phospy.api.results import KinasePredictionResult
from phospy.datasets.builders.executor import DatasetBuildExecutor
from phospy.datasets.builders.interpreter import DatasetBuildRequestInterpreter
from tests.support.transformation_states import supported_linear_state


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [2.0, 1.0],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )


def _site_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
            ],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )


@dataclass(slots=True)
class _CopyCounts:
    dataframe_deep: int = 0


@contextmanager
def _count_dataframe_deep_copies() -> Iterator[_CopyCounts]:
    counts = _CopyCounts()
    original_copy = pd.DataFrame.copy

    def wrapped_copy(self, *args, **kwargs):  # type: ignore[no-untyped-def]
        deep = kwargs.get("deep", args[0] if args else True)
        if bool(deep):
            counts.dataframe_deep += 1
        return original_copy(self, *args, **kwargs)

    pd.DataFrame.copy = wrapped_copy
    try:
        yield counts
    finally:
        pd.DataFrame.copy = original_copy


def test_public_dataset_isolated_from_caller_mutation() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()
    dataset = AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        transformation_state=supported_linear_state(has_total_matrix=False),
    )

    phospho.iloc[0, 0] = 999.0
    site_metadata.iloc[0, 0] = "CHANGED"

    assert float(dataset.phospho.iloc[0, 0]) == 1.0
    assert str(dataset.site_metadata.iloc[0, 0]) == "MAPK14"


def test_builder_result_isolated_from_caller_mutation_after_build() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )

    phospho.iloc[1, 1] = 777.0
    site_metadata.iloc[1, 0] = "CHANGED"

    assert float(built.phospho.iloc[1, 1]) == 1.0
    assert str(built.site_metadata.iloc[1, 0]) == "GSK3B"


def test_builder_stage_handoff_transfers_owned_frames_without_recopies() -> None:
    request = DatasetBuildRequest(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
    )
    interpreted = DatasetBuildRequestInterpreter().run(request)
    built = DatasetBuildExecutor().run(interpreted)

    assert built.phospho is interpreted.phospho
    assert built.site_metadata is interpreted.site_metadata


def test_builder_dataframe_copy_churn_regression_budget() -> None:
    request = DatasetBuildRequest(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
    )

    with _count_dataframe_deep_copies() as counts:
        AnalysisReadyDatasetBuilder().run(request)

    assert counts.dataframe_deep == 2


def test_internal_activity_inputs_alias_owned_frames() -> None:
    pred_mat = _phospho()
    phospho_matrix = _phospho()
    overlap_summary = PredMatOverlapSummary(
        overlap_count=2,
        pred_mat_rows=2,
        phospho_rows=2,
    )

    inputs = KinaseActivityInputs(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=0.5,
        min_substrates=1,
        top_n_substrates=2,
        overlap_summary=overlap_summary,
    )

    assert inputs.pred_mat is pred_mat
    assert inputs.phospho_matrix is phospho_matrix


def test_prediction_result_boundary_copy_and_owned_transfer_modes() -> None:
    pred_mat = _phospho()
    substrate_list = pd.DataFrame(
        {
            "kinase": ["MAP2K6"],
            "substrate_site": ["MAPK14;Y182;"],
            "score": [0.75],
            "rank": [1],
        }
    )

    copied_result = KinasePredictionResult(
        pred_mat=pred_mat,
        substrate_list=substrate_list,
    )
    owned_result = KinasePredictionResult._from_owned(
        pred_mat=pred_mat,
        substrate_list=substrate_list,
    )

    assert copied_result.pred_mat is not pred_mat
    assert copied_result.substrate_list is not substrate_list
    assert owned_result.pred_mat is pred_mat
    assert owned_result.substrate_list is substrate_list

    pred_mat.iloc[0, 0] = 999.0
    substrate_list.iloc[0, 0] = "CHANGED"
    assert float(copied_result.pred_mat.iloc[0, 0]) == 1.0
    assert str(copied_result.substrate_list.iloc[0, 0]) == "MAP2K6"
