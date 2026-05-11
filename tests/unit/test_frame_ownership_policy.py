from __future__ import annotations

import inspect
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    KinaseWorkflow,
    SignalomeWorkflow,
)
from phospy.activities.models import (
    KinaseActivityInputs,
    KinaseActivityResult,
    PredMatOverlapSummary,
)
from phospy.api import (
    DatasetBuildRequest,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    SignalomeWorkflowRequest,
)
from phospy.api.results import (
    KinasePredictionResult,
    KinaseWorkflowResult,
    SignalomeWorkflowResult,
)
from phospy.datasets.builders.executor import DatasetBuildExecutor
from phospy.datasets.builders.interpreter import DatasetBuildRequestInterpreter
from phospy.datasets.models import DatasetPreprocessingReport
from phospy.datasets.preprocessing.report_schema import (
    ComparisonGroupStatsRow,
    ComparisonPairStatsRow,
    DuplicateSiteResolutionRow,
    MetadataConflictRow,
    PreprocessingOperationRow,
    PreprocessingRowAuditRow,
    PreprocessingRowCountRow,
)
from phospy.prediction.models import KinaseScoringResult
from phospy.provenance.hashing import fingerprint_table
from phospy.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeModules,
)
from phospy.tables.base import TableSchema
from phospy.tables.datasets import PhosphoIntensityMatrix
from phospy.workflows.signalome.interpreter import SignalomeWorkflowInterpreter
from phospy.workflows.signalome.validator import SignalomeWorkflowValidator
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.signalome_config import build_signalome_config


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
            "protein_id": ["MAPK14", "GSK3B"],
            "localisation_confidence": [0.95, 0.9],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )


def _references() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6"],
                "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                    "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
                ]
            },
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
        ),
    )


def _kinase_result():
    site_ids = ["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"]
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 4.0],
            "sample_b": [2.0, 4.0, 1.0],
        },
        index=site_ids,
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1"],
            "site": ["Y182", "S9", "T308"],
            "site_sequence": ["A" * 31, "B" * 31, "C" * 31],
            "protein_id": ["MAPK14", "GSK3B", "AKT1"],
            "localisation_confidence": [0.95, 0.9, 0.92],
        },
        index=site_ids,
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6", "AKT1", "AKT1"],
                "substrate_site": [
                    "MAPK14;Y182;",
                    "GSK3B;S9;",
                    "GSK3B;S9;",
                    "AKT1;T308;",
                ],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31, "B" * 31, "C" * 31]},
            index=pd.Index(site_ids, name="site_id"),
        ),
    )
    return KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=AnalysisReadyPhosphoDataset(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                intensity_scale_state=supported_linear_intensity_scale_state(
                    has_total_matrix=False
                ),
                processing_state=supported_linear_processing_state(
                    has_total_matrix=False
                ),
            ),
            references=references,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=1,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
            ),
            activity_config=None,
        )
    )


def _signalome_request_with_borrowed_frames() -> SignalomeWorkflowRequest:
    dataset = AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )
    prediction_matrix = pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.8],
            "AKT1": [0.1, 0.2],
        },
        index=dataset._borrow_phospho_frame().index.copy(),
    )
    score_matrix = pd.DataFrame(
        {
            "MAP2K6": [1.5, 1.2],
            "AKT1": [0.6, 0.7],
        },
        index=dataset._borrow_phospho_frame().index.copy(),
    )
    return SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=dataset,
            references=_references(),
            scoring_result=KinaseScoringResult(
                profile_scores=score_matrix,
                rank_weighted_fusion_scores=score_matrix,
            ),
            prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
            activity_result=None,
        ),
        config=build_signalome_config(substrate_support_cutoff=0.5),
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


def _mutate_first_frame_cell(frame: pd.DataFrame) -> None:
    if pd.api.types.is_numeric_dtype(frame.dtypes.iloc[0]):
        frame.iloc[0, 0] = float(frame.iloc[0, 0]) + 1.0
    else:
        frame.iloc[0, 0] = f"{frame.iloc[0, 0]}_changed"


def _assert_dataframe_getter_defensive_snapshot(
    getter: Callable[[], pd.DataFrame],
) -> None:
    exported = getter()
    _mutate_first_frame_cell(exported)
    reread = getter()
    assert exported is not reread
    assert exported.iloc[0, 0] != reread.iloc[0, 0]


def _assert_optional_dataframe_getter_defensive_snapshot(
    getter: Callable[[], pd.DataFrame | None],
) -> None:
    exported = getter()
    assert exported is not None
    _mutate_first_frame_cell(exported)
    reread = getter()
    assert reread is not None
    assert exported is not reread
    assert exported.iloc[0, 0] != reread.iloc[0, 0]


def _assert_optional_dataframe_getter_defensive_snapshot_when_present(
    getter: Callable[[], pd.DataFrame | None],
) -> None:
    exported = getter()
    if exported is None:
        return

    _mutate_first_frame_cell(exported)
    reread = getter()
    assert reread is not None
    assert exported is not reread
    assert exported.iloc[0, 0] != reread.iloc[0, 0]


def _assert_copy_keyword_rejected(
    export: Callable[..., object],
) -> None:
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        export(copy=False)  # type: ignore[call-arg]


@pytest.mark.parametrize(
    (
        "builder",
        "mutated_phospho",
        "mutated_site_metadata",
        "expected_phospho",
        "expected_gene",
    ),
    [
        pytest.param(
            lambda p, s: AnalysisReadyPhosphoDataset(
                phospho=p,
                site_metadata=s,
                organism=Organism.RAT,
                intensity_scale_state=supported_linear_intensity_scale_state(
                    has_total_matrix=False
                ),
                processing_state=supported_linear_processing_state(
                    has_total_matrix=False
                ),
            ),
            (0, 0, 999.0),
            (0, 0, "CHANGED"),
            1.0,
            "MAPK14",
            id="public-constructor-copies-caller-inputs",
        ),
        pytest.param(
            lambda p, s: AnalysisReadyDatasetBuilder().run(
                DatasetBuildRequest(
                    phospho=p,
                    site_metadata=s,
                    organism=Organism.RAT,
                )
            ),
            (1, 1, 777.0),
            (1, 0, "CHANGED"),
            1.0,
            "GSK3B",
            id="builder-result-copies-caller-inputs",
        ),
    ],
)
def test_public_constructor_copy_contract_matrix(
    builder: Callable[[pd.DataFrame, pd.DataFrame], AnalysisReadyPhosphoDataset],
    mutated_phospho: tuple[int, int, float],
    mutated_site_metadata: tuple[int, int, str],
    expected_phospho: float,
    expected_gene: str,
) -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()
    built = builder(phospho, site_metadata)

    phospho.iloc[mutated_phospho[0], mutated_phospho[1]] = mutated_phospho[2]
    site_metadata.iloc[mutated_site_metadata[0], mutated_site_metadata[1]] = (
        mutated_site_metadata[2]
    )

    assert float(built.phospho.iloc[mutated_phospho[0], mutated_phospho[1]]) == (
        expected_phospho
    )
    assert (
        str(
            built.site_metadata.iloc[mutated_site_metadata[0], mutated_site_metadata[1]]
        )
        == expected_gene
    )


def test_builder_stage_handoff_transfers_owned_frames_without_recopies() -> None:
    request = DatasetBuildRequest(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
    )
    interpreted = DatasetBuildRequestInterpreter().run(request)
    built = DatasetBuildExecutor().run(interpreted)

    assert built._phospho is interpreted.phospho
    assert built._site_metadata is interpreted.site_metadata


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
    pred_mat = pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.8],
            "AKT1": [0.2, 0.1],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )
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
    pred_mat = pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.8],
            "AKT1": [0.2, 0.1],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )
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

    assert copied_result._pred_mat is not pred_mat
    assert copied_result._substrate_list is not substrate_list
    assert owned_result._pred_mat is pred_mat
    assert owned_result._substrate_list is substrate_list

    pred_mat.iloc[0, 0] = 999.0
    substrate_list.iloc[0, 0] = "CHANGED"
    assert float(copied_result.pred_mat.iloc[0, 0]) == 0.9
    assert str(copied_result.substrate_list.iloc[0, 0]) == "MAP2K6"


def test_public_dataframe_accessors_do_not_accept_copy_keyword() -> None:
    public_accessors = (
        (AnalysisReadyPhosphoDataset, "to_dataframe"),
        (AnalysisReadyPhosphoDataset, "site_metadata_dataframe"),
        (AnalysisReadyPhosphoDataset, "sample_metadata_dataframe"),
        (AnalysisReadyPhosphoDataset, "total_dataframe"),
        (AnalysisReadyPhosphoDataset, "comparisons_dataframe"),
        (DatasetPreprocessingReport, "row_counts_dataframe"),
        (DatasetPreprocessingReport, "operations_dataframe"),
        (DatasetPreprocessingReport, "row_audit_dataframe"),
        (DatasetPreprocessingReport, "duplicate_site_resolution_dataframe"),
        (DatasetPreprocessingReport, "metadata_conflicts_dataframe"),
        (DatasetPreprocessingReport, "comparison_group_stats_dataframe"),
        (DatasetPreprocessingReport, "comparison_pair_stats_dataframe"),
        (KinaseScoringResult, "to_dataframe"),
        (KinaseScoringResult, "motif_scores_dataframe"),
        (KinaseScoringResult, "rank_weighted_fusion_scores_dataframe"),
        (KinaseScoringResult, "score_fusion_weights_dataframe"),
        (KinasePredictionResult, "to_dataframe"),
        (KinasePredictionResult, "substrate_list_dataframe"),
        (KinaseActivityResult, "to_dataframe"),
        (KinaseActivityResult, "thresholded_substrate_mean_activity_dataframe"),
        (KinaseActivityResult, "target_table_dataframe"),
        (ReferenceBundle, "kinase_substrate_map_dataframe"),
        (ReferenceBundle, "site_sequences_dataframe"),
        (SignalomeWorkflowResult, "to_dataframe"),
        (SignalomeWorkflowResult, "site_membership_dataframe"),
        (SignalomeWorkflowResult, "protein_site_context_dataframe"),
        (SignalomeAssignments, "to_pandas"),
        (SignalomeModules, "to_pandas"),
        (KinaseNetwork, "to_pandas"),
        (KinaseNetwork, "nodes_dataframe"),
        (KinaseNetwork, "candidate_correlations_dataframe"),
        (TableSchema, "to_pandas"),
    )

    for owner, method_name in public_accessors:
        signature = inspect.signature(getattr(owner, method_name))
        assert "copy" not in signature.parameters


def test_internal_borrowed_dataset_access_is_mutation_isolated_without_deep_copy() -> (
    None
):
    dataset = AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )

    with _count_dataframe_deep_copies() as counts:
        borrowed = dataset._borrow_phospho_frame()
        borrowed.iloc[0, 0] = 999.0
        borrowed.loc[:, "borrowed_only"] = [1.0, 2.0]

    assert borrowed is not dataset._phospho
    assert counts.dataframe_deep == 0
    assert not hasattr(dataset, "borrow_phospho_frame")
    assert float(dataset._phospho.iloc[0, 0]) == 1.0
    assert "borrowed_only" not in dataset._phospho.columns


def test_internal_borrowed_prediction_and_scoring_access_is_mutation_isolated() -> None:
    pred_mat = pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.8],
            "AKT1": [0.2, 0.1],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )
    profile_scores = pd.DataFrame(
        {"MAP2K6": [0.8, 0.2]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"]),
    )
    rank_weighted = pd.DataFrame(
        {"MAP2K6": [0.75, 0.15]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"]),
    )
    prediction_result = KinasePredictionResult._from_owned(pred_mat=pred_mat)
    scoring_result = KinaseScoringResult._from_owned(
        profile_scores=profile_scores,
        rank_weighted_fusion_scores=rank_weighted,
    )

    with _count_dataframe_deep_copies() as counts:
        borrowed_pred = prediction_result._borrow_pred_mat_frame()
        borrowed_profile = scoring_result._borrow_profile_scores_frame()
        borrowed_rank_weighted = (
            scoring_result._borrow_rank_weighted_fusion_scores_frame()
        )
        borrowed_pred.iloc[0, 0] = 99.0
        borrowed_profile.iloc[0, 0] = 88.0
        assert borrowed_rank_weighted is not None
        borrowed_rank_weighted.iloc[0, 0] = 77.0

    assert counts.dataframe_deep == 0
    assert borrowed_pred is not prediction_result._pred_mat
    assert borrowed_profile is not scoring_result._profile_scores
    assert borrowed_rank_weighted is not scoring_result._rank_weighted_fusion_scores
    assert float(prediction_result._pred_mat.iloc[0, 0]) == 0.9
    assert float(scoring_result._profile_scores.iloc[0, 0]) == 0.8
    assert scoring_result._rank_weighted_fusion_scores is not None
    assert float(scoring_result._rank_weighted_fusion_scores.iloc[0, 0]) == 0.75
    assert not hasattr(prediction_result, "borrow_pred_mat_frame")
    assert not hasattr(scoring_result, "borrow_profile_scores_frame")


def test_signalome_validator_borrowed_reads_do_not_mutate_internal_frames() -> None:
    request = _signalome_request_with_borrowed_frames()
    dataset = request.kinase_result.dataset
    prediction_result = request.kinase_result.prediction_result
    scoring_result = request.kinase_result.scoring_result

    dataset_phospho_before = fingerprint_table(
        dataset._borrow_phospho_frame(),
        name="dataset.phospho",
    )
    dataset_site_metadata_before = fingerprint_table(
        dataset._borrow_site_metadata_frame(),
        name="dataset.site_metadata",
    )
    prediction_before = fingerprint_table(
        prediction_result._borrow_pred_mat_frame(),
        name="prediction_result.pred_mat",
    )
    score_before = fingerprint_table(
        scoring_result._borrow_profile_scores_frame(),
        name="scoring_result.profile_scores",
    )
    assert scoring_result._borrow_rank_weighted_fusion_scores_frame() is not None
    rank_weighted_before = fingerprint_table(
        scoring_result._borrow_rank_weighted_fusion_scores_frame(),
        name="scoring_result.rank_weighted_fusion_scores",
    )

    validated = SignalomeWorkflowValidator().run(request)
    assert validated is request

    dataset_phospho_after = fingerprint_table(
        dataset._borrow_phospho_frame(),
        name="dataset.phospho",
    )
    dataset_site_metadata_after = fingerprint_table(
        dataset._borrow_site_metadata_frame(),
        name="dataset.site_metadata",
    )
    prediction_after = fingerprint_table(
        prediction_result._borrow_pred_mat_frame(),
        name="prediction_result.pred_mat",
    )
    score_after = fingerprint_table(
        scoring_result._borrow_profile_scores_frame(),
        name="scoring_result.profile_scores",
    )
    assert scoring_result._borrow_rank_weighted_fusion_scores_frame() is not None
    rank_weighted_after = fingerprint_table(
        scoring_result._borrow_rank_weighted_fusion_scores_frame(),
        name="scoring_result.rank_weighted_fusion_scores",
    )

    assert dataset_phospho_before.hash_value == dataset_phospho_after.hash_value
    assert (
        dataset_site_metadata_before.hash_value
        == dataset_site_metadata_after.hash_value
    )
    assert prediction_before.hash_value == prediction_after.hash_value
    assert score_before.hash_value == score_after.hash_value
    assert rank_weighted_before.hash_value == rank_weighted_after.hash_value


def test_signalome_interpreter_read_path_does_not_mutate_borrowed_dataset_frames() -> (
    None
):
    request = _signalome_request_with_borrowed_frames()
    dataset = request.kinase_result.dataset

    phospho_before = fingerprint_table(
        dataset._borrow_phospho_frame(),
        name="dataset.phospho",
    )
    site_metadata_before = fingerprint_table(
        dataset._borrow_site_metadata_frame(),
        name="dataset.site_metadata",
    )

    SignalomeWorkflowInterpreter().run(request)

    phospho_after = fingerprint_table(
        dataset._borrow_phospho_frame(),
        name="dataset.phospho",
    )
    site_metadata_after = fingerprint_table(
        dataset._borrow_site_metadata_frame(),
        name="dataset.site_metadata",
    )

    assert phospho_before.hash_value == phospho_after.hash_value
    assert site_metadata_before.hash_value == site_metadata_after.hash_value


def test_owned_construction_frames_can_be_mutated_after_owned_transfer() -> None:
    phospho = _phospho()
    site_metadata = _site_metadata()
    dataset = AnalysisReadyPhosphoDataset._from_owned(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )

    phospho.iloc[0, 0] = 321.0
    site_metadata.iloc[0, 0] = "UPDATED_GENE"
    assert float(dataset._borrow_phospho_frame().iloc[0, 0]) == 321.0
    assert str(dataset._borrow_site_metadata_frame().iloc[0, 0]) == "UPDATED_GENE"

    public_snapshot = dataset.phospho
    public_snapshot.iloc[0, 0] = 123.0
    assert float(dataset._borrow_phospho_frame().iloc[0, 0]) == 321.0


def test_internal_borrowed_accessors_are_not_public_api_exports() -> None:
    import phospy
    import phospy._frame_ownership as frame_ownership
    import phospy.datasets as datasets

    assert not any("borrow" in name for name in phospy.__all__)
    assert not any("borrow" in name for name in datasets.__all__)
    assert not any("borrow" in name for name in frame_ownership.__all__)
    assert not hasattr(AnalysisReadyPhosphoDataset, "borrow_phospho_frame")
    assert not hasattr(AnalysisReadyPhosphoDataset, "borrow_site_metadata_frame")
    assert not hasattr(DatasetPreprocessingReport, "borrow_row_counts_frame")


def test_prediction_result_public_export_copy_default_is_safe() -> None:
    pred_mat = pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.8],
            "AKT1": [0.2, 0.1],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )
    result = KinasePredictionResult._from_owned(pred_mat=pred_mat)

    exported = result.to_dataframe()
    exported.iloc[0, 0] = 0.0

    assert float(result.pred_mat.iloc[0, 0]) == 0.9


def test_safe_public_export_does_not_change_owned_provenance_state() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            organism=Organism.RAT,
        )
    )
    fingerprint_before = fingerprint_table(built.phospho, name="dataset.phospho")

    safe_copy = built.to_dataframe()
    safe_copy.iloc[0, 0] = 999.0

    fingerprint_after = fingerprint_table(built.phospho, name="dataset.phospho")
    assert fingerprint_before.hash_value == fingerprint_after.hash_value


@pytest.mark.parametrize(
    ("export_factory",),
    [
        pytest.param(
            lambda: (
                AnalysisReadyPhosphoDataset(
                    phospho=_phospho(),
                    site_metadata=_site_metadata(),
                    organism=Organism.RAT,
                    intensity_scale_state=supported_linear_intensity_scale_state(
                        has_total_matrix=False
                    ),
                    processing_state=supported_linear_processing_state(
                        has_total_matrix=False
                    ),
                ).to_dataframe
            ),
            id="dataset-export-rejects-copy-keyword",
        ),
        pytest.param(
            lambda: (
                AnalysisReadyDatasetBuilder()
                .run(
                    DatasetBuildRequest(
                        phospho=_phospho(),
                        site_metadata=_site_metadata(),
                        organism=Organism.RAT,
                    )
                )
                .to_dataframe
            ),
            id="builder-output-export-rejects-copy-keyword",
        ),
        pytest.param(
            lambda: (
                KinasePredictionResult._from_owned(
                    pred_mat=pd.DataFrame(
                        {
                            "MAP2K6": [0.9, 0.8],
                            "AKT1": [0.2, 0.1],
                        },
                        index=["MAPK14;Y182;", "GSK3B;S9;"],
                    )
                ).to_dataframe
            ),
            id="prediction-export-rejects-copy-keyword",
        ),
        pytest.param(
            lambda: (
                SignalomeWorkflow()
                .run(
                    SignalomeWorkflowRequest(
                        kinase_result=_kinase_result(),
                        config=build_signalome_config(substrate_support_cutoff=0.5),
                    )
                )
                .to_dataframe
            ),
            id="signalome-result-export-rejects-copy-keyword",
        ),
        pytest.param(
            lambda: (
                SignalomeWorkflow()
                .run(
                    SignalomeWorkflowRequest(
                        kinase_result=_kinase_result(),
                        config=build_signalome_config(substrate_support_cutoff=0.5),
                    )
                )
                .module_assignments.to_pandas
            ),
            id="signalome-module-assignments-export-rejects-copy-keyword",
        ),
        pytest.param(
            lambda: (
                SignalomeWorkflow()
                .run(
                    SignalomeWorkflowRequest(
                        kinase_result=_kinase_result(),
                        config=build_signalome_config(substrate_support_cutoff=0.5),
                    )
                )
                .signalome_modules.to_pandas
            ),
            id="signalome-modules-export-rejects-copy-keyword",
        ),
        pytest.param(
            lambda: (
                SignalomeWorkflow()
                .run(
                    SignalomeWorkflowRequest(
                        kinase_result=_kinase_result(),
                        config=build_signalome_config(substrate_support_cutoff=0.5),
                    )
                )
                .kinase_network.to_pandas
            ),
            id="signalome-network-export-rejects-copy-keyword",
        ),
        pytest.param(
            lambda: PhosphoIntensityMatrix(frame=_phospho()).to_pandas,
            id="table-schema-export-rejects-copy-keyword",
        ),
    ],
)
def test_public_export_copy_keyword_rejection_contract_matrix(
    export_factory: Callable[[], Callable[..., object]],
) -> None:
    _assert_copy_keyword_rejected(export_factory())


def test_public_signalome_exports_isolated_from_mutation() -> None:
    result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=_kinase_result(),
            config=build_signalome_config(substrate_support_cutoff=0.5),
        )
    )

    # Public exports must be defensive snapshots across representative result types.
    for getter in (
        result.module_assignments.to_pandas,
        result.signalome_modules.to_pandas,
        result.kinase_network.to_pandas,
    ):
        _assert_dataframe_getter_defensive_snapshot(getter)


@pytest.mark.parametrize(
    ("getter_factory",),
    [
        pytest.param(
            lambda: (
                AnalysisReadyPhosphoDataset(
                    phospho=_phospho(),
                    site_metadata=_site_metadata(),
                    organism=Organism.RAT,
                    intensity_scale_state=supported_linear_intensity_scale_state(
                        has_total_matrix=False
                    ),
                    processing_state=supported_linear_processing_state(
                        has_total_matrix=False
                    ),
                ).to_dataframe
            ),
            id="dataset-export-snapshot",
        ),
        pytest.param(
            lambda: (
                KinasePredictionResult._from_owned(
                    pred_mat=pd.DataFrame(
                        {
                            "MAP2K6": [0.9, 0.8],
                            "AKT1": [0.2, 0.1],
                        },
                        index=["MAPK14;Y182;", "GSK3B;S9;"],
                    )
                ).to_dataframe
            ),
            id="prediction-export-snapshot",
        ),
        pytest.param(
            lambda: _references().kinase_substrate_map_dataframe,
            id="reference-export-snapshot",
        ),
    ],
)
def test_public_export_snapshot_contract_matrix(
    getter_factory: Callable[[], Callable[[], pd.DataFrame]],
) -> None:
    _assert_dataframe_getter_defensive_snapshot(getter_factory())


def test_dataset_dataframe_properties_are_defensive_snapshots() -> None:
    sample_metadata = pd.DataFrame(
        {"batch": [1, 2]},
        index=pd.Index(["sample_a", "sample_b"]),
    )
    dataset = AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        sample_metadata=sample_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )

    _assert_dataframe_getter_defensive_snapshot(lambda: dataset.phospho)
    _assert_dataframe_getter_defensive_snapshot(lambda: dataset.site_metadata)
    _assert_optional_dataframe_getter_defensive_snapshot(
        lambda: dataset.sample_metadata
    )


def test_preprocessing_report_dataframe_properties_are_defensive_snapshots() -> None:
    report = DatasetPreprocessingReport.from_rows(
        row_count_rows=(
            PreprocessingRowCountRow(
                stage="missing_data",
                input_rows=2,
                output_rows=2,
                dropped_rows=0,
            ),
        ),
        operation_rows=(
            PreprocessingOperationRow(
                step_order=1,
                stage="missing_data",
                operation="forbid",
                parameters={},
                input_rows=2,
                output_rows=2,
                notes=None,
            ),
        ),
        row_audit_rows=(
            PreprocessingRowAuditRow(
                stage="missing_data",
                action="retained",
                reason="complete",
                source_row_id="MAPK14;Y182;",
                site_id="MAPK14;Y182;",
                retained=True,
                retained_row_id="MAPK14;Y182;",
                source_rows=None,
                retained_row=None,
                parameter_snapshot={},
            ),
        ),
        duplicate_site_resolution_rows=(
            DuplicateSiteResolutionRow(
                site_id="MAPK14;Y182;",
                source_row_id="MAPK14;Y182;",
                retained=True,
                resolution_policy="max_mean_signal",
                aggregation_method="max_mean_signal",
                missing_value_policy=None,
                metadata_resolution_policy="retain_row_ranked_by_observed_values_then_mean_signal_then_input_order",
                retained_reason=None,
                dropped_reason=None,
                observed_values=None,
                mean_signal=1.0,
                n_source_rows=1,
                n_aggregated_rows=1,
                source_protein_id="MAPK14",
                source_gene_symbol="MAPK14",
                source_site="Y182",
                source_site_sequence="A",
                metadata_conflict_detected=False,
            ),
        ),
        metadata_conflict_rows=(
            MetadataConflictRow(
                site_id="MAPK14;Y182;",
                field="protein_id",
                values=["MAPK14"],
                n_distinct_values=1,
                source_row_ids=["MAPK14;Y182;"],
            ),
        ),
        comparison_group_stats_rows=(
            ComparisonGroupStatsRow(
                site_id="MAPK14;Y182;",
                group="group_a",
                n=1,
                mean=1.0,
                sd=None,
                sem=None,
                median=1.0,
                min=1.0,
                max=1.0,
                sample_ids=["sample_a"],
            ),
        ),
        comparison_pair_stats_rows=(
            ComparisonPairStatsRow(
                site_id="MAPK14;Y182;",
                comparison="p_group_a_group_b",
                left_group="group_a",
                right_group="group_b",
                left_n=1,
                right_n=1,
                left_mean=1.0,
                right_mean=2.0,
                left_sd=None,
                right_sd=None,
                left_sem=None,
                right_sem=None,
                effect_size=-1.0,
                left_median=1.0,
                right_median=2.0,
                left_min=1.0,
                right_min=2.0,
                left_max=1.0,
                right_max=2.0,
            ),
        ),
    )

    for getter in (
        lambda: report.row_counts,
        lambda: report.operations,
        lambda: report.row_audit,
    ):
        _assert_dataframe_getter_defensive_snapshot(getter)

    for optional_getter in (
        lambda: report.duplicate_site_resolution,
        lambda: report.metadata_conflicts,
        lambda: report.comparison_group_stats,
        lambda: report.comparison_pair_stats,
    ):
        _assert_optional_dataframe_getter_defensive_snapshot(optional_getter)


def test_kinase_result_table_properties_are_defensive_snapshots() -> None:
    scoring_result = KinaseScoringResult._from_owned(
        profile_scores=pd.DataFrame(
            {"MAP2K6": [0.8, 0.2]},
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"]),
        ),
        motif_scores=pd.DataFrame(
            {"MAP2K6": [0.7, 0.1]},
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"]),
        ),
        rank_weighted_fusion_scores=pd.DataFrame(
            {"MAP2K6": [0.75, 0.15]},
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"]),
        ),
        score_fusion_weights=pd.DataFrame(
            {"MAP2K6": [1.0, 1.0]},
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"]),
        ),
    )
    _assert_dataframe_getter_defensive_snapshot(lambda: scoring_result.profile_scores)
    for optional_getter in (
        lambda: scoring_result.motif_scores,
        lambda: scoring_result.rank_weighted_fusion_scores,
        lambda: scoring_result.score_fusion_weights,
    ):
        _assert_optional_dataframe_getter_defensive_snapshot(optional_getter)

    prediction_result = KinasePredictionResult._from_owned(
        pred_mat=pd.DataFrame(
            {"MAP2K6": [0.9, 0.8]},
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"]),
        ),
        substrate_list=pd.DataFrame(
            {
                "kinase": ["MAP2K6"],
                "substrate_site": ["MAPK14;Y182;"],
                "score": [0.9],
                "rank": [1],
            }
        ),
    )
    _assert_dataframe_getter_defensive_snapshot(lambda: prediction_result.pred_mat)
    _assert_optional_dataframe_getter_defensive_snapshot(
        lambda: prediction_result.substrate_list
    )

    activity_result = KinaseActivityResult._from_owned(
        weighted_activity=pd.DataFrame(
            {"MAP2K6": [1.0, 2.0]},
            index=pd.Index(["sample_a", "sample_b"]),
        ),
        thresholded_substrate_mean_activity=pd.DataFrame(
            {"MAP2K6": [0.5, 1.5]},
            index=pd.Index(["sample_a", "sample_b"]),
        ),
        thresholded_substrate_counts=pd.Series(
            [2, 2],
            index=pd.Index(["MAP2K6", "AKT1"]),
            name="n_substrates",
        ),
        target_counts=pd.Series(
            [1, 1],
            index=pd.Index(["MAP2K6", "AKT1"]),
            name="n_targets",
        ),
        target_table=pd.DataFrame(
            {
                "site_id": ["MAPK14;Y182;"],
                "kinase": ["MAP2K6"],
                "score": [0.9],
            }
        ),
    )
    for getter in (
        lambda: activity_result.weighted_activity,
        lambda: activity_result.activity_scores,
        lambda: activity_result.thresholded_substrate_mean_activity,
        lambda: activity_result.target_table,
    ):
        _assert_dataframe_getter_defensive_snapshot(getter)


def test_kinase_activity_result_series_properties_are_defensive_snapshots() -> None:
    activity_result = KinaseActivityResult._from_owned(
        weighted_activity=pd.DataFrame(
            {"MAP2K6": [1.0, 2.0]},
            index=pd.Index(["sample_a", "sample_b"]),
        ),
        thresholded_substrate_mean_activity=pd.DataFrame(
            {"MAP2K6": [0.5, 1.5]},
            index=pd.Index(["sample_a", "sample_b"]),
        ),
        thresholded_substrate_counts=pd.Series(
            [2, 2],
            index=pd.Index(["MAP2K6", "AKT1"]),
            name="n_substrates",
        ),
        target_counts=pd.Series(
            [1, 1],
            index=pd.Index(["MAP2K6", "AKT1"]),
            name="n_targets",
        ),
        target_table=pd.DataFrame(
            {
                "site_id": ["MAPK14;Y182;"],
                "kinase": ["MAP2K6"],
                "score": [0.9],
            }
        ),
    )

    assert hasattr(activity_result, "thresholded_substrate_counts")
    assert hasattr(activity_result, "target_counts")

    thresholded_before = fingerprint_table(
        activity_result.thresholded_substrate_counts.to_frame(name="n_substrates"),
        name="outputs.activity.thresholded_substrate_counts",
    )
    exported_thresholded = activity_result.thresholded_substrate_counts
    exported_thresholded.iloc[0] = 999
    reread_thresholded = activity_result.thresholded_substrate_counts
    assert exported_thresholded is not reread_thresholded
    assert reread_thresholded.to_dict() == {"MAP2K6": 2, "AKT1": 2}
    thresholded_after = fingerprint_table(
        activity_result.thresholded_substrate_counts.to_frame(name="n_substrates"),
        name="outputs.activity.thresholded_substrate_counts",
    )
    assert thresholded_before.hash_value == thresholded_after.hash_value

    target_before = fingerprint_table(
        activity_result.target_counts.to_frame(name="n_targets"),
        name="outputs.activity.target_counts",
    )
    exported_target = activity_result.target_counts
    exported_target.iloc[0] = 999
    reread_target = activity_result.target_counts
    assert exported_target is not reread_target
    assert reread_target.to_dict() == {"MAP2K6": 1, "AKT1": 1}
    target_after = fingerprint_table(
        activity_result.target_counts.to_frame(name="n_targets"),
        name="outputs.activity.target_counts",
    )
    assert target_before.hash_value == target_after.hash_value


def test_signalome_result_table_properties_are_defensive_snapshots() -> None:
    signalome_result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=_kinase_result(),
            config=build_signalome_config(substrate_support_cutoff=0.5),
        )
    )

    for getter in (
        lambda: signalome_result.module_assignments.table,
        lambda: signalome_result.signalome_modules.table,
        lambda: signalome_result.kinase_network.edges,
    ):
        _assert_dataframe_getter_defensive_snapshot(getter)

    for optional_getter in (
        lambda: signalome_result.kinase_network.nodes,
        lambda: signalome_result.kinase_network.candidate_correlations,
        lambda: signalome_result.expanded_signalome,
        lambda: signalome_result.site_membership,
        lambda: signalome_result.protein_site_context,
    ):
        _assert_optional_dataframe_getter_defensive_snapshot_when_present(
            optional_getter
        )
