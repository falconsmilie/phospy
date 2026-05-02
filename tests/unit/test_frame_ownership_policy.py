from __future__ import annotations

import inspect
from collections.abc import Iterator
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
from phospy.api.results import KinasePredictionResult, SignalomeWorkflowResult
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
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
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


def test_dataset_public_export_copy_default_is_safe() -> None:
    dataset = AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )

    exported = dataset.to_dataframe()
    exported.iloc[0, 0] = 999.0

    assert float(dataset.phospho.iloc[0, 0]) == 1.0


def test_dataset_public_export_rejects_legacy_copy_keyword() -> None:
    dataset = AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )

    with pytest.raises(TypeError, match="unexpected keyword argument"):
        dataset.to_dataframe(copy=False)  # type: ignore[call-arg]


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


def test_prediction_result_public_export_rejects_legacy_copy_keyword() -> None:
    pred_mat = pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.8],
            "AKT1": [0.2, 0.1],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )
    result = KinasePredictionResult._from_owned(pred_mat=pred_mat)

    with pytest.raises(TypeError, match="unexpected keyword argument"):
        result.to_dataframe(copy=False)  # type: ignore[call-arg]


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


def test_public_export_rejects_legacy_copy_keyword_for_dataset_builder_output() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            organism=Organism.RAT,
        )
    )
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        built.to_dataframe(copy=False)  # type: ignore[call-arg]


def test_public_reference_export_isolated_from_mutation() -> None:
    references = _references()
    exported = references.kinase_substrate_map_dataframe()
    exported.iloc[0, 0] = "CHANGED"

    assert str(references.kinase_substrate_map.iloc[0, 0]) == "MAP2K6"


def test_public_signalome_exports_isolated_from_mutation() -> None:
    result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=_kinase_result(),
            config=build_signalome_config(substrate_support_cutoff=0.5),
        )
    )

    assignments_export = result.module_assignments.to_pandas()
    assignments_export.iloc[0, 0] = "CHANGED"
    assert assignments_export.iloc[0, 0] != result.module_assignments.table.iloc[0, 0]

    modules_export = result.signalome_modules.to_pandas()
    modules_export.iloc[0, 0] = float(modules_export.iloc[0, 0]) + 1.0
    assert modules_export.iloc[0, 0] != result.signalome_modules.table.iloc[0, 0]

    network_export = result.kinase_network.to_pandas()
    network_export.iloc[0, 0] = "CHANGED"
    assert network_export.iloc[0, 0] != result.kinase_network.edges.iloc[0, 0]


def test_public_signalome_and_table_exports_reject_legacy_copy_keyword() -> None:
    result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=_kinase_result(),
            config=build_signalome_config(substrate_support_cutoff=0.5),
        )
    )

    with pytest.raises(TypeError, match="unexpected keyword argument"):
        result.to_dataframe(copy=False)  # type: ignore[call-arg]
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        result.module_assignments.to_pandas(copy=False)  # type: ignore[call-arg]
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        result.signalome_modules.to_pandas(copy=False)  # type: ignore[call-arg]
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        result.kinase_network.to_pandas(copy=False)  # type: ignore[call-arg]

    table = PhosphoIntensityMatrix(frame=_phospho())
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        table.to_pandas(copy=False)  # type: ignore[call-arg]


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

    exported_phospho = dataset.phospho
    exported_phospho.iloc[0, 0] = 999.0
    assert float(dataset.phospho.iloc[0, 0]) == 1.0

    exported_site_metadata = dataset.site_metadata
    exported_site_metadata.iloc[0, 0] = "CHANGED"
    assert str(dataset.site_metadata.iloc[0, 0]) == "MAPK14"

    exported_sample_metadata = dataset.sample_metadata
    assert exported_sample_metadata is not None
    exported_sample_metadata.iloc[0, 0] = 999
    current_sample_metadata = dataset.sample_metadata
    assert current_sample_metadata is not None
    assert int(current_sample_metadata.iloc[0, 0]) == 1


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

    row_counts = report.row_counts
    row_counts.iloc[0, 1] = 999
    assert int(report.row_counts.iloc[0, 1]) == 2

    operations = report.operations
    operations.iloc[0, 2] = "changed"
    assert str(report.operations.iloc[0, 2]) == "forbid"

    row_audit = report.row_audit
    row_audit.iloc[0, 1] = "changed"
    assert str(report.row_audit.iloc[0, 1]) == "retained"

    duplicate_site_resolution = report.duplicate_site_resolution
    assert duplicate_site_resolution is not None
    duplicate_site_resolution.iloc[0, 0] = "changed"
    reread_duplicate = report.duplicate_site_resolution
    assert reread_duplicate is not None
    assert str(reread_duplicate.iloc[0, 0]) == "MAPK14;Y182;"

    metadata_conflicts = report.metadata_conflicts
    assert metadata_conflicts is not None
    metadata_conflicts.iloc[0, 0] = "changed"
    reread_conflicts = report.metadata_conflicts
    assert reread_conflicts is not None
    assert str(reread_conflicts.iloc[0, 0]) == "MAPK14;Y182;"

    comparison_group_stats = report.comparison_group_stats
    assert comparison_group_stats is not None
    comparison_group_stats.iloc[0, 2] = 999
    reread_group_stats = report.comparison_group_stats
    assert reread_group_stats is not None
    assert int(reread_group_stats.iloc[0, 2]) == 1

    comparison_pair_stats = report.comparison_pair_stats
    assert comparison_pair_stats is not None
    comparison_pair_stats.iloc[0, 1] = "changed"
    reread_pair_stats = report.comparison_pair_stats
    assert reread_pair_stats is not None
    assert str(reread_pair_stats.iloc[0, 1]) == "p_group_a_group_b"


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
    profile_scores = scoring_result.profile_scores
    profile_scores.iloc[0, 0] = 999.0
    assert float(scoring_result.profile_scores.iloc[0, 0]) == 0.8

    motif_scores = scoring_result.motif_scores
    assert motif_scores is not None
    motif_scores.iloc[0, 0] = 999.0
    reread_motif = scoring_result.motif_scores
    assert reread_motif is not None
    assert float(reread_motif.iloc[0, 0]) == 0.7

    rank_weighted = scoring_result.rank_weighted_fusion_scores
    assert rank_weighted is not None
    rank_weighted.iloc[0, 0] = 999.0
    reread_rank_weighted = scoring_result.rank_weighted_fusion_scores
    assert reread_rank_weighted is not None
    assert float(reread_rank_weighted.iloc[0, 0]) == 0.75

    fusion_weights = scoring_result.score_fusion_weights
    assert fusion_weights is not None
    fusion_weights.iloc[0, 0] = 999.0
    reread_weights = scoring_result.score_fusion_weights
    assert reread_weights is not None
    assert float(reread_weights.iloc[0, 0]) == 1.0

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
    pred_mat = prediction_result.pred_mat
    pred_mat.iloc[0, 0] = 0.0
    assert float(prediction_result.pred_mat.iloc[0, 0]) == 0.9

    substrate_list = prediction_result.substrate_list
    assert substrate_list is not None
    substrate_list.iloc[0, 0] = "changed"
    reread_substrate_list = prediction_result.substrate_list
    assert reread_substrate_list is not None
    assert str(reread_substrate_list.iloc[0, 0]) == "MAP2K6"

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
    weighted_activity = activity_result.weighted_activity
    weighted_activity.iloc[0, 0] = 999.0
    assert float(activity_result.weighted_activity.iloc[0, 0]) == 1.0

    thresholded_mean = activity_result.thresholded_substrate_mean_activity
    thresholded_mean.iloc[0, 0] = 999.0
    assert float(activity_result.thresholded_substrate_mean_activity.iloc[0, 0]) == 0.5

    target_table = activity_result.target_table
    target_table.iloc[0, 1] = "changed"
    assert str(activity_result.target_table.iloc[0, 1]) == "MAP2K6"


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

    assignments = signalome_result.module_assignments.table
    assignments.iloc[0, 0] = "changed"
    assert (
        assignments.iloc[0, 0] != signalome_result.module_assignments.table.iloc[0, 0]
    )

    modules = signalome_result.signalome_modules.table
    modules.iloc[0, 0] = float(modules.iloc[0, 0]) + 1.0
    assert modules.iloc[0, 0] != signalome_result.signalome_modules.table.iloc[0, 0]

    edges = signalome_result.kinase_network.edges
    edges.iloc[0, 0] = "changed"
    assert edges.iloc[0, 0] != signalome_result.kinase_network.edges.iloc[0, 0]

    nodes = signalome_result.kinase_network.nodes
    if nodes is not None:
        if pd.api.types.is_numeric_dtype(nodes.dtypes.iloc[0]):
            nodes.iloc[0, 0] = float(nodes.iloc[0, 0]) + 1.0
        else:
            nodes.iloc[0, 0] = f"{nodes.iloc[0, 0]}_changed"
        reread_nodes = signalome_result.kinase_network.nodes
        assert reread_nodes is not None
        assert nodes.iloc[0, 0] != reread_nodes.iloc[0, 0]

    candidate_correlations = signalome_result.kinase_network.candidate_correlations
    if candidate_correlations is not None:
        if pd.api.types.is_numeric_dtype(candidate_correlations.dtypes.iloc[0]):
            candidate_correlations.iloc[0, 0] = (
                float(candidate_correlations.iloc[0, 0]) + 1.0
            )
        else:
            candidate_correlations.iloc[0, 0] = (
                f"{candidate_correlations.iloc[0, 0]}_changed"
            )
        reread_candidate = signalome_result.kinase_network.candidate_correlations
        assert reread_candidate is not None
        assert candidate_correlations.iloc[0, 0] != reread_candidate.iloc[0, 0]

    expanded_signalome = signalome_result.expanded_signalome
    assert expanded_signalome is not None
    expanded_signalome.iloc[0, 0] = "changed"
    reread_expanded = signalome_result.expanded_signalome
    assert reread_expanded is not None
    assert expanded_signalome.iloc[0, 0] != reread_expanded.iloc[0, 0]

    site_membership = signalome_result.site_membership
    assert site_membership is not None
    site_membership.iloc[0, 0] = "changed"
    reread_site_membership = signalome_result.site_membership
    assert reread_site_membership is not None
    assert site_membership.iloc[0, 0] != reread_site_membership.iloc[0, 0]

    protein_site_context = signalome_result.protein_site_context
    assert protein_site_context is not None
    protein_site_context.iloc[0, 0] = "changed"
    reread_protein_site_context = signalome_result.protein_site_context
    assert reread_protein_site_context is not None
    assert protein_site_context.iloc[0, 0] != reread_protein_site_context.iloc[0, 0]
