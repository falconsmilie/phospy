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

    assert copied_result.pred_mat is not pred_mat
    assert copied_result.substrate_list is not substrate_list
    assert owned_result.pred_mat is pred_mat
    assert owned_result.substrate_list is substrate_list

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
