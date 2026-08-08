from __future__ import annotations

from collections.abc import Callable

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.advanced.configs import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api import (
    Contrast,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    Organism,
    ReferencePreset,
    SampleDesignRecord,
    SignalomeWorkflowRequest,
)
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
)
from phospy.errors import (
    DatasetValidationError,
    ReferenceCompatibilityError,
    WorkflowValidationError,
)
from phospy.science.references.models import ReferenceBundle
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
from phospy.workflows.signalome.validator import SignalomeWorkflowValidator
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.signalome_config import build_signalome_config
from tests.support.site_keys import (
    protein_site_key_index,
    site_key_context_columns,
    site_key_index_from_display_ids,
)
from tests.support.unsafe_dataset_states import (
    unsafe_corrupt_dataset_to_display_index,
    unsafe_drop_dataset_site_metadata_columns,
    unsafe_reverse_dataset_site_metadata_index,
)

_WorkflowRunner = Callable[[AnalysisReadyPhosphoDataset], object]
_DatasetFactory = Callable[[], AnalysisReadyPhosphoDataset]


def _site_keys(display_ids: list[str]) -> pd.Index:
    return site_key_index_from_display_ids(
        display_ids,
        protein_namespace="gene_symbol",
    )


def _differential_dataset(
    *, allow_opaque_site_values: bool = False
) -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;", "AKT1;T308;"]
    site_ids = _site_keys(display_ids)
    metadata_sites = ["Y182", "T308"]
    sequence_values = [
        ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
        for site in metadata_sites
    ]
    if allow_opaque_site_values:
        display_ids = ["MAPK14;FOO;", "AKT1;BAR;"]
        site_ids = _site_keys(["MAPK14;Y182;", "AKT1;T308;"])
        metadata_sites = ["FOO", "BAR"]
        sequence_values = [("A" * 15) + "Y" + ("A" * 15)] * 2
    return trusted_analysis_ready_dataset_from_tables(
        phospho=pd.DataFrame(
            {
                "A_1": [1.0, 2.0],
                "A_2": [1.1, 2.1],
                "B_1": [2.1, 2.0],
                "B_2": [2.0, 2.2],
            },
            index=site_ids.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_ids.tolist(),
                "display_id": display_ids,
                **site_key_context_columns(site_ids),
                "gene_symbol": ["MAPK14", "AKT1"],
                "site": metadata_sites,
                "site_sequence": sequence_values,
            },
            index=site_ids.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
        allow_opaque_site_values=allow_opaque_site_values,
    )


def _differential_request(
    *, dataset: AnalysisReadyPhosphoDataset | None = None
) -> DifferentialAnalysisRequest:
    resolved = _differential_dataset() if dataset is None else dataset
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A"),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B"),
            SampleDesignRecord(sample_id="B_2", condition="B"),
        )
    )
    return DifferentialAnalysisRequest(
        dataset=resolved,
        design=design,
        contrasts=(
            Contrast(name="B_vs_A", numerator_condition="B", denominator_condition="A"),
        ),
    )


def _kinase_dataset() -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;", "AKT1;T308;"]
    site_ids = _site_keys(display_ids)
    return trusted_analysis_ready_dataset_from_tables(
        phospho=pd.DataFrame(
            {"sample_a": [1.0, 2.0], "sample_b": [1.1, 2.1]},
            index=site_ids.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_ids.tolist(),
                "display_id": display_ids,
                **site_key_context_columns(site_ids),
                "gene_symbol": ["MAPK14", "AKT1"],
                "site": ["Y182", "T308"],
                "site_sequence": ["AAAAAAAYAAAAAAA", "AAAAAAATAAAAAAA"],
                "protein_id": ["P28482", "P31749"],
                "protein_accession": ["P28482-1", "P31749-1"],
            },
            index=site_ids.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _references() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "AKT1"],
                "substrate_site": ["MAPK14;Y182;", "AKT1;T308;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["AAAAAAAYAAAAAAA", "AAAAAAATAAAAAAA"]},
            index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
        ),
    )


def _signalome_request(
    *, dataset: AnalysisReadyPhosphoDataset
) -> SignalomeWorkflowRequest:
    site_ids = dataset.phospho.index.copy()
    score_matrix = pd.DataFrame({"MAP2K6": [1.0] * len(site_ids)}, index=site_ids)
    prediction_matrix = pd.DataFrame({"MAP2K6": [0.8] * len(site_ids)}, index=site_ids)
    return SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=dataset,
            references=_references(),
            scoring_result=KinaseScoringResult._from_owned(
                profile_scores=score_matrix,
                rank_weighted_fusion_scores=score_matrix,
            ),
            prediction_result=KinasePredictionResult._from_owned(
                pred_mat=prediction_matrix
            ),
            activity_result=None,
        ),
        config=build_signalome_config(
            substrate_support_cutoff=0.5,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
    )


def _kinase_request(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    references: ReferencePreset | ReferenceBundle | None = None,
) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=dataset,
        references=_references() if references is None else references,
        scoring_config=KinaseScoringConfig(
            reliability_profile="custom",
            min_substrates=2,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=5,
            deterministic_max_selected_kinases=5,
            adaptive_ensemble_runs=5,
        ),
        activity_config=None,
    )


def _run_differential_validator(dataset: AnalysisReadyPhosphoDataset) -> object:
    return DifferentialAnalysisValidator().run(_differential_request(dataset=dataset))


def _run_kinase_validator(dataset: AnalysisReadyPhosphoDataset) -> object:
    return KinaseWorkflowValidator().run(_kinase_request(dataset=dataset))


def _run_signalome_validator(dataset: AnalysisReadyPhosphoDataset) -> object:
    return SignalomeWorkflowValidator().run(_signalome_request(dataset=dataset))


def _workflow_cases() -> tuple[tuple[str, _WorkflowRunner, _DatasetFactory], ...]:
    return (
        ("differential", _run_differential_validator, _differential_dataset),
        ("kinase", _run_kinase_validator, _kinase_dataset),
        ("signalome", _run_signalome_validator, _kinase_dataset),
    )


def _make_display_indexed(dataset: AnalysisReadyPhosphoDataset) -> None:
    unsafe_corrupt_dataset_to_display_index(dataset)


def _drop_site_metadata_column(
    dataset: AnalysisReadyPhosphoDataset, column_name: str
) -> None:
    unsafe_drop_dataset_site_metadata_columns(dataset, column_name)


def _reverse_site_metadata_index(dataset: AnalysisReadyPhosphoDataset) -> None:
    unsafe_reverse_dataset_site_metadata_index(dataset)


def _duplicate_display_differential_dataset() -> AnalysisReadyPhosphoDataset:
    site_ids = protein_site_key_index(
        protein_identifiers=["P28482", "Q99999"],
        sites=["Y182", "Y182"],
    )
    return trusted_analysis_ready_dataset_from_tables(
        phospho=pd.DataFrame(
            {
                "A_1": [1.0, 2.0],
                "A_2": [1.1, 2.1],
                "B_1": [2.1, 2.0],
                "B_2": [2.0, 2.2],
            },
            index=site_ids.copy(),
        ),
        site_metadata=_duplicate_display_site_metadata(site_ids),
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


def _duplicate_display_kinase_dataset() -> AnalysisReadyPhosphoDataset:
    site_ids = protein_site_key_index(
        protein_identifiers=["P28482", "Q99999"],
        sites=["Y182", "Y182"],
    )
    return trusted_analysis_ready_dataset_from_tables(
        phospho=pd.DataFrame(
            {"sample_a": [1.0, 2.0], "sample_b": [1.1, 2.1]},
            index=site_ids.copy(),
        ),
        site_metadata=_duplicate_display_site_metadata(site_ids),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _duplicate_display_site_metadata(site_ids: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_key": site_ids.astype(str).tolist(),
            "display_id": ["MAPK14;Y182;", "MAPK14;Y182;"],
            **site_key_context_columns(site_ids),
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["AAAAAAAYAAAAAAA", "AAAAAAAYAAAAAAA"],
            "protein_id": ["P28482", "Q99999"],
            "protein_accession": ["P28482-1", "Q99999-1"],
            "localisation_confidence": [0.95, 0.95],
        },
        index=site_ids.copy(),
    )


def test_differential_validator_accepts_valid_site_key_dataset() -> None:
    validated = DifferentialAnalysisValidator().run(
        _differential_request(dataset=_differential_dataset())
    )
    site_metadata = validated.dataset._borrow_site_metadata_frame()
    assert site_metadata.index.name == "site_key"
    assert site_metadata.index.tolist() == site_metadata.loc[:, "site_key"].tolist()


def test_kinase_validator_accepts_valid_site_key_dataset() -> None:
    request = _kinase_request(dataset=_kinase_dataset())
    validated = KinaseWorkflowValidator().run(request)
    assert validated.request is request


def test_signalome_validator_accepts_valid_site_key_dataset() -> None:
    request = _signalome_request(dataset=_kinase_dataset())
    validated = SignalomeWorkflowValidator().run(request)
    assert validated.request is request


@pytest.mark.parametrize(
    ("workflow_name", "run_workflow", "dataset_factory"),
    _workflow_cases(),
)
def test_workflow_validators_reject_display_indexed_datasets(
    workflow_name: str,
    run_workflow: _WorkflowRunner,
    dataset_factory: _DatasetFactory,
) -> None:
    dataset = dataset_factory()
    _make_display_indexed(dataset)

    with pytest.raises(WorkflowValidationError, match="display-indexed"):
        run_workflow(dataset)


@pytest.mark.parametrize(
    ("workflow_name", "run_workflow", "dataset_factory"),
    _workflow_cases(),
)
@pytest.mark.parametrize(
    ("column_name", "pattern"),
    (
        ("display_id", "missing required columns: display_id"),
        ("site_key", "missing required columns: site_key"),
    ),
)
def test_workflow_validators_require_site_key_and_display_id_columns(
    workflow_name: str,
    run_workflow: _WorkflowRunner,
    dataset_factory: _DatasetFactory,
    column_name: str,
    pattern: str,
) -> None:
    dataset = dataset_factory()
    _drop_site_metadata_column(dataset, column_name)

    with pytest.raises(WorkflowValidationError, match=pattern):
        run_workflow(dataset)


@pytest.mark.parametrize(
    "column_name",
    ("organism", "protein_namespace", "protein_identifier"),
)
def test_differential_validator_requires_protein_scoped_identity_context(
    column_name: str,
) -> None:
    dataset = _differential_dataset()
    _drop_site_metadata_column(dataset, column_name)

    with pytest.raises(
        WorkflowValidationError,
        match=f"missing required columns: {column_name}",
    ):
        _run_differential_validator(dataset)


@pytest.mark.parametrize("column_name", ("site_key", "display_id"))
def test_signalome_site_key_and_display_id_fail_as_identity_contract(
    column_name: str,
) -> None:
    dataset = _kinase_dataset()
    _drop_site_metadata_column(dataset, column_name)

    with pytest.raises(WorkflowValidationError) as exc_info:
        _run_signalome_validator(dataset)

    message = str(exc_info.value)
    assert (
        "identity requirement failed (contract=protein_scoped_site_identity)" in message
    )
    assert f"missing required columns: {column_name}" in message
    assert "protein grouping metadata requirement failed" not in message


@pytest.mark.parametrize(
    ("workflow_name", "run_workflow", "dataset_factory"),
    _workflow_cases(),
)
def test_workflow_validators_reject_site_key_index_mismatch(
    workflow_name: str,
    run_workflow: _WorkflowRunner,
    dataset_factory: _DatasetFactory,
) -> None:
    dataset = dataset_factory()
    _reverse_site_metadata_index(dataset)

    with pytest.raises(WorkflowValidationError, match="site_key must exactly match"):
        run_workflow(dataset)


@pytest.mark.parametrize(
    ("workflow_name", "run_workflow", "dataset_factory"),
    (
        (
            "differential",
            _run_differential_validator,
            _duplicate_display_differential_dataset,
        ),
        ("kinase", _run_kinase_validator, _duplicate_display_kinase_dataset),
        ("signalome", _run_signalome_validator, _duplicate_display_kinase_dataset),
    ),
)
def test_workflow_validators_allow_duplicate_display_ids_when_site_keys_differ(
    workflow_name: str,
    run_workflow: _WorkflowRunner,
    dataset_factory: _DatasetFactory,
) -> None:
    dataset = dataset_factory()
    site_metadata = dataset._borrow_site_metadata_frame()
    assert site_metadata.loc[:, "display_id"].duplicated().any()
    assert site_metadata.loc[:, "site_key"].is_unique

    run_workflow(dataset)


def test_signalome_validator_requires_prediction_matrix_site_key_alignment() -> None:
    request = _signalome_request(dataset=_kinase_dataset())
    reversed_index = pd.Index(
        list(reversed(request.kinase_result.dataset.phospho.index.astype(str))),
        name="site_key",
    )
    request.kinase_result.prediction_result._pred_mat.index = reversed_index

    with pytest.raises(
        WorkflowValidationError,
        match="prediction_result\\.pred_mat\\.index must exactly match",
    ):
        SignalomeWorkflowValidator().run(request)


def test_signalome_validator_requires_score_matrix_site_key_alignment() -> None:
    request = _signalome_request(dataset=_kinase_dataset())
    reversed_index = pd.Index(
        list(reversed(request.kinase_result.dataset.phospho.index.astype(str))),
        name="site_key",
    )
    request.kinase_result.scoring_result._rank_weighted_fusion_scores.index = (
        reversed_index
    )

    with pytest.raises(
        WorkflowValidationError,
        match="rank_weighted_fusion_scores\\.index must exactly match",
    ):
        SignalomeWorkflowValidator().run(request)


def test_differential_identity_contract_rejects_opaque_sites_even_with_explicit_opt_in() -> (
    None
):
    with pytest.raises(DatasetValidationError, match="strict 'S/T/Y<position>' tokens"):
        _differential_dataset(allow_opaque_site_values=True)


def test_signalome_validator_requires_signalome_protein_grouping_metadata() -> None:
    source = _kinase_dataset()
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=source.phospho,
        site_metadata=source.site_metadata.drop(columns=["protein_id"]),
        sample_metadata=source.sample_metadata,
        total=source.total,
        organism=source.organism,
        intensity_scale_state=source.intensity_scale_state,
        processing_state=source.processing_state,
    )
    request = _signalome_request(dataset=dataset)

    assert "protein_id" not in dataset.site_metadata.columns
    with pytest.raises(WorkflowValidationError) as exc_info:
        SignalomeWorkflowValidator().run(request)

    message = str(exc_info.value)
    assert "signalome protein grouping metadata requirement failed" in message
    assert "protein grouping metadata" in message
    assert "dataset.site_metadata.protein_group_id" in message
    assert "legacy dataset.site_metadata.protein_id" in message
    assert "gene_symbol" in message
    assert "display_id" in message
    assert "identity requirement failed" not in message


def test_signalome_validator_accepts_legacy_protein_id_grouping_alias() -> None:
    request = _signalome_request(dataset=_kinase_dataset())

    validated = SignalomeWorkflowValidator().run(request)

    assert validated.request is request


def test_signalome_validator_rejects_conflicting_grouping_alias_values() -> None:
    source = _kinase_dataset()
    site_metadata = source.site_metadata.copy(deep=True)
    site_metadata.loc[:, "protein_group_id"] = ["GROUP_A", "GROUP_B"]
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=source.phospho,
        site_metadata=site_metadata,
        sample_metadata=source.sample_metadata,
        total=source.total,
        organism=source.organism,
        intensity_scale_state=source.intensity_scale_state,
        processing_state=source.processing_state,
    )

    with pytest.raises(WorkflowValidationError) as exc_info:
        _run_signalome_validator(dataset)

    message = str(exc_info.value)
    assert "Conflicting signalome protein grouping metadata" in message
    assert "protein_group_id" in message
    assert "legacy alias" in message
    assert "protein_id" in message
    assert "protein_identifier" in message


def test_kinase_workflow_rejects_reference_organism_mismatch_where_applicable() -> None:
    request = _kinase_request(
        dataset=_kinase_dataset(),
        references=ReferencePreset.HUMAN,
    )
    with pytest.raises(ReferenceCompatibilityError):
        KinaseWorkflow().run(request)
