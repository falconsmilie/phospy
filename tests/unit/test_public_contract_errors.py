from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import phospy
import phospy.errors as public_errors
from phospy import (
    AnalysisReadyPhosphoDataset,
    DatasetBuildRequest,
    KinaseWorkflowResult,
    Organism,
    PhosPyInputError,
    PhosPyValidationError,
    ReferenceBundle,
    UnsupportedInputFormatError,
)
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    SignalomeWorkflowResult,
)
from phospy.datasets.builders.validator import DatasetBuildRequestValidator
from phospy.errors import (
    DatasetValidationError,
    ReferenceValidationError,
    WorkflowValidationError,
)
from phospy.io.readers.tables import read_table
from phospy.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeModules,
)
from tests.support.transformation_states import supported_linear_state

TOP_LEVEL_ERROR_FACADE = {
    "PhosPyError",
    "PhosPyInputError",
    "UnsupportedInputFormatError",
    "PhosPyBuildError",
    "PhosPyValidationError",
    "PhosPyReferenceError",
    "UnsupportedOrganismError",
    "PhosPyTransformationError",
    "PhosPyWorkflowError",
    "WorkflowBoundaryError",
}

NON_FACADE_ERROR_TYPES = {
    "DatasetBuildError",
    "DatasetValidationError",
    "ReferenceValidationError",
    "TransformationValidationError",
    "WorkflowValidationError",
    "ReferenceCompatibilityError",
    "ReferenceResolutionError",
    "InvalidTransformationStateError",
    "TransformationStateEstablishmentError",
    "TransformerExecutionError",
    "WorkflowStageError",
}


def _dataset() -> AnalysisReadyPhosphoDataset:
    index = pd.Index(["MAPK14;Y182;"], name="site_id")
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=index),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            },
            index=index,
        ),
        organism=Organism.RAT,
        transformation_state=supported_linear_state(has_total_matrix=False),
    )


def _references() -> ReferenceBundle:
    index = pd.Index(["MAPK14;Y182;"], name="site_id")
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["MAP2K6"], "substrate_site": ["MAPK14;Y182;"]}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"]},
            index=index,
        ),
    )


def _scoring_result() -> KinaseScoringResult:
    return KinaseScoringResult(
        profile_scores=pd.DataFrame(
            {"MAP2K6": [1.0]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        )
    )


def _prediction_result() -> KinasePredictionResult:
    return KinasePredictionResult(
        pred_mat=pd.DataFrame(
            {"MAP2K6": [0.8]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        )
    )


def _kinase_result() -> KinaseWorkflowResult:
    return KinaseWorkflowResult(
        dataset=_dataset(),
        references=_references(),
        scoring_result=_scoring_result(),
        prediction_result=_prediction_result(),
        activity_result=None,
    )


def test_top_level_exception_exports_match_curated_facade() -> None:
    assert TOP_LEVEL_ERROR_FACADE.issubset(set(public_errors.__all__))
    assert TOP_LEVEL_ERROR_FACADE.issubset(set(phospy.__all__))
    for exported in TOP_LEVEL_ERROR_FACADE:
        assert getattr(phospy, exported) is getattr(public_errors, exported)
    for exported in NON_FACADE_ERROR_TYPES:
        assert exported in public_errors.__all__
        assert exported not in phospy.__all__


def test_dataset_build_request_uses_phospy_exception_for_invalid_sources() -> None:
    request = DatasetBuildRequest(
        phospho=object(),
        site_metadata=object(),
    )
    with pytest.raises(
        UnsupportedInputFormatError,
        match="dataset build request phospho must be a pandas DataFrame or a file path",
    ):
        DatasetBuildRequestValidator().run(request)


def test_dataset_constructor_rejects_non_dataframe_with_dataset_validation_error() -> (
    None
):
    with pytest.raises(
        DatasetValidationError,
        match="dataset.phospho must be a pandas DataFrame",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=object(),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "site": ["Y182"],
                    "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                },
                index=["MAPK14;Y182;"],
            ),
            organism=Organism.RAT,
            transformation_state=supported_linear_state(has_total_matrix=False),
        )


def test_dataset_constructor_rejects_blank_gene_symbol_values() -> None:
    with pytest.raises(
        DatasetValidationError,
        match="dataset.site_metadata.gene_symbol must contain non-empty string values",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"]),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["  "],
                    "site": ["Y182"],
                    "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                },
                index=["MAPK14;Y182;"],
            ),
            organism=Organism.RAT,
            transformation_state=supported_linear_state(has_total_matrix=False),
        )


def test_dataset_constructor_rejects_blank_site_values() -> None:
    with pytest.raises(
        DatasetValidationError,
        match="dataset.site_metadata.site must contain non-empty string values",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"]),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "site": ["\t"],
                    "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                },
                index=["MAPK14;Y182;"],
            ),
            organism=Organism.RAT,
            transformation_state=supported_linear_state(has_total_matrix=False),
        )


def test_dataset_constructor_allows_missing_site_sequence_column() -> None:
    dataset = AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"]),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
            },
            index=["MAPK14;Y182;"],
        ),
        organism=Organism.RAT,
        transformation_state=supported_linear_state(has_total_matrix=False),
    )
    assert "site_sequence" not in dataset.site_metadata.columns


def test_dataset_constructor_rejects_blank_site_sequence_values() -> None:
    with pytest.raises(
        DatasetValidationError,
        match="dataset.site_metadata.site_sequence must contain non-empty string values",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"]),
            site_metadata=pd.DataFrame(
                {
                    "gene_symbol": ["MAPK14"],
                    "site": ["Y182"],
                    "site_sequence": [""],
                },
                index=["MAPK14;Y182;"],
            ),
            organism=Organism.RAT,
            transformation_state=supported_linear_state(has_total_matrix=False),
        )


def test_reference_bundle_constructor_rejects_non_dataframe_with_reference_validation_error() -> (
    None
):
    with pytest.raises(
        ReferenceValidationError,
        match="references.kinase_substrate_map must be a pandas DataFrame",
    ):
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=object(),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
        )


def test_result_containers_use_phospy_validation_errors_for_dataframe_fields() -> None:
    with pytest.raises(
        PhosPyValidationError,
        match="scoring_result.profile_scores must be a pandas DataFrame",
    ):
        KinaseScoringResult(profile_scores=object())


def test_workflow_results_are_typed_containers_not_nested_type_validators() -> None:
    dataset = object()
    references = object()
    scoring_result = object()
    prediction_result = object()
    result = KinaseWorkflowResult(
        dataset=dataset,
        references=references,
        scoring_result=scoring_result,
        prediction_result=prediction_result,
    )
    assert result.dataset is dataset
    assert result.references is references
    assert result.scoring_result is scoring_result
    assert result.prediction_result is prediction_result


def test_signalome_result_validates_expanded_signalome_field_type() -> None:
    kinase_result = _kinase_result()
    with pytest.raises(
        WorkflowValidationError,
        match="signalome_result.expanded_signalome must be a pandas DataFrame",
    ):
        SignalomeWorkflowResult(
            dataset=kinase_result.dataset,
            kinase_result=kinase_result,
            module_assignments=SignalomeAssignments(
                table=pd.DataFrame({"site_id": ["MAPK14;Y182;"], "module": [1]})
            ),
            signalome_modules=SignalomeModules(
                table=pd.DataFrame({"module": [1], "size": [1]})
            ),
            kinase_network=KinaseNetwork(
                edges=pd.DataFrame(
                    {"source": ["MAP2K6"], "target": ["MAP2K6"], "weight": [1.0]}
                )
            ),
            expanded_signalome=[],
        )


def test_table_reader_translates_missing_path_to_phospy_input_error(
    tmp_path: Path,
) -> None:
    with pytest.raises(
        PhosPyInputError,
        match="input file does not exist:",
    ):
        read_table(tmp_path / "missing.csv")
