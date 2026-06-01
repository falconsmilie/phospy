from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import phospy
import phospy.api as public_api
import phospy.errors as public_errors
from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
)
from phospy.api import (
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
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
from phospy.errors import (
    DatasetValidationError,
    ReferenceValidationError,
    WorkflowValidationError,
)
from phospy.io.readers.tables import read_table
from phospy.science.datasets.builders.validator import DatasetBuildRequestValidator
from phospy.science.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeModules,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import (
    site_key_from_display_id,
    site_key_index_from_display_ids,
)

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
    "SignalomeScaleError",
}


def _supported_dataset_state(*, has_total_matrix: bool) -> dict[str, object]:
    return {
        "intensity_scale_state": supported_linear_intensity_scale_state(
            has_total_matrix=has_total_matrix
        ),
        "processing_state": supported_linear_processing_state(
            has_total_matrix=has_total_matrix
        ),
    }


def _dataset() -> AnalysisReadyPhosphoDataset:
    display_id = "MAPK14;Y182;"
    index = site_key_index_from_display_ids([display_id])
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame({"sample_a": [1.0]}, index=index),
        site_metadata=pd.DataFrame(
            {
                "site_key": index.astype(str).tolist(),
                "display_id": [display_id],
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                "localisation_confidence": [0.95],
            },
            index=index,
        ),
        organism=Organism.RAT,
        **_supported_dataset_state(has_total_matrix=False),
    )


def _coherent_site_identity_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    display_ids = ["MAPK14;Y182;", "AKT1;T308;"]
    index = site_key_index_from_display_ids(display_ids)
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [3.0, 4.0],
        },
        index=index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": index.astype(str).tolist(),
            "display_id": display_ids,
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "T308"]
            ],
            "localisation_confidence": [0.95, 0.9],
        },
        index=index,
    )
    return phospho, site_metadata


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
    site_key = site_key_from_display_id("MAPK14;Y182;")
    return KinaseScoringResult(
        profile_scores=pd.DataFrame(
            {"MAP2K6": [1.0]},
            index=pd.Index([site_key], name="site_key"),
        )
    )


def _prediction_result() -> KinasePredictionResult:
    site_key = site_key_from_display_id("MAPK14;Y182;")
    return KinasePredictionResult(
        pred_mat=pd.DataFrame(
            {"MAP2K6": [0.8]},
            index=pd.Index([site_key], name="site_key"),
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


def _valid_signalome_assignments_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "protein_id": ["MAPK14"],
            "module_id": [1],
            "top_kinase": ["MAP2K6"],
            "top_score": [0.8],
            "top_kinase_candidates": [("MAP2K6",)],
            "top_kinase_weights": [(("MAP2K6", 1.0),)],
            "top_kinase_tie_count": [1],
            "top_kinase_is_ambiguous": [False],
            "top_kinase_selection_policy": ["max_score_then_lexicographic_tiebreak"],
            "module_top_kinase": ["MAP2K6"],
            "module_top_kinase_candidates": [("MAP2K6",)],
            "module_top_kinase_tie_count": [1],
            "module_top_kinase_is_ambiguous": [False],
            "module_top_kinase_selection_policy": [
                "max_score_then_lexicographic_tiebreak"
            ],
        },
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )


def _valid_signalome_modules_table() -> pd.DataFrame:
    return pd.DataFrame(
        {"MAP2K6": [100.0]},
        index=pd.Index([1], name="module_id"),
    )


def _valid_kinase_network_edges_table() -> pd.DataFrame:
    return pd.DataFrame(
        {"source_kinase": ["MAP2K6"], "target_kinase": ["MAP2K6"], "correlation": [1.0]}
    )


def test_top_level_exception_exports_match_curated_facade() -> None:
    assert TOP_LEVEL_ERROR_FACADE.issubset(set(public_errors.__all__))
    assert TOP_LEVEL_ERROR_FACADE.issubset(set(public_api.__all__))
    for exported in TOP_LEVEL_ERROR_FACADE:
        assert getattr(public_api, exported) is getattr(public_errors, exported)
        assert exported not in phospy.__all__
    for exported in NON_FACADE_ERROR_TYPES:
        assert exported in public_errors.__all__
        assert exported in public_api.__all__
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


def test_builder_site_matrix_reports_no_retained_rows_when_all_sequence_support_is_missing() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=pd.Index(["FAKE1;S1;", "FAKE2;T2;"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["FAKE1", "FAKE2"],
            "site": ["S1", "T2"],
            "localisation_confidence": [0.95, 0.9],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "site-matrix construction produced no retained rows after filtering; "
            "input_rows=2, dropped_missing_sequence=2"
        ),
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                preprocessing_config=DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
                ),
            )
        )


@pytest.mark.parametrize(
    "missing_data_policy",
    ("retain_missing", "require_min_observed_values"),
)
def test_builder_rejects_internal_only_site_matrix_missing_modes_at_public_boundary(
    missing_data_policy: str,
) -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [float("nan")]},
        index=pd.Index(["ROW_1"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            "localisation_confidence": [0.95],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            f"missing_data_policy='{missing_data_policy}' is not supported for strict "
            "AnalysisReadyPhosphoDataset construction"
        ),
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                preprocessing_config=DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(
                        policy="build_from_metadata",
                        missing_data_policy=missing_data_policy,  # type: ignore[arg-type]
                    )
                ),
            )
        )


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
            **_supported_dataset_state(has_total_matrix=False),
        )


def test_dataset_constructor_rejects_blank_gene_symbol_values() -> None:
    site_key = site_key_from_display_id("MAPK14;Y182;")
    index = pd.Index([site_key], name="site_key")
    with pytest.raises(
        DatasetValidationError,
        match="dataset.site_metadata.gene_symbol must contain non-empty string values",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame({"sample_a": [1.0]}, index=index),
            site_metadata=pd.DataFrame(
                {
                    "site_key": [site_key],
                    "display_id": ["MAPK14;Y182;"],
                    "gene_symbol": ["  "],
                    "site": ["Y182"],
                    "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                },
                index=index.copy(),
            ),
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


def test_dataset_constructor_rejects_blank_site_values() -> None:
    site_key = site_key_from_display_id("MAPK14;Y182;")
    index = pd.Index([site_key], name="site_key")
    with pytest.raises(
        DatasetValidationError,
        match="dataset.site_metadata.site must contain non-empty string values",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame({"sample_a": [1.0]}, index=index),
            site_metadata=pd.DataFrame(
                {
                    "site_key": [site_key],
                    "display_id": ["MAPK14;Y182;"],
                    "gene_symbol": ["MAPK14"],
                    "site": ["\t"],
                    "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                },
                index=index.copy(),
            ),
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


def test_dataset_constructor_accepts_site_identity_coherence() -> None:
    phospho, site_metadata = _coherent_site_identity_inputs()
    dataset = AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        **_supported_dataset_state(has_total_matrix=False),
    )
    assert list(dataset.phospho.index) == list(phospho.index)


def test_dataset_constructor_rejects_site_identity_gene_symbol_mismatch() -> None:
    phospho, site_metadata = _coherent_site_identity_inputs()
    site_metadata.loc[site_metadata.index[0], "gene_symbol"] = "MAPK1"
    with pytest.raises(
        DatasetValidationError,
        match="site-identity coherence failed",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


def test_dataset_constructor_rejects_site_identity_site_mismatch() -> None:
    phospho, site_metadata = _coherent_site_identity_inputs()
    site_metadata.loc[site_metadata.index[0], "site"] = "T185"
    with pytest.raises(
        DatasetValidationError,
        match="site_sequence central residue must agree with site/residue metadata",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


def test_dataset_constructor_rejects_when_one_row_has_site_identity_mismatch() -> None:
    phospho, site_metadata = _coherent_site_identity_inputs()
    akt1_site_key = site_metadata.index[1]
    site_metadata.loc[akt1_site_key, "site"] = "S473"
    with pytest.raises(
        DatasetValidationError,
        match="site_sequence central residue must agree with site/residue metadata",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


def test_dataset_constructor_rejects_missing_site_sequence_column() -> None:
    site_key = site_key_from_display_id("MAPK14;Y182;")
    index = pd.Index([site_key], name="site_key")
    with pytest.raises(
        DatasetValidationError,
        match="dataset.site_metadata is missing required columns: site_sequence",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame({"sample_a": [1.0]}, index=index),
            site_metadata=pd.DataFrame(
                {
                    "site_key": [site_key],
                    "display_id": ["MAPK14;Y182;"],
                    "gene_symbol": ["MAPK14"],
                    "site": ["Y182"],
                },
                index=index.copy(),
            ),
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


def test_dataset_constructor_rejects_blank_site_sequence_values() -> None:
    site_key = site_key_from_display_id("MAPK14;Y182;")
    index = pd.Index([site_key], name="site_key")
    with pytest.raises(
        DatasetValidationError,
        match="dataset.site_metadata.site_sequence must contain non-empty string values",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame({"sample_a": [1.0]}, index=index),
            site_metadata=pd.DataFrame(
                {
                    "site_key": [site_key],
                    "display_id": ["MAPK14;Y182;"],
                    "gene_symbol": ["MAPK14"],
                    "site": ["Y182"],
                    "site_sequence": [""],
                },
                index=index.copy(),
            ),
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
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
                table=_valid_signalome_assignments_table()
            ),
            signalome_modules=SignalomeModules(table=_valid_signalome_modules_table()),
            kinase_network=KinaseNetwork(edges=_valid_kinase_network_edges_table()),
            expanded_signalome=[],
        )


def test_signalome_result_validates_site_membership_field_type() -> None:
    kinase_result = _kinase_result()
    with pytest.raises(
        WorkflowValidationError,
        match="signalome_result.site_membership must be a pandas DataFrame",
    ):
        SignalomeWorkflowResult(
            dataset=kinase_result.dataset,
            kinase_result=kinase_result,
            module_assignments=SignalomeAssignments(
                table=_valid_signalome_assignments_table()
            ),
            signalome_modules=SignalomeModules(table=_valid_signalome_modules_table()),
            kinase_network=KinaseNetwork(edges=_valid_kinase_network_edges_table()),
            site_membership=[],
        )


def test_signalome_result_validates_protein_site_context_field_type() -> None:
    kinase_result = _kinase_result()
    with pytest.raises(
        WorkflowValidationError,
        match="signalome_result.protein_site_context must be a pandas DataFrame",
    ):
        SignalomeWorkflowResult(
            dataset=kinase_result.dataset,
            kinase_result=kinase_result,
            module_assignments=SignalomeAssignments(
                table=_valid_signalome_assignments_table()
            ),
            signalome_modules=SignalomeModules(table=_valid_signalome_modules_table()),
            kinase_network=KinaseNetwork(edges=_valid_kinase_network_edges_table()),
            protein_site_context=[],
        )


def test_table_reader_translates_missing_path_to_phospy_input_error(
    tmp_path: Path,
) -> None:
    with pytest.raises(
        PhosPyInputError,
        match="input file does not exist:",
    ):
        read_table(tmp_path / "missing.csv")
