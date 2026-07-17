from __future__ import annotations

import inspect
import re
from collections.abc import Mapping
from pathlib import Path

import pandas as pd
import pytest

import phospy.io as phospy_io
from phospy.api import AnalysisReadyDatasetBuilder
from phospy.api.configs import (
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
)
from phospy.api.requests import DatasetBuildRequest
from phospy.api.results import KinasePredictionResult
from phospy.errors import (
    DatasetValidationError,
    PhosPyInputError,
    PhosPyValidationError,
    ReferenceCompatibilityError,
    ReferenceResolutionError,
    ReferenceValidationError,
    TransformationValidationError,
    UnsupportedInputFormatError,
    UnsupportedOrganismError,
)
from phospy.provenance import (
    TrustedDatasetConstructionAssertions,
    TrustedDatasetConstructionEvidence,
)
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.references.models import Organism, ReferenceBundle, ReferencePreset
from phospy.science.references.resolution import ReferenceResolver
from phospy.science.transformations.models import (
    MatrixIntensityScaleState,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import protein_site_key, site_key_context_columns

ROOT = Path(__file__).resolve().parents[2]
_SITE_KEY = protein_site_key(protein_identifier="MAPK14", site="Y182")
_SITE_INDEX = pd.Index([_SITE_KEY], name="site_key")
_DISPLAY_ID = "MAPK14;Y182;"
_AKT1_T308_KEY = protein_site_key(protein_identifier="AKT1", site="T308")
_REFERENCE_IDENTIFIER_GUARD_TARGETS = (
    "src/phospy/workflows/kinase/contracts.py",
    "src/phospy/workflows/kinase/scoring_runner.py",
    "src/phospy/workflows/kinase/prediction_runner.py",
    "src/phospy/workflows/kinase/activity_runner.py",
    "src/phospy/science/prediction/scoring.py",
    "src/phospy/science/activities/scoring.py",
)
_FORBIDDEN_IDENTIFIER_CLEANUP = re.compile(r"\.(?:upper|lower|strip)\(")
_REFERENCE_IDENTIFIER_COLUMN_HINT = re.compile(
    r"['\"](?:kinase|substrate_site)['\"]",
    re.IGNORECASE,
)
_REFERENCE_IDENTIFIER_BOUNDARY_OWNER = "src/phospy/science/references/identifiers.py"
_PRODUCTION_REFERENCE_TABLES = "src/phospy/tables/references.py"
_FORBIDDEN_PROTEIN_ACCESSION_BOUNDARY_LEAK_TARGETS = (
    "src/phospy/science/sequences",
    "src/phospy/science/scoring",
    "src/phospy/workflows",
    "src/phospy/science/datasets/preprocessing",
)
_FORBIDDEN_ACCESSION_CASE_NORMALISATION = re.compile(
    r"accession[^\n]{0,120}\.(?:upper|lower)\(",
    re.IGNORECASE,
)
_WORKFLOW_BORROW_METHOD_PATTERN = re.compile(r"\._borrow_[A-Za-z0-9_]*\s*\(")
_WORKFLOW_BORROW_GUARD_DIRS = (
    "src/phospy/workflows",
    "src/phospy/validation/workflows",
)
_DATASET_INTERNAL_VIEW_FRAME_PROPERTIES = {
    "comparisons",
    "phospho",
    "sample_metadata",
    "site_metadata",
    "total",
}
_DATASET_INTERNAL_VIEW_PUBLIC_MEMBERS = _DATASET_INTERNAL_VIEW_FRAME_PROPERTIES | {
    "aggregate_imputation_observation_mask",
    "imputation_observation_summary",
}


def _phospho() -> pd.DataFrame:
    return pd.DataFrame({"sample_a": [1.0]}, index=_SITE_INDEX.copy())


def _site_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_key": [_SITE_KEY],
            "display_id": [_DISPLAY_ID],
            **site_key_context_columns(_SITE_INDEX),
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
            "protein_id": ["MAPK14"],
            "localisation_confidence": [0.95],
        },
        index=_SITE_INDEX.copy(),
    )


def _site_key_for_display_id(site_metadata: pd.DataFrame, display_id: str) -> str:
    matches = site_metadata.index[
        site_metadata.loc[:, "display_id"].astype(str) == display_id
    ].astype(str)
    assert len(matches) == 1
    return str(matches[0])


def _as_builder_input(
    frame: pd.DataFrame,
    *,
    use_file_path: bool,
    tmp_path,
    filename: str,
):
    if not use_file_path:
        return frame
    path = tmp_path / filename
    frame.to_csv(path)
    return path


def _valid_bundle(organism: Organism = Organism.RAT) -> ReferenceBundle:
    return ReferenceBundle(
        organism=organism,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["MAP2K6"], "substrate_site": ["MAPK14;Y182;"]}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
    )


def _supported_dataset_state(*, has_total_matrix: bool) -> dict[str, object]:
    return {
        "intensity_scale_state": supported_linear_intensity_scale_state(
            has_total_matrix=has_total_matrix
        ),
        "processing_state": supported_linear_processing_state(
            has_total_matrix=has_total_matrix
        ),
    }


def _trusted_assertions() -> TrustedDatasetConstructionAssertions:
    return TrustedDatasetConstructionAssertions(
        identity=TrustedDatasetConstructionEvidence.evidence(
            source="unit-test site_key fixtures"
        ),
        intensity_scale=TrustedDatasetConstructionEvidence.evidence(
            source="unit-test intensity scale state",
            policy="require_established_intensity_scale_state",
        ),
        quantitative_meaning=TrustedDatasetConstructionEvidence.evidence(
            source="unit-test intensity matrix",
            policy="linear analysis-ready matrix",
        ),
        aligned_structure=TrustedDatasetConstructionEvidence.evidence(
            source="unit-test aligned phospho/site metadata fixtures",
            policy="require_matching_site_key_index",
        ),
        localisation=TrustedDatasetConstructionEvidence.evidence(
            source="localisation_confidence column",
            policy="require_threshold",
            threshold=0.75,
        ),
        sequence=TrustedDatasetConstructionEvidence.evidence(
            source="unit-test site_sequence fixtures"
        ),
        reference_context=TrustedDatasetConstructionEvidence.waiver(
            reason="unit-test fixtures do not carry reference context"
        ),
        asserted_by="unit-test",
        assertion_source="test_domain_boundaries",
    )


def test_analysis_ready_from_trusted_tables_enforces_site_sequence() -> None:
    bad_site_metadata = _site_metadata().drop(columns=["site_sequence"])
    with pytest.raises(
        DatasetValidationError,
        match="dataset.site_metadata is missing required columns: site_sequence",
    ):
        AnalysisReadyPhosphoDataset.from_trusted_tables(
            phospho=_phospho(),
            site_metadata=bad_site_metadata,
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


def test_analysis_ready_from_trusted_tables_matches_constructor_validation() -> None:
    bad_sample_metadata = pd.DataFrame(
        {"condition": ["control"]},
        index=pd.Index(["wrong_sample"], name="sample_id"),
    )

    with pytest.raises(DatasetValidationError) as factory_exc:
        AnalysisReadyPhosphoDataset.from_trusted_tables(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            sample_metadata=bad_sample_metadata,
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )
    with pytest.raises(DatasetValidationError) as constructor_exc:
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            sample_metadata=bad_sample_metadata,
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )

    assert type(factory_exc.value) is type(constructor_exc.value)
    assert str(factory_exc.value) == str(constructor_exc.value)


def test_analysis_ready_from_trusted_tables_records_trusted_construction_marker() -> (
    None
):
    dataset = AnalysisReadyPhosphoDataset.from_trusted_tables(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        trusted_construction_assertions=_trusted_assertions(),
        **_supported_dataset_state(has_total_matrix=False),
    )

    provenance = dataset.provenance

    assert provenance is not None
    assert provenance.workflow_name == "analysis_ready_dataset_direct_construction"
    construction = provenance.workflow_parameters["construction"]
    assert isinstance(construction, Mapping)
    assert construction["method"] == "AnalysisReadyPhosphoDataset.__init__"
    assert construction["source"] == "direct_trusted_construction"
    assert construction["builder_used"] is False
    assert construction["warning"] == (
        "Direct construction cannot prove biological correctness of "
        "caller-provided analysis-ready state."
    )
    assert construction["trusted_assertion_metadata_provided"] is True
    assert construction["missing_trusted_assertions"] == []


def test_direct_construction_records_provenance_marker() -> None:
    dataset = AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        **_supported_dataset_state(has_total_matrix=False),
    )

    provenance = dataset.provenance

    assert provenance is not None
    assert provenance.workflow_name == "analysis_ready_dataset_direct_construction"
    assert provenance.workflow_parameters["construction"]["source"] == (
        "direct_trusted_construction"
    )


def test_trusted_factory_records_same_direct_construction_marker() -> None:
    direct_dataset = AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        **_supported_dataset_state(has_total_matrix=False),
    )
    trusted_dataset = AnalysisReadyPhosphoDataset.from_trusted_tables(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        trusted_construction_assertions=_trusted_assertions(),
        **_supported_dataset_state(has_total_matrix=False),
    )

    assert direct_dataset.provenance is not None
    assert trusted_dataset.provenance is not None
    assert trusted_dataset.provenance.workflow_name == (
        direct_dataset.provenance.workflow_name
    )
    trusted_construction = trusted_dataset.provenance.workflow_parameters[
        "construction"
    ]
    direct_construction = direct_dataset.provenance.workflow_parameters["construction"]
    assert trusted_construction["method"] == direct_construction["method"]
    assert trusted_construction["source"] == direct_construction["source"]
    assert trusted_construction["builder_used"] == direct_construction["builder_used"]


def test_dataset_rejects_missing_site_sequence_column() -> None:
    bad_site_metadata = _site_metadata().drop(columns=["site_sequence"])
    with pytest.raises(
        DatasetValidationError,
        match="dataset.site_metadata is missing required columns: site_sequence",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=bad_site_metadata,
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


def test_analysis_ready_dataset_still_exposes_pandas_frames_after_schema_wrappers() -> (
    None
):
    dataset = AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        **_supported_dataset_state(has_total_matrix=False),
    )
    assert isinstance(dataset.phospho, pd.DataFrame)
    assert isinstance(dataset.site_metadata, pd.DataFrame)


def test_kinase_prediction_result_rejects_malformed_pred_mat_immediately() -> None:
    with pytest.raises(PhosPyValidationError, match="between 0.0 and 1.0"):
        KinasePredictionResult(
            pred_mat=pd.DataFrame(
                {"MAP2K6": [1.5]},
                index=_SITE_INDEX.copy(),
            )
        )


def test_dataset_exposes_optional_preprocessing_report_field() -> None:
    dataset = AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        **_supported_dataset_state(has_total_matrix=False),
    )

    assert hasattr(dataset, "preprocessing_report")
    assert dataset.preprocessing_report is None


def test_dataset_rejects_blank_gene_symbol_values() -> None:
    bad_site_metadata = _site_metadata().copy(deep=True)
    bad_site_metadata.loc[:, "gene_symbol"] = ["   "]
    with pytest.raises(
        DatasetValidationError,
        match="dataset.site_metadata.gene_symbol must contain non-empty string values",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=bad_site_metadata,
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


def test_dataset_rejects_blank_site_values() -> None:
    bad_site_metadata = _site_metadata().copy(deep=True)
    bad_site_metadata.loc[:, "site"] = [""]
    with pytest.raises(
        DatasetValidationError,
        match="dataset.site_metadata.site must contain non-empty string values",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=bad_site_metadata,
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


def test_dataset_rejects_blank_site_sequence_values() -> None:
    bad_site_metadata = _site_metadata().copy(deep=True)
    bad_site_metadata.loc[:, "site_sequence"] = [" \t "]
    with pytest.raises(
        DatasetValidationError,
        match=(
            "dataset.site_metadata.site_sequence must contain non-empty string values"
        ),
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=bad_site_metadata,
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


def test_dataset_rejects_null_site_sequence_values() -> None:
    bad_site_metadata = _site_metadata().copy(deep=True)
    bad_site_metadata.loc[:, "site_sequence"] = [None]
    with pytest.raises(
        DatasetValidationError,
        match="dataset.site_metadata.site_sequence must not contain missing values",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=bad_site_metadata,
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


def test_dataset_rejects_empty_phospho_matrix() -> None:
    with pytest.raises(DatasetValidationError):
        AnalysisReadyPhosphoDataset(
            phospho=pd.DataFrame(),
            site_metadata=pd.DataFrame(),
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


def test_dataset_rejects_nan_in_phospho_matrix() -> None:
    phospho = _phospho()
    phospho.loc[_SITE_KEY, "sample_a"] = float("nan")
    with pytest.raises(DatasetValidationError, match="must not contain missing values"):
        AnalysisReadyPhosphoDataset(
            phospho=phospho,
            site_metadata=_site_metadata(),
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


@pytest.mark.parametrize("invalid_value", [float("inf"), float("-inf")])
def test_dataset_rejects_infinite_values_in_phospho_matrix(
    invalid_value: float,
) -> None:
    phospho = _phospho()
    phospho.loc[_SITE_KEY, "sample_a"] = invalid_value
    with pytest.raises(
        DatasetValidationError, match="must contain finite numeric values"
    ):
        AnalysisReadyPhosphoDataset(
            phospho=phospho,
            site_metadata=_site_metadata(),
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


@pytest.mark.parametrize("invalid_value", [float("inf"), float("-inf")])
def test_dataset_rejects_infinite_values_in_comparisons_matrix(
    invalid_value: float,
) -> None:
    with pytest.raises(DatasetValidationError) as exc_info:
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            comparisons=pd.DataFrame(
                {"p_group1_group4": [invalid_value]},
                index=_SITE_INDEX.copy(),
            ),
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )
    message = str(exc_info.value)
    assert "dataset.comparisons must contain finite numeric values" in message
    assert f"({_SITE_KEY!r}, 'p_group1_group4')" in message


def test_dataset_rejects_nan_in_comparisons_matrix() -> None:
    with pytest.raises(
        DatasetValidationError,
        match="dataset.comparisons must not contain missing values",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            comparisons=pd.DataFrame(
                {"p_group1_group4": [float("nan")]},
                index=_SITE_INDEX.copy(),
            ),
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


def test_dataset_rejects_boolean_columns_in_total_matrix() -> None:
    with pytest.raises(DatasetValidationError) as exc_info:
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            total=pd.DataFrame({"sample_a": [True]}, index=["MAPK14"]),
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=True),
        )
    message = str(exc_info.value)
    assert "dataset.total must contain only scientific numeric columns" in message
    assert "boolean columns are invalid: sample_a" in message


def test_dataset_rejects_boolean_columns_in_comparisons_matrix() -> None:
    with pytest.raises(DatasetValidationError) as exc_info:
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            comparisons=pd.DataFrame(
                {"p_group1_group4": [True]},
                index=_SITE_INDEX.copy(),
            ),
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )
    message = str(exc_info.value)
    assert "dataset.comparisons must contain only scientific numeric columns" in message
    assert "boolean columns are invalid: p_group1_group4" in message


def test_builder_rejects_sparse_missingness_in_phospho_matrix() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan")],
            "sample_b": [2.0, 3.0],
        },
        index=["MAPK14;Y182;", "AKT1;T308;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "RPHFPQFSYSASGTA",
            ],
            "protein_id": ["MAPK14", "AKT1"],
            "localisation_confidence": [0.95, 0.9],
        },
        index=phospho.index,
    )
    with pytest.raises(
        PhosPyInputError,
        match="dataset preprocessing stage 'missing_data' rejected missing values",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                input_intensity_scale="linear",
            )
        )


def test_dataset_validates_intensity_scale_state_for_total_matrix() -> None:
    with pytest.raises(TransformationValidationError):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            total=pd.DataFrame({"sample_a": [2.0]}, index=["GENEA"]),
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


def test_dataset_rejects_mixed_intensity_scale_kind_between_phospho_and_total() -> None:
    supported_state = supported_linear_intensity_scale_state(has_total_matrix=True)
    object.__setattr__(supported_state, "total", MatrixIntensityScaleState.log2())
    with pytest.raises(TransformationValidationError):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            total=pd.DataFrame({"sample_a": [2.0]}, index=["GENEA"]),
            organism=Organism.RAT,
            intensity_scale_state=supported_state,
            processing_state=supported_linear_processing_state(has_total_matrix=True),
        )


def test_dataset_rejects_inf_in_total_matrix() -> None:
    with pytest.raises(
        DatasetValidationError, match="dataset.total must contain finite numeric values"
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            total=pd.DataFrame({"sample_a": [float("-inf")]}, index=["MAPK14"]),
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=True),
        )


def test_dataset_rejects_total_matrix_sample_mismatch() -> None:
    with pytest.raises(
        DatasetValidationError,
        match="dataset.total.columns must exactly match dataset.phospho.columns",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            total=pd.DataFrame({"sample_b": [2.0]}, index=["MAPK14"]),
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=True),
        )


def test_dataset_accepts_aligned_numeric_comparisons() -> None:
    dataset = AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        comparisons=pd.DataFrame(
            {"p_group1_group4": [3.0]},
            index=_SITE_INDEX.copy(),
        ),
        organism=Organism.RAT,
        **_supported_dataset_state(has_total_matrix=False),
    )
    assert dataset.comparisons is not None
    assert dataset.comparisons.loc[_SITE_KEY, "p_group1_group4"] == 3.0


def test_dataset_rejects_comparisons_index_mismatch() -> None:
    with pytest.raises(
        DatasetValidationError,
        match="dataset.comparisons.index must exactly match dataset.phospho.index",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            comparisons=pd.DataFrame(
                {"p_group1_group4": [3.0]},
                index=pd.Index([_AKT1_T308_KEY], name="site_key"),
            ),
            organism=Organism.RAT,
            **_supported_dataset_state(has_total_matrix=False),
        )


def test_reference_bundle_requires_non_empty_resources() -> None:
    with pytest.raises(ReferenceValidationError):
        ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(),
            site_sequences=pd.DataFrame(
                {"site_sequence": ["AAAAAAA"]},
                index=pd.Index(["GENEA;S1;"]),
            ),
        )


def test_reference_resolver_auto_uses_dataset_organism() -> None:
    resolved = ReferenceResolver().run(
        ReferencePreset.AUTO,
        dataset_organism=Organism.RAT,
    )
    assert resolved.organism is Organism.RAT
    assert not resolved.kinase_substrate_map.empty
    assert not resolved.site_sequences.empty


def test_reference_resolver_rejects_preset_dataset_mismatch() -> None:
    with pytest.raises(ReferenceCompatibilityError):
        ReferenceResolver().run(
            ReferencePreset.HUMAN,
            dataset_organism=Organism.RAT,
        )


def test_reference_resolver_rejects_unsupported_human_bundled_preset() -> None:
    with pytest.raises(UnsupportedOrganismError):
        ReferenceResolver().run(
            ReferencePreset.HUMAN,
            dataset_organism=Organism.HUMAN,
        )


def test_reference_resolver_auto_requires_dataset_organism() -> None:
    with pytest.raises(ReferenceResolutionError):
        ReferenceResolver().run(
            ReferencePreset.AUTO,
            dataset_organism=None,
        )


def test_builder_rejects_non_dataframe_phospho_input() -> None:
    with pytest.raises(UnsupportedInputFormatError):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho="not-a-dataframe",
                site_metadata=_site_metadata(),
                organism=Organism.RAT,
                input_intensity_scale="linear",
            )
        )


def test_builder_supports_file_path_inputs(tmp_path) -> None:
    phospho_path = tmp_path / "phospho.csv"
    site_metadata_path = tmp_path / "site_metadata.csv"
    _phospho().to_csv(phospho_path)
    _site_metadata().to_csv(site_metadata_path)

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho_path,
            site_metadata=site_metadata_path,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )
    assert built.phospho.index.name == "site_key"
    site_key = _site_key_for_display_id(built.site_metadata, _DISPLAY_ID)
    assert list(built.phospho.index.astype(str)) == [site_key]
    assert built.site_metadata.loc[site_key, "site_sequence"]


def test_builder_accepts_string_dtype_identity_columns() -> None:
    site_metadata = _site_metadata().astype(
        {
            "site_key": "string",
            "display_id": "string",
            "gene_symbol": "string",
            "site": "string",
            "site_sequence": "string",
            "protein_id": "string",
        }
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    assert built.phospho.index.name == "site_key"
    assert built.site_metadata.loc[:, "display_id"].astype(str).tolist() == [
        _DISPLAY_ID
    ]


def test_builder_rejects_duplicate_sample_metadata_columns_before_returning_dataset() -> (
    None
):
    sample_metadata = pd.DataFrame(
        [["condition_a", "duplicate_condition"]],
        columns=["condition", "condition"],
        index=pd.Index(["sample_a"], name="sample_id"),
    )

    with pytest.raises(
        DatasetValidationError,
        match="dataset.sample_metadata.columns must be unique",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=_phospho(),
                site_metadata=_site_metadata(),
                sample_metadata=sample_metadata,
                organism=Organism.RAT,
                input_intensity_scale="linear",
            )
        )


def test_builder_rejects_preprocessing_threshold_above_sample_count() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="cannot exceed the number of phospho samples",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=_phospho(),
                site_metadata=_site_metadata(),
                organism=Organism.RAT,
                input_intensity_scale="linear",
                preprocessing_config=DatasetPreprocessingConfig(
                    missing_data=DatasetMissingDataConfig(
                        policy="impute_row_median",
                        min_observed_values=3,
                    ),
                ),
            )
        )


def test_builder_derives_site_sequence_for_supported_organism() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata().drop(columns=["site_sequence"]),
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )
    site_key = _site_key_for_display_id(built.site_metadata, _DISPLAY_ID)
    assert built.site_metadata.loc[site_key, "site_sequence"]


def test_builder_derives_gene_symbol_and_site_from_supported_index_convention() -> None:
    display_phospho = pd.DataFrame({"sample_a": [1.0]}, index=[_DISPLAY_ID])
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=display_phospho,
            site_metadata=pd.DataFrame(
                {
                    "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                    "protein_id": ["MAPK14"],
                    "localisation_confidence": [0.95],
                },
                index=[_DISPLAY_ID],
            ),
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )
    site_key = _site_key_for_display_id(built.site_metadata, _DISPLAY_ID)
    assert built.site_metadata.loc[site_key, "gene_symbol"] == "MAPK14"
    assert built.site_metadata.loc[site_key, "site"] == "Y182"


def test_builder_rejects_ambiguous_sequence_column_without_explicit_convention() -> (
    None
):
    with pytest.raises(
        UnsupportedInputFormatError,
        match="column 'sequence' is unsupported",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=_phospho(),
                site_metadata=pd.DataFrame(
                    {
                        "gene_symbol": ["MAPK14"],
                        "site": ["Y182"],
                        "sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                    },
                    index=["MAPK14;Y182;"],
                ),
                organism=Organism.RAT,
                input_intensity_scale="linear",
            )
        )


def test_builder_rejects_unsupported_historical_gene_alias() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match="column 'gene' is unsupported",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=_phospho(),
                site_metadata=pd.DataFrame(
                    {
                        "gene": ["MAPK14"],
                        "site": ["Y182"],
                        "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                    },
                    index=["MAPK14;Y182;"],
                ),
                organism=Organism.RAT,
                input_intensity_scale="linear",
            )
        )


@pytest.mark.parametrize("use_file_path", [False, True], ids=["dataframe", "file_path"])
@pytest.mark.parametrize(
    ("column_name", "blank_value", "error_type", "message"),
    [
        (
            "gene_symbol",
            "   ",
            DatasetValidationError,
            "dataset.site_metadata.gene_symbol must contain non-empty string values",
        ),
        (
            "site",
            "   ",
            DatasetValidationError,
            "dataset.site_metadata.site must contain non-empty string values",
        ),
        (
            "site_sequence",
            "   ",
            UnsupportedInputFormatError,
            "dataset build request site_metadata.site_sequence must contain non-empty string values",
        ),
    ],
)
def test_builder_rejects_blank_required_site_metadata_fields_across_input_routes(
    use_file_path: bool,
    column_name: str,
    blank_value: str,
    error_type,
    message: str,
    tmp_path,
) -> None:
    site_metadata = _site_metadata().copy(deep=True)
    site_metadata.loc[:, column_name] = [blank_value]
    with pytest.raises(
        error_type,
        match=message,
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=_as_builder_input(
                    _phospho(),
                    use_file_path=use_file_path,
                    tmp_path=tmp_path,
                    filename="phospho.csv",
                ),
                site_metadata=_as_builder_input(
                    site_metadata,
                    use_file_path=use_file_path,
                    tmp_path=tmp_path,
                    filename="site_metadata.csv",
                ),
                organism=Organism.RAT,
                input_intensity_scale="linear",
            )
        )


@pytest.mark.parametrize("use_file_path", [False, True], ids=["dataframe", "file_path"])
def test_builder_rejects_unsupported_historical_site_alias_across_input_routes(
    use_file_path: bool,
    tmp_path,
) -> None:
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "residue": ["Y182"],
            "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
        },
        index=["MAPK14;Y182;"],
    )
    with pytest.raises(
        UnsupportedInputFormatError,
        match="column 'residue' is unsupported",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=_as_builder_input(
                    _phospho(),
                    use_file_path=use_file_path,
                    tmp_path=tmp_path,
                    filename="phospho.csv",
                ),
                site_metadata=_as_builder_input(
                    site_metadata,
                    use_file_path=use_file_path,
                    tmp_path=tmp_path,
                    filename="site_metadata.csv",
                ),
                organism=Organism.RAT,
                input_intensity_scale="linear",
            )
        )


def test_builder_fails_when_missing_gene_or_site_cannot_be_derived_from_index() -> None:
    with pytest.raises(
        UnsupportedInputFormatError,
        match="site_metadata is missing required metadata columns",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=_phospho(),
                site_metadata=pd.DataFrame(
                    {
                        "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
                        "localisation_confidence": [0.95],
                    },
                    index=["MAPK14"],
                ),
                organism=Organism.RAT,
                input_intensity_scale="linear",
            )
        )


def test_builder_rejects_missing_site_sequence_when_not_provided_or_derivable() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="cannot construct AnalysisReadyPhosphoDataset",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=_phospho(),
                site_metadata=_site_metadata().drop(columns=["site_sequence"]),
                organism=None,
            )
        )


def test_builder_site_matrix_excludes_only_unresolved_rows_in_mixed_support_case() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [1.5, 2.5, 3.5],
        },
        index=pd.Index(["MAPK14;Y182;", "FAKE1;S1;", "GSK3B;S9;"], name="input_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "FAKE1", "GSK3B"],
            "site": ["Y182", "S1", "S9"],
            "protein_id": ["MAPK14", "FAKE1", "GSK3B"],
            "localisation_confidence": [0.95, 0.9, 0.92],
        },
        index=phospho.index.copy(),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
            ),
        )
    )

    assert built.phospho.index.name == "site_key"
    assert built.site_metadata.loc[:, "display_id"].astype(str).tolist() == [
        "GSK3B;S9;",
        "MAPK14;Y182;",
    ]
    assert built.preprocessing_report is not None
    assert built.preprocessing_report.row_audit is not None
    dropped = built.preprocessing_report.row_audit.loc[
        (built.preprocessing_report.row_audit.loc[:, "stage"] == "site_matrix")
        & (built.preprocessing_report.row_audit.loc[:, "action"] == "dropped")
    ]
    assert set(dropped.loc[:, "source_row_id"].astype(str)) == {"FAKE1;S1;"}


def test_builder_site_matrix_reports_no_retained_rows_when_all_rows_lack_sequence_support() -> (
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
            "protein_id": ["FAKE1", "FAKE2"],
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
                input_intensity_scale="linear",
                preprocessing_config=DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
                ),
            )
        )


def test_io_namespace_does_not_expose_parallel_dataset_file_builder() -> None:
    assert not hasattr(phospy_io, "build_dataset_from_files")


def test_reference_resolver_accepts_explicit_bundle() -> None:
    bundle = _valid_bundle(organism=Organism.MOUSE)
    resolved = ReferenceResolver().run(
        bundle,
        dataset_organism=Organism.MOUSE,
    )
    assert resolved is bundle


def test_builder_establishes_intensity_scale_state_with_supported_path() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )
    assert built.intensity_scale_state.label == "linear"
    assert built.intensity_scale_state.is_established
    assert built.intensity_scale_state.established_via is not None
    assert (
        built.intensity_scale_state.phospho.established_by
        == "phospy.science.datasets.builders.executor.input_intensity_scale"
    )


def test_workflows_do_not_call_borrow_methods_directly() -> None:
    violations: list[str] = []
    for relative_dir in _WORKFLOW_BORROW_GUARD_DIRS:
        for path in (ROOT / relative_dir).rglob("*.py"):
            source = path.read_text(encoding="utf-8")
            lines = source.splitlines()
            for match in _WORKFLOW_BORROW_METHOD_PATTERN.finditer(source):
                line_number = source.count("\n", 0, match.start()) + 1
                line = lines[line_number - 1].strip()
                relative_path = path.relative_to(ROOT).as_posix()
                violations.append(f"{relative_path}:{line_number}: {line}")

    assert not violations, (
        "workflow modules must use narrow domain-owned internal views, "
        "not direct dataset _borrow_* methods:\n" + "\n".join(violations)
    )


def test_dataset_internal_view_exposes_only_required_frames() -> None:
    public_members = {
        name for name in dir(DatasetInternalView) if not name.startswith("_")
    }
    view_properties = {
        name
        for name, value in inspect.getmembers(DatasetInternalView)
        if isinstance(value, property)
    }

    assert public_members == _DATASET_INTERNAL_VIEW_PUBLIC_MEMBERS
    assert view_properties == _DATASET_INTERNAL_VIEW_FRAME_PROPERTIES
    assert hasattr(DatasetInternalView, "imputation_observation_summary")
    assert not hasattr(DatasetInternalView, "imputation_observed_mask")
    assert not hasattr(
        DatasetInternalView(
            AnalysisReadyPhosphoDataset(
                phospho=_phospho(),
                site_metadata=_site_metadata(),
                organism=Organism.RAT,
                **_supported_dataset_state(has_total_matrix=False),
            )
        ),
        "dataset",
    )


def test_dataset_internal_view_returns_defensive_frame_snapshots() -> None:
    dataset = AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        sample_metadata=pd.DataFrame(
            {"condition": ["treated"]},
            index=pd.Index(["sample_a"], name="sample_id"),
        ),
        total=pd.DataFrame({"sample_a": [2.0]}, index=pd.Index(["MAPK14"])),
        comparisons=pd.DataFrame({"p_group1_group4": [0.05]}, index=_SITE_INDEX.copy()),
        imputation_observation_mask=pd.DataFrame(
            {"sample_a": [True]},
            index=_SITE_INDEX.copy(),
        ),
        organism=Organism.RAT,
        **_supported_dataset_state(has_total_matrix=True),
    )
    view = DatasetInternalView(dataset)

    phospho_snapshot = view.phospho
    site_metadata_snapshot = view.site_metadata
    sample_metadata_snapshot = view.sample_metadata
    total_snapshot = view.total
    comparisons_snapshot = view.comparisons
    imputation_summary = view.imputation_observation_summary(
        feature_ids=_SITE_INDEX,
        sample_ids=["sample_a"],
    )

    def _mutate_snapshot(mutator) -> None:
        try:
            mutator()
        except (TypeError, ValueError):
            pass

    _mutate_snapshot(lambda: phospho_snapshot.iloc.__setitem__((0, 0), 999.0))
    _mutate_snapshot(
        lambda: site_metadata_snapshot.loc.__setitem__(
            (_SITE_KEY, "gene_symbol"),
            "MUTATED",
        )
    )
    assert sample_metadata_snapshot is not None
    _mutate_snapshot(
        lambda: sample_metadata_snapshot.loc.__setitem__(
            ("sample_a", "condition"),
            "mutated",
        )
    )
    assert total_snapshot is not None
    _mutate_snapshot(lambda: total_snapshot.iloc.__setitem__((0, 0), 999.0))
    assert comparisons_snapshot is not None
    _mutate_snapshot(lambda: comparisons_snapshot.iloc.__setitem__((0, 0), 0.99))
    assert imputation_summary is not None
    imputation_summary.loc[_SITE_KEY, "observed_cell_count"] = 0

    assert phospho_snapshot is not dataset._phospho
    assert site_metadata_snapshot is not dataset._site_metadata
    assert float(dataset.phospho.iloc[0, 0]) == 1.0
    assert str(dataset.site_metadata.loc[_SITE_KEY, "gene_symbol"]) == "MAPK14"
    assert dataset.sample_metadata is not None
    assert str(dataset.sample_metadata.loc["sample_a", "condition"]) == "treated"
    assert dataset.total is not None
    assert float(dataset.total.iloc[0, 0]) == 2.0
    assert dataset.comparisons is not None
    assert float(dataset.comparisons.iloc[0, 0]) == 0.05
    reread_summary = dataset.imputation_observation_summary_dataframe(
        feature_ids=_SITE_INDEX,
        sample_ids=["sample_a"],
    )
    assert reread_summary is not None
    assert int(reread_summary.loc[_SITE_KEY, "observed_cell_count"]) == 1


def test_kinase_reference_identifier_cleanup_is_owned_by_reference_ingestion_boundary() -> (
    None
):
    violations: list[str] = []
    for relative_path in _REFERENCE_IDENTIFIER_GUARD_TARGETS:
        path = ROOT / relative_path
        source = path.read_text(encoding="utf-8")
        lines = source.splitlines()
        for match in _FORBIDDEN_IDENTIFIER_CLEANUP.finditer(source):
            line_number = source.count("\n", 0, match.start()) + 1
            window_start = max(1, line_number - 2)
            window_end = min(len(lines), line_number + 2)
            window_text = "\n".join(lines[window_start - 1 : window_end])
            if _REFERENCE_IDENTIFIER_COLUMN_HINT.search(window_text) is None:
                continue
            line = lines[line_number - 1].strip()
            violations.append(f"{relative_path}:{line_number}: {line}")

    assert not violations, (
        "identifier normalisation must stay inside reference ingestion; "
        "forbidden workflow-local cleanup detected:\n" + "\n".join(violations)
    )


def test_production_reference_table_owns_protein_accession_normalisation() -> None:
    source = (ROOT / _PRODUCTION_REFERENCE_TABLES).read_text(encoding="utf-8")
    assert "class ProteinAccessionReference" in source, (
        "protein accession normalisation must be owned by an explicit production "
        "reference ingestion boundary in src/phospy/tables/references.py"
    )
    assert "normalise_reference_protein_accession" in source, (
        "production protein accession ingestion must call the reference identifier "
        f"boundary helper in {_REFERENCE_IDENTIFIER_BOUNDARY_OWNER}"
    )


def test_non_reference_domains_do_not_call_reference_protein_accession_normaliser() -> (
    None
):
    violations: list[str] = []
    for relative_dir in _FORBIDDEN_PROTEIN_ACCESSION_BOUNDARY_LEAK_TARGETS:
        directory = ROOT / relative_dir
        for path in directory.rglob("*.py"):
            source = path.read_text(encoding="utf-8")
            if "normalise_reference_protein_accession(" not in source:
                continue
            relative_path = path.relative_to(ROOT).as_posix()
            violations.append(relative_path)

    assert not violations, (
        "protein accession normalisation belongs to reference identifier boundary; "
        "forbidden calls found outside reference ingestion:\n" + "\n".join(violations)
    )


def test_sequence_domain_does_not_apply_accession_case_cleanup() -> None:
    violations: list[str] = []
    for path in (ROOT / "src/phospy/science/sequences").rglob("*.py"):
        source = path.read_text(encoding="utf-8")
        for match in _FORBIDDEN_ACCESSION_CASE_NORMALISATION.finditer(source):
            line_number = source.count("\n", 0, match.start()) + 1
            line = source.splitlines()[line_number - 1].strip()
            relative_path = path.relative_to(ROOT).as_posix()
            violations.append(f"{relative_path}:{line_number}: {line}")

    assert not violations, (
        "sequence-domain code must not apply reference-style accession casing "
        "normalisation:\n" + "\n".join(violations)
    )
