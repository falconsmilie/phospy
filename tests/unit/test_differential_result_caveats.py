from __future__ import annotations

from collections.abc import Mapping

import pandas as pd
import pytest

from phospy.advanced import DifferentialAnalysisConfig
from phospy.api import (
    Contrast,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.provenance.immutability import thaw_json_mapping
from phospy.science.datasets.builders.preprocessing import (
    build_dataset_processing_state,
)
from phospy.science.datasets.builders.transformation_resolver import (
    DatasetIntensityScaleResolver,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT,
)
from phospy.science.transformations.models import (
    IntensityScaleEstablishmentMode,
    IntensityScaleEvidenceLevel,
    IntensityScaleKind,
    IntensityScaleState,
    MatrixIntensityScaleState,
)
from phospy.science.transformations.transformers import (
    IdentityTransformer,
    Log2Transformer,
)
from phospy.workflows.differential.caveats import (
    DIFFERENTIAL_DECLARED_SCALE_OVERRIDE_CAVEAT_CODE,
    DIFFERENTIAL_DIRECT_TRUSTED_DATASET_CAVEAT_CODE,
    DIFFERENTIAL_IMPUTATION_WITHHOLDING_POLICY_CAVEAT_CODE,
    DIFFERENTIAL_WITHHELD_FEATURES_CAVEAT_CODE,
)
from phospy.workflows.intensity_scale_evidence import (
    INPUT_INTENSITY_SCALE_DECLARED_CAVEAT_CODE,
)
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import supported_log2_intensity_scale_state
from tests.support.processing_state import (
    imputed_processing_state as valid_imputed_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

_GENES = ["MAPK14", "AKT1", "GSK3B", "RPS6"]
_SITES = ["Y182", "T308", "S9", "S235"]
_SAMPLES = ("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")


def _site_index() -> pd.Index:
    return protein_site_key_index(protein_identifiers=_GENES, sites=_SITES)


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_key": index.astype(str).tolist(),
            "display_id": [
                f"{gene};{site};" for gene, site in zip(_GENES, _SITES, strict=True)
            ],
            **site_key_context_columns(index),
            "gene_symbol": _GENES,
            "site": _SITES,
            "site_sequence": [("A" * 15) + site[0] + ("A" * 15) for site in _SITES],
            "protein_id": _GENES,
        },
        index=index.copy(),
    )


def _phospho(index: pd.Index, *, constant_first: bool = False) -> pd.DataFrame:
    first_row = [10.0] * 6 if constant_first else [10.0, 10.1, 9.9, 15.0, 15.1, 14.9]
    return pd.DataFrame(
        {
            "A_1": [first_row[0], 20.0, 5.0, 30.0],
            "A_2": [first_row[1], 20.2, 5.1, 30.2],
            "A_3": [first_row[2], 20.1, 5.2, 30.1],
            "B_1": [first_row[3], 24.0, 6.0, 29.0],
            "B_2": [first_row[4], 24.1, 6.2, 29.2],
            "B_3": [first_row[5], 24.2, 6.1, 29.1],
        },
        index=index.copy(),
    )


def _observed_mask(index: pd.Index) -> pd.DataFrame:
    mask = pd.DataFrame(True, index=index.copy(), columns=pd.Index(_SAMPLES))
    mask.loc[index[1], "B_3"] = False
    mask.loc[index[2], ["A_3", "B_3"]] = False
    mask.loc[index[3], ["B_2", "B_3"]] = False
    return mask


def _dataset(
    phospho: pd.DataFrame,
    *,
    intensity_scale_state: IntensityScaleState | None = None,
    imputed: bool = False,
    imputation_observation_mask: pd.DataFrame | None = None,
) -> AnalysisReadyPhosphoDataset:
    if intensity_scale_state is None:
        intensity_scale_state = supported_log2_intensity_scale_state(
            has_total_matrix=False
        )
    processing_state = build_dataset_processing_state(
        plan=PreprocessingPlan.default(),
        intensity_scale_state=intensity_scale_state,
    )
    intensity_scale_state = processing_state.intensity_scale
    if imputed:
        processing_state = valid_imputed_processing_state(processing_state)
        intensity_scale_state = processing_state.intensity_scale
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        imputation_observation_mask=imputation_observation_mask,
        organism=Organism.RAT,
        intensity_scale_state=intensity_scale_state,
        processing_state=processing_state,
    )


def _design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=tuple(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=sample_id.split("_", maxsplit=1)[0],
                biological_replicate_id=sample_id,
            )
            for sample_id in _SAMPLES
        )
    )


def _request(
    dataset: AnalysisReadyPhosphoDataset,
    *,
    config: DifferentialAnalysisConfig | None = None,
) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=dataset,
        design=_design(),
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
            ),
        ),
        config=DifferentialAnalysisConfig() if config is None else config,
    )


def _declared_log2_intensity_scale_state(
    phospho: pd.DataFrame,
) -> IntensityScaleState:
    declared_state = IntensityScaleState(
        phospho=MatrixIntensityScaleState.log2(established_by="test.declaration"),
        total=None,
    )
    return (
        DatasetIntensityScaleResolver(transformer=IdentityTransformer())
        .run(
            phospho=phospho,
            total=None,
            declared_input_scale_state=declared_state,
            declared_input_establishment_mode=IntensityScaleEstablishmentMode.DECLARED,
            input_declaration_source="tests.unit.test_differential_result_caveats",
        )
        .intensity_scale_state
    )


def _observed_log2_dataset(index: pd.Index) -> AnalysisReadyPhosphoDataset:
    target_log2_values = _phospho(index)
    linear_values = 2.0**target_log2_values
    resolved = DatasetIntensityScaleResolver(
        transformer=Log2Transformer(pseudocount=0.0)
    ).run(
        phospho=linear_values,
        total=None,
        expected_scale_kind=IntensityScaleKind.LOG2,
        scale_establishment_evidence_level=(
            IntensityScaleEvidenceLevel.OBSERVED_TRANSFORMATION
        ),
        scale_establishment_parameters={"operation": "log2", "pseudocount": 0.0},
    )
    return _dataset(
        resolved.phospho,
        intensity_scale_state=resolved.intensity_scale_state,
    )


def _caveat_by_code(result, code: str):
    matches = [caveat for caveat in result.caveats if caveat.code == code]
    assert len(matches) == 1
    return matches[0]


def test_differential_result_caveats_include_trusted_table_reconstruction() -> None:
    index = _site_index()
    result = DifferentialAnalysisWorkflow().run(_request(_dataset(_phospho(index))))

    caveat = _caveat_by_code(result, DIFFERENTIAL_DIRECT_TRUSTED_DATASET_CAVEAT_CODE)

    assert caveat.severity == "info"
    assert caveat.details["construction_source"] == "direct_trusted_construction"
    assert caveat.details["trusted_assertion_metadata_provided"] is True
    assert caveat.details["missing_trusted_assertions"] == ()
    assert "input_tables" not in caveat.details


def test_differential_result_caveats_include_withheld_features() -> None:
    index = _site_index()
    result = DifferentialAnalysisWorkflow().run(
        _request(_dataset(_phospho(index, constant_first=True)))
    )

    caveat = _caveat_by_code(result, DIFFERENTIAL_WITHHELD_FEATURES_CAVEAT_CODE)

    assert caveat.severity == "warning"
    assert caveat.details["withheld_feature_count"] == 1
    status_counts = caveat.details["status_counts"]
    assert isinstance(status_counts, Mapping)
    assert status_counts[DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT] == 1
    assert caveat.details["result_status_column"] == "result_status"
    assert caveat.details["result_status_reason_column"] == "result_status_reason"


def test_differential_result_caveats_include_imputation_policy_when_used() -> None:
    index = _site_index()
    dataset = _dataset(
        _phospho(index),
        imputed=True,
        imputation_observation_mask=_observed_mask(index),
    )

    result = DifferentialAnalysisWorkflow().run(
        _request(
            dataset,
            config=DifferentialAnalysisConfig(
                imputed_value_policy="withhold_imputed_features",
                imputed_value_max_fraction=0.20,
            ),
        )
    )

    caveat = _caveat_by_code(
        result,
        DIFFERENTIAL_IMPUTATION_WITHHOLDING_POLICY_CAVEAT_CODE,
    )

    assert caveat.severity == "warning"
    assert caveat.details["policy"] == "withhold_imputed_features"
    assert caveat.details["imputed_value_max_fraction"] == pytest.approx(0.20)
    assert caveat.details["withheld_feature_count"] == 2
    assert caveat.details["testable_feature_count"] == 2


def test_differential_result_caveats_include_declared_scale_override_when_used() -> (
    None
):
    index = _site_index()
    suspicious_matrix = _phospho(index) * 10000.0
    state = _declared_log2_intensity_scale_state(suspicious_matrix)

    result = DifferentialAnalysisWorkflow().run(
        _request(
            _dataset(suspicious_matrix, intensity_scale_state=state),
            config=DifferentialAnalysisConfig(
                allow_suspicious_declared_input_scale=True
            ),
        )
    )

    caveat = _caveat_by_code(result, DIFFERENTIAL_DECLARED_SCALE_OVERRIDE_CAVEAT_CODE)

    assert caveat.severity == "warning"
    assert caveat.details["establishment_mode"] == "declared"
    assert caveat.details["diagnostic_warning_count"] >= 1
    assert caveat.details["override_config"] == "allow_suspicious_declared_input_scale"


def test_differential_result_provenance_includes_observed_transformation_evidence() -> (
    None
):
    index = _site_index()
    result = DifferentialAnalysisWorkflow().run(_request(_observed_log2_dataset(index)))

    assert result.workflow_provenance is not None
    assert result.workflow_provenance["input_intensity_scale"] == "log2"
    assert (
        result.workflow_provenance["input_intensity_scale_evidence_level"]
        == "observed_transformation"
    )
    assert (
        result.workflow_provenance["input_intensity_scale_source"]
        == "transformed_by_phospy"
    )
    assert result.policy_provenance is not None
    testing_policy = result.policy_provenance.statistical_testing
    assert testing_policy.input_intensity_scale_evidence_level == (
        "observed_transformation"
    )
    assert testing_policy.input_intensity_scale_source == "transformed_by_phospy"
    payload = result.to_payload()
    assert payload["workflow_provenance"] == thaw_json_mapping(
        result.workflow_provenance,
        field_name="differential_result.workflow_provenance",
    )
    assert all(
        caveat.code != INPUT_INTENSITY_SCALE_DECLARED_CAVEAT_CODE
        for caveat in result.caveats
    )


def test_differential_result_provenance_includes_declared_scale_evidence_and_caveat() -> (
    None
):
    index = _site_index()
    result = DifferentialAnalysisWorkflow().run(_request(_dataset(_phospho(index))))

    assert result.workflow_provenance is not None
    assert result.workflow_provenance["input_intensity_scale"] == "log2"
    assert (
        result.workflow_provenance["input_intensity_scale_evidence_level"]
        == "declared_by_user"
    )
    assert (
        result.workflow_provenance["input_intensity_scale_source"] == "declared_by_user"
    )
    assert result.policy_provenance is not None
    testing_policy = result.policy_provenance.statistical_testing
    assert testing_policy.input_intensity_scale_evidence_level == "declared_by_user"
    assert testing_policy.input_intensity_scale_source == "declared_by_user"

    caveat = _caveat_by_code(result, INPUT_INTENSITY_SCALE_DECLARED_CAVEAT_CODE)
    assert caveat.severity == "warning"
    assert caveat.details["input_intensity_scale"] == "log2"
    assert caveat.details["input_intensity_scale_evidence_level"] == (
        "declared_by_user"
    )
    assert caveat.details["workflow_scope"] == "differential"
