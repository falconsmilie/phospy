from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetTotalProteinCorrectionConfig,
)
from phospy.errors import DatasetValidationError, PhosPyInputError
from phospy.errors.transformations import TransformationStateEstablishmentError
from phospy.science.transformations.models import IntensityScaleEstablishmentSource

pytestmark = pytest.mark.integration


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    gene_symbols: list[str] = []
    sites: list[str] = []
    for site_id in index.astype(str):
        parts = site_id.split(";")
        gene_symbols.append(parts[0])
        sites.append(parts[1] if len(parts) > 1 else "S1")
    rows = len(gene_symbols)
    return pd.DataFrame(
        {
            "gene_symbol": gene_symbols,
            "site": sites,
            "protein_id": gene_symbols,
            "organism": ["rat"] * rows,
            "site_sequence": ["SEQ_A"] * rows,
            "localisation_confidence": [0.95] * rows,
        },
        index=index.copy(),
    )


def _workflow_establishment_payload(dataset) -> dict[str, object]:
    assert dataset.provenance is not None
    payload = dataset.provenance.workflow_parameters.get(
        "intensity_scale_establishment"
    )
    assert isinstance(payload, dict)
    return payload


def _workflow_meaning_payload(dataset) -> dict[str, object]:
    assert dataset.provenance is not None
    payload = dataset.provenance.workflow_parameters.get(
        "quantitative_meaning_provenance"
    )
    assert isinstance(payload, dict)
    return payload


def _final_stage_establishment_payload(dataset) -> dict[str, object]:
    assert dataset.preprocessing_report is not None
    final_stage = dataset.preprocessing_report.operations.loc[
        dataset.preprocessing_report.operations.loc[:, "stage"]
        == "final_dataset_construction"
    ].iloc[0]
    parameters = final_stage["parameters"]
    assert isinstance(parameters, dict)
    payload = parameters.get("intensity_scale_establishment")
    assert isinstance(payload, dict)
    return payload


def _final_stage_meaning_payload(dataset) -> dict[str, object]:
    parameters = _final_stage_parameters(dataset)
    payload = parameters.get("quantitative_meaning_provenance")
    assert isinstance(payload, dict)
    return payload


def _final_stage_parameters(dataset) -> dict[str, object]:
    assert dataset.preprocessing_report is not None
    final_stage = dataset.preprocessing_report.operations.loc[
        dataset.preprocessing_report.operations.loc[:, "stage"]
        == "final_dataset_construction"
    ].iloc[0]
    parameters = final_stage["parameters"]
    assert isinstance(parameters, dict)
    return parameters


def test_builder_declared_log2_records_declared_establishment_mode_in_provenance() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [1.1, 2.0], "sample_b": [2.0, 3.1]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(policy="identity")
            ),
        )
    )

    workflow_payload = _workflow_establishment_payload(built)
    final_stage_payload = _final_stage_establishment_payload(built)
    assert workflow_payload["establishment_mode"] == "declared"
    assert workflow_payload["evidence_level"] == "declared_by_user"
    assert (
        workflow_payload["establishment_source"]
        == IntensityScaleEstablishmentSource.DECLARED_BY_USER.value
    )
    assert workflow_payload["input_declaration_source"] == (
        "dataset_build_request.input_intensity_scale"
    )
    assert final_stage_payload["establishment_mode"] == "declared"
    meaning_payload = _workflow_meaning_payload(built)
    assert meaning_payload["evidence_mode"] == "inferred_from_scale_contract"
    assert meaning_payload["source_quantity"] is None
    assert meaning_payload["target_quantity"] == "phosphosite_log_abundance"
    assert _final_stage_meaning_payload(built) == meaning_payload


def test_builder_caller_declared_base_meaning_records_declaration_provenance() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [100.0], "sample_b": [200.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            input_intensity_scale="linear",
            quantitative_meaning="phosphosite_abundance",
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(policy="identity")
            ),
        )
    )

    payload = _workflow_meaning_payload(built)
    assert payload["evidence_mode"] == "declared_by_caller"
    assert payload["target_quantity"] == "phosphosite_abundance"
    assert payload["operation_id"] == (
        "phospy.dataset_builder.quantitative_meaning.declaration"
    )
    assert payload["diagnostic_caveat_codes"] == ["quantitative_meaning_user_declared"]
    assert built.provenance is not None
    assert built.provenance.workflow_parameters[
        "quantitative_meaning_caveat_codes"
    ] == ["quantitative_meaning_user_declared"]


def test_builder_rejects_operation_derived_public_quantitative_meaning() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.1, 2.0], "sample_b": [2.0, 3.1]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )

    with pytest.raises(PhosPyInputError, match="may only declare direct input"):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho.index),
                input_intensity_scale="log2",
                quantitative_meaning="contrast_log2_fold_change",
                preprocessing_config=DatasetPreprocessingConfig(
                    intensity_transform=DatasetIntensityTransformConfig(
                        policy="identity"
                    )
                ),
            )
        )


def test_builder_rejects_incompatible_base_meaning_declaration() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.1, 2.0], "sample_b": [2.0, 3.1]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )

    with pytest.raises(
        PhosPyInputError,
        match="phosphosite_abundance.*linear intensity scale",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho.index),
                input_intensity_scale="log2",
                quantitative_meaning="phosphosite_abundance",
                preprocessing_config=DatasetPreprocessingConfig(
                    intensity_transform=DatasetIntensityTransformConfig(
                        policy="identity"
                    )
                ),
            )
        )


def test_total_protein_correction_records_derived_meaning_transition() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [3.0, 7.0], "sample_b": [15.0, 31.0]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    total = pd.DataFrame(
        {"sample_a": [2.0, 5.0], "sample_b": [8.0, 16.0]},
        index=pd.Index(["MAPK14", "GSK3B"], name="protein_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            total=total,
            site_metadata=_site_metadata(phospho.index),
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                ),
                total_protein_correction=DatasetTotalProteinCorrectionConfig(
                    policy="subtract_log_total"
                ),
            ),
        )
    )

    meaning_payload = _workflow_meaning_payload(built)
    assert meaning_payload["evidence_mode"] == "derived_by_phospy_operation"
    assert meaning_payload["source_quantity"] == "phosphosite_log_abundance"
    assert meaning_payload["target_quantity"] == "phospho_total_log_ratio"
    assert meaning_payload["operation_id"] == (
        "phospy.dataset_preprocessing.total_protein_correction.subtract_log_total"
    )
    assert str(meaning_payload["trace_id"]).startswith(
        "total_protein_correction:subtract_log_total:"
    )
    input_names = {
        str(item["name"])
        for item in meaning_payload["input_table_fingerprints"]
        if isinstance(item, dict)
    }
    assert {"dataset.phospho", "dataset.total"}.issubset(input_names)
    output_fingerprint = meaning_payload["output_table_fingerprint"]
    assert isinstance(output_fingerprint, dict)
    assert output_fingerprint["name"] == "dataset.phospho"
    parameters = meaning_payload["parameters"]
    assert isinstance(parameters, dict)
    assert parameters["total_protein_correction_policy"] == "subtract_log_total"
    assert parameters["diagnostics.formula"] == "log2_phospho - log2_total"
    assert built.intensity_scale_state.establishment_provenance is not None
    assert (
        built.intensity_scale_state.establishment_provenance.to_payload()
        == _workflow_establishment_payload(built)
    )


def test_builder_log2_transformation_records_transformed_mode_in_provenance() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [3.0, 7.0], "sample_b": [15.0, 31.0]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                )
            ),
        )
    )

    workflow_payload = _workflow_establishment_payload(built)
    assert workflow_payload["establishment_mode"] == "transformed"
    assert workflow_payload["evidence_level"] == "observed_transformation"
    assert (
        workflow_payload["establishment_source"]
        == IntensityScaleEstablishmentSource.TRANSFORMED_BY_PHOSPY.value
    )
    assert workflow_payload["transformer_name"] == (
        "phospy.science.transformations.transformers.log2.Log2Transformer"
    )
    assert workflow_payload["parameters"]["operation"] == "log2"
    assert workflow_payload["parameters"]["pseudocount"] == 1.0
    assert workflow_payload["parameters"]["affected_matrices"] == ["phospho"]

    assert built.provenance is not None
    stage = next(
        item
        for item in built.provenance.preprocessing_stages
        if item.stage == "intensity_transform"
    )
    diagnostics = stage.diagnostics or {}
    assert diagnostics["pseudocount"] == 1.0
    assert diagnostics["affected_matrices"] == ["phospho"]
    assert isinstance(diagnostics.get("input_phospho_hash"), str)
    assert isinstance(diagnostics.get("output_phospho_hash"), str)


def test_builder_intensity_scale_establishment_payloads_are_fresh() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [3.0, 7.0], "sample_b": [15.0, 31.0]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                )
            ),
        )
    )

    workflow_payload = _workflow_establishment_payload(built)
    workflow_parameters = workflow_payload["parameters"]
    assert isinstance(workflow_parameters, dict)
    workflow_affected = workflow_parameters["affected_matrices"]
    assert isinstance(workflow_affected, list)
    workflow_affected.append("payload-only")

    final_stage_payload = _final_stage_establishment_payload(built)
    final_parameters = final_stage_payload["parameters"]
    assert isinstance(final_parameters, dict)
    final_affected = final_parameters["affected_matrices"]
    assert isinstance(final_affected, list)
    final_affected.append("report-only")

    assert _workflow_establishment_payload(built)["parameters"][
        "affected_matrices"
    ] == ["phospho"]
    assert _final_stage_establishment_payload(built)["parameters"][
        "affected_matrices"
    ] == ["phospho"]


def test_builder_identity_pass_through_without_declared_scale_fails_establishment() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [100.0], "sample_b": [200.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    with pytest.raises(
        TransformationStateEstablishmentError,
        match="pass-through/identity transformer cannot establish scientific input scale",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho.index),
            )
        )


def test_builder_identity_pass_through_with_declared_linear_records_declared_mode() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [100.0], "sample_b": [200.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(policy="identity")
            ),
        )
    )

    workflow_payload = _workflow_establishment_payload(built)
    assert workflow_payload["establishment_mode"] == "declared"
    assert workflow_payload["evidence_level"] == "declared_by_user"
    assert (
        workflow_payload["establishment_source"]
        == IntensityScaleEstablishmentSource.DECLARED_BY_USER.value
    )
    assert workflow_payload["transformer_name"] is None


def test_builder_rejects_suspicious_declared_log2_by_default() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [12000.0, 14000.0], "sample_b": [18000.0, 22000.0]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )

    with pytest.raises(
        TransformationStateEstablishmentError,
        match="declared log2 intensity scale produced",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho.index),
                input_intensity_scale="log2",
                preprocessing_config=DatasetPreprocessingConfig(
                    intensity_transform=DatasetIntensityTransformConfig(
                        policy="identity"
                    )
                ),
            )
        )


def test_builder_records_suspicious_declared_log2_override_in_provenance() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [12000.0, 14000.0], "sample_b": [18000.0, 22000.0]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            input_intensity_scale="log2",
            allow_suspicious_declared_input_intensity_scale=True,
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(policy="identity")
            ),
        )
    )

    workflow_payload = _workflow_establishment_payload(built)
    warnings = workflow_payload["diagnostic_warnings"]
    assert isinstance(warnings, list)
    assert any("declared log2 scale is suspicious" in warning for warning in warnings)
    assert workflow_payload["establishment_mode"] == "declared"
    assert workflow_payload["input_declaration_source"] == (
        "dataset_build_request.input_intensity_scale"
    )
    assert built.provenance is not None
    assert (
        built.provenance.workflow_parameters[
            "allow_suspicious_declared_input_intensity_scale"
        ]
        is True
    )
    assert (
        built.provenance.workflow_parameters[
            "effective_declared_input_intensity_scale_diagnostic_policy"
        ]
        == "warn"
    )
    final_parameters = _final_stage_parameters(built)
    assert final_parameters["allow_suspicious_declared_input_intensity_scale"] is True
    assert (
        final_parameters["effective_declared_input_intensity_scale_diagnostic_policy"]
        == "warn"
    )
    final_payload = _final_stage_establishment_payload(built)
    final_warnings = final_payload["diagnostic_warnings"]
    assert isinstance(final_warnings, list)
    assert any(
        "declared log2 scale is suspicious" in warning for warning in final_warnings
    )


def test_builder_rejects_declared_linear_negative_values_at_final_boundary() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [-1.0, 3.0], "sample_b": [2.0, 4.0]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )

    with pytest.raises(
        DatasetValidationError,
        match="linear phosphosite_abundance must be non-negative",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=_site_metadata(phospho.index),
                input_intensity_scale="linear",
                preprocessing_config=DatasetPreprocessingConfig(
                    intensity_transform=DatasetIntensityTransformConfig(
                        policy="identity"
                    )
                ),
            )
        )


def test_builder_plausible_declared_scale_has_no_warning() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [0.2, 1.0], "sample_b": [2.2, 4.1]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(policy="identity")
            ),
        )
    )

    workflow_payload = _workflow_establishment_payload(built)
    warnings = workflow_payload["diagnostic_warnings"]
    assert isinstance(warnings, list)
    assert warnings == []


def test_final_report_records_intensity_transformation_state_before_and_after() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [3.0, 7.0], "sample_b": [15.0, 31.0]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                )
            ),
        )
    )

    parameters = _final_stage_parameters(built)
    transformation_state = parameters.get("intensity_transformation_state")
    assert isinstance(transformation_state, dict)
    assert transformation_state["before_preprocessing"] == "linear_or_unknown"
    assert transformation_state["after_preprocessing"] == "log2"


def test_preprocessing_operations_include_execution_summary_for_imputation() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 2.0], "sample_b": [2.0, float("nan")]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=1,
                )
            ),
        )
    )

    assert built.preprocessing_report is not None
    operations = built.preprocessing_report.operations
    missing_data_row = operations.loc[
        operations.loc[:, "stage"] == "missing_data"
    ].iloc[0]
    parameters = missing_data_row["parameters"]
    assert isinstance(parameters, dict)
    summary = parameters.get("execution_summary")
    assert isinstance(summary, dict)
    assert summary["imputed_cell_count"] == 1
    assert summary["imputed_row_count"] == 1
    assert summary["imputation_scope"] == "per_row"
