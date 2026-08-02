from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import replace
from typing import cast

import numpy as np
import pandas as pd
import pytest

from phospy.errors import DatasetProcessingStateError, DatasetValidationError
from phospy.errors.input import PhosPyInputError
from phospy.provenance.hashing import hash_json_payload
from phospy.provenance.immutability import FrozenJsonMapping
from phospy.science.datasets.preprocessing.policy_models import (
    MissingDataPolicy,
    TotalProteinCorrectionPolicy,
)
from phospy.science.datasets.processing_state import (
    JsonValue,
    MissingDataDiagnostics,
    MissingDataDiagnosticsV1,
    MissingDataState,
    NormalisationState,
    SiteSequenceResolutionRowDiagnostic,
    SiteSequenceResolutionState,
    TotalProteinCorrectionDiagnostics,
    TotalProteinCorrectionDiagnosticsV1,
    TotalProteinCorrectionState,
)
from phospy.science.references.models import Organism
from phospy.science.transformations.models import QuantitativeMeaning
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
    with_restored_quantitative_meaning_for_tests,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns


class _DuplicateKeyMapping(Mapping[str, object]):
    def __iter__(self) -> Iterator[str]:
        return iter(("duplicate", "duplicate"))

    def __len__(self) -> int:
        return 2

    def __getitem__(self, key: str) -> object:
        if key == "duplicate":
            return "value"
        raise KeyError(key)


def _missing_data_diagnostics() -> MissingDataDiagnosticsV1:
    return MissingDataDiagnosticsV1(
        missing_data_policy="impute_row_median",
        imputation_method_id="row_median",
        imputation_method_family="deterministic_row_statistic",
        input_missing_cell_count=1,
        output_missing_cell_count=0,
        imputed_cell_count=1,
        affected_row_count=1,
        affected_column_count=1,
        affected_row_ids=("row_a",),
        affected_column_ids=("sample_a",),
        imputed_row_ids=("row_a",),
        imputed_column_ids=("sample_a",),
        dropped_row_ids=(),
        method_parameters={
            "min_observed_values": 1,
            "input_scale": "linear",
            "imputation_operation_order": "no_intensity_transform",
        },
        imputation_input_scale="linear",
        imputation_input_scale_source="caller_selected",
        imputation_operation_order="no_intensity_transform",
        stage_order=("missing_data",),
        missingness_mask_hash="missingness-hash",
        imputation_mask_hash="imputation-hash",
        rows_not_imputable=(),
    )


def _total_correction_diagnostics(
    policy: str = "subtract_log_total",
) -> TotalProteinCorrectionDiagnosticsV1:
    return TotalProteinCorrectionDiagnosticsV1(
        policy=policy,
        requested_policy=policy,
        resolved_policy=policy,
        formula=("log2_phospho - log2_total" if policy != "none" else None),
        requires_log_scale=(True if policy != "none" else None),
        input_scale=("log2" if policy != "none" else None),
        output_scale=("log2_ratio" if policy != "none" else None),
        quantitative_meaning="phospho_total_log_ratio",
    )


def _applied_total_correction_state(**overrides: object) -> TotalProteinCorrectionState:
    payload = {
        "policy": "subtract_log_total",
        "applied": True,
        "formula": "log2_phospho - log2_total",
        "requires_log_scale": True,
        "input_scale": "log2",
        "output_scale": "log2_ratio",
        "quantitative_meaning": QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO,
        "diagnostics": _total_correction_diagnostics(),
    }
    payload.update(overrides)
    return TotalProteinCorrectionState(**payload)


def test_missing_data_diagnostics_constructor_inputs_are_recursively_detached() -> None:
    method_parameters = {"nested": {"items": [1, {"score": 2.0}]}}
    row_medians_used = {"row_a": 1.25}
    per_column_distribution_parameters = {"sample_a": {"items": [1, {"score": 3.0}]}}

    diagnostics = MissingDataDiagnosticsV1(
        missing_data_policy="impute_row_median",
        imputation_method_id="row_median",
        imputation_method_family="deterministic_row_statistic",
        input_missing_cell_count=1,
        output_missing_cell_count=0,
        imputed_cell_count=1,
        affected_row_count=1,
        affected_column_count=1,
        affected_row_ids=("row_a",),
        affected_column_ids=("sample_a",),
        imputed_row_ids=("row_a",),
        imputed_column_ids=("sample_a",),
        dropped_row_ids=(),
        method_parameters=method_parameters,
        imputation_input_scale="linear",
        imputation_input_scale_source="caller_selected",
        imputation_operation_order="no_intensity_transform",
        stage_order=("missing_data",),
        missingness_mask_hash="missingness-hash",
        imputation_mask_hash="imputation-hash",
        rows_not_imputable=(),
        row_medians_used=row_medians_used,
        per_column_distribution_parameters=per_column_distribution_parameters,
    )

    method_parameters["nested"]["items"].append("constructor-only")
    method_parameters["nested"]["items"][1]["score"] = 9.0
    row_medians_used["row_a"] = 9.0
    per_column_distribution_parameters["sample_a"]["items"].append("constructor-only")
    per_column_distribution_parameters["sample_a"]["items"][1]["score"] = 9.0

    assert isinstance(diagnostics.method_parameters, FrozenJsonMapping)
    assert isinstance(diagnostics.row_medians_used, FrozenJsonMapping)
    assert isinstance(diagnostics.per_column_distribution_parameters, FrozenJsonMapping)
    assert diagnostics.method_parameters["nested"]["items"] == (1, {"score": 2.0})
    assert diagnostics.row_medians_used["row_a"] == 1.25
    assert diagnostics.per_column_distribution_parameters["sample_a"]["items"] == (
        1,
        {"score": 3.0},
    )
    assert diagnostics.to_payload()["method_parameters"]["nested"]["items"] == [
        1,
        {"score": 2.0},
    ]
    assert diagnostics.to_payload()["imputation_input_scale"] == "linear"


def test_missing_data_state_derives_imputation_scale_order_from_diagnostics() -> None:
    state = MissingDataState(
        policy="impute_row_median",
        min_observed_values=1,
        complete_matrix=True,
        imputed=True,
        diagnostics=_missing_data_diagnostics(),
    )

    assert state.imputation_input_scale == "linear"
    assert state.imputation_operation_order == "no_intensity_transform"


def test_missing_data_diagnostics_payload_is_fresh_and_nested_detached() -> None:
    diagnostics = MissingDataDiagnosticsV1(
        missing_data_policy="impute_minprob",
        imputation_method_id="minprob",
        imputation_method_family="left_censored_random",
        input_missing_cell_count=1,
        output_missing_cell_count=0,
        imputed_cell_count=1,
        affected_row_count=1,
        affected_column_count=1,
        affected_row_ids=("row_a",),
        affected_column_ids=("sample_a",),
        imputed_row_ids=("row_a",),
        imputed_column_ids=("sample_a",),
        dropped_row_ids=(),
        method_parameters={"nested": {"items": [1, {"score": 2.0}]}},
        stage_order=("intensity_transform", "missing_data"),
        missingness_mask_hash="missingness-hash",
        imputation_mask_hash="imputation-hash",
        rows_not_imputable=(),
        per_column_distribution_parameters={"sample_a": {"items": [1, {"score": 3.0}]}},
    )

    payload = diagnostics.to_payload()
    payload["method_parameters"]["nested"]["items"].append("payload-only")
    payload["method_parameters"]["nested"]["items"][1]["score"] = 9.0
    payload["per_column_distribution_parameters"]["sample_a"]["items"].append(
        "payload-only"
    )
    payload["per_column_distribution_parameters"]["sample_a"]["items"][1]["score"] = 9.0
    mapping_value = diagnostics["method_parameters"]
    mapping_value["nested"]["items"].append("mapping-only")

    fresh_payload = diagnostics.to_payload()
    assert fresh_payload["method_parameters"]["nested"]["items"] == [
        1,
        {"score": 2.0},
    ]
    assert fresh_payload["per_column_distribution_parameters"]["sample_a"]["items"] == [
        1,
        {"score": 3.0},
    ]
    assert diagnostics["method_parameters"]["nested"]["items"] == [
        1,
        {"score": 2.0},
    ]


def test_total_protein_diagnostics_attributes_and_payload_cannot_diverge() -> None:
    corrected_mapping = {"row_a": "tp_a"}
    uncorrected_reasons = {"row_b": "no_total_match"}
    diagnostics = TotalProteinCorrectionDiagnosticsV1(
        policy="subtract_log_total",
        requested_policy="subtract_log_total",
        resolved_policy="subtract_log_total",
        formula="log2_phospho - log2_total",
        requires_log_scale=True,
        input_scale="log2",
        output_scale="log2_ratio",
        quantitative_meaning="mixed_phospho_total_log_ratio_and_phosphosite_log_abundance",
        corrected_row_count=1,
        uncorrected_row_count=1,
        corrected_phosphosite_row_ids=("row_a",),
        corrected_phosphosite_to_total_protein_row_id=corrected_mapping,
        unmatched_phosphosite_row_ids=("row_b",),
        uncorrected_phosphosite_row_reasons=uncorrected_reasons,
    )

    corrected_mapping["row_a"] = "constructor-only"
    uncorrected_reasons["row_b"] = "constructor-only"
    assert isinstance(
        diagnostics.corrected_phosphosite_to_total_protein_row_id, FrozenJsonMapping
    )
    assert isinstance(
        diagnostics.uncorrected_phosphosite_row_reasons, FrozenJsonMapping
    )
    assert diagnostics.corrected_phosphosite_to_total_protein_row_id == {
        "row_a": "tp_a"
    }
    assert diagnostics.to_payload()[
        "corrected_phosphosite_to_total_protein_row_id"
    ] == {"row_a": "tp_a"}

    payload = diagnostics.to_payload()
    payload["corrected_phosphosite_to_total_protein_row_id"]["row_a"] = "payload-only"
    mapping_value = diagnostics["uncorrected_phosphosite_row_reasons"]
    mapping_value["row_b"] = "mapping-only"

    assert diagnostics.to_payload()[
        "corrected_phosphosite_to_total_protein_row_id"
    ] == {"row_a": "tp_a"}
    assert diagnostics.to_payload()["uncorrected_phosphosite_row_reasons"] == {
        "row_b": "no_total_match"
    }


def test_processing_diagnostics_round_trip_and_hash_are_stable() -> None:
    diagnostics_items = (
        (
            _missing_data_diagnostics(),
            MissingDataDiagnostics,
            "dataset.metadata.processing_state.missing_data.diagnostics",
        ),
        (
            _total_correction_diagnostics(),
            TotalProteinCorrectionDiagnostics,
            "dataset.metadata.processing_state.total_protein_correction.diagnostics",
        ),
    )

    for diagnostics, diagnostics_type, field_name in diagnostics_items:
        first_payload = diagnostics.to_payload()
        second_payload = diagnostics.to_payload()
        restored = diagnostics_type.from_payload(first_payload, field_name=field_name)

        assert first_payload == second_payload
        assert first_payload is not second_payload
        assert restored.to_payload() == first_payload
        assert hash_json_payload(cast(JsonValue, first_payload)) == hash_json_payload(
            cast(JsonValue, second_payload)
        )


def test_processing_diagnostics_reject_invalid_json_keys_and_values() -> None:
    missing_data_kwargs = {
        **_missing_data_diagnostics().to_payload(),
        "imputation_method_id": "forbid",
        "imputation_mask_hash": None,
    }
    missing_data_kwargs.pop("row_medians_used", None)
    total_kwargs = {
        **_total_correction_diagnostics().to_payload(),
        "corrected_phosphosite_row_ids": ("row_a",),
    }

    with pytest.raises(PhosPyInputError, match="keys must be strings"):
        MissingDataDiagnosticsV1(
            **{**missing_data_kwargs, "method_parameters": {1: "bad"}}
        )
    with pytest.raises(PhosPyInputError, match="duplicate JSON object key"):
        MissingDataDiagnosticsV1(
            **{**missing_data_kwargs, "method_parameters": _DuplicateKeyMapping()}
        )
    with pytest.raises(PhosPyInputError, match="finite JSON number"):
        MissingDataDiagnosticsV1(
            **{**missing_data_kwargs, "method_parameters": {"bad": float("nan")}}
        )
    with pytest.raises(PhosPyInputError, match="JSON-compatible"):
        MissingDataDiagnosticsV1(
            **{**missing_data_kwargs, "method_parameters": {"bad": object()}}
        )
    with pytest.raises(PhosPyInputError, match="finite JSON number"):
        MissingDataDiagnosticsV1(
            **{**missing_data_kwargs, "row_medians_used": {"row_a": float("inf")}}
        )
    with pytest.raises(PhosPyInputError, match="keys must be strings"):
        TotalProteinCorrectionDiagnosticsV1(
            **{
                **total_kwargs,
                "corrected_phosphosite_to_total_protein_row_id": {1: "tp_a"},
            }
        )
    with pytest.raises(PhosPyInputError, match="duplicate JSON object key"):
        TotalProteinCorrectionDiagnosticsV1(
            **{
                **total_kwargs,
                "corrected_phosphosite_to_total_protein_row_id": _DuplicateKeyMapping(),
            }
        )
    with pytest.raises(PhosPyInputError, match="finite JSON number"):
        TotalProteinCorrectionDiagnosticsV1(
            **{
                **total_kwargs,
                "corrected_phosphosite_to_total_protein_row_id": {
                    "row_a": float("nan")
                },
            }
        )
    with pytest.raises(PhosPyInputError, match="JSON-compatible"):
        TotalProteinCorrectionDiagnosticsV1(
            **{
                **total_kwargs,
                "corrected_phosphosite_to_total_protein_row_id": {"row_a": object()},
            }
        )


def test_processing_diagnostics_optional_mapping_fields_preserve_presence() -> None:
    missing_absent = _missing_data_diagnostics().to_payload()
    assert missing_absent["row_medians_used"] == {}
    assert "per_column_distribution_parameters" not in missing_absent

    missing_present = MissingDataDiagnosticsV1(
        **{
            **_missing_data_diagnostics().to_payload(),
            "per_column_distribution_parameters": {"sample_a": {"observed_count": 2}},
        }
    )
    assert missing_present.to_payload()["per_column_distribution_parameters"] == {
        "sample_a": {"observed_count": 2}
    }

    total_absent = _total_correction_diagnostics().to_payload()
    assert "corrected_phosphosite_to_total_protein_row_id" not in total_absent
    assert "uncorrected_phosphosite_row_reasons" not in total_absent

    total_present = TotalProteinCorrectionDiagnosticsV1(
        **{
            **total_absent,
            "corrected_phosphosite_to_total_protein_row_id": {"row_a": "tp_a"},
            "uncorrected_phosphosite_row_reasons": {"row_b": "no_total_match"},
        }
    )
    assert total_present.to_payload()[
        "corrected_phosphosite_to_total_protein_row_id"
    ] == {"row_a": "tp_a"}
    assert total_present.to_payload()["uncorrected_phosphosite_row_reasons"] == {
        "row_b": "no_total_match"
    }


def test_missing_data_diagnostics_existing_scientific_failures_still_specific() -> None:
    with pytest.raises(PhosPyInputError, match="imputed_row_count must match"):
        MissingDataDiagnosticsV1(
            **{**_missing_data_diagnostics().to_payload(), "imputed_row_count": 2}
        )
    with pytest.raises(PhosPyInputError, match="imputation_mask_hash is required"):
        MissingDataDiagnosticsV1(
            **{**_missing_data_diagnostics().to_payload(), "imputation_mask_hash": None}
        )
    with pytest.raises(PhosPyInputError, match="neighbour_count must be >= 1"):
        MissingDataDiagnosticsV1(
            **{**_missing_data_diagnostics().to_payload(), "neighbour_count": 0}
        )


def _conflict_row(
    *,
    action: str = "preserve_existing",
    conflict_policy: str = "preserve_existing",
) -> SiteSequenceResolutionRowDiagnostic:
    return SiteSequenceResolutionRowDiagnostic(
        row_index=0,
        row_id="MAPK14;Y182;",
        site_id="MAPK14;Y182;",
        status="existing_sequence_conflict",
        existing_site_sequence="XXXXX",
        fasta_site_sequence="AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
        resolved_site_sequence="XXXXX",
        action=action,
        reason="existing site_sequence conflicts with FASTA-derived sequence",
        conflict_policy=conflict_policy,
        resolver_version="resolver.v1",
        fasta_source_path="proteome.fasta",
        fasta_sha256="sha256",
    )


def _site_sequence_state(**overrides: object) -> SiteSequenceResolutionState:
    payload = {
        "configured": True,
        "mode": "fill_missing_only",
        "flank_size": 7,
        "fasta_source_path": "proteome.fasta",
        "fasta_source_label": "test reference",
        "fasta_sha256": "sha256",
        "resolver_version": "resolver.v1",
        "resolved_site_count": 0,
        "unresolved_site_count": 0,
        "unresolved_counts_by_reason": {},
        "filled_missing_count": 0,
        "replaced_existing_count": 0,
        "preserved_existing_count": 0,
        "existing_sequence_conflict_count": 0,
        "conflict_policy": "preserve_existing",
        "row_diagnostics": (),
    }
    payload.update(overrides)
    return SiteSequenceResolutionState(**payload)


def _dataset_frames() -> tuple[pd.DataFrame, pd.DataFrame]:
    index = protein_site_key_index(
        protein_identifiers=["MAPK14", "AKT1"],
        sites=["Y182", "T308"],
    )
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 2.0], "sample_b": [1.5, 2.5]},
        index=index.copy(),
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": index.astype(str).tolist(),
            "display_id": ["MAPK14;Y182;", "AKT1;T308;"],
            **site_key_context_columns(index),
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [
                "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
            ],
            "protein_id": ["MAPK14", "AKT1"],
        },
        index=index.copy(),
    )
    return phospho, site_metadata


def test_missing_data_state_rejects_complete_matrix_with_missing_flags() -> None:
    with pytest.raises(DatasetProcessingStateError, match="complete_matrix"):
        MissingDataState(
            policy=MissingDataPolicy.FORBID,
            min_observed_values=None,
            complete_matrix=True,
            imputed=False,
            has_missing_values=True,
        )


def test_missing_data_state_rejects_complete_matrix_with_missing_count() -> None:
    with pytest.raises(DatasetProcessingStateError, match="missing_value_count"):
        MissingDataState(
            policy=MissingDataPolicy.FORBID,
            min_observed_values=None,
            complete_matrix=True,
            imputed=False,
            missing_value_count=1,
        )


def test_missing_data_state_rejects_no_missing_flag_with_missing_count() -> None:
    with pytest.raises(DatasetProcessingStateError, match="has_missing_values"):
        MissingDataState(
            policy=MissingDataPolicy.FORBID,
            min_observed_values=None,
            complete_matrix=False,
            imputed=False,
            has_missing_values=False,
            missing_value_count=1,
        )


def test_missing_data_state_rejects_negative_counts() -> None:
    with pytest.raises(DatasetProcessingStateError, match="missing_value_count"):
        MissingDataState(
            policy=MissingDataPolicy.FORBID,
            min_observed_values=None,
            complete_matrix=False,
            imputed=False,
            missing_value_count=-1,
        )

    with pytest.raises(DatasetProcessingStateError, match="min_observed_values"):
        MissingDataState(
            policy=MissingDataPolicy.FORBID,
            min_observed_values=-1,
            complete_matrix=True,
            imputed=False,
        )


def test_missing_data_state_rejects_imputed_without_provenance() -> None:
    with pytest.raises(DatasetProcessingStateError, match="imputation provenance"):
        MissingDataState(
            policy=MissingDataPolicy.IMPUTE_ROW_MEDIAN,
            min_observed_values=1,
            complete_matrix=True,
            imputed=True,
        )


def test_missing_data_state_rejects_diagnostics_contradicting_missing_flag() -> None:
    with pytest.raises(DatasetProcessingStateError, match="has_missing_values"):
        MissingDataState(
            policy=MissingDataPolicy.IMPUTE_ROW_MEDIAN,
            min_observed_values=1,
            complete_matrix=True,
            imputed=True,
            diagnostics=_missing_data_diagnostics(),
            has_missing_values=True,
        )


def test_normalisation_state_rejects_empty_or_unknown_method() -> None:
    with pytest.raises(DatasetProcessingStateError, match="normalisation.policy"):
        NormalisationState(policy="")

    with pytest.raises(DatasetProcessingStateError, match="must be one of"):
        NormalisationState(policy="not_supported")


def test_site_sequence_state_rejects_empty_sequence_source() -> None:
    with pytest.raises(DatasetProcessingStateError, match="source identifier"):
        _site_sequence_state(
            fasta_source_path=None,
            fasta_source_label=None,
            fasta_sha256=None,
        )


def test_site_sequence_state_rejects_reference_resolution_without_reference_id() -> (
    None
):
    with pytest.raises(DatasetProcessingStateError, match="reference-resolved"):
        _site_sequence_state(
            fasta_sha256=None,
            resolver_version=None,
            resolved_site_count=1,
        )


def test_site_sequence_state_rejects_conflict_rows_with_zero_count() -> None:
    with pytest.raises(DatasetProcessingStateError, match="conflict row diagnostics"):
        _site_sequence_state(row_diagnostics=(_conflict_row(),))


def test_site_sequence_state_rejects_negative_counts() -> None:
    with pytest.raises(DatasetProcessingStateError, match="conflict_count"):
        _site_sequence_state(existing_sequence_conflict_count=-1)

    with pytest.raises(DatasetProcessingStateError, match="unresolved_site_count"):
        _site_sequence_state(unresolved_site_count=-1)


def test_site_sequence_state_rejects_empty_or_unknown_conflict_policy() -> None:
    with pytest.raises(DatasetProcessingStateError, match="conflict_policy"):
        _site_sequence_state(conflict_policy="")

    with pytest.raises(DatasetProcessingStateError, match="must be one of"):
        _site_sequence_state(conflict_policy="prefer_reference")


def test_site_sequence_state_rejects_error_policy_with_tolerated_conflicts() -> None:
    with pytest.raises(DatasetProcessingStateError, match="cannot be 'error'"):
        _site_sequence_state(
            existing_sequence_conflict_count=1,
            conflict_policy="error",
            row_diagnostics=(
                _conflict_row(action="preserve_existing", conflict_policy="error"),
            ),
        )


def test_total_correction_state_rejects_applied_without_method() -> None:
    with pytest.raises(DatasetProcessingStateError, match="formula"):
        _applied_total_correction_state(formula=None)


def test_total_correction_state_rejects_not_applied_with_correction_provenance() -> (
    None
):
    with pytest.raises(DatasetProcessingStateError, match="correction provenance"):
        TotalProteinCorrectionState(
            policy=TotalProteinCorrectionPolicy.NONE,
            applied=False,
            quantitative_meaning=QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
            diagnostics=_total_correction_diagnostics("subtract_log_total"),
        )


def test_total_correction_state_rejects_empty_method() -> None:
    with pytest.raises(DatasetProcessingStateError, match="formula"):
        _applied_total_correction_state(formula="")


def test_dataset_construction_rejects_total_correction_state_without_total_matrix() -> (
    None
):
    phospho, site_metadata = _dataset_frames()
    scale = with_restored_quantitative_meaning_for_tests(
        supported_log2_intensity_scale_state(has_total_matrix=True),
        QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO,
    )
    base_state = supported_log2_processing_state(has_total_matrix=True)
    processing_state = replace(
        base_state,
        intensity_scale=scale,
        total_protein_correction=_applied_total_correction_state(),
    )

    with pytest.raises(DatasetValidationError, match="requires dataset.total"):
        trusted_analysis_ready_dataset_from_tables(
            phospho=phospho,
            site_metadata=site_metadata,
            intensity_scale_state=scale,
            processing_state=processing_state,
            organism=Organism.RAT,
            total=None,
        )


def test_dataset_rejects_no_missing_state_when_matrix_contains_missing_values() -> None:
    phospho, site_metadata = _dataset_frames()
    phospho.iloc[0, 0] = np.nan

    with pytest.raises(DatasetValidationError, match="claims no missing values"):
        trusted_analysis_ready_dataset_from_tables(
            phospho=phospho,
            site_metadata=site_metadata,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
            organism=Organism.RAT,
        )


def test_dataset_rejects_observation_mask_with_mismatched_labels() -> None:
    phospho, site_metadata = _dataset_frames()
    mask = pd.DataFrame(
        True,
        index=pd.Index(["wrong_a", "wrong_b"], name=phospho.index.name),
        columns=phospho.columns.copy(),
    )

    with pytest.raises(
        DatasetValidationError, match="imputation_observation_mask.index"
    ):
        trusted_analysis_ready_dataset_from_tables(
            phospho=phospho,
            site_metadata=site_metadata,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
            organism=Organism.RAT,
            imputation_observation_mask=mask,
        )


def test_valid_processing_state_constructs_dataset_successfully() -> None:
    phospho, site_metadata = _dataset_frames()

    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
        organism=Organism.RAT,
    )

    assert dataset.processing_state.missing_data.complete_matrix is True
    assert dataset.processing_state.missing_data.missing_value_count == 0
