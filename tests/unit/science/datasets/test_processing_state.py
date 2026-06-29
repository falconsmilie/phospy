from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from phospy.errors import DatasetProcessingStateError, DatasetValidationError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.datasets.preprocessing.policy_models import (
    MissingDataPolicy,
    TotalProteinCorrectionPolicy,
)
from phospy.science.datasets.processing_state import (
    MissingDataDiagnosticsV1,
    MissingDataState,
    NormalisationState,
    SiteSequenceResolutionRowDiagnostic,
    SiteSequenceResolutionState,
    TotalProteinCorrectionDiagnosticsV1,
    TotalProteinCorrectionState,
)
from phospy.science.references.models import Organism
from phospy.science.transformations.models import QuantitativeMeaning
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns


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
        method_parameters={"min_observed_values": 1},
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
    scale = supported_log2_intensity_scale_state(
        has_total_matrix=True
    ).with_quantitative_meaning(QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO)
    base_state = supported_log2_processing_state(has_total_matrix=True)
    processing_state = replace(
        base_state,
        intensity_scale=scale,
        total_protein_correction=_applied_total_correction_state(),
    )

    with pytest.raises(DatasetValidationError, match="requires dataset.total"):
        AnalysisReadyPhosphoDataset(
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
        AnalysisReadyPhosphoDataset(
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
        AnalysisReadyPhosphoDataset(
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

    dataset = AnalysisReadyPhosphoDataset(
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
