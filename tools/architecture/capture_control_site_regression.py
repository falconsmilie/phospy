"""Emit deterministic SPS and control-site regression data."""

from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from typing import Any, cast

import pandas as pd

from phospy.advanced import (
    ControlSiteSet,
    CorrectionMissingnessPolicy,
    ObservationMask,
    OriginallyMissingCellTracking,
    SpsDiscoveryConfig,
    SpsDiscoveryRequest,
    SpsDiscoveryResult,
    SpsDiscoveryValidationError,
    SpsDiscoveryWorkflow,
    SpsReferenceDataset,
)
from phospy.preprocessing.batch_correction import (
    BatchCorrectionWorkflow,
    BatchCorrectionWorkflowRequest,
)
from phospy.preprocessing.batch_correction.contracts.request import (
    InternalBatchCorrectionControlSiteMode,
    InternalBatchCorrectionControlSiteSource,
    InternalBatchCorrectionImputationPolicy,
    InternalBatchCorrectionMethod,
    InternalBatchCorrectionMissingValuePolicy,
    InternalBatchCorrectionRequest,
    InternalBatchCorrectionStageOrder,
)
from phospy.provenance.hashing import hash_json_payload
from phospy.provenance.models import JsonValue
from phospy.provenance.serialization.batch_correction import (
    batch_correction_provenance_to_payload,
)


def _site_keys() -> pd.Index:
    return pd.Index(
        tuple(
            "phospy:v1|organism=rat|protein_namespace=protein_id|"
            f"protein_identifier=P{position}|residue=S|position={position}"
            for position in range(1, 5)
        ),
        name="site_key",
    )


def _reference(
    dataset_id: str,
    values: dict[str, tuple[float, ...]],
    *,
    index: pd.Index | None = None,
) -> SpsReferenceDataset:
    return SpsReferenceDataset.from_condition_relative_log2(
        dataset_id=dataset_id,
        intensities=pd.DataFrame(
            values, index=_site_keys() if index is None else index
        ),
        condition_by_sample={
            "a_1": "a",
            "a_2": "a",
            "b_1": "b",
            "b_2": "b",
        },
        log2_scale_established_by=f"{dataset_id} governed log2 transform",
        baseline_centering_established_by=f"{dataset_id} governed centring",
        organism="rat",
        baseline_context="condition a",
        reference_context="control-site regression synthetic reference",
        source_name=f"control-site-regression-{dataset_id}",
        source_version="1",
        source_uri=f"https://example.test/control-site-regression/{dataset_id}",
    )


def _references() -> tuple[SpsReferenceDataset, ...]:
    return (
        _reference(
            "reference_a",
            {
                "a_1": (0.0, 0.0, 0.0, 0.0),
                "a_2": (0.0, 0.0, 0.0, 0.0),
                "b_1": (0.1, 0.2, 0.2, 2.0),
                "b_2": (0.1, 0.2, 0.2, 2.0),
            },
        ),
        _reference(
            "reference_b",
            {
                "a_1": (0.0, 0.0, 0.0, 0.0),
                "a_2": (0.0, 0.0, 0.0, 0.0),
                "b_1": (0.2, 0.1, 0.1, 3.0),
                "b_2": (0.2, 0.1, 0.1, 3.0),
            },
        ),
        _reference(
            "reference_c",
            {
                "a_1": (0.0, 0.0, 0.0, 0.0),
                "a_2": (0.0, 0.0, 0.0, 0.0),
                "b_1": (0.3, 0.2, float("nan"), 4.0),
                "b_2": (0.3, 0.2, float("nan"), 4.0),
            },
        ),
    )


def _run(
    references: tuple[SpsReferenceDataset, ...],
    *,
    minimum_datasets_per_site: int,
) -> SpsDiscoveryResult:
    return SpsDiscoveryWorkflow().run(
        SpsDiscoveryRequest(
            reference_datasets=references,
            config=SpsDiscoveryConfig(
                top_n=3,
                minimum_reference_datasets=2,
                minimum_datasets_per_site=minimum_datasets_per_site,
                minimum_shared_sites=3,
            ),
        )
    )


def _legacy_payload(result: SpsDiscoveryResult) -> dict[str, object]:
    legacy = replace(
        result,
        provenance=replace(result.provenance, algorithm_version="1.0.0"),
    )
    payload = legacy.to_payload()
    payload["identity_schema"] = "phospy-sps-discovery-result-v1"
    for record in cast(list[dict[str, Any]], payload["site_ranking"]):
        for statistic in cast(list[dict[str, Any]], record["dataset_statistics"]):
            statistic.pop("rank_quantile")
    identity_payload = {
        key: value for key, value in payload.items() if key != "discovery_identity"
    }
    payload["discovery_identity"] = (
        f"sha256-stable-json-v1:{hash_json_payload(cast(JsonValue, identity_payload))}"
    )
    return payload


def _invalid_validation_payload() -> dict[str, object]:
    matrix = pd.DataFrame(
        {"a_1": (0.0,), "a_2": (0.0,), "b_1": (0.1,), "b_2": (0.1,)},
        index=_site_keys()[:1],
    )
    references = tuple(
        SpsReferenceDataset.from_condition_relative_log2(
            dataset_id="duplicate",
            intensities=matrix,
            condition_by_sample={
                "a_1": "a",
                "a_2": "a",
                "b_1": "b",
                "b_2": "b",
            },
            log2_scale_established_by="invalid fixture scale",
            baseline_centering_established_by="invalid fixture baseline",
        )
        for _ in range(2)
    )
    return (
        SpsDiscoveryWorkflow()
        .validate(
            SpsDiscoveryRequest(
                reference_datasets=references,
                config=SpsDiscoveryConfig(top_n=1, minimum_shared_sites=2),
            )
        )
        .to_payload()
    )


def _invalid_quantitative_state_payload() -> dict[str, object]:
    """Capture rejection of unestablished SPS state without mutating its input."""

    matrix = pd.DataFrame(
        {
            "a_1": (0.0, 0.0, 0.0, 0.0),
            "a_2": (0.0, 0.0, 0.0, 0.0),
            "b_1": (0.1, 0.2, 0.3, 2.0),
            "b_2": (0.1, 0.2, 0.3, 2.0),
        },
        index=_site_keys(),
    )
    governed = SpsReferenceDataset.from_condition_relative_log2(
        dataset_id="invalid-quantitative-state-template",
        intensities=matrix,
        condition_by_sample={
            "a_1": "a",
            "a_2": "a",
            "b_1": "b",
            "b_2": "b",
        },
        log2_scale_established_by="invalid-state fixture scale",
        baseline_centering_established_by="invalid-state fixture baseline",
        organism="rat",
        baseline_context="condition a",
        reference_context="control-site invalid quantitative-state fixture",
        source_name="control-site-invalid-state-template",
        source_version="1",
        source_uri="https://example.test/control-site/invalid-state-template",
    )
    unestablished_state = replace(governed.intensity_scale_state)
    invalid = SpsReferenceDataset(
        dataset_id="invalid-quantitative-state",
        intensities=matrix,
        condition_by_sample=governed.condition_by_sample,
        intensity_scale_state=unestablished_state,
        organism="rat",
        baseline_context="condition a",
        reference_context="control-site invalid quantitative-state fixture",
        source_name="control-site-invalid-state",
        source_version="1",
        source_uri="https://example.test/control-site/invalid-state",
    )
    matrix_before = _frame_payload(invalid.intensities)
    state_before = {
        "kind": unestablished_state.kind.value,
        "quantity": (
            None
            if unestablished_state.quantity is None
            else unestablished_state.quantity.value
        ),
        "is_established": unestablished_state.is_established,
        "established_via": unestablished_state.established_via,
        "has_scale_establishment": (
            unestablished_state.establishment_provenance is not None
        ),
        "has_meaning_establishment": (
            unestablished_state.quantitative_meaning_provenance is not None
        ),
    }
    request = SpsDiscoveryRequest(
        reference_datasets=(invalid, _references()[0]),
        config=SpsDiscoveryConfig(top_n=1, minimum_shared_sites=2),
    )
    try:
        SpsDiscoveryWorkflow().run(request)
    except SpsDiscoveryValidationError as exc:
        validation = exc.validation_result.to_payload()
    else:
        raise AssertionError("invalid SPS quantitative state reached ranking")
    matrix_after = _frame_payload(invalid.intensities)
    state_after = {
        "kind": invalid.intensity_scale_state.kind.value,
        "quantity": (
            None
            if invalid.intensity_scale_state.quantity is None
            else invalid.intensity_scale_state.quantity.value
        ),
        "is_established": invalid.intensity_scale_state.is_established,
        "established_via": invalid.intensity_scale_state.established_via,
        "has_scale_establishment": (
            invalid.intensity_scale_state.establishment_provenance is not None
        ),
        "has_meaning_establishment": (
            invalid.intensity_scale_state.quantitative_meaning_provenance is not None
        ),
    }
    if matrix_after != matrix_before or state_after != state_before:
        raise AssertionError("invalid SPS rejection mutated the submitted matrix/state")
    if invalid.intensity_scale_state is not unestablished_state:
        raise AssertionError("invalid SPS rejection replaced the submitted state")
    return {
        "validation_error": validation,
        "matrix_before": matrix_before,
        "matrix_after": matrix_after,
        "quantitative_state_before": state_before,
        "quantitative_state_after": state_after,
        "quantitative_state_identity_preserved": True,
        "rejected_before_ranking": True,
    }


def _frame_payload(frame: pd.DataFrame) -> dict[str, object]:
    return {
        "index": [str(value) for value in frame.index],
        "index_name": frame.index.name,
        "columns": [str(value) for value in frame.columns],
        "column_name": frame.columns.name,
        "dtypes": [str(dtype) for dtype in frame.dtypes],
        "values": frame.to_numpy().tolist(),
    }


def _correction_payload(controls: ControlSiteSet) -> dict[str, object]:
    sites = _site_keys()
    samples = ("sample_1", "sample_2", "sample_3", "sample_4")
    phospho = pd.DataFrame(
        {
            "sample_1": (1.0, 2.0, 5.0, 8.0),
            "sample_2": (2.0, 4.0, 5.5, 7.0),
            "sample_3": (3.0, 1.0, 6.0, 9.0),
            "sample_4": (4.0, 3.0, 7.0, 6.0),
        },
        index=sites,
    )
    config = InternalBatchCorrectionRequest(
        method=InternalBatchCorrectionMethod.SPS_RUV_STYLE,
        batch_column="batch",
        condition_columns=("condition",),
        replicate_column=None,
        control_site_source=(InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED),
        control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
        missing_value_policy=InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING,
        imputation_policy=InternalBatchCorrectionImputationPolicy.NONE,
        n_unwanted_factors=1,
        stage_order=(
            InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM
        ),
        diagnostics_enabled=True,
    )
    result = BatchCorrectionWorkflow().run(
        BatchCorrectionWorkflowRequest(
            phospho=phospho,
            config=config,
            sample_metadata=pd.DataFrame(
                {
                    "batch": ("batch_1", "batch_1", "batch_2", "batch_2"),
                    "condition": ("a", "b", "a", "b"),
                },
                index=samples,
            ),
            control_site_set=controls,
            missingness_policy=CorrectionMissingnessPolicy(
                originally_missing_cells_tracked_by=(
                    OriginallyMissingCellTracking.OBSERVATION_MASK
                ),
                observation_mask=ObservationMask(
                    feature_ids=tuple(str(value) for value in sites),
                    sample_ids=samples,
                    originally_missing_cells=(),
                ),
            ),
        )
    )
    applied = result.corrected_preprocessing_output
    if applied.output_observation_mask is None or applied.corrected_cell_status is None:
        raise AssertionError("correction omitted its mask or cell classification")
    return {
        "corrected_matrix": _frame_payload(result.corrected_matrix),
        "diagnostics": dict(result.diagnostics),
        "warnings": list(result.warnings),
        "output_observation_mask": _frame_payload(applied.output_observation_mask),
        "corrected_cell_status": _frame_payload(applied.corrected_cell_status),
        "provenance": batch_correction_provenance_to_payload(result.provenance),
    }


def build_comparison_payload() -> dict[str, object]:
    references = _references()
    partial = _run(references, minimum_datasets_per_site=2)
    permuted = _run(tuple(reversed(references)), minimum_datasets_per_site=2)
    complete = _run(references[:2], minimum_datasets_per_site=2)
    suspicious = _reference(
        "suspicious",
        {
            "a_1": (1000.0, 2000.0, 3000.0, 4000.0),
            "a_2": (1000.0, 2000.0, 3000.0, 4000.0),
            "b_1": (1000.0, 2000.0, 3000.0, 4000.0),
            "b_2": (1000.0, 2000.0, 3000.0, 4000.0),
        },
    )
    suspicious_provenance = suspicious.intensity_scale_state.establishment_provenance
    if suspicious_provenance is None:
        raise AssertionError("SPS assertion omitted scale-establishment provenance")
    legacy_payload = _legacy_payload(complete)
    legacy_restored = SpsDiscoveryResult.from_payload(legacy_payload)
    legacy_controls_payload = legacy_restored.control_site_set.to_payload()
    legacy_controls = ControlSiteSet.from_payload(legacy_controls_payload)
    if legacy_controls.to_payload() != legacy_controls_payload:
        raise AssertionError("legacy SPS controls changed during restoration")
    current_restored = SpsDiscoveryResult.from_payload(partial.to_payload())
    controls = ControlSiteSet.from_payload(partial.control_site_set.to_payload())
    cases: dict[str, object] = {
        "complete": complete.to_payload(),
        "partial": partial.to_payload(),
        "permuted": permuted.to_payload(),
        "current_restored": current_restored.to_payload(),
        "legacy_input": legacy_payload,
        "legacy_restored": legacy_restored.to_payload(),
        "legacy_controls": legacy_controls_payload,
        "legacy_controls_restored": legacy_controls.to_payload(),
        "controls": controls.to_payload(),
        "correction_from_restored_controls": _correction_payload(controls),
        "suspicious_warnings": list(suspicious_provenance.diagnostic_warnings),
        "invalid_metadata": _invalid_validation_payload(),
        "invalid_quantitative_state": _invalid_quantitative_state_payload(),
    }
    canonical = json.dumps(
        cases,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
    ).encode("utf-8")
    return {
        "schema": "phospy-control-site-regression-v1",
        "cases": cases,
        "cases_sha256": hashlib.sha256(canonical).hexdigest(),
    }


def main() -> None:
    print(
        json.dumps(
            build_comparison_payload(),
            ensure_ascii=False,
            allow_nan=False,
            separators=(",", ":"),
        )
    )


if __name__ == "__main__":
    main()
