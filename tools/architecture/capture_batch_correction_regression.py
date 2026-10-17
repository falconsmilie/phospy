"""Emit deterministic batch-correction regression data."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from typing import Any, cast
from unittest.mock import patch

import numpy as np
import pandas as pd

from phospy.advanced import (
    CorrectionMissingnessPolicy,
    ObservationMask,
    OriginallyMissingCellTracking,
    TemporaryImputationMethod,
    TemporaryImputationPolicy,
)
from phospy.errors import PhosPyInputError
from phospy.preparation_records.batch_correction import (
    BATCH_CORRECTION_CONFOUNDING_NOT_APPLICABLE,
    BATCH_CORRECTION_STATUS_DISABLED,
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
)
from phospy.preprocessing.batch_correction import (
    BatchCorrectionWorkflow,
    BatchCorrectionWorkflowRequest,
    BatchCorrectionWorkflowResult,
)
from phospy.preprocessing.batch_correction.contracts.configuration import (
    DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
    DatasetBatchCorrectionConfig,
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
from phospy.preprocessing.batch_correction.execution.linear import (
    LinearResidualizeBatchCorrectionEngine,
)
from phospy.preprocessing.batch_correction.execution.ruv_iii_executor import (
    _prepare_matrix,
)
from phospy.preprocessing.control_sites.contracts import (
    ControlSiteSet,
    ControlSiteSourceMetadata,
)
from phospy.preprocessing.provenance.models import PreprocessingStageProvenance
from phospy.provenance.models import EnvironmentProvenance
from phospy.provenance.serialization.batch_correction import (
    batch_correction_provenance_from_payload,
    batch_correction_provenance_to_payload,
)
from phospy.run_lineage.models import RunProvenance
from phospy.run_lineage.serialization import from_payload, to_payload

_CORRECTION_ENVIRONMENT_COLLECTOR = (
    "phospy.preprocessing.batch_correction.provenance.assembly."
    "collect_batch_correction_environment_provenance"
)


def _json_scalar(value: object) -> object:
    if pd.isna(cast(Any, value)):
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value


def _frame_payload(frame: pd.DataFrame) -> dict[str, object]:
    return {
        "index": [str(value) for value in frame.index],
        "index_name": frame.index.name,
        "columns": [str(value) for value in frame.columns],
        "column_name": frame.columns.name,
        "dtypes": [str(dtype) for dtype in frame.dtypes],
        "values": [
            [_json_scalar(value) for value in row]
            for row in frame.to_numpy(dtype="object").tolist()
        ],
    }


def _is_numerical_frame_payload(value: object) -> bool:
    if not isinstance(value, dict):
        return False
    dtypes = value.get("dtypes")
    return (
        isinstance(dtypes, list)
        and bool(dtypes)
        and all(
            isinstance(dtype, str)
            and (dtype.startswith("float") or dtype.startswith("int"))
            for dtype in dtypes
        )
        and isinstance(value.get("values"), list)
    )


def _normalise_numerical_frames(
    value: object,
    *,
    field_name: str | None = None,
) -> object:
    """Normalise runtime-sensitive numerical data for the invariant hash.

    Matrix values remain covered by the separate tolerance-aware fixture.  Codec
    bytes are checked exactly against their payloads, but cannot be frozen across
    numerical-library lanes because those payloads contain numerical fingerprints.
    """

    if _is_numerical_frame_payload(value):
        assert isinstance(value, dict)
        return {
            key: (
                "compared-separately-with-rtol-1e-10-atol-1e-12"
                if key == "values"
                else _normalise_numerical_frames(item, field_name=str(key))
            )
            for key, item in value.items()
        }
    if isinstance(value, dict):
        normalised: dict[str, object] = {}
        for key, item in value.items():
            key_text = str(key)
            if key_text in {"direct_provenance_bytes", "payload_bytes"}:
                normalised[key_text] = "checked-exactly-against-runtime-payload"
            elif key_text in {
                "exact_hash_value",
                "tolerance_hash_value",
                "input_hash",
                "output_hash",
                "phospho_input_hash",
                "phospho_output_hash",
            }:
                normalised[key_text] = "runtime-numerical-fingerprint"
            else:
                normalised[key_text] = _normalise_numerical_frames(
                    item,
                    field_name=key_text,
                )
        return normalised
    if isinstance(value, list):
        return [
            _normalise_numerical_frames(item, field_name=field_name) for item in value
        ]
    if isinstance(value, float):
        rounded = round(value, 10)
        return 0.0 if rounded == 0.0 else rounded
    if (
        isinstance(value, str)
        and field_name in {"dtype", "dtypes"}
        and (value in {"object", "str"} or value.startswith("string"))
    ):
        return "runtime-string-dtype"
    return value


def _canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
    ).encode("utf-8")


def _control_site_set(site_keys: tuple[str, ...]) -> ControlSiteSet:
    return ControlSiteSet.from_site_keys(
        site_keys,
        source_metadata=ControlSiteSourceMetadata(
            organism="rat",
            identifier_namespace="site_key",
            source_name="batch-correction-regression-controls",
            source_version="regression-v1",
            license="synthetic fixture",
            redistribution="permitted",
        ),
    )


def _site_key(position: int) -> str:
    return (
        "phospy:v1|organism=rat|protein_namespace=protein_id|"
        f"protein_identifier=P{position:04d}|residue=S|position={position}"
    )


def _internal_request(
    *,
    method: InternalBatchCorrectionMethod,
    replicate_column: str | None,
    allow_temporary_completion: bool = False,
) -> InternalBatchCorrectionRequest:
    return InternalBatchCorrectionRequest(
        method=method,
        batch_column="batch",
        condition_columns=("condition",),
        replicate_column=replicate_column,
        control_site_source=InternalBatchCorrectionControlSiteSource.CALLER_SUPPLIED,
        control_site_mode=InternalBatchCorrectionControlSiteMode.SITE_KEY_LIST,
        missing_value_policy=(
            InternalBatchCorrectionMissingValuePolicy.ALLOW_TEMPORARY_IMPUTATION
            if allow_temporary_completion
            else InternalBatchCorrectionMissingValuePolicy.REJECT_MISSING
        ),
        imputation_policy=(
            InternalBatchCorrectionImputationPolicy.ROW_MEDIAN_TEMPORARY
            if allow_temporary_completion
            else InternalBatchCorrectionImputationPolicy.NONE
        ),
        n_unwanted_factors=1,
        stage_order=(
            InternalBatchCorrectionStageOrder.AFTER_MISSING_DATA_BEFORE_DOWNSTREAM
        ),
        diagnostics_enabled=True,
    )


def _complete_missingness_policy(
    matrix: pd.DataFrame,
) -> CorrectionMissingnessPolicy:
    return CorrectionMissingnessPolicy(
        originally_missing_cells_tracked_by=(
            OriginallyMissingCellTracking.OBSERVATION_MASK
        ),
        observation_mask=ObservationMask(
            feature_ids=tuple(str(value) for value in matrix.index),
            sample_ids=tuple(str(value) for value in matrix.columns),
            originally_missing_cells=(),
        ),
    )


def _run_workflow(
    request: BatchCorrectionWorkflowRequest,
) -> BatchCorrectionWorkflowResult:
    with patch(
        _CORRECTION_ENVIRONMENT_COLLECTOR,
        return_value=_stable_environment_provenance(),
    ):
        return BatchCorrectionWorkflow().run(request)


def _result_payload(result: BatchCorrectionWorkflowResult) -> dict[str, object]:
    corrected = cast(pd.DataFrame, result.corrected_matrix)
    applied = result.corrected_preprocessing_output
    output_mask = cast(pd.DataFrame, applied.output_observation_mask)
    cell_status = cast(pd.DataFrame, applied.corrected_cell_status)
    provenance = result.provenance
    direct_payload = batch_correction_provenance_to_payload(provenance)
    restored = batch_correction_provenance_from_payload(direct_payload)
    if batch_correction_provenance_to_payload(restored) != direct_payload:
        raise AssertionError("direct batch-correction provenance round trip changed")

    defensive = applied.corrected
    defensive.iloc[0, 0] = -999_999.0
    defensive_access_preserved = not bool(applied.corrected.iloc[0, 0] == -999_999.0)
    if not defensive_access_preserved:
        raise AssertionError("corrected output exposed its owned dataframe")
    return {
        "corrected_matrix": _frame_payload(corrected),
        "output_observation_mask": _frame_payload(output_mask),
        "corrected_cell_status": _frame_payload(cell_status),
        "report": applied.batch_correction_report.to_payload(),
        "diagnostics": dict(result.diagnostics),
        "warnings": list(result.warnings),
        "caveats": [item.to_payload() for item in result.caveats],
        "stage_order": list(applied.stage_order),
        "consumed_by_downstream": applied.consumed_by_downstream,
        "defensive_dataframe_access_preserved": defensive_access_preserved,
        "direct_provenance": direct_payload,
        "direct_provenance_bytes": _canonical_bytes(direct_payload).decode("utf-8"),
    }


def _sps_case(*, permuted: bool) -> dict[str, object]:
    samples = ("sample_1", "sample_2", "sample_3", "sample_4")
    matrix = pd.DataFrame(
        {
            "sample_1": (1.0, 2.0, 5.0, 8.0),
            "sample_2": (2.0, 4.0, 5.5, 7.0),
            "sample_3": (3.0, 1.0, 6.0, 9.0),
            "sample_4": (4.0, 3.0, 7.0, 6.0),
        },
        index=pd.Index(("site_a", "site_b", "site_c", "site_d"), name="site_key"),
    )
    metadata = pd.DataFrame(
        {
            "batch": ("run_1", "run_1", "run_2", "run_2"),
            "condition": ("control", "treated", "control", "treated"),
        },
        index=pd.Index(samples, name="sample_id"),
    )
    if permuted:
        order = ("sample_4", "sample_2", "sample_3", "sample_1")
        matrix = matrix.reindex(
            index=list(("site_d", "site_b", "site_a", "site_c")),
            columns=list(order),
        )
        metadata = metadata.reindex(index=list(order))
    before = _frame_payload(matrix)
    result = _run_workflow(
        BatchCorrectionWorkflowRequest(
            phospho=matrix,
            config=_internal_request(
                method=InternalBatchCorrectionMethod.SPS_RUV_STYLE,
                replicate_column=None,
            ),
            sample_metadata=metadata,
            control_site_set=_control_site_set(("site_a", "site_b")),
            missingness_policy=_complete_missingness_policy(matrix),
        )
    )
    if _frame_payload(matrix) != before:
        raise AssertionError("SPS/RUV workflow mutated its input matrix")
    payload = _result_payload(result)
    payload["input_non_mutation"] = True
    return payload


def _ruv_iii_case() -> dict[str, object]:
    samples = tuple(f"sample_{position}" for position in range(1, 7))
    sites = tuple(_site_key(position) for position in range(1, 6))
    unwanted = np.asarray((-2.0, 2.0, -1.0, 1.0, -3.0, 3.0))
    biological = np.asarray((-4.0, -4.0, 1.0, 1.0, 5.0, 5.0))
    matrix = pd.DataFrame(
        np.vstack(
            (
                15.0 + unwanted,
                12.0 - (2.0 * unwanted),
                18.0 + biological + (0.5 * unwanted),
                12.0 + biological,
                np.repeat(7.0, len(samples)),
            )
        ),
        index=pd.Index(sites, name="site_key"),
        columns=pd.Index(samples, name="sample_id"),
    )
    metadata = pd.DataFrame(
        {
            "batch": ("run_1", "run_2", "run_1", "run_2", "run_1", "run_2"),
            "condition": ("a", "a", "b", "b", "c", "c"),
            "replicate": ("r1", "r1", "r2", "r2", "r3", "r3"),
        },
        index=matrix.columns.copy(),
    )
    before = _frame_payload(matrix)
    result = _run_workflow(
        BatchCorrectionWorkflowRequest(
            phospho=matrix,
            config=_internal_request(
                method=InternalBatchCorrectionMethod.RUV_III_STYLE,
                replicate_column="replicate",
            ),
            sample_metadata=metadata,
            control_site_set=_control_site_set(sites[:2]),
            missingness_policy=_complete_missingness_policy(matrix),
        )
    )
    if _frame_payload(matrix) != before:
        raise AssertionError("RUV-III workflow mutated its input matrix")
    payload = _result_payload(result)
    payload["input_non_mutation"] = True
    return payload


@dataclass(frozen=True)
class _PreparationPlan:
    observation_mask: ObservationMask
    temporary_imputation_policy: TemporaryImputationPolicy


def _missingness_preparation_cases() -> dict[str, object]:
    matrix = pd.DataFrame(
        {
            "sample_1": (10.0, 1.0, 5.0),
            "sample_2": (12.0, float("nan"), 5.0),
            "sample_3": (14.0, 3.0, 5.0),
            "sample_4": (16.0, 5.0, 5.0),
        },
        index=pd.Index(("site_a", "site_b", "site_c"), name="site_key"),
    )
    policy = TemporaryImputationPolicy(
        allowed=True,
        method=TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY,
        method_parameters=(("min_observed_values", 2),),
    )
    actual_plan = _PreparationPlan(
        observation_mask=ObservationMask(
            feature_ids=tuple(str(value) for value in matrix.index),
            sample_ids=tuple(str(value) for value in matrix.columns),
            originally_missing_cells=(("site_b", "sample_2"),),
        ),
        temporary_imputation_policy=policy,
    )
    before = _frame_payload(matrix)
    prepared = _prepare_matrix(phospho=matrix, plan=actual_plan)
    restored = prepared.working.mask(prepared.actual_missing, np.nan)
    if _frame_payload(matrix) != before:
        raise AssertionError("missingness preparation mutated its input matrix")

    upstream_matrix = matrix.fillna(3.0)
    upstream_plan = _PreparationPlan(
        observation_mask=actual_plan.observation_mask,
        temporary_imputation_policy=policy,
    )
    upstream = _prepare_matrix(phospho=upstream_matrix, plan=upstream_plan)
    return {
        "actual_missing": {
            "input": before,
            "temporary_completion": _frame_payload(prepared.working),
            "actual_missing": _frame_payload(prepared.actual_missing),
            "governed_missing": _frame_payload(prepared.originally_missing),
            "restored": _frame_payload(restored),
            "missing_cells": [list(cell) for cell in prepared.missing_cells],
            "temporary_completion_applied": prepared.temporary_completion_applied,
            "input_non_mutation": True,
        },
        "upstream_imputed": {
            "working": _frame_payload(upstream.working),
            "actual_missing": _frame_payload(upstream.actual_missing),
            "governed_missing": _frame_payload(upstream.originally_missing),
            "temporary_completion_applied": upstream.temporary_completion_applied,
        },
    }


def _linear_case() -> dict[str, object]:
    matrix = pd.DataFrame(
        {
            "sample_1": (10.0, 8.0, 12.0),
            "sample_2": (14.0, 9.0, 15.0),
            "sample_3": (12.0, 7.0, 14.0),
            "sample_4": (16.0, 8.0, 17.0),
            "sample_5": (11.0, 9.0, 13.0),
            "sample_6": (15.0, 10.0, 16.0),
        },
        index=pd.Index(("site_a", "site_b", "site_c"), name="site_key"),
    )
    before = _frame_payload(matrix)
    result = LinearResidualizeBatchCorrectionEngine().run(
        phospho=matrix,
        batch_labels=("run_1", "run_1", "run_2", "run_2", "run_3", "run_3"),
        condition_labels=("control", "treated") * 3,
        config=DatasetBatchCorrectionConfig(
            method=DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
            batch_column="batch",
            condition_column="condition",
        ),
    )
    if _frame_payload(matrix) != before:
        raise AssertionError("linear residualisation mutated its input matrix")
    corrected = result.corrected_matrix
    exported = result.corrected_matrix
    exported.iloc[0, 0] = -999_999.0
    defensive_access = not bool(result.corrected_matrix.iloc[0, 0] == -999_999.0)
    if not defensive_access:
        raise AssertionError("linear result exposed its owned dataframe")
    return {
        "input": before,
        "corrected": _frame_payload(corrected),
        "report": result.report.to_payload(),
        "diagnostics": dict(result.diagnostics),
        "input_non_mutation": True,
        "defensive_dataframe_access_preserved": defensive_access,
    }


def _disabled_report_case() -> dict[str, object]:
    report = BatchCorrectionReport(
        status=BATCH_CORRECTION_STATUS_DISABLED,
        policy=BatchCorrectionPolicy(
            method="none",
            batch_column="batch",
            condition_column="condition",
        ),
        diagnostics=BatchCorrectionDiagnostics(
            confounding_check_status=BATCH_CORRECTION_CONFOUNDING_NOT_APPLICABLE,
            matrix_shape_before=(3, 6),
            matrix_shape_after=(3, 6),
            limitations=("batch correction disabled by preprocessing configuration",),
        ),
    )
    return report.to_payload()


def _error_payload(callable_: object) -> dict[str, str]:
    try:
        callable_()  # type: ignore[operator]
    except PhosPyInputError as exc:
        return {"type": type(exc).__name__, "message": str(exc)}
    raise AssertionError("expected PhosPyInputError")


def _supported_errors() -> dict[str, object]:
    matrix = pd.DataFrame(
        {"sample_1": (1.0, 2.0), "sample_2": (float("nan"), 3.0)},
        index=pd.Index(("site_a", "site_b"), name="site_key"),
    )
    plan = _PreparationPlan(
        observation_mask=ObservationMask(
            feature_ids=("site_a", "site_b"),
            sample_ids=("sample_1", "sample_2"),
            originally_missing_cells=(),
        ),
        temporary_imputation_policy=TemporaryImputationPolicy(),
    )
    malformed = _sps_case(permuted=False)["direct_provenance"]
    assert isinstance(malformed, dict)
    malformed = dict(malformed)
    malformed["schema_version"] = 999
    return {
        "ungoverned_missing_cell": _error_payload(
            lambda: _prepare_matrix(phospho=matrix, plan=plan)
        ),
        "malformed_historical_payload": _error_payload(
            lambda: batch_correction_provenance_from_payload(malformed)
        ),
    }


def _nested_provenance_case(
    direct_payload: dict[str, object],
) -> dict[str, object]:
    correction = batch_correction_provenance_from_payload(direct_payload)
    input_fingerprint = correction.input_matrix_fingerprint
    output_fingerprint = correction.output_matrix_fingerprint
    if output_fingerprint is None:
        raise AssertionError("comparison correction omitted output fingerprint")
    stage = PreprocessingStageProvenance(
        stage="batch_correction",
        operation=correction.requested_method,
        parameters={"method": correction.requested_method},
        input_shape=(input_fingerprint.rows, input_fingerprint.columns),
        output_shape=(output_fingerprint.rows, output_fingerprint.columns),
        input_hash=input_fingerprint.exact_hash_value,
        output_hash=output_fingerprint.exact_hash_value,
        phospho_input_hash=input_fingerprint.exact_hash_value,
        phospho_output_hash=output_fingerprint.exact_hash_value,
        dropped_row_ids=(),
        dropped_row_count=0,
        batch_correction_provenance=correction,
    )
    run = RunProvenance(
        environment=_stable_environment_provenance(),
        input_tables=(input_fingerprint,),
        preprocessing_stages=(stage,),
        reference=None,
        workflow_name="batch_correction_regression",
        workflow_parameters={"method": correction.requested_method},
        random_state=None,
        random_seed_policy=None,
        output_tables=(output_fingerprint,),
    )
    payload = to_payload(run)
    restored = from_payload(payload)
    restored_payload = to_payload(restored)
    if restored_payload != payload:
        raise AssertionError("nested RunProvenance round trip changed")
    return {
        "payload": payload,
        "payload_bytes": _canonical_bytes(payload).decode("utf-8"),
        "round_trip": restored_payload,
        "runtime_identity_preserved": (
            restored.preprocessing_stages[0].batch_correction_provenance.__class__
            is correction.__class__
        ),
    }


def _stable_environment_provenance() -> EnvironmentProvenance:
    """Return the controlled collaborator used by the regression fixture."""

    return EnvironmentProvenance(
        package_name="phospy",
        package_version="1.7.6",
        python_version="controlled-regression-runtime",
        dependency_versions={
            "numpy": "controlled",
            "pandas": "controlled",
            "scipy": "controlled",
            "scikit-learn": "controlled",
            "pyarrow": None,
            "openpyxl": None,
        },
        platform={
            "platform": "controlled",
            "system": "controlled",
            "release": "controlled",
            "version": "controlled",
            "machine": "controlled",
            "processor": "controlled",
            "python_implementation": "CPython",
        },
        blas_lapack={"source": "controlled"},
        thread_environment={
            "OMP_NUM_THREADS": None,
            "OPENBLAS_NUM_THREADS": None,
            "MKL_NUM_THREADS": None,
            "NUMEXPR_NUM_THREADS": None,
        },
        timezone="controlled",
        locale={
            "language_code": "controlled",
            "encoding": "UTF-8",
            "lc_all": "controlled",
            "preferred_encoding": "UTF-8",
        },
        constraints_fingerprint={
            "algorithm": "controlled-fixture-v1",
            "value": "batch-correction-regression",
            "sources": "constraints/ci.txt,constraints/minimum.txt",
        },
    )


def build_comparison_payload() -> dict[str, object]:
    sps = _sps_case(permuted=False)
    direct = cast(dict[str, object], sps["direct_provenance"])
    cases: dict[str, object] = {
        "disabled_report": _disabled_report_case(),
        "linear_residualisation": _linear_case(),
        "sps_ruv_style": sps,
        "sps_ruv_style_permuted": _sps_case(permuted=True),
        "ruv_iii_style": _ruv_iii_case(),
        "missingness_preparation": _missingness_preparation_cases(),
        "nested_provenance": _nested_provenance_case(direct),
        "supported_errors": _supported_errors(),
    }
    invariant_cases = _normalise_numerical_frames(cases)
    canonical = _canonical_bytes(invariant_cases)
    return {
        "schema": "phospy-batch-correction-regression-v1",
        "cases": cases,
        "cases_sha256": hashlib.sha256(canonical).hexdigest(),
    }


def main() -> None:
    print(_canonical_bytes(build_comparison_payload()).decode("utf-8"))


if __name__ == "__main__":
    main()
