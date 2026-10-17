from __future__ import annotations

import hashlib
import json
import sys
from collections.abc import Mapping
from dataclasses import fields, is_dataclass
from enum import Enum
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd
import pytest

from phospy import SignalomeWorkflow
from phospy.api import (
    AnalysisReadyPhosphoDataset,
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
    Organism,
    ReferenceBundle,
    SignalomeWorkflowRequest,
)
from phospy.run_lineage.serialization import to_payload as provenance_to_payload
from phospy.signalome.contracts import ResolvedSignalomeWorkflowRequest
from phospy.signalome.executor import SignalomeWorkflowExecutor
from phospy.signalome.results import SignalomeWorkflowResult
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.signalome_config import build_signalome_config
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)

pytestmark = pytest.mark.integration

_BASELINE_SHA = "177558170b205a63da4881192878caf3113568b3"
_GOLDEN_PATH = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "signalome"
    / "score_source_regression_baseline.json"
)
_SOURCES = (
    "profile_scores",
    "rank_weighted_fusion_scores",
    "kinase_library_motif_scores",
    "combined_profile_motif_scores",
)
_SUPPORTED_PYTHON_MINORS = frozenset({(3, 11), (3, 12)})
requires_supported_python = pytest.mark.skipif(
    sys.version_info[:2] not in _SUPPORTED_PYTHON_MINORS,
    reason=(
        "the exact Signalome regression baseline is defined for supported "
        "Python 3.11 and 3.12 runtimes; running "
        f"{sys.version_info.major}.{sys.version_info.minor}"
    ),
)


def _fixture_inputs() -> tuple[
    AnalysisReadyPhosphoDataset,
    ReferenceBundle,
    pd.DataFrame,
    dict[str, pd.DataFrame],
]:
    display_ids = [f"P{row // 2 + 1};S{row + 1};" for row in range(12)]
    site_index = site_key_index_from_display_ids(display_ids)
    protein_ids = [f"PROT{row // 2 + 1}" for row in range(12)]
    phospho = pd.DataFrame(
        {
            "sample_a": np.linspace(1.0, 2.1, num=12),
            "sample_b": np.linspace(1.5, 2.6, num=12),
        },
        index=site_index.copy(),
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": display_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": [f"P{row // 2 + 1}" for row in range(12)],
            "site": [f"S{row + 1}" for row in range(12)],
            "site_sequence": ["A" * 15 + "S" + "A" * 15] * 12,
            "protein_id": protein_ids,
        },
        index=site_index.copy(),
    )
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1", "K2", "K3"],
                "substrate_site": display_ids[:3],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 15 + "S" + "A" * 15] * 12},
            index=pd.Index(display_ids, name="site_id"),
        ),
    )
    columns = pd.Index(["K1", "K2", "K3"], name="kinase")
    prediction = pd.DataFrame(
        [
            [0.90, 0.15, 0.05],
            [0.85, 0.20, 0.10],
            [0.75, 0.30, 0.15],
            [0.65, 0.40, 0.20],
            [0.55, 0.50, 0.25],
            [0.45, 0.60, 0.30],
            [0.35, 0.70, 0.35],
            [0.25, 0.80, 0.40],
            [0.15, 0.90, 0.45],
            [0.10, 0.75, 0.60],
            [0.20, 0.55, 0.80],
            [0.30, 0.35, 0.95],
        ],
        index=site_index.copy(),
        columns=columns.copy(),
    )
    lane_values = {
        "profile_scores": np.asarray(
            [
                [0.95, 0.10, 0.05],
                [0.90, 0.15, 0.10],
                [0.80, 0.25, 0.15],
                [0.70, 0.35, 0.20],
                [0.60, 0.45, 0.25],
                [0.50, 0.55, 0.30],
                [0.40, 0.65, 0.35],
                [0.30, 0.75, 0.40],
                [0.20, 0.85, 0.45],
                [0.15, 0.70, 0.65],
                [0.25, 0.50, 0.85],
                [0.35, 0.30, 0.90],
            ]
        ),
        "rank_weighted_fusion_scores": np.asarray(
            [
                [0.15, 0.90, 0.05],
                [0.20, 0.85, 0.10],
                [0.30, 0.75, 0.15],
                [0.40, 0.65, 0.20],
                [0.50, 0.55, 0.25],
                [0.60, 0.45, 0.30],
                [0.70, 0.35, 0.35],
                [0.80, 0.25, 0.40],
                [0.90, 0.15, 0.45],
                [0.75, 0.10, 0.60],
                [0.55, 0.20, 0.80],
                [0.35, 0.30, 0.95],
            ]
        ),
        "kinase_library_motif_scores": np.asarray(
            [
                [0.05, 0.15, 0.95],
                [0.10, 0.20, 0.90],
                [0.15, 0.30, 0.80],
                [0.20, 0.40, 0.70],
                [0.25, 0.50, 0.60],
                [0.30, 0.60, 0.50],
                [0.35, 0.70, 0.40],
                [0.40, 0.80, 0.30],
                [0.45, 0.90, 0.20],
                [0.60, 0.75, 0.15],
                [0.80, 0.55, 0.25],
                [0.95, 0.35, 0.30],
            ]
        ),
        "combined_profile_motif_scores": np.asarray(
            [
                [0.90, 0.05, 0.15],
                [0.80, 0.10, 0.25],
                [0.70, 0.15, 0.35],
                [0.60, 0.20, 0.45],
                [0.50, 0.25, 0.55],
                [0.40, 0.30, 0.65],
                [0.30, 0.35, 0.75],
                [0.20, 0.40, 0.85],
                [0.10, 0.45, 0.95],
                [0.25, 0.60, 0.75],
                [0.45, 0.80, 0.55],
                [0.65, 0.90, 0.35],
            ]
        ),
    }
    lanes: dict[str, pd.DataFrame] = {}
    for source, values in lane_values.items():
        values = values.copy()
        values[4, :] = np.nan
        lanes[source] = pd.DataFrame(
            values,
            index=site_index.copy(),
            columns=columns.copy(),
        )
    return dataset, references, prediction, lanes


def _kinase_result_for_source(source: str) -> KinaseWorkflowResult:
    dataset, references, prediction, lanes = _fixture_inputs()
    optional_lane: dict[str, pd.DataFrame] = {}
    if source != "profile_scores":
        optional_lane[source] = lanes[source]
    return KinaseWorkflowResult(
        dataset=dataset,
        references=references,
        scoring_result=KinaseScoringResult(
            profile_scores=lanes["profile_scores"],
            score_source=source,
            **cast(Any, optional_lane),
        ),
        prediction_result=KinasePredictionResult(pred_mat=prediction),
        activity_result=None,
    )


def _normalize(value: Any) -> Any:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, Enum):
        return _normalize(value.value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        if pd.isna(value):
            return None
        return round(float(value), 10)
    if isinstance(value, Mapping):
        return {str(key): _normalize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_normalize(item) for item in value]
    if is_dataclass(value):
        return {
            field.name: _normalize(getattr(value, field.name))
            for field in fields(value)
        }
    if pd.isna(value):
        return None
    return str(value)


def _frame_evidence(frame: pd.DataFrame | None) -> dict[str, object] | None:
    if frame is None:
        return None
    payload = {
        "index_name": _normalize(frame.index.name),
        "index": _normalize(frame.index.tolist()),
        "columns_name": _normalize(frame.columns.name),
        "columns": _normalize(frame.columns.tolist()),
        "dtypes": [str(dtype) for dtype in frame.dtypes.tolist()],
        "values": _normalize(frame.to_numpy(dtype=object).tolist()),
    }
    encoded = json.dumps(
        payload,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return {
        "shape": [int(frame.shape[0]), int(frame.shape[1])],
        "missing_values": int(frame.isna().to_numpy().sum()),
        "index": payload["index"],
        "columns": payload["columns"],
        "sha256": hashlib.sha256(encoded).hexdigest(),
    }


def _payload_digest(value: object) -> str:
    encoded = json.dumps(
        _normalize(value),
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _result_evidence(
    *,
    source: str,
    selected_matrix: pd.DataFrame,
    result: SignalomeWorkflowResult,
) -> dict[str, object]:
    provenance = result.provenance
    assert provenance is not None
    provenance_payload = provenance_to_payload(provenance)
    provenance_payload.pop("environment", None)
    workflow_parameters = provenance_payload["workflow_parameters"]
    assert isinstance(workflow_parameters, dict)
    scientific_policies = cast(
        list[dict[str, object]], provenance_payload["scientific_policies"]
    )
    output_tables = cast(list[dict[str, object]], provenance_payload["output_tables"])
    return {
        "selected_source": source,
        "selected_matrix": _frame_evidence(selected_matrix),
        "score_preconditioning_diagnostics": _normalize(
            result.score_preconditioning_diagnostics
        ),
        "alignment_diagnostics": _normalize(result.alignment_diagnostics),
        "row_attrition": _normalize(workflow_parameters["row_attrition"]),
        "tables": {
            "module_assignments": _frame_evidence(result.module_assignments.table),
            "signalome_modules": _frame_evidence(result.signalome_modules.table),
            "network_edges": _frame_evidence(result.kinase_network.edges),
            "network_nodes": _frame_evidence(result.kinase_network.nodes),
            "candidate_correlations": _frame_evidence(
                result.kinase_network.candidate_correlations
            ),
            "expanded_signalome": _frame_evidence(result.expanded_signalome),
            "site_membership": _frame_evidence(result.site_membership),
            "protein_site_context": _frame_evidence(result.protein_site_context),
        },
        "caveats": [_normalize(caveat.to_payload()) for caveat in result.caveats],
        "provenance": {
            "workflow_name": provenance_payload["workflow_name"],
            "signalome_score_semantics": _normalize(
                workflow_parameters["signalome_score_semantics"]
            ),
            "row_attrition": _normalize(workflow_parameters["row_attrition"]),
            "score_preconditioning_diagnostics": _normalize(
                workflow_parameters["score_preconditioning_diagnostics"]
            ),
            "alignment_diagnostics": _normalize(
                workflow_parameters["alignment_diagnostics"]
            ),
            "scientific_policies": _normalize(
                {
                    "count": len(scientific_policies),
                    "ids": [policy["id"] for policy in scientific_policies],
                    "sha256": _payload_digest(scientific_policies),
                }
            ),
            "output_tables": _normalize(
                {
                    "count": len(output_tables),
                    "names": [table["name"] for table in output_tables],
                    "sha256": _payload_digest(output_tables),
                }
            ),
        },
    }


def build_four_source_evidence() -> dict[str, object]:
    captured: list[tuple[str, pd.DataFrame]] = []
    original_run = SignalomeWorkflowExecutor.run

    def _capture_and_run(
        executor: SignalomeWorkflowExecutor,
        request: ResolvedSignalomeWorkflowRequest,
    ) -> SignalomeWorkflowResult:
        captured.append(
            (
                str(request.downstream_score_source),
                request.downstream_score_matrix.copy(deep=True),
            )
        )
        return original_run(executor, request)

    SignalomeWorkflowExecutor.run = _capture_and_run  # type: ignore[method-assign]
    try:
        cases: dict[str, object] = {}
        for source in _SOURCES:
            result = SignalomeWorkflow().run(
                SignalomeWorkflowRequest(
                    kinase_result=_kinase_result_for_source(source),
                    config=build_signalome_config(
                        substrate_support_cutoff=0.5,
                        score_preconditioning_policy="allow_and_report",
                    ),
                )
            )
            observed_source, selected_matrix = captured[-1]
            cases[source] = _result_evidence(
                source=observed_source,
                selected_matrix=selected_matrix,
                result=result,
            )
    finally:
        SignalomeWorkflowExecutor.run = original_run  # type: ignore[method-assign]
    return {"baseline_sha": _BASELINE_SHA, "cases": cases}


@requires_supported_python
def test_all_four_authoritative_sources_match_signalome_regression_baseline() -> None:
    expected = json.loads(_GOLDEN_PATH.read_text(encoding="utf-8"))
    assert expected["baseline_sha"] == _BASELINE_SHA
    assert set(expected["cases"]) == set(_SOURCES)

    observed = build_four_source_evidence()

    assert observed == expected
