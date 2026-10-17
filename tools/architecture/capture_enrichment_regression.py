"""Emit deterministic enrichment regression payloads through the public API."""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import fields
from enum import Enum
from typing import Any
from unittest.mock import patch

import pandas as pd

from phospy.api import (
    EnrichmentConfig,
    EnrichmentIdentifierSetProvenance,
    EnrichmentIdentifierSetSourceType,
    EnrichmentSetCollection,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    EnrichmentWorkflowResult,
    GeneSetCollection,
    PtmSetCollection,
    WorkflowValidationError,
)
from phospy.provenance.models.environment import EnvironmentProvenance
from phospy.run_lineage.serialization import to_payload as run_to_payload

GENE = "gene_symbol"
SITE = "site_key"


def _fixed_environment() -> EnvironmentProvenance:
    return EnvironmentProvenance(
        package_name="phospy",
        package_version="1.7.6",
        python_version="3.11.15",
        dependency_versions={
            "numpy": "2.3.5",
            "pandas": "2.3.3",
            "scipy": "1.16.3",
        },
        platform={
            "implementation": "CPython",
            "system": "Windows",
            "machine": "AMD64",
        },
        blas_lapack={"provider": "controlled-enrichment-regression"},
        thread_environment={
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
        },
        timezone="UTC",
        locale={"preferred_encoding": "UTF-8"},
        constraints_fingerprint={
            "algorithm": "sha256",
            "value": "enrichment-regression-constraints",
        },
    )


def _provenance_patch_target() -> str:
    return "phospy.enrichment.provenance.collect_environment_provenance"


def _identifier_provenance(
    *, source_label: str, identifier_count: int
) -> EnrichmentIdentifierSetProvenance:
    return EnrichmentIdentifierSetProvenance(
        source_type=EnrichmentIdentifierSetSourceType.RAW_IDENTIFIER_LIST,
        source_label=source_label,
        identifier_count=identifier_count,
        upstream_workflow_id="enrichment-regression-control",
        upstream_result_id=f"{source_label}-v1",
    )


def _gene_request() -> EnrichmentWorkflowRequest:
    return EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=GENE,
        set_collection=GeneSetCollection(
            sets={
                "SIGNALING": ("AKT1", "MAPK1", "MTOR"),
                "CELL_CYCLE": ("CDK1", "CDK2", "MAPK1"),
                "NO_BACKGROUND_MEMBERS": ("OUTSIDE_A", "OUTSIDE_B"),
                "TOO_LARGE": ("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
            },
            identifier_kind=GENE,
            term_names={
                "SIGNALING": "Signaling",
                "CELL_CYCLE": "Cell cycle",
                "NO_BACKGROUND_MEMBERS": "No background members",
                "TOO_LARGE": "Too large",
            },
            source_name="synthetic_gene_sets",
            source_version="regression-v1",
        ),
        selected_identifiers=("AKT1", "UNKNOWN", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=EnrichmentConfig(
            min_set_size=1,
            max_set_size=3,
            selected_outside_background_policy="drop",
        ),
        selected_identifier_provenance=_identifier_provenance(
            source_label="gene-selected", identifier_count=3
        ),
        background_identifier_provenance=_identifier_provenance(
            source_label="gene-background", identifier_count=5
        ),
    )


def _ptm_request() -> EnrichmentWorkflowRequest:
    return EnrichmentWorkflowRequest(
        identifier_column="site_key",
        identifier_kind=SITE,
        set_collection=PtmSetCollection(
            sets={
                "MOTIF_A": ("rat|P12345|S10", "rat|P12345|T20"),
                "MOTIF_B": ("rat|P12345|Y30", "rat|P99999|S5"),
                "MOTIF_EMPTY": ("rat|OUTSIDE|S1",),
            },
            identifier_kind=SITE,
            source_name="synthetic_ptm_sets",
            source_version="regression-v1",
        ),
        input_table=pd.DataFrame({"site_key": ["rat|P12345|S10", " rat|P12345|T20 "]}),
        background_universe=(
            "rat|P12345|S10",
            "rat|P12345|T20",
            "rat|P12345|Y30",
            "rat|P99999|S5",
        ),
        selected_identifier_provenance=_identifier_provenance(
            source_label="ptm-selected", identifier_count=2
        ),
        background_identifier_provenance=_identifier_provenance(
            source_label="ptm-background", identifier_count=4
        ),
    )


def _result_payload(result: EnrichmentWorkflowResult) -> dict[str, object]:
    provenance = result.provenance
    if provenance is None:
        raise AssertionError("comparison result omitted RunProvenance")
    table = result.table
    payload = {
        "identifier_kind": result.identifier_kind,
        "set_collection": _collection_payload(result.set_collection),
        "config": {
            field.name: getattr(result.config, field.name)
            for field in fields(result.config)
        },
        "records": [
            {field.name: getattr(record, field.name) for field in fields(record)}
            for record in result.records
        ],
        "table": {
            "columns": list(table.columns),
            "records": table.to_dict(orient="records"),
        },
        "unmatched_identifiers": result.unmatched_identifiers,
        "warnings": result.warnings,
        "caveats": [caveat.to_payload() for caveat in result.caveats],
        "diagnostics": result.diagnostics,
        "method_metadata": result.method_metadata,
        "background_summary": result.background_summary,
        "set_collection_summary": result.set_collection_summary,
        "selected_identifier_provenance": (
            None
            if result.selected_identifier_provenance is None
            else result.selected_identifier_provenance.to_payload()
        ),
        "background_identifier_provenance": (
            None
            if result.background_identifier_provenance is None
            else result.background_identifier_provenance.to_payload()
        ),
        "run_provenance": run_to_payload(provenance),
    }
    return _json_value(payload)


def _collection_payload(collection: EnrichmentSetCollection) -> dict[str, object]:
    return {
        "identifier_kind": collection.identifier_kind,
        "collection_kind": collection.collection_kind,
        "source_name": collection.source_name,
        "source_version": collection.source_version,
        "sets": [
            {
                field.name: getattr(enrichment_set, field.name)
                for field in fields(enrichment_set)
            }
            for enrichment_set in collection.enrichment_sets
        ],
    }


def _json_value(value: object) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, str | bytes | bytearray):
        return [_json_value(item) for item in value]
    if isinstance(value, Enum):
        return _json_value(value.value)
    if isinstance(value, float) and math.isnan(value):
        return {"special_float": "nan"}
    if value is None or isinstance(value, str | int | float | bool):
        return value
    raise TypeError(f"unsupported comparison payload value: {type(value).__name__}")


def _payload_bytes(payload: object) -> bytes:
    return json.dumps(
        payload,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
    ).encode("utf-8")


def build_comparison_payload() -> dict[str, object]:
    """Return full deterministic Gene/PTM and exact failure comparisons."""

    with patch(_provenance_patch_target(), return_value=_fixed_environment()):
        workflow = EnrichmentWorkflow()
        cases = {
            "gene": _result_payload(workflow.run(_gene_request())),
            "ptm": _result_payload(workflow.run(_ptm_request())),
        }
        failure_request = EnrichmentWorkflowRequest(
            identifier_column="gene_symbol",
            identifier_kind=GENE,
            set_collection=_gene_request().set_collection,
            selected_identifiers=("UNKNOWN_A", "UNKNOWN_B"),
            background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
            config=EnrichmentConfig(selected_outside_background_policy="drop"),
        )
        try:
            workflow.run(failure_request)
        except WorkflowValidationError as exc:
            failure = {"type": type(exc).__name__, "message": str(exc)}
        else:
            raise AssertionError("expected empty-foreground comparison failure")

    case_sha256 = {
        name: hashlib.sha256(_payload_bytes(payload)).hexdigest()
        for name, payload in cases.items()
    }
    return {
        "schema": "phospy-enrichment-regression-v1",
        "cases": cases,
        "case_sha256": case_sha256,
        "failure": failure,
    }


def main() -> None:
    print(_payload_bytes(build_comparison_payload()).decode("utf-8"))


if __name__ == "__main__":
    main()
