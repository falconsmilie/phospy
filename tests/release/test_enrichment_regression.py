from __future__ import annotations

import hashlib
from typing import Any, cast

import pytest

from tools.architecture.capture_enrichment_regression import (
    _payload_bytes,
    build_comparison_payload,
)

pytestmark = pytest.mark.reproducibility


def test_enrichment_matches_regression_payload_bytes() -> None:
    payload = build_comparison_payload()

    assert hashlib.sha256(_payload_bytes(payload)).hexdigest() == (
        "431ad63879d4199e460a93fc43a936f5a91d7cd7d1f7657f5647d6c53fb98c11"
    )
    assert payload["case_sha256"] == {
        "gene": "5ca26a4328008886501bbdc905ced2e5045b1044f80240c4c5f41aa1543f6851",
        "ptm": "e641ae074fffe0e2492387349ad8a3c5d8f01044d79027508e585024d01ffc2c",
    }
    cases = cast(dict[str, dict[str, Any]], payload["cases"])
    assert {
        name: (
            case["run_provenance"]["output_tables"][0]["exact_hash_value"],
            case["run_provenance"]["output_tables"][0]["tolerance_hash_value"],
        )
        for name, case in cases.items()
    } == {
        "gene": (
            "515cc654bff57bb72d133c2a27a9f9d891e08ba813712994854aab1237c0508a",
            "c38127cbe38985d3c506ebc090e1c83f28fa51edff48284dc91aa5a40dcc2e67",
        ),
        "ptm": (
            "7b512e7a2ce44e7fad93024cf4e6ab6d84f642f4a7697d37f2a11238b5414ee5",
            "f044568e67561e5416e4f0db8297cf076141fe6d5f0cd250ef8320921e686bde",
        ),
    }
    assert payload["failure"] == {
        "type": "WorkflowValidationError",
        "message": (
            "enrichment workflow request selected_identifiers are empty after "
            "background filtering; "
            "config.selected_outside_background_policy='drop'; "
            "selected_identifier_count=2; "
            "foreground_identifiers_missing_from_background="
            "('UNKNOWN_A', 'UNKNOWN_B')."
        ),
    }


def test_enrichment_regression_locks_correction_and_caveat_order() -> None:
    cases = cast(dict[str, dict[str, Any]], build_comparison_payload()["cases"])

    for case in cases.values():
        provenance_parameters = case["run_provenance"]["workflow_parameters"]
        assert provenance_parameters["correction_owner"] == "ora_engine"
        assert provenance_parameters["number_of_tests"] == len(case["records"])
        assert [caveat["code"] for caveat in case["caveats"]] == [
            "enrichment_offline_ora_only_scope",
            "enrichment_background_universe_assumption",
            "enrichment_identifier_kind_assumption",
            "enrichment_no_rank_based_semantics",
        ]
        assert [caveat["severity"] for caveat in case["caveats"]] == [
            "info",
            "warning",
            "info",
            "info",
        ]
