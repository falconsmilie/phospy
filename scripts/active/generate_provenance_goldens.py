#!/usr/bin/env python3
"""Regenerate integration provenance golden fixtures.

This script refreshes:
  - kinase_public_predmat_provenance_golden.json
  - signalome_l6_provenance_golden.json

Run it in the same environment used by CI to avoid cross-environment drift.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

PUBLIC_WORKFLOW_REFERENCE = ROOT / "tests" / "fixtures" / "public_workflow_reference"
KINASE_GOLDEN_PATH = (
    PUBLIC_WORKFLOW_REFERENCE / "kinase_public_predmat_provenance_golden.json"
)
SIGNALOME_GOLDEN_PATH = (
    PUBLIC_WORKFLOW_REFERENCE / "signalome_l6_provenance_golden.json"
)


def _fingerprints_by_name(
    fingerprints: tuple[object, ...],
) -> dict[str, dict[str, Any]]:
    return {
        str(item.name): {
            "rows": int(item.rows),
            "columns": int(item.columns),
            "exact_hash_algorithm": str(item.exact_hash_algorithm),
            "exact_hash_value": str(item.exact_hash_value),
            "tolerance_hash_algorithm": str(item.tolerance_hash_algorithm),
            "tolerance_hash_value": str(item.tolerance_hash_value),
        }
        for item in fingerprints
    }


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _refresh_kinase_golden() -> None:
    from phospy import KinaseWorkflow
    from phospy.api import (
        KinasePredictionConfig,
        KinaseScoringConfig,
        KinaseWorkflowRequest,
        ReferenceContextCompatibilityPolicy,
    )
    from tests.support.public_predmat_parity_metrics import (
        _build_public_predmat_dataset,
        _build_public_predmat_references,
    )

    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_build_public_predmat_dataset(),
            references=_build_public_predmat_references(reverse_order=False),
            scoring_config=KinaseScoringConfig(
                reliability_profile="custom",
                min_substrates=2,
                reference_context_compatibility_policy=(
                    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                ),
            ),
            prediction_config=KinasePredictionConfig(
                top_k=4,
                deterministic_max_selected_kinases=3,
                adaptive_ensemble_runs=3,
                mode="adaptive_ensemble",
                adaptive_policy="stable",
                n_iterations=2,
                random_state=17,
            ),
            activity_config=None,
        )
    )
    provenance = result.provenance
    if provenance is None:
        raise RuntimeError("kinase provenance missing; cannot refresh golden fixture")

    golden = _load_json(KINASE_GOLDEN_PATH)
    golden["input_tables"] = _fingerprints_by_name(provenance.input_tables)
    golden["output_tables"] = _fingerprints_by_name(provenance.output_tables)
    if provenance.reference is None:
        raise RuntimeError("kinase reference provenance missing")
    golden["reference"]["table_fingerprints"] = _fingerprints_by_name(
        provenance.reference.table_fingerprints
    )
    _write_json(KINASE_GOLDEN_PATH, golden)


def _refresh_signalome_golden() -> None:
    from phospy import KinaseWorkflow, SignalomeWorkflow
    from phospy.api import (
        KinasePredictionConfig,
        KinaseScoringConfig,
        KinaseWorkflowRequest,
        ReferenceContextCompatibilityPolicy,
        ReferencePreset,
        SignalomeWorkflowRequest,
    )
    from tests.support.rewrite_fixture_data import build_rat_l6_dataset
    from tests.support.signalome_config import build_signalome_config

    dataset = build_rat_l6_dataset(n_sites=260)
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(
                reliability_profile="custom",
                min_substrates=2,
                reference_context_compatibility_policy=(
                    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                ),
            ),
            prediction_config=KinasePredictionConfig(
                top_k=6,
                deterministic_max_selected_kinases=12,
                adaptive_ensemble_runs=12,
            ),
            activity_config=None,
        )
    )
    result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=build_signalome_config(substrate_support_cutoff=0.5),
        )
    )
    provenance = result.provenance
    if provenance is None:
        raise RuntimeError(
            "signalome provenance missing; cannot refresh golden fixture"
        )

    golden = _load_json(SIGNALOME_GOLDEN_PATH)
    golden["input_tables"] = _fingerprints_by_name(provenance.input_tables)
    golden["output_tables"] = _fingerprints_by_name(provenance.output_tables)
    if provenance.reference is None:
        raise RuntimeError("signalome reference provenance missing")
    golden["reference"]["table_fingerprints"] = _fingerprints_by_name(
        provenance.reference.table_fingerprints
    )
    _write_json(SIGNALOME_GOLDEN_PATH, golden)


def main() -> int:
    _refresh_kinase_golden()
    _refresh_signalome_golden()
    print("Updated provenance golden fixtures:")
    print(f"  - {KINASE_GOLDEN_PATH}")
    print(f"  - {SIGNALOME_GOLDEN_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
