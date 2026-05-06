#!/usr/bin/env python3
"""Regenerate integration provenance golden fixtures.

This script refreshes:
  - kinase_public_predmat_provenance_golden.json
  - signalome_l6_provenance_golden.json

Run it in the same environment used by CI to avoid cross-environment drift.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Any

import pandas as pd

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
_PUBLIC_SITE_ID_PATTERN = re.compile(r"^\s*[^;]+\s*;\s*[^;]+\s*;\s*$")


def _canonical_public_site_components(site_id: object) -> tuple[str, str, str]:
    raw_site = str(site_id).strip()
    if _PUBLIC_SITE_ID_PATTERN.fullmatch(raw_site):
        parts = raw_site.split(";")
        gene_symbol = parts[0].strip()
        site = parts[1].strip()
        return f"{gene_symbol};{site};", gene_symbol, site

    gene_symbol = raw_site.split("_", 1)[0].strip()
    site = raw_site
    return f"{gene_symbol};{site};", gene_symbol, site


def _fingerprints_by_name(
    fingerprints: tuple[object, ...],
) -> dict[str, dict[str, Any]]:
    return {
        str(item.name): {
            "rows": int(item.rows),
            "columns": int(item.columns),
            "hash_algorithm": str(item.hash_algorithm),
            "hash_value": str(item.hash_value),
        }
        for item in fingerprints
    }


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _refresh_kinase_golden() -> None:
    from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
    from phospy.api import (
        DatasetBuildRequest,
        KinasePredictionConfig,
        KinaseScoringConfig,
        KinaseWorkflowRequest,
        Organism,
        ReferenceBundle,
    )
    from tests.support.rewrite_fixture_data import (
        load_public_predmat_input_phospho,
        load_public_predmat_input_site_sequences,
        load_public_predmat_input_substrate_map,
    )

    input_phospho = load_public_predmat_input_phospho()
    site_sequences = load_public_predmat_input_site_sequences()
    canonical_components = [
        _canonical_public_site_components(site_id) for site_id in input_phospho.index
    ]
    phospho = input_phospho.copy(deep=True)
    phospho.index = pd.Index(
        [site_id for site_id, _, _ in canonical_components],
        name=input_phospho.index.name,
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": [gene_symbol for _, gene_symbol, _ in canonical_components],
            "site": [site for _, _, site in canonical_components],
            "site_sequence": [
                str(site_sequences[str(site_id).strip()])
                for site_id in input_phospho.index.astype(str)
            ],
        },
        index=phospho.index.copy(),
    )
    substrate_map = load_public_predmat_input_substrate_map()
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            [
                {
                    "kinase": str(kinase),
                    "substrate_site": _canonical_public_site_components(site_id)[0],
                }
                for kinase, site_ids in substrate_map.items()
                for site_id in site_ids
            ]
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    str(sequence) for _, sequence in site_sequences.items()
                ]
            },
            index=pd.Index(
                [
                    _canonical_public_site_components(site_id)[0]
                    for site_id, _ in site_sequences.items()
                ],
                name="site_id",
            ),
        ),
    )
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(min_substrates=2),
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
            scoring_config=KinaseScoringConfig(min_substrates=2),
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
