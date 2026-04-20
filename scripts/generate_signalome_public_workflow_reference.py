#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
DEFAULT_OUTPUT_DIR = ROOT / "tests" / "fixtures" / "public_workflow_reference"

NUMERIC_RTOL = 1e-9
NUMERIC_ATOL = 1e-12


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate rewrite-owned full-output signalome reference fixtures for "
            "the supported L6 public workflow parity lane."
        )
    )
    parser.add_argument(
        "--outdir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory where signalome_rewrite_l6_*.csv fixtures are written.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.outdir)
    contract_path = output_dir / "signalome_rewrite_l6_contract.json"
    output_dir.mkdir(parents=True, exist_ok=True)

    if str(SRC) not in sys.path:
        sys.path.insert(0, str(SRC))
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))

    from phospy import (
        KinasePredictionConfig,
        KinaseScoringConfig,
        KinaseWorkflow,
        KinaseWorkflowRequest,
        ReferencePreset,
        SignalomeConfig,
        SignalomeWorkflow,
        SignalomeWorkflowRequest,
    )
    from tests.support.rewrite_fixture_data import build_rat_l6_dataset

    dataset = build_rat_l6_dataset(n_sites=260)
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(top_k=6, ensemble_size=12),
            activity_config=None,
        )
    )
    signalome = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(substrate_support_cutoff=0.5),
        )
    )

    assignments = signalome.module_assignments.table.copy(deep=True)
    modules = signalome.signalome_modules.table.copy(deep=True)
    nodes = signalome.kinase_network.nodes.copy(deep=True)
    edges = signalome.kinase_network.edges.copy(deep=True)
    expanded = signalome.expanded_signalome.copy(deep=True)

    assignments.to_csv(output_dir / "signalome_rewrite_l6_module_assignments.csv")
    modules.to_csv(output_dir / "signalome_rewrite_l6_modules.csv")
    nodes.to_csv(output_dir / "signalome_rewrite_l6_network_nodes.csv")
    edges.to_csv(output_dir / "signalome_rewrite_l6_network_edges.csv", index=False)
    expanded.to_csv(
        output_dir / "signalome_rewrite_l6_expanded_signalome.csv",
        index=False,
    )

    if contract_path.exists():
        contract = json.loads(contract_path.read_text(encoding="utf-8"))
    else:
        contract = {}
    contract.setdefault("generation_date", "2026-04-20")
    contract.update(
        {
            "fixture_set_id": "signalome_rewrite_l6_supported_lane_v2_full_outputs",
            "donor_source": "tests/fixtures/rewrite_parity/r_reference_l6/l6_phospho_matrix.csv + bundled ReferencePreset.AUTO (rat/l6_native)",
            "generation_path": "scripts/generate_signalome_public_workflow_reference.py",
            "supported_outputs": [
                "signalome_modules",
                "module_assignments.table",
                "kinase_network.nodes",
                "kinase_network.edges",
                "expanded_signalome",
            ],
            "n_assignments": int(assignments.shape[0]),
            "n_modules": int(modules.shape[0]),
            "n_module_kinases": int(modules.shape[1]),
            "n_nodes": int(nodes.shape[0]),
            "n_edges": int(edges.shape[0]),
            "n_expanded_rows": int(expanded.shape[0]),
            "n_expanded_kinases": int(expanded.loc[:, "kinase"].nunique()),
            "expanded_row_kind_counts": expanded.loc[:, "row_kind"]
            .value_counts()
            .sort_index()
            .astype("int64")
            .to_dict(),
            "ambiguous_assignment_count": int(
                assignments.loc[:, "top_kinase_is_ambiguous"].sum()
            ),
            "module_id_counts": assignments.loc[:, "module_id"]
            .value_counts()
            .sort_index()
            .astype("int64")
            .to_dict(),
            "tie_count_distribution": assignments.loc[:, "top_kinase_tie_count"]
            .value_counts()
            .sort_index()
            .astype("int64")
            .to_dict(),
            "positive_edge_count": int((edges.loc[:, "correlation"] > 0.0).sum()),
            "negative_edge_count": int((edges.loc[:, "correlation"] < 0.0).sum()),
            "comparison_rules": {
                "signalome_modules": {
                    "row_order": "fixture_order_semantic",
                    "column_order": "fixture_order_semantic",
                    "comparison": "exact_equality",
                },
                "module_assignments.table": {
                    "row_order": "site_id_sorted",
                    "column_order": "fixture_order_semantic",
                    "comparison": "numeric_tolerance",
                    "numeric_rtol": NUMERIC_RTOL,
                    "numeric_atol": NUMERIC_ATOL,
                    "normalization": "canonicalize_collection_columns",
                },
                "kinase_network.nodes": {
                    "row_order": "kinase_sorted",
                    "column_order": "fixture_order_semantic",
                    "comparison": "exact_equality",
                },
                "kinase_network.edges": {
                    "row_order": "source_target_sorted",
                    "column_order": "fixture_order_semantic",
                    "comparison": "numeric_tolerance",
                    "numeric_rtol": NUMERIC_RTOL,
                    "numeric_atol": NUMERIC_ATOL,
                },
                "expanded_signalome": {
                    "row_order": "kinase_row_kind_site_sorted",
                    "column_order": "fixture_order_semantic",
                    "comparison": "numeric_tolerance",
                    "numeric_rtol": NUMERIC_RTOL,
                    "numeric_atol": NUMERIC_ATOL,
                },
            },
            "normalization_rules": {
                "module_assignments.table": [
                    "cast index to site_id strings",
                    "canonicalize collection-like columns (top_kinase_candidates, top_kinase_weights, module_top_kinase_candidates)",
                    "cast module and tie-count columns to int64",
                    "cast ambiguity flags to bool",
                    "sort rows by site_id",
                ],
                "kinase_network.nodes": [
                    "cast kinase index to strings",
                    "cast degree and n_substrates to int64",
                    "sort rows by kinase",
                ],
                "kinase_network.edges": [
                    "cast source_kinase and target_kinase to strings",
                    "cast correlation to float",
                    "sort rows by source_kinase,target_kinase",
                ],
                "expanded_signalome": [
                    "cast stable schema columns to expected dtypes",
                    "sort rows by kinase,row_kind,site_id,site_order,module_id,protein_id,top_kinase",
                ],
            },
            "stability_rationale": "Outputs are generated from the supported rewrite public lane with fixed fixture input and deterministic workflow settings; parity locks all supported downstream signalome tables.",
            "assignment_policy": "cutoff_binary",
            "network_policy": "signed",
        }
    )
    contract_path.write_text(
        json.dumps(contract, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
