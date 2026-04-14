#!/usr/bin/env python3
"""Advanced native workflow demo.

Use this lane when you want direct control over workflow-shaped inputs such as
substrate_map, motif_sequences, site_sequences, and intermediate scoring or
prediction outputs. For the shorter common end-to-end path, use
examples/simple_workflow_demo.py and SimpleKinaseWorkflow instead.
"""

from __future__ import annotations

import pandas as pd

from phospy.api import KinaseWorkflow, PredictionRunConfig
from phospy.api.workflow_results import KinaseWorkflowResult


def run_demo() -> KinaseWorkflowResult:
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 1.1, 0.9, 1.2, 3.0, 2.9, 3.1, 2.8],
            "sample_2": [2.0, 2.1, 1.9, 2.2, 2.0, 2.1, 1.9, 2.2],
            "sample_3": [3.0, 3.1, 2.9, 3.2, 1.0, 1.1, 0.9, 1.2],
        },
        index=[f"SITE_{i}" for i in range(1, 9)],
    )
    substrate_map = {
        "KINASE_A": ["SITE_1", "SITE_2", "SITE_3", "SITE_4"],
        "KINASE_B": ["SITE_5", "SITE_6", "SITE_7", "SITE_8"],
    }
    site_sequences = {
        "SITE_1": "QQAAAAAYY",
        "SITE_2": "QQAAAAAYY",
        "SITE_3": "QQAAAAAYY",
        "SITE_4": "QQAAAAAYY",
        "SITE_5": "QQTTTTTYY",
        "SITE_6": "QQTTTTTYY",
        "SITE_7": "QQTTTTTYY",
        "SITE_8": "QQTTTTTYY",
    }
    motif_sequences = {
        "KINASE_A": ["QQAAAAAYY", "QQAAAAAYY", "QQAAAAAYY"],
        "KINASE_B": ["QQTTTTTYY", "QQTTTTTYY", "QQTTTTTYY"],
    }

    workflow = KinaseWorkflow(flank_size=2)
    result = workflow.run(
        phospho_matrix=phospho_matrix,
        substrate_map=substrate_map,
        site_sequences=site_sequences,
        motif_sequences=motif_sequences,
        prediction_config=PredictionRunConfig(
            min_substrates=2,
            min_motif_size=2,
            ensemble_size=3,
            top=4,
            score_threshold=0.75,
            inclusion=3,
            n_iterations=2,
            random_state=17,
        ),
    )
    return result


def main() -> None:
    result = run_demo()
    print("Prediction matrix")
    print(result.prediction_result.pred_matrix.round(4))
    print()
    print("Top predicted substrate counts")
    print(
        {
            kinase: len(sites)
            for kinase, sites in result.prediction_result.substrate_list.items()
        }
    )


if __name__ == "__main__":
    main()
