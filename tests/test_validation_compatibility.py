from __future__ import annotations

import pandas as pd
import pytest

from phospy import KinaseActivityAnalyzer
from phospy.validation.errors import InputCompatibilityError
from phospy.workflow import KinaseWorkflow


def test_kinase_activity_analyzer_rejects_zero_overlap() -> None:
    analyzer = KinaseActivityAnalyzer(
        pd.DataFrame(
            {
                "PRKACA": [0.9],
            },
            index=["SITE_A"],
        )
    )
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0],
        },
        index=["SITE_B"],
    )

    with pytest.raises(InputCompatibilityError, match="no overlapping phosphosite IDs"):
        analyzer.analyze(phospho_matrix)


def test_kinase_workflow_rejects_missing_site_sequence_coverage() -> None:
    workflow = KinaseWorkflow()
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 2.0],
        },
        index=["SITE_1", "SITE_2"],
    )

    with pytest.raises(
        InputCompatibilityError, match="site_sequences is missing entries"
    ):
        workflow.run(
            phospho_matrix=phospho_matrix,
            substrate_map={"KINASE_A": ["SITE_1", "SITE_2"]},
            site_sequences={"SITE_1": "QQAAAAAYY"},
            motif_sequences={"KINASE_A": ["QQAAAAAYY", "QQAAAAAYY"]},
            min_substrates=1,
            min_motif_size=1,
            ensemble_size=2,
            top=2,
            score_threshold=0.5,
            inclusion=1,
            n_iterations=1,
        )


def test_kinase_workflow_rejects_zero_substrate_overlap() -> None:
    workflow = KinaseWorkflow()
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0],
        },
        index=["SITE_1"],
    )

    with pytest.raises(
        InputCompatibilityError,
        match="no overlap between substrate_map and phospho_matrix",
    ):
        workflow.run(
            phospho_matrix=phospho_matrix,
            substrate_map={"KINASE_A": ["SITE_9"]},
            motif_sequences=None,
            allow_profile_only_fallback=True,
            ensemble_size=2,
            top=2,
            score_threshold=0.5,
            inclusion=1,
            n_iterations=1,
        )
