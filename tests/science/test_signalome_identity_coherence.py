from __future__ import annotations

import pandas as pd

from phospy.science.signalomes.constants import DISPLAY_ID_COLUMN, SITE_KEY_COLUMN
from phospy.science.signalomes.context import build_site_membership_table
from phospy.science.signalomes.science import build_module_assignments
from tests.support.workflow_identity_coherence import (
    DUPLICATE_DISPLAY_ID,
    duplicate_display_site_index,
    duplicate_display_site_metadata,
)


def test_signalome_science_tables_keep_duplicate_display_ids_separate() -> None:
    site_index = duplicate_display_site_index()
    site_metadata = duplicate_display_site_metadata(site_index)
    prediction_matrix = pd.DataFrame(
        {"K1": [0.8, 0.7]},
        index=site_index.copy(),
    )
    site_to_protein = pd.Series(
        ["P28482", "Q99999"],
        index=site_index.copy(),
        dtype=object,
        name="protein_group_id",
    )

    module_assignments = build_module_assignments(
        prediction_matrix=prediction_matrix,
        site_to_protein_group_id=site_to_protein,
        site_metadata=site_metadata,
    )

    assert module_assignments.index.astype(str).tolist() == (
        site_index.astype(str).tolist()
    )
    assert module_assignments.loc[:, SITE_KEY_COLUMN].is_unique
    assert int(module_assignments.loc[:, SITE_KEY_COLUMN].nunique()) == 2
    assert int(module_assignments.loc[:, DISPLAY_ID_COLUMN].nunique()) == 1
    assert module_assignments.loc[:, DISPLAY_ID_COLUMN].astype(str).tolist() == [
        DUPLICATE_DISPLAY_ID,
        DUPLICATE_DISPLAY_ID,
    ]

    site_membership = build_site_membership_table(
        module_assignments=module_assignments,
        site_clusters=pd.Series(
            [1, 2],
            index=site_index.copy(),
            dtype="int64",
            name="site_cluster",
        ),
        site_metadata=site_metadata,
        prediction_matrix=prediction_matrix,
        kinase_substrates={"K1": tuple(site_index.astype(str).tolist())},
        substrate_support_cutoff=0.5,
        assignment_policy="cutoff_binary",
    )

    assert site_membership.loc[:, SITE_KEY_COLUMN].astype(str).tolist() == (
        site_index.astype(str).tolist()
    )
    assert site_membership.loc[:, SITE_KEY_COLUMN].is_unique
    assert int(site_membership.loc[:, SITE_KEY_COLUMN].nunique()) == 2
    assert int(site_membership.loc[:, DISPLAY_ID_COLUMN].nunique()) == 1
