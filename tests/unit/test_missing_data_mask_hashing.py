from __future__ import annotations

import pandas as pd

from phospy.science.datasets.preprocessing.missing_data_mask_hashing import (
    hash_group_aware_imputation_mask,
    hash_group_aware_mechanism_mask,
)


def test_group_aware_mask_hashes_ignore_inferred_string_axis_dtype() -> None:
    values = [[True, False], [False, True]]
    object_axis_mask = pd.DataFrame(
        values,
        index=pd.Index(["site_a", "site_b"], dtype="object", name="site"),
        columns=pd.Index(["sample_a", "sample_b"], dtype="object", name="sample"),
    )
    string_axis_mask = pd.DataFrame(
        values,
        index=pd.Index(["site_a", "site_b"], dtype="string", name="site"),
        columns=pd.Index(["sample_a", "sample_b"], dtype="string", name="sample"),
    )

    assert hash_group_aware_imputation_mask(
        object_axis_mask
    ) == hash_group_aware_imputation_mask(string_axis_mask)
    assert hash_group_aware_mechanism_mask(
        object_axis_mask,
        mechanism="knn",
        mask_kind="target",
    ) == hash_group_aware_mechanism_mask(
        string_axis_mask,
        mechanism="knn",
        mask_kind="target",
    )
